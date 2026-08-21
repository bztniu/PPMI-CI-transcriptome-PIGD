# Biomarker LMM sweep: blood/CSF biomarkers x PC1 group (main effect + interaction)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lme4); library(lmerTest)})
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
bm <- read.csv(file.path(OUT,"data/biomarkers_long.csv"))
vmap <- c(BL=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# collapse duplicate measurements (same PATNO/visit/biomarker) by mean
bm <- bm %>% group_by(PATNO, EVENT_ID, biomarker) %>%
  summarise(value=mean(value, na.rm=TRUE), .groups="drop")

# keep biomarkers with adequate longitudinal coverage
qc <- bm %>% group_by(biomarker) %>%
  summarise(n_obs=n(), n_subj=n_distinct(PATNO),
            n_visits=n_distinct(EVENT_ID), .groups="drop") %>%
  filter(n_visits>=3, n_obs>=80, n_subj>=40)
cat(sprintf("Biomarkers passing QC: %d\n", nrow(qc)))

res <- data.frame()
for(b in qc$biomarker) {
  d <- bm %>% filter(biomarker==b) %>% inner_join(pc1, by="PATNO")
  if(nrow(d)<80 || n_distinct(d$group)<2) next
  # log transform if all positive & skewed
  if(all(d$value>0, na.rm=TRUE)) d$value <- log(d$value)
  d <- d %>% mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
                    Time=vmap[EVENT_ID],
                    VISIT=factor(EVENT_ID, levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & is.finite(d$value), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT)<3 || nrow(d)<80) next
  lastv <- tail(levels(d$VISIT),1)

  # continuous-time: group main + time:group
  grpp<-NA; grpe<-NA; intp<-NA; inte<-NA
  mc <- tryCatch(lmer(value ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(value ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){ cm<-coef(summary(mc)); gr<-grep("^groupLow$",rownames(cm)); ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0){grpp<-cm[gr[1],"Pr(>|t|)"]; grpe<-cm[gr[1],"Estimate"]}
    if(length(ir)>0){intp<-cm[ir[1],"Pr(>|t|)"]; inte<-cm[ir[1],"Estimate"]} }
  # factor-time last visit interaction
  lastp<-NA
  mf <- tryCatch(lmer(value ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)){ cm<-coef(summary(mf)); lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
    if(length(lr)>0) lastp<-cm[lr[1],"Pr(>|t|)"] }
  res <- rbind(res, data.frame(Biomarker=b, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)),
    LastVisit=lastv, Grp_est=grpe, Grp_p=grpp, Int_est=inte, Int_p=intp, LastV_int_p=lastp, stringsAsFactors=FALSE))
}
res$Grp_padj <- p.adjust(res$Grp_p,"BH")
res$Int_padj <- p.adjust(res$Int_p,"BH")
res <- res[order(res$Grp_p),]

cat(sprintf("\n========== BIOMARKER LMM SWEEP (%d biomarkers) ==========\n", nrow(res)))
cat("\n=== GROUP MAIN EFFECT: nominal p<0.05 (low vs high CI overall) ===\n")
print(res[!is.na(res$Grp_p)&res$Grp_p<0.05, c("Biomarker","N_subj","Grp_est","Grp_p","Grp_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (group main padj<0.05): %s\n", paste(res$Biomarker[!is.na(res$Grp_padj)&res$Grp_padj<0.05],collapse=", ")))
cat("\n=== TIME x GROUP INTERACTION: nominal p<0.05 ===\n")
print(res[!is.na(res$Int_p)&res$Int_p<0.05, c("Biomarker","Int_est","Int_p","Int_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (interaction padj<0.05): %s\n", paste(res$Biomarker[!is.na(res$Int_padj)&res$Int_padj<0.05],collapse=", ")))
write.csv(res, file.path(OUT,"biomarker_LMM_sweep.csv"), row.names=FALSE)
cat(sprintf("\nSaved biomarker_LMM_sweep.csv (%d biomarkers)\n", nrow(res)))
