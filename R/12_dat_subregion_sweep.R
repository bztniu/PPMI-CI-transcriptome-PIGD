# DAT-SBR sub-region longitudinal LMM sweep x PC1 group
suppressPackageStartupMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(emmeans)})
setwd("E:/PPMI帕金森数据库专用")
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
PACKAGE_TABLES <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
vmap <- c(SC=0,BL=0,V04=1,V06=2,V08=3,V10=4,V12=5)

dat <- suppressWarnings(read_csv("运动症状数据/Xing_Core_Lab_-_Quant_SBR_23Feb2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% c("SC","BL","V04","V06","V08","V10","V12"))

# all summary sub-regions (bilateral mean)
regions <- c("STRIATUM_REF_CWM","CAUDATE_REF_CWM","PUTAMEN_REF_CWM","PRECAUDATE_REF_CWM",
  "POSCAUDATE_REF_CWM","PRECOMMISSURAL_PUTAMEN_REF_CWM","POSCOMMISSURAL_PUTAMEN_REF_CWM",
  "PREDORSALPUTAMEN_REF_CWM","PREVENTRALPUTAMEN_REF_CWM","POSDORSALPUTAMEN_REF_CWM","POSVENTRALPUTAMEN_REF_CWM")
# asymmetry indices (|L-R|/mean) for caudate & putamen
mkasym <- function(L,R) ifelse((L+R)>0, abs(L-R)/((L+R)/2), NA)
dat$PUTAMEN_asym  <- mkasym(as.numeric(dat$PUTAMEN_L_REF_CWM),  as.numeric(dat$PUTAMEN_R_REF_CWM))
dat$CAUDATE_asym  <- mkasym(as.numeric(dat$CAUDATE_L_REF_CWM),  as.numeric(dat$CAUDATE_R_REF_CWM))
dat$STRIATUM_asym <- mkasym(as.numeric(dat$STRIATUM_L_REF_CWM), as.numeric(dat$STRIATUM_R_REF_CWM))
# lowest-side putamen (most affected hemisphere) — clinically key
dat$PUTAMEN_min <- pmin(as.numeric(dat$PUTAMEN_L_REF_CWM), as.numeric(dat$PUTAMEN_R_REF_CWM), na.rm=TRUE)

outcomes <- c(regions, "PUTAMEN_asym","CAUDATE_asym","STRIATUM_asym","PUTAMEN_min")

res <- data.frame()
for(rg in outcomes) {
  d <- dat %>% mutate(Score=as.numeric(.data[[rg]])) %>%
    dplyr::select(PATNO, EVENT_ID, Score) %>% filter(!is.na(Score)) %>%
    inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT)<3 || nrow(d)<80) { cat(sprintf("SKIP %s\n", rg)); next }
  lastv <- tail(levels(d$VISIT),1)
  grpp<-NA;grpe<-NA;intp<-NA;inte<-NA;lastp<-NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){cm<-coef(summary(mc));gr<-grep("^groupLow$",rownames(cm));ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0){grpe<-cm[gr[1],"Estimate"];grpp<-cm[gr[1],"Pr(>|t|)"]}
    if(length(ir)>0){inte<-cm[ir[1],"Estimate"];intp<-cm[ir[1],"Pr(>|t|)"]}}
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)){cm<-coef(summary(mf));lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
    if(length(lr)>0) lastp<-cm[lr[1],"Pr(>|t|)"]}
  res <- rbind(res, data.frame(Region=rg, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)), LastVisit=lastv,
    Grp_est=grpe, Grp_p=grpp, Int_est=inte, Int_p=intp, LastV_int_p=lastp, stringsAsFactors=FALSE))
}
res$Grp_padj <- p.adjust(res$Grp_p,"BH"); res$Int_padj <- p.adjust(res$Int_p,"BH")
res <- res[order(res$Grp_p),]
cat(sprintf("\n===== DAT-SBR SUB-REGION SWEEP (%d regions, follow-up to %s) =====\n", nrow(res), paste(unique(res$LastVisit),collapse="/")))
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)
cat("\n=== Group main effect nominal p<0.05 ===\n")
print(res[!is.na(res$Grp_p)&res$Grp_p<0.05,c("Region","Grp_est","Grp_p","Grp_padj")], row.names=FALSE)
cat(sprintf("SURVIVE FDR (group): %s\n", paste(res$Region[!is.na(res$Grp_padj)&res$Grp_padj<0.05],collapse=", ")))
cat("\n=== Interaction nominal p<0.05 ===\n")
print(res[!is.na(res$Int_p)&res$Int_p<0.05,c("Region","Int_est","Int_p","Int_padj")], row.names=FALSE)
cat(sprintf("SURVIVE FDR (interaction): %s\n", paste(res$Region[!is.na(res$Int_padj)&res$Int_padj<0.05],collapse=", ")))
write.csv(res, file.path(OUT,"DAT_subregion_sweep.csv"), row.names=FALSE)
write.csv(res, file.path(PACKAGE_TABLES,"DAT_subregion_sweep.csv"), row.names=FALSE)
cat("\nSaved DAT_subregion_sweep.csv\n")

# Model-estimated trajectory for the clinically representative posterior dorsal putamen.
# The full 15-region screen remains in DAT_subregion_sweep.csv.
representative_region <- "POSDORSALPUTAMEN_REF_CWM"
d_traj <- dat %>%
  mutate(Score=as.numeric(.data[[representative_region]])) %>%
  dplyr::select(PATNO, EVENT_ID, Score) %>%
  filter(!is.na(Score)) %>%
  inner_join(pc1, by="PATNO") %>%
  mutate(
    PATNO=factor(PATNO),
    group=factor(group, levels=c("High","Low")),
    VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID),
      levels=c("BL","V04","V06","V08","V10","V12"))
  ) %>%
  filter(!is.na(VISIT)) %>%
  droplevels()

m_traj <- lmer(Score ~ VISIT*group + (1|PATNO), data=d_traj, REML=FALSE)
emm_traj <- as.data.frame(emmeans(m_traj, ~ group | VISIT))
n_traj <- d_traj %>%
  count(VISIT, group, name="N")
screen_row <- res %>%
  filter(Region == representative_region)

emm_traj <- emm_traj %>%
  left_join(n_traj, by=c("VISIT","group")) %>%
  mutate(
    Region=representative_region,
    Baseline_group_p=screen_row$Grp_p[[1]],
    Linear_group_by_year_p=screen_row$Int_p[[1]],
    Last_visit_interaction_p=screen_row$LastV_int_p[[1]],
    N_subjects=nlevels(d_traj$PATNO)
  )

write.csv(emm_traj,
  file.path(PACKAGE_TABLES,"DAT_posterior_dorsal_putamen_trajectory.csv"),
  row.names=FALSE)
cat("Saved DAT_posterior_dorsal_putamen_trajectory.csv\n")
