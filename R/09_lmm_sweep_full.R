# ============================================================
# COMPREHENSIVE LMM sweep: ~43 longitudinal indicators x PC1 group
# Totals + DAT-SBR imaging + MoCA subdomains + UPDRS subdomains/items + SCOPA domains
# Factor-time (last-visit interaction) + continuous-time, FDR across ALL
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})
setwd("E:/PPMI帕金森数据库专用")
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
VIS <- c("SC","BL","V04","V06","V08","V10","V12")
vmap <- c(BL=0,SC=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# Loaders (cached reads)
rd <- function(f) suppressWarnings(read_csv(f, show_col_types=FALSE)) %>% filter(EVENT_ID %in% VIS)
p1c <- rd("运动症状数据/MDS-UPDRS_Part_I_31Jan2026.csv")
p1p <- rd("运动症状数据/MDS-UPDRS_Part_I_Patient_Questionnaire_31Jan2026.csv")
p2  <- rd("运动症状数据/MDS_UPDRS_Part_II__Patient_Questionnaire_31Jan2026.csv")
p3  <- rd("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv") %>% arrange(PATNO,EVENT_ID) %>% distinct(PATNO,EVENT_ID,.keep_all=TRUE)
p4  <- rd("运动症状数据/MDS-UPDRS_Part_IV__Motor_Complications_29Jan2026.csv")
moca<- rd("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv")
rb  <- rd("REM_Sleep_Behavior_Disorder_Questionnaire_03Feb2026.csv")
sca <- rd("SCOPA-AUT_07Feb2026.csv")
dat <- rd("运动症状数据/Xing_Core_Lab_-_Quant_SBR_23Feb2026.csv")

SS <- function(d, cols) rowSums(sapply(intersect(cols,colnames(d)), function(c) as.numeric(d[[c]])), na.rm=TRUE)
mk <- function(d, score, nm) data.frame(PATNO=d$PATNO, EVENT_ID=d$EVENT_ID, Score=score, Outcome=nm)

L <- list()
# A. Totals
L[["UPDRS_I_clin"]] <- mk(p1c, p1c$NP1RTOT, "UPDRS_I_clin")
L[["UPDRS_I_pat"]]  <- mk(p1p, p1p$NP1PTOT, "UPDRS_I_pat")
L[["UPDRS_II"]]     <- mk(p2,  p2$NP2PTOT,  "UPDRS_II")
L[["UPDRS_III"]]    <- mk(p3,  p3$NP3TOT,   "UPDRS_III")
L[["UPDRS_IV"]]     <- mk(p4,  p4$NP4TOT,   "UPDRS_IV")
L[["MoCA"]]         <- mk(moca,moca$MCATOT, "MoCA")
L[["NHY"]]          <- mk(p3,  p3$NHY,      "NHY")
L[["PIGD"]]         <- mk(p3, SS(p3,c("NP3GAIT","NP3PSTBL","NP3FRZGT")), "PIGD")
L[["Tremor"]]       <- mk(p3, SS(p3,c("NP3PTRMR","NP3PTRML","NP3KTRMR","NP3KTRML","NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL")), "Tremor")
L[["RestTremor"]]   <- mk(p3, SS(p3,c("NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL")), "RestTremor")
it12 <- c("DRMVIVID","DRMAGRAC","DRMNOCTB","SLPLMBMV","SLPINJUR","DRMVERBL","DRMFIGHT","DRMUMV","DRMOBJFL","MVAWAKEN","DRMREMEM","SLPDSTRB")
neuro <- c("STROKE","HETRA","PARKISM","RLS","NARCLPSY","DEPRS","EPILEPSY","BRNINFM","CNSOTH")
L[["RBDSQ"]] <- mk(rb, SS(rb,it12) + as.integer(SS(rb,neuro)>0), "RBDSQ")
sca_na <- sca; for(c in paste0("SCAU",1:21)) if(c %in% colnames(sca_na)) sca_na[[c]][sca_na[[c]]==9] <- NA
L[["SCOPA_total"]] <- mk(sca_na, SS(sca_na, paste0("SCAU",1:21)), "SCOPA_total")

# B. DAT-SBR imaging (lower=worse; SC/V04/V06/V10, no V12)
L[["DAT_striatum"]] <- mk(dat, as.numeric(dat$STRIATUM_REF_CWM), "DAT_striatum")
L[["DAT_caudate"]]  <- mk(dat, as.numeric(dat$CAUDATE_REF_CWM),  "DAT_caudate")
L[["DAT_putamen"]]  <- mk(dat, as.numeric(dat$PUTAMEN_REF_CWM),  "DAT_putamen")

# C. MoCA subdomains
L[["MoCA_visuospatial"]] <- mk(moca, SS(moca,c("MCAALTTM","MCACUBE","MCACLCKC","MCACLCKN","MCACLCKH")), "MoCA_visuospatial")
L[["MoCA_naming"]]       <- mk(moca, SS(moca,c("MCALION","MCARHINO","MCACAMEL")), "MoCA_naming")
L[["MoCA_attention"]]    <- mk(moca, SS(moca,c("MCAFDS","MCABDS","MCAVIGIL","MCASER7")), "MoCA_attention")
L[["MoCA_language"]]     <- mk(moca, SS(moca,c("MCASNTNC","MCAVF")), "MoCA_language")
L[["MoCA_abstraction"]]  <- mk(moca, SS(moca,c("MCAABSTR")), "MoCA_abstraction")
L[["MoCA_recall"]]       <- mk(moca, SS(moca,c("MCAREC1","MCAREC2","MCAREC3","MCAREC4","MCAREC5")), "MoCA_recall")
L[["MoCA_orientation"]]  <- mk(moca, SS(moca,c("MCADATE","MCAMONTH","MCAYR","MCADAY","MCAPLACE","MCACITY")), "MoCA_orientation")

# D. UPDRS III motor subdomains
L[["UPDRS3_rigidity"]]    <- mk(p3, SS(p3,c("NP3RIGN","NP3RIGRU","NP3RIGLU","NP3RIGRL","NP3RIGLL")), "UPDRS3_rigidity")
L[["UPDRS3_bradykinesia"]]<- mk(p3, SS(p3,c("NP3FTAPR","NP3FTAPL","NP3HMOVR","NP3HMOVL","NP3PRSPR","NP3PRSPL","NP3TTAPR","NP3TTAPL","NP3LGAGR","NP3LGAGL")), "UPDRS3_bradykinesia")
L[["UPDRS3_axial"]]       <- mk(p3, SS(p3,c("NP3SPCH","NP3GAIT","NP3FRZGT","NP3PSTBL")), "UPDRS3_axial")

# E. UPDRS I non-motor individual items
for(it in c("NP1COG","NP1HALL","NP1DPRS","NP1ANXS","NP1APAT","NP1DDS")) L[[it]] <- mk(p1c, as.numeric(p1c[[it]]), it)
for(it in c("NP1SLPN","NP1SLPD","NP1PAIN","NP1URIN","NP1CNST","NP1LTHD","NP1FATG")) L[[it]] <- mk(p1p, as.numeric(p1p[[it]]), it)

# F. SCOPA subdomains
L[["SCOPA_GI"]]       <- mk(sca_na, SS(sca_na, paste0("SCAU",1:7)),   "SCOPA_GI")
L[["SCOPA_urinary"]]  <- mk(sca_na, SS(sca_na, paste0("SCAU",8:13)),  "SCOPA_urinary")
L[["SCOPA_cardio"]]   <- mk(sca_na, SS(sca_na, paste0("SCAU",14:16)), "SCOPA_cardio")
L[["SCOPA_thermo"]]   <- mk(sca_na, SS(sca_na, c("SCAU17","SCAU18","SCAU20","SCAU21")), "SCOPA_thermo")
L[["SCOPA_pupillo"]]  <- mk(sca_na, SS(sca_na, c("SCAU19")), "SCOPA_pupillo")

cat(sprintf("Total indicators to test: %d\n", length(L)))

# ---- Run LMM ----
res <- data.frame()
for(nm in names(L)) {
  d <- L[[nm]] %>% inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & !is.na(d$Score), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT) < 3 || nrow(d) < 50) { cat(sprintf("SKIP %s\n", nm)); next }
  lastv <- tail(levels(d$VISIT), 1)

  lastp<-NA; lastest<-NA; minp<-NA
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)) { cm<-coef(summary(mf)); ir<-grep(":groupLow$", rownames(cm))
    if(length(ir)>0){ minp<-min(cm[ir,"Pr(>|t|)"],na.rm=TRUE)
      lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
      if(length(lr)>0){lastp<-cm[lr[1],"Pr(>|t|)"];lastest<-cm[lr[1],"Estimate"]} } }
  grpp<-NA; intp<-NA; intest<-NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){ cm<-coef(summary(mc)); gr<-grep("^groupLow$",rownames(cm)); ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0) grpp<-cm[gr[1],"Pr(>|t|)"]
    if(length(ir)>0){intp<-cm[ir[1],"Pr(>|t|)"];intest<-cm[ir[1],"Estimate"]} }
  res <- rbind(res, data.frame(Outcome=nm, N_obs=nrow(d), LastVisit=lastv,
    LastV_int_est=lastest, LastV_int_p=lastp, MinVisit_int_p=minp,
    Cont_grp_p=grpp, Cont_int_est=intest, Cont_int_p=intp, stringsAsFactors=FALSE))
}
res$LastV_int_padj <- p.adjust(res$LastV_int_p, "BH")
res$Cont_int_padj  <- p.adjust(res$Cont_int_p, "BH")
res$Cont_grp_padj  <- p.adjust(res$Cont_grp_p, "BH")
res <- res[order(res$LastV_int_p),]

cat(sprintf("\n========== COMPREHENSIVE LMM SWEEP (%d indicators) ==========\n", nrow(res)))
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)
cat("\n=== Last-visit interaction nominal p<0.05 ===\n")
print(res[!is.na(res$LastV_int_p)&res$LastV_int_p<0.05, c("Outcome","LastVisit","LastV_int_est","LastV_int_p","LastV_int_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (last-visit padj<0.05): %s\n", paste(res$Outcome[!is.na(res$LastV_int_padj)&res$LastV_int_padj<0.05],collapse=", ")))
cat("\n=== Continuous interaction nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_int_p)&res$Cont_int_p<0.05, c("Outcome","Cont_int_est","Cont_int_p","Cont_int_padj")], row.names=FALSE)
cat("\n=== Continuous group main effect nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_grp_p)&res$Cont_grp_p<0.05, c("Outcome","Cont_grp_p","Cont_grp_padj")], row.names=FALSE)
write.csv(res, file.path(OUT,"LMM_sweep_COMPREHENSIVE.csv"), row.names=FALSE)
cat(sprintf("\nSaved LMM_sweep_COMPREHENSIVE.csv (%d indicators)\n", nrow(res)))
