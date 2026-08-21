# ============================================================
# LMM sweep: PC1 group x time interaction across ALL longitudinal scales
# Factor-time (V12 interaction) + continuous-time, FDR corrected
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

# PC1 groups (PATNO -> group)
pc1 <- read.csv(file.path(FD,"pc1_scores.csv"))[,c("PATNO","group")]
VIS <- c("SC","BL","V04","V06","V08","V10","V12")
vmap <- c(BL=0,SC=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# ---- helper to load a longitudinal total ----
get_long <- function(file, scorefun, name) {
  d <- suppressWarnings(read_csv(file, show_col_types=FALSE))
  d <- d %>% filter(EVENT_ID %in% VIS)
  d$Score <- scorefun(d)
  d %>% dplyr::select(PATNO, EVENT_ID, Score) %>% filter(!is.na(Score)) %>%
    mutate(Outcome=name)
}

scales <- list()
scales[["UPDRS_I_clin"]]  <- get_long("运动症状数据/MDS-UPDRS_Part_I_31Jan2026.csv", function(d) d$NP1RTOT, "UPDRS_I_clin")
scales[["UPDRS_I_pat"]]   <- get_long("运动症状数据/MDS-UPDRS_Part_I_Patient_Questionnaire_31Jan2026.csv", function(d) d$NP1PTOT, "UPDRS_I_pat")
scales[["UPDRS_II"]]      <- get_long("运动症状数据/MDS_UPDRS_Part_II__Patient_Questionnaire_31Jan2026.csv", function(d) d$NP2PTOT, "UPDRS_II")
scales[["UPDRS_IV"]]      <- get_long("运动症状数据/MDS-UPDRS_Part_IV__Motor_Complications_29Jan2026.csv", function(d) d$NP4TOT, "UPDRS_IV")
scales[["MoCA"]]          <- get_long("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", function(d) d$MCATOT, "MoCA")

# Part III: total, NHY, PIGD, Tremor, RestTremor
p3 <- suppressWarnings(read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% VIS) %>% arrange(PATNO, EVENT_ID) %>% distinct(PATNO, EVENT_ID, .keep_all=TRUE)
mkp3 <- function(col) p3 %>% mutate(Score=.data[[col]]) %>% dplyr::select(PATNO,EVENT_ID,Score) %>% filter(!is.na(Score)) %>% mutate(Outcome=col)
scales[["UPDRS_III"]] <- mkp3("NP3TOT"); scales[["UPDRS_III"]]$Outcome <- "UPDRS_III"
scales[["NHY"]]       <- mkp3("NHY");    scales[["NHY"]]$Outcome <- "NHY"
p3 <- p3 %>% mutate(PIGD=NP3GAIT+NP3PSTBL+NP3FRZGT,
                    Tremor=NP3PTRMR+NP3PTRML+NP3KTRMR+NP3KTRML+NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL,
                    RestTremor=NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL)
scales[["PIGD"]]       <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=PIGD) %>% filter(!is.na(Score)) %>% mutate(Outcome="PIGD")
scales[["Tremor"]]     <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=Tremor) %>% filter(!is.na(Score)) %>% mutate(Outcome="Tremor")
scales[["RestTremor"]] <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=RestTremor) %>% filter(!is.na(Score)) %>% mutate(Outcome="RestTremor")

# RBDSQ total = sum(12 items) + (any neuro disorder ? 1 : 0)
rb <- suppressWarnings(read_csv("REM_Sleep_Behavior_Disorder_Questionnaire_03Feb2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% VIS)
it12 <- c("DRMVIVID","DRMAGRAC","DRMNOCTB","SLPLMBMV","SLPINJUR","DRMVERBL","DRMFIGHT","DRMUMV","DRMOBJFL","MVAWAKEN","DRMREMEM","SLPDSTRB")
neuro <- c("STROKE","HETRA","PARKISM","RLS","NARCLPSY","DEPRS","EPILEPSY","BRNINFM","CNSOTH")
rb$Score <- rowSums(rb[,it12], na.rm=TRUE) + as.integer(rowSums(rb[,neuro], na.rm=TRUE) > 0)
scales[["RBDSQ"]] <- rb %>% dplyr::select(PATNO,EVENT_ID,Score) %>% mutate(Outcome="RBDSQ")

# SCOPA-AUT total = sum SCAU1..SCAU21 (proxy)
sc <- suppressWarnings(read_csv("SCOPA-AUT_07Feb2026.csv", show_col_types=FALSE)) %>% filter(EVENT_ID %in% VIS)
scau_cols <- paste0("SCAU", 1:21); scau_cols <- intersect(scau_cols, colnames(sc))
sc_num <- sc[,scau_cols]; sc_num[sc_num==9] <- NA   # 9 = N/A
sc$Score <- rowSums(sapply(sc_num, as.numeric), na.rm=TRUE)
scales[["SCOPA_AUT"]] <- sc %>% dplyr::select(PATNO,EVENT_ID,Score) %>% mutate(Outcome="SCOPA_AUT")

# ---- Run LMM for each scale ----
res <- data.frame()
for(nm in names(scales)) {
  d <- scales[[nm]] %>% inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & !is.na(d$Score), ]
  if(nlevels(droplevels(d$VISIT)) < 3 || nrow(d) < 50) { cat(sprintf("SKIP %s (insufficient)\n", nm)); next }

  # Factor-time RI-LMM: V12 interaction + min interaction p across visits
  v12p <- NA; v12est <- NA; minp <- NA
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)) {
    cm <- coef(summary(mf))
    int_rows <- grep(":groupLow$", rownames(cm))
    if(length(int_rows)>0) {
      ips <- cm[int_rows,"Pr(>|t|)"]; minp <- min(ips, na.rm=TRUE)
      v12r <- grep("VISITV12:groupLow", rownames(cm), fixed=TRUE)
      if(length(v12r)>0) { v12p <- cm[v12r[1],"Pr(>|t|)"]; v12est <- cm[v12r[1],"Estimate"] }
    }
  }
  # Continuous-time RS-LMM: group main + time:group interaction
  grp_p <- NA; int_p <- NA; int_est <- NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)) {
    cm <- coef(summary(mc))
    gr <- grep("^groupLow$", rownames(cm)); ir <- grep("Time:groupLow", rownames(cm), fixed=TRUE)
    if(length(gr)>0) grp_p <- cm[gr[1],"Pr(>|t|)"]
    if(length(ir)>0) { int_p <- cm[ir[1],"Pr(>|t|)"]; int_est <- cm[ir[1],"Estimate"] }
  }
  res <- rbind(res, data.frame(Outcome=nm, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)),
    V12_int_est=v12est, V12_int_p=v12p, MinVisit_int_p=minp,
    Cont_grp_p=grp_p, Cont_int_est=int_est, Cont_int_p=int_p, stringsAsFactors=FALSE))
}

# FDR correction across outcomes (on V12 interaction p)
res$V12_int_padj <- p.adjust(res$V12_int_p, "BH")
res$Cont_int_padj <- p.adjust(res$Cont_int_p, "BH")
res <- res[order(res$V12_int_p),]

cat("\n========== LMM SWEEP RESULTS (sorted by V12 interaction p) ==========\n")
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)

cat("\n=== V12 interaction: nominal p<0.05 ===\n")
print(res[!is.na(res$V12_int_p) & res$V12_int_p<0.05, c("Outcome","V12_int_est","V12_int_p","V12_int_padj")], row.names=FALSE)
cat(sprintf("\nSurvive FDR (V12 padj<0.05): %s\n",
            paste(res$Outcome[!is.na(res$V12_int_padj) & res$V12_int_padj<0.05], collapse=", ")))
cat("\n=== Continuous group main effect: nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_grp_p) & res$Cont_grp_p<0.05, c("Outcome","Cont_grp_p")], row.names=FALSE)

write.csv(res, file.path(OUT,"LMM_sweep_all_outcomes.csv"), row.names=FALSE)
cat(sprintf("\nSaved: %s/LMM_sweep_all_outcomes.csv\n", OUT))
