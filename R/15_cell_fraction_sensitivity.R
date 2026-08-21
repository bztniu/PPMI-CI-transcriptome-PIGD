# Part 2 (standalone, no dds): PIGD year-5 LMM cell-fraction sensitivity
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(dplyr); library(readr); library(tidyr) })
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
TAB <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"

pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))
imf  <- read.csv(file.path(FD,"immune_fractions.csv"), check.names=FALSE)
cf <- imf[, c("SAMPLE_ID","Neutrophils","T cells CD8","T cells CD4 memory resting",
              "T cells CD4 naive","B cells naive","NK cells resting","Monocytes")]
colnames(cf) <- c("SAMPLE_ID","Neu","CD8","CD4mr","CD4n","Bn","NKr","Mono")
cf2 <- merge(cf, pc1s[,c("SAMPLE_ID","PATNO")], by="SAMPLE_ID")

m3 <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("BL","V04","V06","V08","V10","V12")) %>%
  mutate(PIGD = NP3GAIT+NP3PSTBL+NP3FRZGT) %>%
  dplyr::select(PATNO, EVENT_ID, PIGD) %>% filter(!is.na(PIGD)) %>%
  group_by(PATNO, EVENT_ID) %>% summarise(PIGD=mean(PIGD), .groups="drop")   # dedupe

ld <- m3 %>% inner_join(pc1s[,c("PATNO","PC1","group")], by="PATNO") %>%
  left_join(cf2[,c("PATNO","Neu","CD8","CD4mr","CD4n","Bn","NKr","Mono")], by="PATNO")
ld$visit <- factor(ld$EVENT_ID, levels=c("BL","V04","V06","V08","V10","V12"))
ld$group <- factor(ld$group, levels=c("Low","High"))
ld$PC1z  <- scale(ld$PC1)[,1]
ld <- ld[complete.cases(ld[,c("Neu","CD8","CD4mr")]),]
cat(sprintf("LMM long data: %d rows, %d patients\n", nrow(ld), length(unique(ld$PATNO))))

getV12 <- function(m){ s<-summary(m)$coefficients
  r <- rownames(s)[grepl("V12", rownames(s)) & grepl(":", rownames(s))]
  s[r, c("Estimate","Pr(>|t|)"), drop=FALSE] }

getV12_full <- function(m, model, adjustment, direction) {
  s <- summary(m)$coefficients
  r <- rownames(s)[grepl("V12", rownames(s)) & grepl(":", rownames(s))]
  out <- as.data.frame(s[r, , drop=FALSE], check.names=FALSE)
  out$Term <- rownames(out)
  out$Model <- model
  out$Adjustment <- adjustment
  out$Outcome <- "PIGD"
  out$Visit <- "V12"
  out$Estimate_lowCI_direction <- direction * out$Estimate
  out$N_obs <- nobs(m)
  out$N_patients <- length(unique(ld$PATNO))
  out[, c("Outcome","Visit","Model","Adjustment","Term","Estimate",
          "Estimate_lowCI_direction","Std. Error","df","t value","Pr(>|t|)",
          "N_obs","N_patients")]
}

m_bin_raw <- lmer(PIGD~group*visit+(1|PATNO),ld)
m_bin_adj <- lmer(PIGD~group*visit+Neu+CD8+CD4mr+CD4n+Bn+NKr+Mono+(1|PATNO),ld)
m_pc1_raw <- lmer(PIGD~PC1z*visit+(1|PATNO),ld)
m_pc1_adj <- lmer(PIGD~PC1z*visit+Neu+CD8+CD4mr+CD4n+Bn+NKr+Mono+(1|PATNO),ld)

cat("\n(a) binary group×visit UNADJUSTED:\n");      print(round(getV12(m_bin_raw),4))
cat("\n(b) binary group×visit CELL-ADJUSTED:\n");   print(round(getV12(m_bin_adj),4))
cat("\n(c) continuous PC1z×visit UNADJUSTED:\n");    print(round(getV12(m_pc1_raw),4))
cat("\n(d) continuous PC1z×visit CELL-ADJUSTED:\n"); print(round(getV12(m_pc1_adj),4))

cell_sens <- bind_rows(
  getV12_full(m_bin_raw, "Binary PC1 median group x visit", "Unadjusted", -1),
  getV12_full(m_bin_adj, "Binary PC1 median group x visit", "Cell-fraction adjusted", -1),
  getV12_full(m_pc1_raw, "Continuous PC1_z x visit", "Unadjusted", -1),
  getV12_full(m_pc1_adj, "Continuous PC1_z x visit", "Cell-fraction adjusted", -1)
) %>%
  mutate(
    Effect_interpretation = ifelse(
      Model == "Continuous PC1_z x visit",
      "Additional year-5 delta PIGD per 1-SD lower PC1",
      "Additional year-5 delta PIGD in low-CI vs high-CI"
    ),
    Cell_covariates = ifelse(
      Adjustment == "Cell-fraction adjusted",
      "Neutrophils, CD8 T cells, CD4 memory resting T cells, CD4 naive T cells, naive B cells, resting NK cells, monocytes",
      "None"
    )
  )
write.csv(cell_sens, file.path(TAB, "PIGD_cell_composition_sensitivity.csv"), row.names=FALSE)
cat(sprintf("Saved %s\n", file.path(TAB, "PIGD_cell_composition_sensitivity.csv")))
