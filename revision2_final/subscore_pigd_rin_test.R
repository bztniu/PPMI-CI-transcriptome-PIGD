# Decisive test: nuclear-59 vs mtDNA-7 subscores
# (a) which subscore drives the year-5 PIGD interaction
# (b) is the mtDNA signal a RIN-degradation artifact? (subscore ~ RIN)
suppressPackageStartupMessages({library(lme4); library(lmerTest)})

FD   <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
TAB  <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"
OUT  <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订"

# subscores (per PATNO, from morning run)
sub_sc <- read.csv(file.path(OUT, "nuclear_mtDNA_subscores.csv"), stringsAsFactors = FALSE)

# longitudinal PIGD (deduplicated, same as manuscript primary)
long <- read.csv(file.path(TAB, "longitudinal_PIGD_dedup.csv"), stringsAsFactors = FALSE)
long$VISIT <- factor(long$VISIT, levels = c("BL","V04","V06","V08","V10","V12"))
m <- merge(long, sub_sc[, c("PATNO","nuclear59_meanZ","mtDNA7_meanZ")], by = "PATNO", all.x = TRUE)
m$nuc_z <- as.numeric(scale(m$nuclear59_meanZ))
m$mt_z  <- as.numeric(scale(m$mtDNA7_meanZ))
pc1 <- read.csv(file.path(FD, "pc1_scores.csv"), stringsAsFactors = FALSE)
pc1$PC1z <- as.numeric(scale(pc1$PC1))
m$PC1z <- pc1$PC1z[match(m$PATNO, pc1$PATNO)]

cat("merged rows:", nrow(m), " | subjects:", length(unique(m$PATNO)), "\n")

d <- m[!is.na(m$PIGD), ]
cat("PIGD rows:", nrow(d), "\n\n")

# --- validation: reproduce primary model with 66-gene groups ---
m0 <- lmer(PIGD ~ VISIT * group + (1 | PATNO), data = d, REML = TRUE)
cm <- coef(summary(m0))
i <- grep("VISITV12:groupLow", rownames(cm))
cat("=== VALIDATION: 66-gene group model ===\n")
cat(sprintf("V12:groupLow  beta=%.4f SE=%.4f p=%.5f\n", cm[i,"Estimate"], cm[i,"Std. Error"], cm[i,"Pr(>|t|)"]))

# --- continuous subscore models (V12 interaction) ---
fit <- function(form, dat, lab) {
  mm <- lmer(as.formula(form), data = dat, REML = TRUE)
  cm <- coef(summary(mm))
  r <- grep("VISITV12:", rownames(cm), value = TRUE)
  for (term in r) {
    cat(sprintf("%-12s %-28s beta=%.4f SE=%.4f p=%.5f\n", lab, term,
                cm[term,"Estimate"], cm[term,"Std. Error"], cm[term,"Pr(>|t|)"]))
  }
  invisible(mm)
}
cat("\n=== V12 interaction by score (continuous) ===\n")
fit("PIGD ~ VISIT * nuc_z + (1 | PATNO)", d, "nuclear59")
fit("PIGD ~ VISIT * mt_z  + (1 | PATNO)", d, "mtDNA7")
fit("PIGD ~ VISIT * PC1z  + (1 | PATNO)", d, "PC1_66")
# both in one model (competing slopes)
cat("\n=== both subscores in one model ===\n")
fit("PIGD ~ VISIT * nuc_z + VISIT * mt_z + (1 | PATNO)", d, "both")

# --- RIN artifact test ---
meta <- read.csv("E:/PPMI帕金森数据库专用/metaDataIR3.csv", check.names = FALSE)
rin <- data.frame(SAMPLE_ID = as.character(meta[["Specimen Bar Code"]]),
                  RIN = as.numeric(meta[["RIN Value"]]), stringsAsFactors = FALSE)
pc1s <- read.csv(file.path(FD, "pc1_scores.csv"), stringsAsFactors = FALSE)
pc1s$SAMPLE_ID <- as.character(pc1s$SAMPLE_ID)
sm <- merge(pc1s, rin, by = "SAMPLE_ID")
r_pc1_rin <- with(sm, cor.test(PC1, RIN))
cat(sprintf("r(PC1, RIN)            = %.3f  p=%.3g\n", r_pc1_rin$estimate, r_pc1_rin$p.value))
sm <- merge(sm, sub_sc[, c("PATNO","nuclear59_meanZ","mtDNA7_meanZ")], by = "PATNO", all.x = TRUE)
sm <- sm[!is.na(sm$RIN), ]
cat("\n=== RIN artifact test (n =", nrow(sm), "with RIN) ===\n")
cat(sprintf("r(mtDNA7_meanZ, RIN)   = %.3f  p=%.3g\n", with(sm, cor.test(mtDNA7_meanZ, RIN))$estimate, with(sm, cor.test(mtDNA7_meanZ, RIN))$p.value))
cat(sprintf("r(nuclear59_meanZ, RIN)= %.3f  p=%.3g\n", with(sm, cor.test(nuclear59_meanZ, RIN))$estimate, with(sm, cor.test(nuclear59_meanZ, RIN))$p.value))
cat(sprintf("r(PC1, RIN)            = %.3f  p=%.3g\n", with(sm, cor.test(PC1, RIN))$estimate, with(sm, cor.test(PC1, RIN))$p.value))
# partial: mtDNA ~ RIN residualized on nuclear?
sm$mt_resid <- residuals(lm(mtDNA7_meanZ ~ nuclear59_meanZ, data = sm, na.action = na.exclude))
cat(sprintf("r(mtDNA resid on nuclear, RIN) = %.3f  p=%.3g\n", with(sm, cor.test(mt_resid, RIN))$estimate, with(sm, cor.test(mt_resid, RIN))$p.value))
cat("DONE\n")
