# Reviewer item #4: validate deconvolution estimates against measured CBC differential
# Measured: PPMI hematology BL labs (Neutrophils/Lymphocytes/Monocytes %)
# Estimated: CIBERSORT fractions (immune_fractions.csv, 393 PD)
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

ROOT <- "E:/PPMI帕金森数据库专用"
OUT  <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订")
FD   <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data")

## ---- measured CBC (screening) ----
lab <- read.csv(file.path(ROOT, "知识库文件/Blood_Chemistry___Hematology-Archived_05Apr2026.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cbc <- lab[lab$LTSTNAME %in% c("Neutrophils (%)","Lymphocytes (%)","Monocytes (%)") & lab$EVENT_ID == "SC", ]
cbc$VAL <- suppressWarnings(as.numeric(cbc$LSIRES))
cbc <- cbc[!is.na(cbc$VAL), ]
cbc_w <- aggregate(VAL ~ PATNO + LTSTNAME, data = cbc, FUN = mean)
cbc_w <- tidyr::pivot_wider(cbc_w, names_from = LTSTNAME, values_from = VAL)
cat("cbc_w rows:", nrow(cbc_w), " | PATNO class:", class(cbc_w$PATNO), " | first:", head(cbc_w$PATNO, 3), "\n")
cat("measured CBC BL: n subjects =", nrow(cbc_w), "\n")
print(summary(cbc_w[, 2:4]))

## ---- estimated fractions ----
ifr <- read.csv(file.path(FD, "immune_fractions.csv"), check.names = FALSE)
pd_data <- read.csv(file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/data/PD_all_clustering_methods.csv"), stringsAsFactors = FALSE)
ifr$PATNO <- pd_data$PATNO[match(ifr$SAMPLE_ID, pd_data$SAMPLE_ID)]
ifr$Neut <- ifr[["Neutrophils"]]
ifr$Monoc <- ifr[["Monocytes"]]
lymph_cols <- grep("T cells|B cells|NK cells|Plasma", names(ifr), value = TRUE)
ifr$Lymph <- rowSums(ifr[, lymph_cols])

m <- merge(cbc_w, ifr[, c("PATNO","Neut","Lymph","Monoc")], by = "PATNO")
cat("after merge:", nrow(m), " | ifr PATNO head:", head(ifr$PATNO, 3), "\n")
m <- m[complete.cases(m[, c("Neutrophils (%)","Lymphocytes (%)","Monocytes (%)","Neut","Lymph","Monoc")]), ]
cat("\nmatched subjects with measured + estimated:", nrow(m), "\n")

comp <- function(m_est, m_meas, lab) {
  pct <- m_meas / 100
  r_p <- cor.test(m_est, pct, method = "pearson")
  r_s <- cor.test(m_est, pct, method = "spearman", exact = FALSE)
  cat(sprintf("%-11s est mean=%.3f meas mean=%.3f | Pearson r=%.3f (p=%.2g) | Spearman rho=%.3f (p=%.2g) | n=%d\n",
              lab, mean(m_est), mean(pct), r_p$estimate, r_p$p.value, r_s$estimate, r_s$p.value, length(m_est)))
}

cat("\n=== deconvolution vs measured CBC (screening) ===\n")
comp(m$Neut,  m[["Neutrophils (%)"]], "Neutrophil")
comp(m$Lymph, m[["Lymphocytes (%)"]], "Lymphocyte")
comp(m$Monoc, m[["Monocytes (%)"]],   "Monocyte")

## ---- also by CI group (is the low-CI neutrophil signal visible in measured CBC?) ----
grp <- ifr$group[match(m$PATNO, ifr$PATNO)]
cat("\n=== measured Neutrophils (%) by CI group (only matched subjects) ===\n")
print(tapply(m[["Neutrophils (%)"]]/100, grp, mean))
print(tapply(m[["Neutrophils (%)"]]/100, grp, sd))
w <- wilcox.test((m[["Neutrophils (%)"]]/100) ~ grp)
cat(sprintf("Wilcoxon p (measured neutrophils, High vs Low) = %.4g  n=%d\n", w$p.value, sum(!is.na(grp))))
w2 <- wilcox.test(m$Neut ~ grp)
cat(sprintf("Wilcoxon p (estimated neutrophil fraction)    = %.4g\n", w2$p.value))

write.csv(m, file.path(OUT, "cbc_vs_deconvolution.csv"), row.names = FALSE)
cat("DONE\n")
