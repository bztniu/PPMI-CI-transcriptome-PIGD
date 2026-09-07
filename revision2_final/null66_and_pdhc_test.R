# Reviewer items #3 (random gene-set null) and #5 (PD vs HC score comparison)
# Null: 1000 random 66-gene sets -> same score pipeline (log2 norm counts, z per gene,
#       PC1 oriented by colSums) -> same primary LMM (PIGD ~ VISIT*score_z + (1|PATNO))
#       -> empirical p for real PC1's V12 interaction (beta=-0.1831, p=0.015).
suppressPackageStartupMessages({library(DESeq2); library(lme4); library(lmerTest)})

ROOT <- "E:/PPMI帕金森数据库专用"
OUT  <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订")
TAB  <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables")
set.seed(20260831)
N_ITER <- 1000

## ---- shared inputs ----
dds <- readRDS(file.path(ROOT, "dds_pd_control_BL_object.rds"))
pd_data <- read.csv(file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/data/PD_all_clustering_methods.csv"), stringsAsFactors = FALSE)
pd_samples <- intersect(colnames(dds), pd_data$SAMPLE_ID)
cat("dds samples:", ncol(dds), " | PD samples used:", length(pd_samples), "\n")

map_df <- read.csv(file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/de_full.csv"), stringsAsFactors = FALSE)[, c("ENSG","Symbol")]
map_df$ENSG <- sub("\\..*", "", map_df$ENSG)
adj <- read.csv(file.path(TAB, "DE_covariate_adjusted_full.csv"), check.names = FALSE)
names(adj)[1] <- "ENSG"; adj$ENSG <- sub("\\..*", "", adj$ENSG)
ld <- read.csv(file.path(OUT, "PC1_loadings_66genes.csv"), stringsAsFactors = FALSE)
universe <- setdiff(adj$ENSG[!is.na(adj$stat)], map_df$ENSG[map_df$Symbol %in% ld$Gene])
cat("random-sampling universe (non-scoring, DE-filtered):", length(universe), "\n")

# normalized counts, log2 — once
cnt <- log2(counts(dds, normalized = TRUE) + 1)
ensg_all <- sub("\\..*", "", rownames(cnt))
# index lookup
idx_univ <- which(ensg_all %in% universe)
cat("universe rows found in matrix:", length(idx_univ), "\n")

# real PC1 validation numbers (from published score)
pc1s <- read.csv(file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/pc1_scores.csv"), stringsAsFactors = FALSE)
samp2pat <- pd_data$PATNO[match(pd_samples, pd_data$SAMPLE_ID)]

# longitudinal data
long <- read.csv(file.path(TAB, "longitudinal_PIGD_dedup.csv"), stringsAsFactors = FALSE)
long$VISIT <- factor(long$VISIT, levels = c("BL","V04","V06","V08","V10","V12"))

fit_v12 <- function(score_by_sid) {
  d <- long
  d$score <- score_by_sid[match(d$PATNO, samp2pat)]
  d <- d[!is.na(d$score) & !is.na(d$PIGD), ]
  d$score_z <- as.numeric(scale(d$score))
  m <- try(lmer(PIGD ~ VISIT * score_z + (1|PATNO), data = d, REML = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(c(NA, NA))
  cm <- coef(summary(m)); i <- grep("VISITV12:score_z", rownames(cm))
  c(cm[i,"Estimate"], cm[i,"Pr(>|t|)"])
}

## ---- real PC1 (recompute to confirm identical pipeline) ----
ci_df <- read.csv(file.path(ROOT, "complex/CI_genes_converted.csv"), stringsAsFactors = FALSE)
ci_ensg <- intersect(ci_df$original_id, ensg_all)
M <- cnt[match(ci_ensg, ensg_all), pd_samples, drop = FALSE]
z <- t(scale(t(M)))
pr <- prcomp(t(z), center = TRUE, scale. = TRUE)
real_pc1 <- pr$x[, "PC1"]
if (cor(real_pc1, colSums(M)) < 0) real_pc1 <- -real_pc1
cat(sprintf("recomputed real PC1 vs published: r = %.4f\n",
            cor(real_pc1, pc1s$PC1[match(pd_samples, pc1s$SAMPLE_ID)])))
real_fit <- fit_v12(setNames(real_pc1, pd_samples))
cat(sprintf("REAL PC1 V12 interaction: beta = %.4f, p = %.5f\n", real_fit[1], real_fit[2]))
BETA_REAL <- real_fit[1]

## ---- null loop ----
res <- matrix(NA_real_, N_ITER, 2, dimnames = list(NULL, c("beta","p")))
t0 <- Sys.time()
for (i in seq_len(N_ITER)) {
  pick <- sample(idx_univ, 66)
  M <- cnt[pick, pd_samples, drop = FALSE]
  z <- t(scale(t(M)))
  pr <- prcomp(t(z), center = TRUE, scale. = TRUE)
  s <- pr$x[, 1]
  if (cor(s, colSums(M)) < 0) s <- -s
  res[i, ] <- fit_v12(setNames(s, pd_samples))
  if (i %% 200 == 0) cat(sprintf("  iter %d/%d  elapsed %.1f min\n", i, N_ITER,
                                 as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
ok <- complete.cases(res)
cat("\nnull runs completed:", sum(ok), "/", N_ITER, "\n")
cat(sprintf("random |beta| quantiles: 50%%=%.4f  90%%=%.4f  95%%=%.4f  99%%=%.4f  max=%.4f\n",
            quantile(abs(res[ok,1]), .5), quantile(abs(res[ok,1]), .9),
            quantile(abs(res[ok,1]), .95), quantile(abs(res[ok,1]), .99), max(abs(res[ok,1]))))
n_beyond <- sum(abs(res[ok,1]) >= abs(BETA_REAL))
emp_p <- (n_beyond + 1) / (sum(ok) + 1)
cat(sprintf("random sets with |beta| >= |%.4f|: %d  -> empirical p = %.4f\n", BETA_REAL, n_beyond, emp_p))
cat(sprintf("random sets with nominal p<0.05: %d (%.1f%%, expect ~5%%)\n",
            sum(res[ok,2] < 0.05, na.rm = TRUE), 100*mean(res[ok,2] < 0.05, na.rm = TRUE)))
write.csv(data.frame(iter = which(ok), res[ok, ]), file.path(OUT, "random66_null_results.csv"), row.names = FALSE)

## ---- item 5: PD vs HC (bl matrix, same loadings, z across all 216) ----
bl <- read.csv(file.path(ROOT, "gene_expression_data_bl.csv"), check.names = FALSE)
names(bl) <- sub("\\..*$", "", names(bl))
meta <- read.csv(file.path(ROOT, "metaDataIR3.csv"), check.names = FALSE, stringsAsFactors = FALSE)
diag <- meta[!duplicated(meta$PATNO), c("PATNO","DIAGNOSIS")]
bl$DIAG <- diag$DIAGNOSIS[match(bl$PATNO, diag$PATNO)]
ldm <- read.csv(file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/de_full.csv"), stringsAsFactors = FALSE)[, c("ENSG","Symbol")]
ldm$ENSG <- sub("\\..*", "", ldm$ENSG)
ld$ENSG <- ldm$ENSG[match(ld$Gene, ldm$Symbol)]
gcols <- ld$Gene[match(intersect(ld$ENSG, names(bl)), ld$ENSG)]
gcols <- ld$Gene[ld$ENSG %in% names(bl)]
M2 <- bl[, ld$ENSG[ld$ENSG %in% names(bl)]]
names(M2) <- gcols
z2 <- sweep(sweep(M2, 2, colMeans(M2), "-"), 2, apply(M2, 2, sd), "/")
sc <- as.numeric(as.matrix(z2) %*% ld$PC1_loading[match(gcols, ld$Gene)])
gp <- bl$DIAG
pd_s <- sc[gp == "PD"]; hc_s <- sc[gp == "Control"]
wt <- wilcox.test(pd_s, hc_s)
d_pool <- (mean(pd_s) - mean(hc_s)) / sqrt((var(pd_s)*(length(pd_s)-1) + var(hc_s)*(length(hc_s)-1)) / (length(pd_s)+length(hc_s)-2))
cat("\n=== PD vs HC CI score (bl matrix; z across all 216 samples) ===\n")
cat(sprintf("PD  (n=%d): mean = %.3f, sd = %.3f\n", length(pd_s), mean(pd_s), sd(pd_s)))
cat(sprintf("HC  (n=%d): mean = %.3f, sd = %.3f\n", length(hc_s), mean(hc_s), sd(hc_s)))
cat(sprintf("difference = %.3f, Cohen's d = %.3f, Wilcoxon p = %.4g\n", mean(pd_s)-mean(hc_s), d_pool, wt$p.value))
cat("DONE\n")
