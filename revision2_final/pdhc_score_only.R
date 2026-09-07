# PD vs HC CI-score comparison (item 5), NA-safe v2
ROOT <- "E:/PPMI帕金森数据库专用"
OUT  <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订")
FD   <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data")
ld  <- read.csv(file.path(OUT, "PC1_loadings_66genes.csv"), stringsAsFactors = FALSE)
map <- read.csv(file.path(FD, "de_full.csv"), stringsAsFactors = FALSE)[, c("ENSG","Symbol")]
map$ENSG <- sub("\\..*", "", map$ENSG)
ld$ENSG <- map$ENSG[match(ld$Gene, map$Symbol)]
bl <- read.csv(file.path(ROOT, "gene_expression_data_bl.csv"), check.names = FALSE)
names(bl) <- sub("\\..*$", "", names(bl))
meta <- read.csv(file.path(ROOT, "metaDataIR3.csv"), check.names = FALSE, stringsAsFactors = FALSE)
diag <- meta[!duplicated(meta$PATNO), c("PATNO","DIAGNOSIS")]
bl$DIAG <- diag$DIAGNOSIS[match(bl$PATNO, diag$PATNO)]
cat("rows with NA DIAG:", sum(is.na(bl$DIAG)), " PATNOs:", bl$PATNO[is.na(bl$DIAG)], "\n")
bl <- bl[!is.na(bl$DIAG), ]
use <- ld$Gene[ld$ENSG %in% names(bl)]
M2 <- bl[, ld$ENSG[ld$ENSG %in% names(bl)]]; names(M2) <- use
sdv <- apply(M2, 2, sd, na.rm = TRUE)
z2 <- sweep(sweep(M2, 2, colMeans(M2, na.rm = TRUE), "-"), 2, sdv, "/")
sc <- as.numeric(as.matrix(z2) %*% ld$PC1_loading[match(use, ld$Gene)])
gp <- bl$DIAG
pd_s <- sc[gp == "PD"]; hc_s <- sc[gp == "Control"]
cat("PD n =", length(pd_s), " HC n =", length(hc_s), "\n")
wt <- wilcox.test(pd_s, hc_s); tt <- t.test(pd_s, hc_s)
d_pool <- (mean(pd_s) - mean(hc_s)) / sqrt((var(pd_s)*(length(pd_s)-1) + var(hc_s)*(length(hc_s)-1)) / (length(pd_s)+length(hc_s)-2))
cat(sprintf("\nPD (n=%d): mean = %.3f, sd = %.3f\n", length(pd_s), mean(pd_s), sd(pd_s)))
cat(sprintf("HC (n=%d): mean = %.3f, sd = %.3f\n", length(hc_s), mean(hc_s), sd(hc_s)))
cat(sprintf("difference = %.3f, Cohen's d = %.3f\n", mean(pd_s)-mean(hc_s), d_pool))
cat(sprintf("Wilcoxon p = %.4g | Welch t p = %.4g\n", wt$p.value, tt$p.value))
cat("positive score = high CI expression (orientation validated vs published PC1)\n")
