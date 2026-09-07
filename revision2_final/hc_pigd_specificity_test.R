# HC specificity test: project PD-derived PC1 loadings onto healthy controls,
# then ask whether high/low CI-score HC differ in baseline PIGD (gait+postural+freezing).
# Cross-sectional only (HC does not progress), covariates: Age + Sex.
suppressPackageStartupMessages({library(dplyr)})

ROOT <- "E:/PPMI帕金森数据库专用"
OUT  <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订")
FD   <- file.path(ROOT, "结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data")

## ---- 1) load pieces ----
ld  <- read.csv(file.path(OUT, "PC1_loadings_66genes.csv"), stringsAsFactors = FALSE)
map <- read.csv(file.path(FD, "de_full.csv"), stringsAsFactors = FALSE)[, c("ENSG","Symbol")]
map$ENSG <- sub("\\..*", "", map$ENSG)
ld$ENSG <- map$ENSG[match(ld$Gene, map$Symbol)]
stopifnot(sum(is.na(ld$ENSG)) == 0)

bl <- read.csv(file.path(ROOT, "gene_expression_data_bl.csv"), check.names = FALSE)
names(bl) <- sub("\\..*$", "", names(bl))  # strip ENSG version suffix
cat("bl matrix:", nrow(bl), "samples x", ncol(bl) - 2, "columns\n")

# 66 gene columns in bl matrix (ENSG)
ensg_cols <- intersect(ld$ENSG, names(bl))
cat("66-gene ENSG columns found in bl matrix:", length(ensg_cols), "\n")
stopifnot(length(ensg_cols) == 66)

# long -> wide by symbol
expr66 <- bl[, c("PATNO","EVENT_ID", ensg_cols)]
names(expr66) <- c("PATNO","EVENT_ID", ld$Gene[match(ensg_cols, ld$ENSG)])

## ---- 2) PD mean/SD centering basis (same matrix, self-consistent) ----
meta <- read.csv(file.path(ROOT, "metaDataIR3.csv"), check.names = FALSE, stringsAsFactors = FALSE)
diag <- meta[!duplicated(meta$PATNO), c("PATNO","DIAGNOSIS","GENDER")]
expr66$DIAG <- diag$DIAGNOSIS[match(expr66$PATNO, diag$PATNO)]
pd_rows <- expr66$DIAG == "PD"
genes <- ld$Gene[match(names(expr66)[3:68], ld$Gene)]  # order loadings to columns
ld_o <- ld[match(names(expr66)[3:68], ld$Gene), ]
mu  <- colMeans(expr66[pd_rows, 3:68], na.rm = TRUE)
sdv <- apply(expr66[pd_rows, 3:68], 2, sd, na.rm = TRUE)
zmat <- sweep(sweep(expr66[, 3:68], 2, mu, "-"), 2, sdv, "/")
expr66$score <- as.numeric(as.matrix(zmat) %*% ld_o$PC1_loading)  # = PC1 projection

## ---- 3) sign validation against published PD PC1 ----
ce <- read.csv(file.path(FD, "ci_expr.csv"), row.names = 1, check.names = FALSE)
ce_z <- sweep(sweep(ce, 2, colMeans(ce), "-"), 2, apply(ce, 2, sd), "/")
pd_proj <- as.numeric(as.matrix(ce_z) %*% ld$PC1_loading)
pc1s <- read.csv(file.path(FD, "pc1_scores.csv"), stringsAsFactors = FALSE)
chk <- merge(pc1s, data.frame(SAMPLE_ID = rownames(ce), proj = pd_proj), by = "SAMPLE_ID")
r_val <- cor(chk$PC1, chk$proj)
cat(sprintf("validation: r(published PC1, projection) = %.4f  (n=%d)\n", r_val, nrow(chk)))
if (r_val < 0) { expr66$score <- -expr66$score; cat("sign flipped to match published PC1\n") }

## ---- 4) HC subset + baseline PIGD ----
hc <- expr66[expr66$DIAG == "Control", c("PATNO","score")]
p3 <- read.csv(file.path(ROOT, "知识库文件/MDS-UPDRS_Part_III_05Apr2026.csv"), check.names = FALSE)
p3bl <- p3[p3$EVENT_ID == "BL", ]
p3bl$PIGD <- p3bl$NP3GAIT + p3bl$NP3PSTBL + p3bl$NP3FRZGT
hc <- merge(hc, p3bl[, c("PATNO","PIGD")], by = "PATNO", all.x = TRUE)
ag <- read.csv(file.path(ROOT, "知识库文件/Age_at_visit_05Apr2026.csv"), check.names = FALSE)
agbl <- ag[ag$EVENT_ID == "BL", c("PATNO","AGE_AT_VISIT")]
hc <- merge(hc, agbl, by = "PATNO", all.x = TRUE)
hc$Sex <- ifelse(diag$GENDER[match(hc$PATNO, diag$PATNO)] == "Female", "F", "M")
hc <- hc[!is.na(hc$PIGD) & !is.na(hc$score), ]
cat("HC with score + PIGD + age:", nrow(hc), "\n")

## ---- 5) stats ----
med <- median(hc$score)
hc$grp <- ifelse(hc$score <= med, "Low-CI", "High-CI")
cat("\nPIGD distribution (HC, BL): median =", median(hc$PIGD),
    " IQR =", quantile(hc$PIGD, .25), "-", quantile(hc$PIGD, .75),
    " max =", max(hc$PIGD), " % >0 =", round(100*mean(hc$PIGD > 0), 1), "\n")
cat("Age range:", min(hc$AGE_AT_VISIT), "-", max(hc$AGE_AT_VISIT),
    " Sex: F =", sum(hc$Sex == "F"), " M =", sum(hc$Sex == "M"), "\n")

cat("\n=== PIGD by CI-score group (HC, cross-sectional) ===\n")
bygrp <- hc %>% group_by(grp) %>%
  summarise(n = n(), mean = mean(PIGD), sd = sd(PIGD), median = median(PIGD))
print(as.data.frame(bygrp))
w <- wilcox.test(PIGD ~ grp, data = hc)
cat(sprintf("Wilcoxon rank-sum p = %.4f\n", w$p.value))

sp <- cor.test(hc$score, hc$PIGD, method = "spearman", exact = FALSE)
cat(sprintf("Spearman rho(score, PIGD) = %.3f  p = %.4f\n", sp$estimate, sp$p.value))

hc$score_z <- as.numeric(scale(hc$score))
hc$Age_z   <- as.numeric(scale(hc$AGE_AT_VISIT))
m1 <- lm(PIGD ~ score_z + Age_z + Sex, data = hc)
s1 <- coef(summary(m1))
cat("\n=== OLS: PIGD ~ score_z + Age + Sex ===\n")
print(s1)

## ---- 6) sensitivity: PD-only check (does score->PIGD exist cross-sectionally at BL?) ----
pd <- expr66[expr66$DIAG == "PD", c("PATNO","score")]
pdp <- merge(pd, p3bl[, c("PATNO","PIGD")], by = "PATNO", all.x = TRUE)
pdp <- pdp[!is.na(pdp$PIGD), ]
spd <- cor.test(pdp$score, pdp$PIGD, method = "spearman", exact = FALSE)
cat(sprintf("\n[reference] PD BL cross-sectional Spearman rho = %.3f  p = %.2g  (n=%d)\n",
            spd$estimate, spd$p.value, nrow(pdp)))

write.csv(hc[, c("PATNO","score","grp","PIGD","AGE_AT_VISIT","Sex")],
          file.path(OUT, "hc_ci_score_pigd.csv"), row.names = FALSE)
cat("DONE\n")
