# Real ssGSEA: CI 66-gene set enrichment on full transcriptome (393 PD samples)
suppressPackageStartupMessages({ library(DESeq2); library(GSVA) })
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"

dds <- readRDS("dds_pd_control_BL_object.rds")
log2_expr <- log2(counts(dds, normalized=TRUE) + 1)            # full genes x 582 samples
rownames(log2_expr) <- sub("\\..*$", "", rownames(log2_expr))  # strip ENSG version

pd <- read.csv("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/data/PD_all_clustering_methods.csv")   # 393 PD SAMPLE_IDs
pd_ids <- intersect(pd$SAMPLE_ID, colnames(log2_expr))
mat <- log2_expr[, pd_ids]                                      # full transcriptome x 393 PD
cat(sprintf("Full PD matrix: %d genes x %d samples\n", nrow(mat), ncol(mat)))

ci <- read.csv("complex/CI_genes_converted.csv", stringsAsFactors=FALSE)
ci_set <- intersect(ci$original_id, rownames(mat))
cat(sprintf("CI gene set matched in transcriptome: %d / %d\n", length(ci_set), nrow(ci)))
gs <- list(CI = ci_set)

# ssGSEA via GSVA (new ssgseaParam API, fallback to legacy)
ss <- tryCatch({
  par <- GSVA::ssgseaParam(as.matrix(mat), gs, normalize=TRUE)
  GSVA::gsva(par)
}, error=function(e){
  cat("new API failed (", conditionMessage(e), "), trying legacy...\n")
  GSVA::gsva(as.matrix(mat), gs, method="ssgsea", ssgsea.norm=TRUE, verbose=FALSE)
})
ssgsea <- as.numeric(ss["CI", ]); names(ssgsea) <- colnames(ss)
cat("ssGSEA computed for", length(ssgsea), "samples\n")

# merge with existing scores
sc <- read.csv(file.path(FD,"three_scores.csv"))               # SAMPLE_ID, PC1, SumZ, RankComposite
sc$ssGSEA <- ssgsea[match(sc$SAMPLE_ID, names(ssgsea))]
# sign align (increase with CI expression)
if(cor(sc$ssGSEA, sc$SumZ, use="complete.obs") < 0) sc$ssGSEA <- -sc$ssGSEA

cat("\n==== 3-SCORE CONCORDANCE (now with REAL ssGSEA) ====\n")
cc <- cor(sc[,c("PC1","SumZ","ssGSEA")], use="complete.obs", method="pearson")
print(round(cc,3))
cat(sprintf("PC1 vs ssGSEA: Pearson r=%.3f, Spearman=%.3f\n",
            cor(sc$PC1, sc$ssGSEA, use="complete.obs"),
            cor(sc$PC1, sc$ssGSEA, use="complete.obs", method="spearman")))
sp <- function(x) x > median(x, na.rm=TRUE)
agr <- mean(sp(sc$PC1)==sp(sc$ssGSEA), na.rm=TRUE)
cat(sprintf("Median-split agreement PC1 vs ssGSEA: %.1f%%\n", 100*agr))

write.csv(sc, file.path(FD,"three_scores.csv"), row.names=FALSE)
cat("\nUpdated three_scores.csv with REAL ssGSEA column\n")
