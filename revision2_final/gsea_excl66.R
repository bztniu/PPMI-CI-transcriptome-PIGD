# Reviewer-driven rerun: KEGG GSEA excluding the 66 CI scoring genes (removes circularity)
# Convention identical to manuscript SECTION 06: rank = -Wald stat, positive NES = low-CI direction
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr)})

TAB <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"
DATA <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
OUT  <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-31-第四版-评审修订"

adj <- read.csv(file.path(TAB, "DE_covariate_adjusted_full.csv"), check.names = FALSE)
names(adj)[1] <- "ENSG"
adj$ENSG <- sub("\\..*", "", adj$ENSG)
map <- read.csv(file.path(DATA, "de_full.csv"))[, c("ENSG", "Symbol")]
map$ENSG <- sub("\\..*", "", map$ENSG)
adj <- merge(adj, map, by = "ENSG", all.x = TRUE)

ld <- read.csv(file.path(OUT, "PC1_loadings_66genes.csv"))
adj2 <- adj[!(adj$Symbol %in% ld$Gene) & !is.na(adj$stat), ]
cat("genes after removing 66 CI genes:", nrow(adj2), "\n")

rnk <- setNames(-adj2$stat, adj2$Symbol)
rnk <- rnk[!duplicated(names(rnk)) & !is.na(rnk) & names(rnk) != "" & !is.na(names(rnk))]
rnk <- sort(rnk, decreasing = TRUE)

t2g <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol) %>% distinct()
set.seed(20260831)
g <- GSEA(geneList = rnk, TERM2GENE = t2g, minGSSize = 15, maxGSSize = 500,
          pvalueCutoff = 1, eps = 0, verbose = FALSE)
o <- as.data.frame(g)
o <- o[order(o$p.adjust, -abs(o$NES)), ]
write.csv(o, file.path(OUT, "GSEA_KEGG_excl66genes_ADJUSTED.csv"), row.names = FALSE)

cat("\n==== OXPHOS row (excl. 66 genes) ====\n")
print(o[grep("OXIDATIVE_PHOSPHORYLATION", o$ID), c("ID", "setSize", "NES", "p.adjust")])
cat("\n==== top low-CI-side pathways ====\n")
print(head(o[o$NES > 0, c("Description", "NES", "p.adjust")], 8))
cat("\n==== top high-CI-side pathways ====\n")
print(head(o[o$NES < 0, c("Description", "NES", "p.adjust")], 8))
cat("DONE\n")
