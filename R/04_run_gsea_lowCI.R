# ============================================================
# GSEA for the low-CI transcriptomic direction
# Positive rank score / NES = enriched toward low CI expression
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
TAB_DIR <- file.path(ROOT, "04_tables")
FALLBACK_DATA_DIR <- normalizePath(file.path(ROOT, "..", "figures_pc1", "data"),
  winslash = "/", mustWork = FALSE)

read_first <- function(paths, ...) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop("None of these input files exists: ", paste(paths, collapse = " | "))
  }
  read.csv(hit, ...)
}

de <- read_first(c(file.path(FALLBACK_DATA_DIR, "de_full.csv")))
de <- de[!is.na(de$Symbol) & !is.na(de$stat), c("Symbol", "stat")]
de$rank_score <- -de$stat
rank_df <- de %>%
  group_by(Symbol) %>%
  slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list <- rank_df$rank_score
names(gene_list) <- rank_df$Symbol
gene_list <- sort(gene_list, decreasing = TRUE)

run_gsea <- function(collection, subcollection, out_name) {
  term2gene <- msigdbr(
    species = "Homo sapiens",
    collection = collection,
    subcollection = subcollection
  ) %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  set.seed(20260804)
  gsea <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    eps = 0,
    verbose = FALSE
  )

  out <- as.data.frame(gsea)
  out <- out[order(out$p.adjust, -abs(out$NES)), ]
  write.csv(out, file.path(TAB_DIR, out_name), row.names = FALSE)
  out
}

go_bp <- run_gsea("C5", "GO:BP", "GSEA_GO_BP_lowCI_ranked.csv")
kegg <- run_gsea("C2", "CP:KEGG_LEGACY", "GSEA_KEGG_lowCI_ranked.csv")

cat(sprintf(
  "Saved GSEA tables. GO BP FDR<0.05: %d low-CI, %d high-CI; KEGG FDR<0.05: %d low-CI, %d high-CI\n",
  sum(go_bp$p.adjust < 0.05 & go_bp$NES > 0, na.rm = TRUE),
  sum(go_bp$p.adjust < 0.05 & go_bp$NES < 0, na.rm = TRUE),
  sum(kegg$p.adjust < 0.05 & kegg$NES > 0, na.rm = TRUE),
  sum(kegg$p.adjust < 0.05 & kegg$NES < 0, na.rm = TRUE)
))
