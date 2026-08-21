# ============================================================
# Figure 2: DE, KEGG-GSEA, GO enrichment, cross-cohort, and immune panels
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(ggtext)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
FIG_DIR <- file.path(ROOT, "02_figures")
TAB_DIR <- file.path(ROOT, "04_tables")
FALLBACK_DATA_DIR <- normalizePath(file.path(ROOT, "..", "figures_pc1", "data"),
  winslash = "/", mustWork = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

read_first <- function(paths, ...) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop("None of these input files exists: ", paste(paths, collapse = " | "))
  }
  read.csv(hit, ...)
}

LOW <- "#D55E00"
HIGH <- "#0072B2"
PPMI_COL <- "#0072B2"
GSE_COL <- "#E69F00"
SIG_COL <- "#D55E00"
NS_COL <- "grey62"

theme_pub <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.tag = element_text(face = "bold", size = 11),
    plot.margin = margin(4, 5, 4, 5)
  )

# ---- 2A: volcano ----
de <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "de_full.csv"),
  file.path(ROOT, "..", "figures_pc1", "data", "de_full.csv")
))
de <- de[!is.na(de$padj), ]
de$Symbol <- ifelse(is.na(de$Symbol), "", as.character(de$Symbol))
de$padj_plot <- pmax(de$padj, .Machine$double.xmin)

pd_genes <- c(
  "SNCA", "LRRK2", "PRKN", "PINK1", "PARK7", "ATP13A2", "FBXO7",
  "VPS35", "CHCHD2", "VPS13C", "PLA2G6", "DNAJC6", "SYNJ1", "EIF4G1",
  "GBA1", "TMEM175", "GCH1", "MAPT", "GPNMB", "MMP16", "SCARB2",
  "ITGA8", "RIT2", "BST1", "STK39", "MCCC1", "SYT11", "CTSB",
  "RAB29", "PM20D1", "GAK", "WNT3", "DDRGK1", "UBAP2L", "SIPA1L2",
  "INPP5F", "BAG3", "CNTNAP2", "ELOVL7", "SCN2A", "RICTOR", "MED13",
  "BCKDK", "CCDC62", "MEX3C", "ASXL3", "STX1B", "UBQLN1", "SATB1",
  "DMPK", "GLT8D1"
)
label_genes <- c("MMP16", "ITGA8", "MAPT", "GPNMB", "RIT2")

de$pass <- de$padj < 0.05 & abs(de$log2FoldChange) > 1
de$is_risk <- de$Symbol %in% pd_genes
de$direction <- ifelse(de$log2FoldChange < 0, "LowUp", "HighUp")
de$group <- ifelse(!de$pass, "nonsig",
  ifelse(de$is_risk,
    ifelse(de$direction == "LowUp", "risk_low", "risk_high"),
    ifelse(de$direction == "LowUp", "other_low", "other_high")
  )
)
de$group <- factor(de$group,
  levels = c("risk_low", "risk_high", "other_low", "other_high", "nonsig"))
de$label <- ifelse(de$Symbol %in% label_genes, de$Symbol, "")
de <- de[order(de$pass), ]

n_low_up <- sum(de$pass & de$log2FoldChange < 0, na.rm = TRUE)
n_high_up <- sum(de$pass & de$log2FoldChange > 0, na.rm = TRUE)

f2a <- ggplot(de, aes(log2FoldChange, -log10(padj_plot))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.25) +
  geom_point(aes(color = group, size = group), alpha = 0.68, stroke = 0) +
  geom_text_repel(
    data = subset(de, label != ""), aes(label = label),
    size = 2.4, fontface = "italic", color = "black",
    box.padding = 0.25, point.padding = 0.18,
    min.segment.length = 0, max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = c(
    risk_low = LOW, risk_high = HIGH,
    other_low = "#F5C9A5", other_high = "#A8CBE5", nonsig = "grey83"
  ), guide = "none") +
  scale_size_manual(values = c(
    risk_low = 1.8, risk_high = 1.8,
    other_low = 0.45, other_high = 0.45, nonsig = 0.25
  ), guide = "none") +
  coord_cartesian(xlim = c(-3.4, 2.4), ylim = c(0, 125), clip = "off") +
  annotate("text", x = -2.45, y = 78,
    label = sprintf("Low-CI expression\ngroup upregulated\n(%s genes)", format(n_low_up, big.mark = ",")),
    color = LOW, fontface = "bold", size = 2.3, hjust = 0.5) +
  annotate("text", x = 1.7, y = 74,
    label = sprintf("High-CI expression\ngroup upregulated\n(%s genes)", format(n_high_up, big.mark = ",")),
    color = HIGH, fontface = "bold", size = 2.3, hjust = 0.5) +
  labs(tag = "a", x = expression(log[2] ~ "fold change (High-CI expr. vs Low-CI expr.)"),
    y = expression(-log[10] ~ "FDR")) +
  theme_pub +
  theme(aspect.ratio = 1)

# ---- 2B: KEGG GSEA dot plot ----
kegg_gsea <- read_first(c(file.path(TAB_DIR, "GSEA_KEGG_lowCI_ranked.csv")))
kegg_keep <- c(
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_AXON_GUIDANCE",
  "KEGG_FOCAL_ADHESION",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_TIGHT_JUNCTION",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY"
)
kegg_panel <- kegg_gsea[kegg_gsea$ID %in% kegg_keep & kegg_gsea$NES > 0, ]
kegg_panel$Term <- gsub("^KEGG_", "", kegg_panel$ID)
kegg_panel$Term <- gsub("_", " ", kegg_panel$Term)
kegg_panel$Term <- tools::toTitleCase(tolower(kegg_panel$Term))
kegg_panel$Term <- gsub("Ecm", "ECM", kegg_panel$Term)
kegg_panel$Term <- gsub("Dna", "DNA", kegg_panel$Term)
kegg_panel$Term <- gsub("Parkinsons Disease", "Parkinson's Disease", kegg_panel$Term)
kegg_panel$Term <- factor(kegg_panel$Term, levels = kegg_panel$Term[order(kegg_panel$NES)])
kegg_panel$mlog10_fdr <- -log10(pmax(kegg_panel$p.adjust, .Machine$double.xmin))

f2b_gsea <- ggplot(kegg_panel, aes(NES, Term)) +
  geom_point(aes(size = mlog10_fdr), color = LOW, alpha = 0.92) +
  scale_size_continuous(range = c(1.7, 5.0), name = expression(-log[10] ~ "FDR")) +
  coord_cartesian(xlim = c(1.55, 3.05), clip = "off") +
  labs(tag = "b", x = "KEGG GSEA NES\n(low-CI expression direction)", y = NULL) +
  theme_pub +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_text(size = 5.7),
    axis.text.x = element_text(size = 6.2),
    axis.title.x = element_text(size = 6.7, margin = margin(t = 2)),
    legend.position = "right",
    legend.text = element_text(size = 5.8),
    legend.title = element_text(size = 6.0),
    legend.key.size = unit(0.16, "in")
  )

# ---- 2C: GO enrichment dot plot ----
go_lo <- read_first(c(file.path(TAB_DIR, "GO_lowCI_up.csv")))
go_lo <- head(go_lo[order(go_lo$p.adjust), ], 10)
go_lo$Description <- ifelse(nchar(go_lo$Description) > 44,
  paste0(substr(go_lo$Description, 1, 41), "..."), go_lo$Description)
go_lo$Description <- factor(go_lo$Description, levels = rev(go_lo$Description))
go_lo$GeneRatio_num <- vapply(strsplit(go_lo$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
}, numeric(1))

f2c_go <- ggplot(go_lo, aes(GeneRatio_num, Description)) +
  geom_point(aes(size = Count), color = "#0072B2", alpha = 0.92) +
  scale_size_continuous(range = c(1.6, 5.8), name = "Gene count") +
  labs(tag = "c", x = "Gene ratio", y = NULL) +
  theme_pub +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_text(size = 6.1),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6.5)
  )

# ---- 2D-E: PD-risk genes with concordant cross-cohort direction ----
risk <- read_first(c(file.path(TAB_DIR, "PD_risk_gene_screen.csv")))
risk$GSE_plot <- -risk$GSE_logFC
risk$concordant <- sign(risk$PPMI_logFC) == sign(risk$GSE_plot)
risk_conc <- risk[risk$concordant, ]
binom_p <- binom.test(nrow(risk_conc), nrow(risk), 0.5)$p.value

make_effect_panel <- function(df, tag, direction_label, xlim, breaks) {
  df <- df[order(abs(df$PPMI_logFC)), ]
  df$Gene <- factor(df$Symbol, levels = df$Symbol)
  df_long <- df %>%
    select(Gene, PPMI = PPMI_logFC, GENEPARK = GSE_plot) %>%
    pivot_longer(c(PPMI, GENEPARK), names_to = "Cohort", values_to = "logFC")
  df_long$Cohort <- factor(df_long$Cohort, levels = c("PPMI", "GENEPARK"))

  ggplot(df_long, aes(logFC, Gene, fill = Cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.22) +
    geom_col(position = position_dodge2(width = 0.62, preserve = "single"),
      width = 0.52, alpha = 0.95) +
    scale_fill_manual(values = c(PPMI = PPMI_COL, GENEPARK = GSE_COL)) +
    scale_x_continuous(breaks = breaks, expand = expansion(mult = c(0.01, 0.04))) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(tag = tag, x = direction_label, y = NULL) +
    theme_pub +
    theme(
      aspect.ratio = 1,
      axis.text.y = element_text(face = "italic", size = 5.4),
      axis.text.x = element_text(size = 6.2),
      axis.title.x = element_text(size = 6.6, margin = margin(t = 2)),
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(0.24, "in")
    )
}

low_panel <- risk_conc[risk_conc$PPMI_logFC < 0, ]
high_panel <- risk_conc[risk_conc$PPMI_logFC > 0, ]

f2d_risk_low <- make_effect_panel(
  low_panel, "d", "logFC (negative = Low-CI expr. up)",
  c(-3.0, 0), c(-3, -2, -1, 0)
)
f2e_risk_high <- make_effect_panel(
  high_panel, "e", "logFC (positive = High-CI expr. up)",
  c(0, 1.05), c(0, 0.5, 1.0)
) +
  theme(legend.position = "none")
f2e_standalone <- f2e_risk_high + theme(legend.position = "top")

cd_title <- sprintf("PD risk genes, concordant direction: %d/%d (%.0f%%), binomial p = %.2g",
  nrow(risk_conc), nrow(risk), 100 * nrow(risk_conc) / nrow(risk), binom_p)
fig2de <- (f2d_risk_low | f2e_risk_high) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = cd_title,
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  )

# ---- 2F-G: immune-cell shifts, biologically ordered across cohorts ----
ppmi_imm <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "immune_cohend.csv"),
  file.path(ROOT, "..", "figures_pc1", "data", "immune_cohend.csv")
))
ppmi_imm <- ppmi_imm[is.finite(ppmi_imm$d), c("Cell", "d", "P", "d_lo", "d_hi", "Padj")]
names(ppmi_imm)[names(ppmi_imm) == "d"] <- "D"

gse_imm <- read_first(c(
  "E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_immune_NNLS.csv"
))
gse_imm <- gse_imm[is.finite(gse_imm$D), c("Cell", "D", "P", "Padj")]

add_ci <- function(df) {
  if (all(c("d_lo", "d_hi") %in% names(df))) {
    return(df)
  }
  z <- abs(qnorm(df$P / 2))
  se <- abs(df$D) / z
  se[!is.finite(se)] <- NA_real_
  df$d_lo <- df$D - 1.96 * se
  df$d_hi <- df$D + 1.96 * se
  df
}
ppmi_imm <- add_ci(ppmi_imm)
gse_imm <- add_ci(gse_imm)

imm <- merge(ppmi_imm, gse_imm, by = "Cell", suffixes = c("_PPMI", "_GSE"))
imm <- imm[is.finite(imm$D_PPMI) & is.finite(imm$D_GSE), ]
imm$concordant <- sign(imm$D_PPMI) == sign(imm$D_GSE)
imm$both_sig <- imm$Padj_PPMI < 0.05 & imm$Padj_GSE < 0.05
imm$either_sig <- imm$Padj_PPMI < 0.05 | imm$Padj_GSE < 0.05
imm$ci_double_positive <- imm$d_lo_PPMI > 0 & imm$d_lo_GSE > 0
imm$ci_double_negative <- imm$d_hi_PPMI < 0 & imm$d_hi_GSE < 0
imm$ci_concordant <- imm$ci_double_positive | imm$ci_double_negative
imm$priority <- dplyr::case_when(
  imm$ci_concordant & imm$both_sig ~ 1,
  imm$ci_concordant ~ 2,
  imm$concordant & imm$either_sig ~ 3,
  imm$either_sig ~ 4,
  TRUE ~ 5
)

bio_order <- c(
  "Neutrophils",
  "Dendritic cells resting", "Dendritic cells activated",
  "Macrophages M0", "Macrophages M2", "Macrophages M1",
  "Monocytes", "Eosinophils", "NK cells resting", "Mast cells resting",
  "T cells CD8", "T cells CD4 memory resting", "T cells CD4 memory activated",
  "T cells regulatory (Tregs)", "T cells follicular helper",
  "B cells naive", "B cells memory", "Plasma cells"
)
imm$bio_rank <- match(imm$Cell, bio_order)
imm$bio_rank[is.na(imm$bio_rank)] <- length(bio_order) + seq_len(sum(is.na(match(imm$Cell, bio_order))))
imm <- imm[order(imm$priority, imm$bio_rank, -abs(imm$D_PPMI)), ]

display_cells <- imm$Cell
label_cells <- ifelse(imm$ci_concordant,
  paste0("<span style='color:", SIG_COL, "'>", display_cells, "</span>"),
  display_cells)
names(label_cells) <- display_cells

prep_immune_panel <- function(imm, cohort) {
  if (cohort == "PPMI") {
    df <- data.frame(
      Cell = imm$Cell,
      Cell_label = label_cells[imm$Cell],
      d = imm$D_PPMI,
      d_lo = imm$d_lo_PPMI,
      d_hi = imm$d_hi_PPMI,
      sig = ifelse(imm$Padj_PPMI < 0.05, "padj<0.05", "ns")
    )
  } else {
    df <- data.frame(
      Cell = imm$Cell,
      Cell_label = label_cells[imm$Cell],
      d = imm$D_GSE,
      d_lo = imm$d_lo_GSE,
      d_hi = imm$d_hi_GSE,
      sig = ifelse(imm$Padj_GSE < 0.05, "padj<0.05", "ns")
    )
  }
  df$Cell_label <- factor(df$Cell_label, levels = rev(unname(label_cells[display_cells])))
  df$sig <- factor(df$sig, levels = c("padj<0.05", "ns"))
  df
}

make_immune_forest <- function(df, tag, cohort, xlim) {
  ggplot(df, aes(d, Cell_label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.28) +
    geom_errorbar(aes(xmin = d_lo, xmax = d_hi), width = 0.3, linewidth = 0.55,
      orientation = "y", na.rm = TRUE) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("padj<0.05" = SIG_COL, "ns" = NS_COL), drop = FALSE) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(tag = tag, x = "Cohen's d", y = cohort, color = NULL) +
    theme_pub +
    theme(
      legend.position = "none",
      axis.text.y = ggtext::element_markdown(size = 5.8),
      axis.text.x = element_text(size = 6.2),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 8.5, angle = 90, margin = margin(r = 4)),
      aspect.ratio = 1
    )
}

imm_ppmi <- prep_immune_panel(imm, "PPMI")
imm_gse <- prep_immune_panel(imm, "GENEPARK")
f2f_imm_ppmi <- make_immune_forest(imm_ppmi, "f", "PPMI", c(-1.15, 1.05))
f2g_imm_gse <- make_immune_forest(imm_gse, "g", "GENEPARK", c(-1.15, 1.15))

fig2fg <- (f2f_imm_ppmi | f2g_imm_gse) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Immune-cell shifts ordered by cross-cohort concordance and cell lineage",
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  )

# ---- 2H: k-means immune robustness ----
immune_frac <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "immune_fractions.csv"),
  file.path(ROOT, "..", "figures_pc1", "data", "immune_fractions.csv")
), check.names = FALSE)
cluster_df <- read_first(c(
  file.path(ROOT, "..", "data", "PD_all_clustering_methods.csv")
))
cluster_df <- cluster_df[, c("SAMPLE_ID", "KM_group")]
km_dat <- merge(immune_frac, cluster_df, by = "SAMPLE_ID")
km_dat <- km_dat[km_dat$KM_group %in% c("Low", "High"), ]

cohen_ci <- function(lo, hi) {
  lo <- lo[is.finite(lo)]
  hi <- hi[is.finite(hi)]
  n_lo <- length(lo)
  n_hi <- length(hi)
  pooled_sd <- sd(c(lo, hi))
  d <- (mean(lo) - mean(hi)) / pooled_sd
  se <- sqrt((n_lo + n_hi) / (n_lo * n_hi) + d^2 / (2 * (n_lo + n_hi - 2)))
  c(d = d, d_lo = d - 1.96 * se, d_hi = d + 1.96 * se,
    p = wilcox.test(lo, hi)$p.value)
}

km_cells <- display_cells[display_cells %in% names(km_dat)]
km_imm <- do.call(rbind, lapply(km_cells, function(cell) {
  est <- cohen_ci(km_dat[km_dat$KM_group == "Low", cell],
    km_dat[km_dat$KM_group == "High", cell])
  data.frame(Cell = cell, t(est), check.names = FALSE)
}))
km_imm$padj <- p.adjust(km_imm$p, "BH")
km_imm$sig <- factor(ifelse(km_imm$padj < 0.05, "padj<0.05", "ns"),
  levels = c("padj<0.05", "ns"))
km_labels <- ifelse(km_imm$padj < 0.05,
  paste0("<span style='color:", SIG_COL, "'>", km_imm$Cell, "</span>"),
  km_imm$Cell)
km_imm$Cell_label <- factor(km_labels, levels = rev(km_labels))
write.csv(km_imm[, c("Cell", "d", "d_lo", "d_hi", "p", "padj")],
  file.path(TAB_DIR, "kmeans_immune_full.csv"), row.names = FALSE)

f2h_kmeans <- ggplot(km_imm, aes(d, Cell_label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.25) +
  geom_errorbar(aes(xmin = d_lo, xmax = d_hi), width = 0.3, linewidth = 0.55,
    orientation = "y", na.rm = TRUE) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("padj<0.05" = SIG_COL, "ns" = NS_COL), drop = FALSE) +
  coord_cartesian(xlim = c(-0.65, 0.65), clip = "off") +
  labs(tag = "h", x = "Cohen's d\n(k-means low-CI cluster)", y = "k-means", color = NULL) +
  theme_pub +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.text.y = ggtext::element_markdown(size = 5.8),
    axis.text.x = element_text(size = 6.2),
    axis.title.x = element_text(size = 6.7),
    axis.title.y = element_text(size = 8.5, angle = 90, margin = margin(r = 4))
  )

# ---- 2I: cell-composition adjusted key-gene associations ----
cell_adj <- read_first(c(
  file.path(TAB_DIR, "celladj_keygenes.csv"),
  file.path(TAB_DIR, "TableS5_cell_adjusted_keygenes.csv")
))
cell_adj <- cell_adj[order(cell_adj$beta_adj), ]
cell_adj$Gene <- factor(cell_adj$Gene, levels = cell_adj$Gene)
cell_adj_l <- cell_adj %>%
  select(Gene, Raw = beta_raw, `Cell adjusted` = beta_adj) %>%
  pivot_longer(c(Raw, `Cell adjusted`), names_to = "Model", values_to = "beta")
cell_adj_l$Model <- factor(cell_adj_l$Model, levels = c("Raw", "Cell adjusted"))

f2i_celladj <- ggplot(cell_adj_l, aes(beta, Gene, color = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.25) +
  geom_line(aes(group = Gene), color = "grey76", linewidth = 0.35) +
  geom_point(size = 1.9) +
  scale_color_manual(values = c(Raw = "grey50", `Cell adjusted` = SIG_COL)) +
  coord_cartesian(xlim = c(-1.0, 0.05), clip = "off") +
  labs(tag = "i", x = "Association beta\n(negative = higher in Low-CI expr.)", y = NULL, color = NULL) +
  theme_pub +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_text(face = "italic", size = 6.0),
    axis.text.x = element_text(size = 6.2),
    axis.title.x = element_text(size = 6.7),
    legend.position = "top",
    legend.text = element_text(size = 6.1),
    legend.key.size = unit(0.18, "in")
  )

fig2 <- patchwork::wrap_plots(
  f2a, f2b_gsea, f2c_go, f2d_risk_low,
  f2e_risk_high, f2f_imm_ppmi, f2g_imm_gse, f2h_kmeans,
  f2i_celladj,
  ncol = 3
)

ggsave(file.path(FIG_DIR, "Figure2B_KEGG_GSEA.pdf"), f2b_gsea, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2B_KEGG_GSEA.png"), f2b_gsea, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2E_highCI_up_square.pdf"), f2e_standalone, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2E_highCI_up_square.png"), f2e_standalone, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2D_highCI_up_square.pdf"), f2e_standalone, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2D_highCI_up_square.png"), f2e_standalone, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2DE_cross_cohort.pdf"), fig2de, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2DE_cross_cohort.png"), fig2de, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2CD_cross_cohort.pdf"), fig2de, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2CD_cross_cohort.png"), fig2de, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2FG_immune_cross_cohort.pdf"), fig2fg, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2FG_immune_cross_cohort.png"), fig2fg, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2EF_immune_cross_cohort.pdf"), fig2fg, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2EF_immune_cross_cohort.png"), fig2fg, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2H_kmeans_immune.pdf"), f2h_kmeans, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2H_kmeans_immune.png"), f2h_kmeans, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2I_cell_adjusted_keygenes.pdf"), f2i_celladj, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2I_cell_adjusted_keygenes.png"), f2i_celladj, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2_combined.pdf"), fig2, width = 12, height = 11)
ggsave(file.path(FIG_DIR, "Figure2_combined.png"), fig2, width = 12, height = 11, dpi = 300)

message(sprintf(
  "Saved Figure 2 outputs to %s; PD-risk concordance: %d/%d; immune top cells: %s; cell-adjusted genes: %d/%d survive",
  FIG_DIR, nrow(risk_conc), nrow(risk), paste(head(display_cells, 6), collapse = ", "),
  sum(cell_adj$survives_adj), nrow(cell_adj)
))
