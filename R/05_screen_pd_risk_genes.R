# ============================================================
# Full PD risk gene screen in PPMI DESeq2 + GENEPARK limma
# ============================================================

# PPMI DESeq2 (High vs Low CI: negative = higher in Low CI)
ppmi <- read.csv("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/de_full.csv")
ppmi <- ppmi[!is.na(ppmi$padj) & !is.na(ppmi$Symbol), c("Symbol","log2FoldChange","pvalue","padj")]
names(ppmi) <- c("Symbol","PPMI_logFC","PPMI_pval","PPMI_padj")

# GENEPARK limma (Low vs High CI: positive = higher in Low CI)
gse <- read.csv("E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_limma_DE_PC1.csv")
gse <- gse[, c("Gene","logFC","adj.P.Val")]
names(gse) <- c("Symbol","GSE_logFC","GSE_padj")

# PD risk genes from Nalls 2019 GWAS + Mendelian genes
pd_genes <- c(
  # Mendelian
  "SNCA","LRRK2","PRKN","PINK1","PARK7","DJ1","ATP13A2","FBXO7",
  "VPS35","CHCHD2","VPS13C","PLA2G6","DNAJC6","SYNJ1","EIF4G1",
  # Strong risk
  "GBA1","GBA","TMEM175","GCH1",
  # GWAS top hits
  "MAPT","GPNMB","MMP16","SCARB2","BST1","STK39","MCCC1",
  "SYT11","ACMSD","FGF20","ITGA8","CTSB","HLA-DRB5","RAB29",
  "PM20D1","GAK","DDRGK1","UBAP2L","SIPA1L2","INPP5F",
  "BAG3","CNTNAP2","ELOVL7","SCN2A","RICTOR","MED13","BCKDK",
  "CCDC62","MEX3C","ASXL3","WNT3","STX1B","UBQLN1",
  "RIT2","SH3GL2","GLT8D1","DMPK","SATB1"
)

# Filter to genes actually in PPMI
df <- merge(ppmi, gse, by="Symbol", all=FALSE)
df <- df[df$Symbol %in% pd_genes, ]

# Label whether meets threshold
df$PPMI_sig <- (df$PPMI_padj < 0.05 & abs(df$PPMI_logFC) > 1)
df$PPMI_nominal <- (df$PPMI_padj < 0.05 & abs(df$PPMI_logFC) <= 1)
df$GSE_sig <- df$GSE_padj < 0.05

# Direction: PPMI negative=higher in Low CI, GSE positive=higher in Low CI
# Concordant if both point same way
df$concordant <- sign(df$PPMI_logFC) * sign(df$GSE_logFC) < 0

# Sort by PPMI abs logFC
df <- df[order(-abs(df$PPMI_logFC)), ]

cat("=== Full PD Risk Gene Screen (PPMI + GENEPARK) ===\n\n")
cat(sprintf("Total PD genes queried: %d\n", length(pd_genes)))
cat(sprintf("Found in both PPMI & GSE: %d\n", nrow(df)))
cat(sprintf("PPMI |logFC|>1 & padj<0.05: %d\n", sum(df$PPMI_sig)))
cat(sprintf("PPMI nominal only: %d\n", sum(df$PPMI_nominal)))
cat(sprintf("GSE padj<0.05: %d\n", sum(df$GSE_sig)))
cat(sprintf("Direction concordant: %d/%d\n\n", sum(df$concordant), nrow(df)))

cat(sprintf("%-14s %8s %10s %6s %8s %9s %6s\n",
            "Gene","PPMI_FC","PPMI_padj","Sig?","GSE_FC","GSE_padj","Dir?"))
cat(strrep("-",80),"\n")
for(i in 1:nrow(df)){
  r <- df[i,]
  sig_lab <- ifelse(r$PPMI_sig, "****",
             ifelse(r$PPMI_nominal, "nom*", "ns"))
  dir_lab <- ifelse(r$concordant, "SAME", "OPP")
  cat(sprintf("%-14s %8.2f %10.2e %6s %8.2f %9.2e %6s\n",
      r$Symbol, r$PPMI_logFC, r$PPMI_padj, sig_lab, r$GSE_logFC, r$GSE_padj, dir_lab))
}

# Also list genes NOT found
miss <- setdiff(pd_genes, df$Symbol)
cat("\nGenes NOT in expression matrix:\n")
cat(paste(miss, collapse=", "),"\n")

write.csv(df, "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/PD_risk_gene_screen.csv", row.names=FALSE)
cat("\nSaved: PD_risk_gene_screen.csv\n")
