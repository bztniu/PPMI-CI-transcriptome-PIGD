# ============================================================
# Figure 1 via PCAtools: screeplot + biplot(ellipse) + loadings + eigencorplot
# ============================================================
suppressPackageStartupMessages({
  library(PCAtools); library(ggplot2); library(patchwork); library(dplyr); library(readr); library(ggplotify)
})
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
FG <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

# ---- expression matrix: genes x samples ----
expr <- read.csv(file.path(FD,"ci_expr.csv"), row.names=1, check.names=FALSE)  # samples x genes
mat <- t(as.matrix(expr))                                                       # genes x samples
pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))                                # SAMPLE_ID,PATNO,PC1,group

# ---- assemble clinical metadata (samples x variables) ----
demo <- read_csv("Demographics_03Feb2026.csv", show_col_types=FALSE) %>% dplyr::select(PATNO,SEX) %>% distinct(PATNO,.keep_all=TRUE)
age  <- read_csv("Age_at_visit_17Mar2025.csv", show_col_types=FALSE) %>% filter(EVENT_ID=="BL") %>% dplyr::select(PATNO,Age=AGE_AT_VISIT) %>% distinct(PATNO,.keep_all=TRUE)
moca <- read_csv("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("SC","V12")) %>% tidyr::pivot_wider(id_cols=PATNO, names_from=EVENT_ID, values_from=MCATOT) %>%
  dplyr::rename(MoCA_BL=SC, MoCA_Y5=V12)
m3 <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("BL","SC","V12")) %>% arrange(PATNO,EVENT_ID) %>%
  mutate(PIGD=NP3GAIT+NP3PSTBL+NP3FRZGT,
         RestTremor=NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL)
m3_bl <- m3 %>% filter(EVENT_ID %in% c("BL","SC")) %>% distinct(PATNO,.keep_all=TRUE) %>%
  dplyr::select(PATNO, UPDRSIII=NP3TOT, PIGD_BL=PIGD, RestTremor_BL=RestTremor)
m3_v12 <- m3 %>% filter(EVENT_ID=="V12") %>% distinct(PATNO,.keep_all=TRUE) %>% dplyr::select(PATNO, PIGD_Y5=PIGD)

meta <- pc1s %>% left_join(age,by="PATNO") %>% left_join(moca,by="PATNO") %>%
  left_join(m3_bl,by="PATNO") %>% left_join(m3_v12,by="PATNO")
meta$CI_group <- factor(meta$group, levels=c("Low","High"))
rownames(meta) <- meta$SAMPLE_ID
meta <- meta[colnames(mat), ]   # align to expression columns

# ---- PCAtools PCA (scale+center, same as primary analysis) ----
p <- pca(mat, metadata=meta, scale=TRUE, center=TRUE)
# enforce sign: low CI expression -> negative PC1
if(cor(p$rotated$PC1, colMeans(mat)) < 0){ p$rotated$PC1 <- -p$rotated$PC1; p$loadings$PC1 <- -p$loadings$PC1 }
v1 <- round(p$variance["PC1"],1); v2 <- round(p$variance["PC2"],1)
LOW <- "#D55E00"; HIGH <- "#0072B2"

base_thm <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.line=element_line(linewidth=0.35), axis.ticks=element_line(linewidth=0.35),
        axis.title=element_text(size=8,face="bold"), axis.text=element_text(size=7),
        plot.title=element_text(size=9,face="bold"), legend.text=element_text(size=7.5),
        legend.title=element_blank(), panel.grid=element_blank(), plot.tag=element_text(size=11,face="bold"))

# ---- A: biplot (samples, color by CI group, 95% ellipse) ----
pA <- biplot(p, colby="CI_group", colkey=c(Low=LOW, High=HIGH),
             ellipse=TRUE, ellipseLevel=0.95, ellipseFill=TRUE, ellipseAlpha=0.10,
             lab=NULL, pointSize=1.1, legendPosition="top",
             xlab=sprintf("PC1 (%.1f%%)",v1), ylab=sprintf("PC2 (%.1f%%)",v2),
             title="CI transcriptomic axis (66 complex I genes)") +
  base_thm + theme(legend.position=c(0.86,0.92)) + labs(tag="A")

# ---- B: screeplot ----
pB <- screeplot(p, components=getComponents(p,1:8), axisLabSize=7, titleLabSize=9,
                colBar="grey70", drawCumulativeSumLine=FALSE) +
  base_thm + labs(title="Variance explained", tag="B") +
  theme(axis.text.x=element_text(angle=0, size=6))

# ---- C: plotloadings (top genes on PC1) ----
pC <- plotloadings(p, components="PC1", rangeRetain=0.18, absolute=FALSE,
                   labSize=2.2, shape=21, col=c("#D95F0E","grey90","#2C7FB8"),
                   drawConnectors=TRUE, title="PC1 gene loadings", subtitle=NULL, caption=NULL,
                   legendPosition="none") +
  base_thm + labs(tag="C")

# ---- D: eigencorplot (PC1-5 vs clinical phenotypes) ----
metavars <- c("Age","UPDRSIII","PIGD_BL","RestTremor_BL","MoCA_BL","PIGD_Y5","MoCA_Y5")
pD <- eigencorplot(p, components=getComponents(p,1:5), metavars=metavars,
                   col=c("#2166AC","#92C5DE","white","#F4A582","#B2182B"),
                   cexCorval=0.7, fontCorval=1, rotLabX=45, scale=TRUE,
                   main="PC vs clinical phenotype", cexMain=1,
                   corFUN="pearson", corUSE="pairwise.complete.obs",
                   signifSymbols=c("***","**","*",""), signifCutpoints=c(0,0.001,0.01,0.05,1),
                   returnPlot=TRUE) 

# eigencorplot returns a lattice/trellis object -> convert to ggplot for patchwork
pD_gg <- ggplotify::as.ggplot(pD) + labs(tag="D")
fig_top <- pA + (pB / pC) + plot_layout(widths=c(1.5,1))
fig <- fig_top / pD_gg + plot_layout(heights=c(2,1))
w <- 183/25.4; h <- 230/25.4

ragg::agg_tiff(file.path(FG,"Figure1_PCA.tiff"), width=w, height=h, units="in", res=600)
print(fig); dev.off()
grDevices::cairo_pdf(file.path(FG,"Figure1_PCA.pdf"), width=w, height=h)
print(fig); dev.off()
ragg::agg_png(file.path(FG,"Figure1_PCA_preview.png"), width=w, height=h, units="in", res=150)
print(fig); dev.off()
cat("Saved Figure1_PCA (.pdf/.tiff/.png). PC1 var =", v1, "%\n")
cat("eigencorplot vars:", paste(metavars,collapse=", "), "\n")
