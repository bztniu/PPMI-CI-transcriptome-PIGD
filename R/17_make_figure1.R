# ============================================================
# Figure 1 (3 square panels): A=PCA biplot(ellipse) | B=scree | C=PC1 distribution
# ============================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
FG <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

pca   <- readRDS(file.path(FD,"pca.rds"))
pc1s  <- read.csv(file.path(FD,"pc1_scores.csv"))
varexp<- read.csv(file.path(FD,"pca_variance.csv"))
pc1s$group <- factor(pc1s$group, levels=c("Low","High"))

sc <- data.frame(SAMPLE_ID=rownames(pca$x), PC2=pca$x[,2])
# use canonical PC1 from pc1_scores (sign-correct: Low CI = lower PC1, as used throughout paper)
sc$PC1   <- pc1s$PC1[match(sc$SAMPLE_ID, pc1s$SAMPLE_ID)]
sc$group <- factor(pc1s$group[match(sc$SAMPLE_ID, pc1s$SAMPLE_ID)], levels=c("Low","High"))
sc <- sc[!is.na(sc$PC1), ]
v1 <- round(summary(pca)$importance[2,1]*100,1)
v2 <- round(summary(pca)$importance[2,2]*100,1)
LOW <- "#D55E00"; HIGH <- "#0072B2"; pal <- c(Low=LOW, High=HIGH)

theme_pub <- theme_classic(base_size=9, base_family="sans") +
  theme(aspect.ratio=1,                                   # <- square panels
        axis.line=element_line(linewidth=0.4,colour="grey20"),
        axis.ticks=element_line(linewidth=0.4,colour="grey20"),
        axis.title=element_text(size=9,face="bold"), axis.text=element_text(size=8,colour="grey20"),
        plot.title=element_text(size=10,face="bold"), legend.title=element_blank(),
        legend.text=element_text(size=8), legend.key.size=unit(3.5,"mm"),
        legend.background=element_blank(), panel.grid=element_blank())

# ---- A: PCA biplot with 95% ellipses (square) ----
fA <- ggplot(sc, aes(PC1, PC2, colour=group, fill=group)) +
  geom_hline(yintercept=0, linewidth=0.25, colour="grey88") +
  geom_vline(xintercept=0, linewidth=0.25, colour="grey88") +
  stat_ellipse(aes(group=group), type="norm", level=0.95, geom="polygon", alpha=0.10, colour=NA) +
  stat_ellipse(aes(group=group), type="norm", level=0.95, linewidth=0.55) +
  geom_point(size=1.2, alpha=0.8, shape=16) +
  scale_colour_manual(values=pal, labels=c("Low CI","High CI")) +
  scale_fill_manual(values=pal, guide="none") +
  labs(title="A  CI transcriptomic axis", x=sprintf("PC1 (%.1f%% variance)",v1),
       y=sprintf("PC2 (%.1f%% variance)",v2)) +
  theme_pub + theme(legend.position=c(0.84,0.93))

# ---- B: scree (square) ----
varexp$PCn <- factor(varexp$PC, levels=varexp$PC); vb <- varexp[1:8,]
fB <- ggplot(vb, aes(PCn, VarExplained*100)) +
  geom_col(fill="grey75", width=0.74) +
  geom_col(data=vb[1,], aes(PCn, VarExplained*100), fill="grey30", width=0.74) +
  geom_text(aes(label=sprintf("%.1f",VarExplained*100)), vjust=-0.35, size=2.4, colour="grey25") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(title="B  Variance explained", x="Principal component", y="% variance") + theme_pub

# ---- C: PC1 score distribution (square) ----
med <- median(pc1s$PC1)
fC <- ggplot(pc1s, aes(PC1, fill=group)) +
  geom_histogram(bins=40, colour="white", alpha=0.9, linewidth=0.2) +
  geom_vline(xintercept=med, linetype="dashed", linewidth=0.6, colour="grey25") +
  annotate("text", x=med, y=Inf, label="median split", vjust=1.8, hjust=-0.06, size=2.5, colour="grey30") +
  scale_fill_manual(values=pal, labels=c("Low CI","High CI")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.06))) +
  labs(title="C  PC1 score distribution", x="PC1 score (CI transcriptomic axis)", y="Patients") +
  theme_pub + theme(legend.position=c(0.83,0.88))

fig <- fA | fB | fC
w <- 190/25.4; h <- 72/25.4
ragg::agg_tiff(file.path(FG,"Figure1_PCA.tiff"), width=w, height=h, units="in", res=600); print(fig); dev.off()
grDevices::cairo_pdf(file.path(FG,"Figure1_PCA.pdf"), width=w, height=h); print(fig); dev.off()
ragg::agg_png(file.path(FG,"Figure1_PCA_preview.png"), width=w, height=h, units="in", res=160); print(fig); dev.off()
cat("Saved Figure1_PCA (.pdf/.tiff/.png) — 3 square panels. PC1 =", v1, "%\n")
