# Supplementary Figure S1 (redesigned): CI axis is a continuous unimodal gradient
# A=PC1 density+dip | B=GMM BIC | C=consensus PAC | D=3-score concordance
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
FG <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
S  <- readRDS(file.path(FD,"structure_assessment.rds"))
pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))
sc <- S$scores; sc$group <- factor(pc1s$group[match(sc$SAMPLE_ID,pc1s$SAMPLE_ID)], levels=c("Low","High"))
LOW<-"#D55E00"; HIGH<-"#0072B2"; pal<-c(Low=LOW,High=HIGH)

theme_pub <- theme_classic(base_size=9, base_family="sans") +
  theme(aspect.ratio=1, axis.line=element_line(linewidth=0.4,colour="grey20"),
        axis.ticks=element_line(linewidth=0.4,colour="grey20"),
        axis.title=element_text(size=9,face="bold"), axis.text=element_text(size=8,colour="grey20"),
        plot.title=element_text(size=10,face="bold"), legend.title=element_blank(),
        legend.text=element_text(size=7.5), legend.key.size=unit(3.2,"mm"),
        legend.background=element_blank(), panel.grid=element_blank())

# A: PC1 density + dip test (unimodal)
med <- median(sc$PC1)
fA <- ggplot(sc, aes(PC1)) +
  geom_histogram(aes(y=after_stat(density)), bins=40, fill="grey80", colour="white", linewidth=0.2) +
  geom_density(linewidth=0.7, colour="grey20") +
  geom_vline(xintercept=med, linetype="dashed", linewidth=0.5, colour=LOW) +
  annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.4, size=2.7, colour="grey15",
           label=sprintf("Dip test: D=%.3f\np=%.2f (unimodal)", S$dip$statistic, S$dip$p.value)) +
  annotate("text", x=med, y=0, label="median split", angle=90, vjust=-0.4, hjust=-0.1, size=2.3, colour=LOW) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title="A  PC1 is a continuous unimodal axis", x="PC1 score (CI transcriptomic axis)", y="Density") + theme_pub

# B: GMM BIC vs #components
bdf <- data.frame(G=as.integer(names(S$bic_by_G)), BIC=as.numeric(S$bic_by_G))
fB <- ggplot(bdf, aes(G, BIC)) +
  geom_line(linewidth=0.6, colour="grey40") + geom_point(size=1.8, colour="grey30") +
  geom_point(data=bdf[which.max(bdf$BIC),], size=3, colour=LOW) +
  annotate("text", x=bdf$G[which.max(bdf$BIC)], y=max(bdf$BIC), label="best: 1 component",
           hjust=-0.12, vjust=0.2, size=2.6, colour=LOW, fontface="bold") +
  scale_x_continuous(breaks=1:5) +
  labs(title="B  Gaussian mixture model (BIC)", x="Number of components", y="BIC (higher = better)") + theme_pub

# C: consensus PAC vs k
pdf_ <- data.frame(k=2:6, PAC=as.numeric(S$pac))
fC <- ggplot(pdf_, aes(k, PAC)) +
  geom_line(linewidth=0.6, colour="grey40") + geom_point(size=1.8, colour="grey30") +
  geom_point(data=pdf_[1,], size=3, colour=LOW) +
  annotate("text", x=2, y=pdf_$PAC[1], label=sprintf("k=2\nPAC=%.2f",pdf_$PAC[1]), hjust=-0.2, vjust=0.4, size=2.6, colour=LOW) +
  scale_x_continuous(breaks=2:6) + ylim(0, max(0.4,max(pdf_$PAC)*1.1)) +
  labs(title="C  Consensus clustering (no stable partition)", x="Number of clusters (k)", y="PAC (ambiguity; lower=cleaner)") + theme_pub

# D: 3-score concordance (PC1 vs mean z-score)
r1 <- cor(sc$PC1, sc$SumZ); r2 <- cor(sc$PC1, sc$RankComposite)
fD <- ggplot(sc, aes(PC1, SumZ, colour=group)) +
  geom_point(size=1.1, alpha=0.8, shape=16) +
  scale_colour_manual(values=pal, labels=c("Low CI","High CI")) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.5, size=2.7, colour="grey15",
           label=sprintf("PC1 vs mean z: r=%.3f\nPC1 vs rank: r=%.3f", r1, r2)) +
  labs(title="D  CI score robust to scoring method", x="PC1 score", y="Mean z-score (equal-weight)") +
  theme_pub + theme(legend.position=c(0.84,0.16))

fig <- (fA | fB) / (fC | fD)
w <- 175/25.4; h <- 165/25.4
ragg::agg_tiff(file.path(FG,"FigS1_structure.tiff"), width=w, height=h, units="in", res=600); print(fig); dev.off()
grDevices::cairo_pdf(file.path(FG,"FigS1_structure.pdf"), width=w, height=h); print(fig); dev.off()
ragg::agg_png(file.path(FG,"FigS1_structure_preview.png"), width=w, height=h, units="in", res=160); print(fig); dev.off()
cat("Saved FigS1_structure (.pdf/.tiff/.png)\n")
