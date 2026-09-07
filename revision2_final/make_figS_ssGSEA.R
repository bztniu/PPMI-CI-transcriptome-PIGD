suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyr) })
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
FG <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

sc <- read.csv(file.path(FD,"three_scores.csv"))
sc <- sc[complete.cases(sc$PC1, sc$ssGSEA),]
sc$PC1g <- factor(ifelse(sc$PC1 < median(sc$PC1), "Low", "High"), levels=c("Low","High"))
sc$SSg  <- factor(ifelse(sc$ssGSEA < median(sc$ssGSEA), "Low", "High"), levels=c("Low","High"))
pct <- round(100*mean(sc$PC1g == sc$SSg))

LOW<-"#D55E00"; HIGH<-"#0072B2"
theme_pub <- theme_classic(base_size=8, base_family="sans") +
  theme(aspect.ratio=1, axis.line=element_line(linewidth=0.35, colour="grey30"),
        axis.ticks=element_line(linewidth=0.35, colour="grey30"),
        axis.title=element_text(size=8,face="bold"), axis.text=element_text(size=7,colour="grey20"),
        plot.title=element_text(size=10,face="bold"), legend.position="none", panel.grid=element_blank())

# A: ssGSEA distribution + median split
fA <- ggplot(sc, aes(ssGSEA, fill=SSg)) +
  geom_histogram(bins=40, colour="white", alpha=0.85, linewidth=0.15, position="identity") +
  geom_vline(xintercept=median(sc$ssGSEA), linetype="dashed", linewidth=0.4, colour="grey20") +
  annotate("text", x=median(sc$ssGSEA), y=Inf, label="median", vjust=2, hjust=-0.1, size=2.5, colour="grey20") +
  scale_fill_manual(values=c(Low=LOW, High=HIGH)) + scale_y_continuous(expand=expansion(mult=c(0,0.06))) +
  labs(title="ssGSEA score distribution", x="ssGSEA enrichment score", y="Patients") + theme_pub

# B: contingency / agreement matrix
tab <- table(PC1=sc$PC1g, ssGSEA=sc$SSg)
tdf <- as.data.frame(as.table(tab)); names(tdf) <- c("PC1","ssGSEA","n")
fB <- ggplot(tdf, aes(PC1, ssGSEA)) +
  geom_tile(aes(fill=n), colour="white", linewidth=0.5) +
  geom_text(aes(label=n), size=4, colour="grey20", fontface="bold") +
  scale_fill_gradient(low="#EEEDFE", high="#534AB7", guide="none") +
  labs(title=sprintf("Median-split agreement\n%.0f%% (%d/%d)", pct, sum(tab[row(tab)==col(tab)]), sum(tab)),
       x="PC1 assignment", y="ssGSEA assignment") +
  theme_pub + theme(aspect.ratio=NULL)

fig <- fA | fB
w<-123/25.4; h<-55/25.4
ggsave(file.path(FG,"figS_ssGSEA_concordance.pdf"), fig, width=w, height=h, device=cairo_pdf)
ggsave(file.path(FG,"figS_ssGSEA_concordance_preview.png"), fig, width=w, height=h, dpi=160)
cat("Saved figS_ssGSEA_concordance. Agreement:", pct, "%\n")
