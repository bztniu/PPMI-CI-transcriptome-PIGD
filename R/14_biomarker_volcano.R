# Biomarker sweep volcano: group effect across 117 biomarkers
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel)})
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
res <- read.csv(file.path(OUT,"biomarker_LMM_sweep.csv"))
res <- res[!is.na(res$Grp_p),]

# classify family
fam <- function(b){
  if(grepl("GL2|GlcCer|Cer$|Cer_|_Cer|_SM$|SM_|Sphing|Lacto|Hex", b, ignore.case=TRUE)) "Glycosphingolipid"
  else if(grepl("NFL|NfL|GFAP|Tau|Alpha_syn|synuclein|TREM2|Abeta|AB42|sTREM", b, ignore.case=TRUE)) "Neurodegen/inflamm"
  else if(grepl("GCase|GBA|Cathepsin|LAMP|Lyso", b, ignore.case=TRUE)) "Lysosomal"
  else "Other"
}
res$Family <- sapply(res$Biomarker, fam)
res$nlogp <- -log10(res$Grp_p)
res$lab <- ifelse(res$Grp_p<0.05 | grepl("NFL_CSF|GFAP_CSF|CSF_Alpha_synuclein|tTau_CSF|NfL_Serum|sTREM2", res$Biomarker),
                  gsub("_"," ",res$Biomarker), "")

p <- ggplot(res, aes(Grp_est, nlogp, color=Family)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey60") +
  geom_vline(xintercept=0, linetype="dotted", color="grey70") +
  geom_point(size=2, alpha=0.8) +
  geom_text_repel(aes(label=lab), size=2.6, max.overlaps=20, min.segment.length=0, show.legend=FALSE) +
  scale_color_manual(values=c("Glycosphingolipid"="#D55E00","Neurodegen/inflamm"="#0072B2",
                              "Lysosomal"="#009E73","Other"="grey60")) +
  annotate("text", x=min(res$Grp_est,na.rm=TRUE), y=-log10(0.05), label="p=0.05", vjust=-0.5, hjust=0, size=3, color="grey40") +
  labs(title="117 blood/CSF biomarkers x PC1 CI group (overall group effect)",
       subtitle="No biomarker survives FDR. Nominal hits = plasma glycosphingolipids + CSF sTREM2. NfL/GFAP/aSyn/tau all NS.",
       x="Group effect estimate (log scale; positive = higher in low CI)", y="-log10(p)",
       color="Biomarker family") +
  theme_classic(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(OUT,"biomarker_sweep_volcano.pdf"), p, width=10, height=7.5)
cat("Saved biomarker_sweep_volcano.pdf\n")
cat(sprintf("Tested %d biomarkers; nominal p<0.05: %d; survive FDR: %d\n",
            nrow(res), sum(res$Grp_p<0.05), sum(res$Grp_padj<0.05,na.rm=TRUE)))
