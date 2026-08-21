# Forest plot of LMM sweep: V12 interaction across all outcomes
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
res <- read.csv(file.path(OUT,"LMM_sweep_all_outcomes.csv"))

res$sig <- ifelse(res$V12_int_padj<0.05, "FDR<0.05",
           ifelse(res$V12_int_p<0.05, "nominal p<0.05", "ns"))
res$SE <- abs(res$V12_int_est) / qnorm(1 - res$V12_int_p/2)  # approx SE from est & p
res$SE[!is.finite(res$SE)] <- NA
res <- res[order(res$V12_int_p),]
res$Outcome <- factor(res$Outcome, levels=rev(res$Outcome))
lab <- c(PIGD="PIGD",RBDSQ="RBDSQ (REM sleep)",MoCA="MoCA (cognition)",
  SCOPA_AUT="SCOPA-AUT (autonomic)",UPDRS_III="MDS-UPDRS III",UPDRS_II="MDS-UPDRS II",
  UPDRS_I_pat="MDS-UPDRS I (patient)",UPDRS_I_clin="MDS-UPDRS I (clinician)",
  UPDRS_IV="MDS-UPDRS IV",Tremor="Total tremor",RestTremor="Rest tremor",NHY="Hoehn & Yahr")
res$lab <- factor(lab[as.character(res$Outcome)], levels=lab[as.character(res$Outcome)])

p <- ggplot(res, aes(V12_int_est, lab, color=sig)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_errorbarh(aes(xmin=V12_int_est-1.96*SE, xmax=V12_int_est+1.96*SE), height=0.25, linewidth=0.6, na.rm=TRUE) +
  geom_point(size=3) +
  geom_text(aes(label=sprintf("p=%.3f%s", V12_int_p, ifelse(V12_int_padj<0.05," *FDR",""))),
            hjust=-0.15, size=3, color="black") +
  scale_color_manual(values=c("FDR<0.05"="#CB181D","nominal p<0.05"="#FD8D3C","ns"="grey60")) +
  labs(title="Year-5 (V12) group x time interaction across 12 longitudinal scales",
       subtitle="PC1 Low vs High CI. Only PIGD survives FDR correction (n=12 outcomes tested).",
       x="V12 x group interaction estimate (Low vs High)", y=NULL, color=NULL) +
  coord_cartesian(xlim=c(-2.2, 2.6)) +
  theme_classic(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(OUT,"LMM_sweep_forest.pdf"), p, width=11, height=6)
cat("Saved LMM_sweep_forest.pdf\n")
