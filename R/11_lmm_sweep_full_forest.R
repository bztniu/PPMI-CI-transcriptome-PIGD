# Comprehensive forest plot: 43-indicator LMM sweep
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
res <- read.csv(file.path(OUT,"LMM_sweep_COMPREHENSIVE.csv"))

# category
cat_map <- function(o){
  if(o %in% c("UPDRS_III","PIGD","Tremor","RestTremor","NHY","UPDRS3_rigidity","UPDRS3_bradykinesia","UPDRS3_axial","UPDRS_II","UPDRS_IV")) "Motor"
  else if(grepl("MoCA",o)) "Cognition"
  else if(grepl("DAT",o)) "DAT imaging"
  else if(grepl("SCOPA",o)) "Autonomic"
  else "Non-motor (UPDRS I / sleep)"
}
res$Category <- sapply(res$Outcome, cat_map)
res$sig <- ifelse(res$LastV_int_padj<0.05,"FDR<0.05", ifelse(res$LastV_int_p<0.05,"nominal","ns"))
res$SE <- abs(res$LastV_int_est)/qnorm(1-res$LastV_int_p/2); res$SE[!is.finite(res$SE)] <- NA
res <- res[order(res$Category, res$LastV_int_p),]
res$Outcome <- factor(res$Outcome, levels=rev(res$Outcome))

p <- ggplot(res, aes(LastV_int_est, Outcome, color=sig)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=LastV_int_est-1.96*SE, xmax=LastV_int_est+1.96*SE), height=0.3, linewidth=0.5, na.rm=TRUE) +
  geom_point(size=2.3) +
  geom_text(aes(label=ifelse(LastV_int_p<0.05, sprintf("p=%.3f%s",LastV_int_p, ifelse(LastV_int_padj<0.05,"*","")), "")),
            hjust=-0.15, size=2.6, color="black", na.rm=TRUE) +
  scale_color_manual(values=c("FDR<0.05"="#CB181D","nominal"="#FD8D3C","ns"="grey65")) +
  facet_grid(Category~., scales="free_y", space="free_y") +
  labs(title="Comprehensive LMM sweep: 43 longitudinal indicators x PC1 CI group",
       subtitle="Last-visit group x time interaction (FDR across all 43). Only PIGD & axial-motor survive FDR (*).",
       x="Interaction estimate (Low vs High CI)", y=NULL, color=NULL) +
  coord_cartesian(xlim=c(-1.6,2.0)) +
  theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold",size=11), legend.position="top",
        strip.text.y=element_text(angle=0, face="bold", size=8))
ggsave(file.path(OUT,"LMM_sweep_COMPREHENSIVE_forest.pdf"), p, width=11, height=11)
cat("Saved. FDR-significant:", paste(res$Outcome[res$LastV_int_padj<0.05 & !is.na(res$LastV_int_padj)], collapse=", "), "\n")
