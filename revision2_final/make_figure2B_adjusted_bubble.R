# Replacement panel 2B: top-8 low-CI-direction KEGG pathways by adjusted P
# (covariate-adjusted GSEA, same data as Supplementary Table S8), bubble style
suppressPackageStartupMessages({ library(ggplot2) })
S8 <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-27-第二版/Supplementary_Tables_updated.xlsx"
stopifnot(requireNamespace("openxlsx", quietly = TRUE))
d <- openxlsx::read.xlsx(S8, sheet = "TableS8", startRow = 3)
d$adj..P <- as.numeric(d$adj..P); d$NES <- as.numeric(d$NES); for (cc in c("NES","adj..P","Set.size")) d[[cc]] <- as.numeric(d[[cc]])
d <- d[!is.na(d$NES) & d$NES > 0 & !is.na(d$Set.size), ]
d <- d[order(d$adj..P), ][1:8, ]
d$log10padj <- -log10(d$adj..P)
d$label <- gsub("KEGG_", "", d$ID)
d$label <- gsub("_", " ", d$label)
d$label <- tools::toTitleCase(tolower(d$label))
d$label <- gsub('Ecm', 'ECM', d$label)
d$label <- gsub('Abc', 'ABC', d$label)
d$label <- factor(d$label, levels = d$label[order(d$log10padj)])
named <- c("Neuroactive Ligand Receptor Interaction", "Olfactory Transduction",
           "Ecm Receptor Interaction", "Axon Guidance",
           "Complement And Coagulation Cascades")
d$col <- ifelse(d$label %in% named, "#C0392B", "black")

p <- ggplot(d, aes(x = log10padj, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(size = Set.size), shape = 21, colour = "black", fill = "white", stroke = 0.9) +
  geom_text(aes(label = label, colour = col), hjust = -0.04, size = 3.4,
            fontface = "bold", nudge_y = 0.12) +
  scale_colour_identity() +
  scale_size_continuous(range = c(3, 9), breaks = c(25, 50, 75, 100)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.80))) +
  labs(x = expression(-log[10](p.adj)), y = NULL, size = "Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        legend.text = element_text(size = 10), legend.title = element_text(size = 11))
ggsave("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/2026-08-27-第二版/Figure2B_adjusted_replacement.png",
       p, width = 7, height = 4.2, dpi = 600)
cat("saved\n")
print(d[, c("label", "NES", "adj..P", "Count")])
