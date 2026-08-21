# Install bioinformatics R packages for figure generation
cat("=== Checking & Installing Packages ===\n")

cran_pkgs <- c("ggpubr","factoextra","pROC","patchwork","cowplot",
               "emmeans","ggeffects","sjPlot","broom.mixed","ggplot2",
               "dplyr","tidyr","RColorBrewer","scales","ggrepel")
bioc_pkgs <- c("ComplexHeatmap","EnhancedVolcano")

installed <- rownames(installed.packages())

cat("\n--- Already installed ---\n")
for(p in c(cran_pkgs, bioc_pkgs)) {
  cat(sprintf("  %-18s %s\n", p, ifelse(p %in% installed, "YES", "-- missing")))
}

# Install missing CRAN
to_install <- setdiff(cran_pkgs, installed)
if(length(to_install) > 0) {
  cat(sprintf("\n--- Installing CRAN: %s ---\n", paste(to_install, collapse=", ")))
  install.packages(to_install, repos="https://cloud.r-project.org", quiet=TRUE)
}

# Install missing Bioc
bioc_missing <- setdiff(bioc_pkgs, installed)
if(length(bioc_missing) > 0) {
  cat(sprintf("\n--- Installing Bioc: %s ---\n", paste(bioc_missing, collapse=", ")))
  if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install(bioc_missing, update=FALSE, ask=FALSE)
}

# Final report
cat("\n=== Final status ===\n")
installed2 <- rownames(installed.packages())
for(p in c(cran_pkgs, bioc_pkgs)) {
  cat(sprintf("  %-18s %s\n", p, ifelse(p %in% installed2, "OK", "FAILED")))
}
cat("\nDone.\n")
