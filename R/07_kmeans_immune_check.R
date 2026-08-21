# k-means immune concordance check vs PC1 immune
suppressPackageStartupMessages({library(dplyr)})
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
imf <- read.csv(file.path(FD,"immune_fractions.csv"), check.names=FALSE)
cl  <- read.csv("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/data/PD_all_clustering_methods.csv")
cl <- cl[,c("SAMPLE_ID","All_CI_score","KM_cluster")]
d <- imf %>% inner_join(cl, by="SAMPLE_ID")

# map KM cluster -> Low/High by mean CI score (lower CI = "Low")
km_means <- tapply(d$All_CI_score, d$KM_cluster, mean)
low_cl <- names(which.min(km_means))
d$KMgroup <- ifelse(d$KM_cluster==as.integer(low_cl) | d$KM_cluster==low_cl, "Low","High")
cat(sprintf("KM clusters: %s ; CI means: %s ; Low=cluster %s\n",
            paste(names(km_means),collapse="/"), paste(round(km_means,2),collapse="/"), low_cl))
cat(sprintf("KM groups: Low=%d High=%d\n", sum(d$KMgroup=="Low"), sum(d$KMgroup=="High")))

cells <- c("Neutrophils","T cells CD8","T cells CD4 memory resting","Dendritic cells resting",
           "Plasma cells","Macrophages M0","Dendritic cells activated",
           "Eosinophils","Macrophages M2","T cells follicular helper")
cat(sprintf("\n%-30s %8s %8s %8s | %s\n","Cell","KM_d","KM_p","KM_padj","dir"))
res <- data.frame()
for(ct in cells){
  if(!ct %in% colnames(d)) next
  lo<-d[d$KMgroup=="Low",ct]; hi<-d[d$KMgroup=="High",ct]
  dd<-(mean(lo)-mean(hi))/sd(c(lo,hi)); p<-wilcox.test(lo,hi)$p.value
  res<-rbind(res,data.frame(Cell=ct,KM_d=dd,KM_p=p))
}
res$KM_padj<-p.adjust(res$KM_p,"BH")
for(i in 1:nrow(res)){r<-res[i,]
  cat(sprintf("%-30s %+8.3f %8.4f %8.4f | %s%s\n",r$Cell,r$KM_d,r$KM_p,r$KM_padj,
              ifelse(r$KM_d>0,"Low>High","High>Low"),ifelse(r$KM_padj<0.05," *","")))}
cat("\n[PC1 reference: Neutrophils d=+0.39*, CD8 d=-0.34*, CD4mem-resting d=-0.29*, DC-resting +0.39*, Plasma +0.32*, M0 +0.24*, DC-activated -0.42*; Eosino/M2/Tfh NS]\n")
write.csv(res, file.path(FD,"kmeans_immune_check.csv"), row.names=FALSE)
