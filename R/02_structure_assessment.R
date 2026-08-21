# ============================================================
# CI axis structure assessment: is it continuous or two clusters?
# Dip test + GMM/BIC (unimodality) + ConsensusClusterPlus (PAC) + 3-score concordance
# ============================================================
suppressPackageStartupMessages({
  library(diptest); library(mclust); library(ConsensusClusterPlus); library(cluster)
})
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
OUT<- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

expr <- read.csv(file.path(FD,"ci_expr.csv"), row.names=1, check.names=FALSE)  # 393 samples x 66 genes (log2)
pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))
Z <- scale(as.matrix(expr))                                                    # z per gene

# ---- 3 scoring methods (all 393 samples) ----
PC1  <- pc1s$PC1[match(rownames(Z), pc1s$SAMPLE_ID)]
SumZ <- rowMeans(Z)                                                            # equal-weight mean z
RankComposite <- rowMeans(apply(Z, 2, rank)) / nrow(Z)                         # rank-based (ssGSEA-style, non-parametric)
# align sign so all increase with CI expression
if(cor(PC1, SumZ) < 0) PC1 <- -PC1
scores <- data.frame(SAMPLE_ID=rownames(Z), PC1=PC1, SumZ=SumZ, RankComposite=RankComposite)

cat("================ 1. UNIMODALITY: is PC1 one continuous axis or two clusters? ================\n")
dt <- dip.test(PC1)
cat(sprintf("Hartigan's dip test on PC1: D=%.4f, p=%.4f  -> %s\n",
            dt$statistic, dt$p.value, ifelse(dt$p.value>0.05,"UNIMODAL (no evidence of 2 modes)","multimodal")))
for(nm in c("SumZ","RankComposite")){ d<-dip.test(scores[[nm]]); cat(sprintf("  dip %s: p=%.3f\n",nm,d$p.value)) }

cat("\n================ 2. GMM / BIC: how many Gaussian components best fit PC1? ================\n")
set.seed(123); mc <- Mclust(PC1, G=1:5, verbose=FALSE)
bic <- mclustBIC(PC1, G=1:5, verbose=FALSE)
bic_by_G <- apply(bic, 1, max, na.rm=TRUE)
cat(sprintf("Best model: G=%d (%s), BIC=%.1f\n", mc$G, mc$modelName, max(mc$BIC,na.rm=TRUE)))
cat("Max BIC by #components (higher=better):\n"); print(round(bic_by_G,1))
cat(sprintf("Delta BIC (G=1 minus G=2): %.1f  (positive => 1 component preferred)\n", bic_by_G["1"]-bic_by_G["2"]))

cat("\n================ 3. CONSENSUS CLUSTERING (ConsensusClusterPlus): cluster stability ================\n")
set.seed(123)
tmpdir <- tempfile(); dir.create(tmpdir)
ccp <- ConsensusClusterPlus(t(Z), maxK=6, reps=100, pItem=0.8, pFeature=1,
                            clusterAlg="km", distance="euclidean", seed=123,
                            plot=NULL, verbose=FALSE)
# PAC (proportion of ambiguous clustering); lower = more stable/real clusters
pac <- sapply(2:6, function(k){ M<-ccp[[k]]$consensusMatrix; v<-M[lower.tri(M)]
  (sum(v<=0.9)-sum(v<0.1))/length(v) })
names(pac) <- paste0("k=",2:6)
cat("PAC by k (lower=cleaner clusters; high & flat => no real cluster structure):\n"); print(round(pac,3))

cat("\n================ 4. SCORE ROBUSTNESS: PC1 vs SumZ vs RankComposite ================\n")
cm_p <- cor(scores[,c("PC1","SumZ","RankComposite")], method="pearson")
cm_s <- cor(scores[,c("PC1","SumZ","RankComposite")], method="spearman")
cat("Pearson:\n"); print(round(cm_p,3))
cat("Spearman:\n"); print(round(cm_s,3))
# median-split concordance (kappa) among the 3 scores
splits <- sapply(scores[,c("PC1","SumZ","RankComposite")], function(x) x > median(x))
ck <- function(a,b){ t<-table(a,b); sum(diag(t))/sum(t) }
cat(sprintf("Median-split agreement: PC1~SumZ=%.1f%%, PC1~Rank=%.1f%%, SumZ~Rank=%.1f%%\n",
            100*ck(splits[,1],splits[,2]), 100*ck(splits[,1],splits[,3]), 100*ck(splits[,2],splits[,3])))

# save results
res <- list(dip_PC1_p=dt$p.value, dip_PC1_D=as.numeric(dt$statistic),
            GMM_bestG=mc$G, GMM_model=mc$modelName, deltaBIC_1minus2=as.numeric(bic_by_G["1"]-bic_by_G["2"]),
            PAC_k2=pac["k=2"], pearson_PC1_Sum=cm_p["PC1","SumZ"], pearson_PC1_Rank=cm_p["PC1","RankComposite"])
saveRDS(list(scores=scores, dip=dt, mclust=mc, bic_by_G=bic_by_G, pac=pac, cor=cm_p),
        file.path(FD,"structure_assessment.rds"))
write.csv(scores, file.path(FD,"three_scores.csv"), row.names=FALSE)
cat("\nSaved structure_assessment.rds + three_scores.csv\n")
cat("\n=== HEADLINE ===\n")
cat(sprintf("dip p=%.3f (%s); GMM best G=%d; PAC(k=2)=%.2f; PC1~Sum r=%.3f, PC1~Rank r=%.3f\n",
            dt$p.value, ifelse(dt$p.value>0.05,"unimodal","multimodal"), mc$G, pac["k=2"],
            cm_p["PC1","SumZ"], cm_p["PC1","RankComposite"]))
