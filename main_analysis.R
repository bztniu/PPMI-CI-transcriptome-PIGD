# =====================================================================
# Peripheral blood mitochondrial complex I transcriptomic axis in
# Parkinson's disease: complete analysis code (single file)
# =====================================================================
# Manuscript: A peripheral-blood mitochondrial complex I (CI)
#   transcriptomic axis and late PIGD progression in Parkinson's disease
# Repository: https://github.com/bztniu/PPMI-CI-transcriptome-PIGD
#
# This single file consolidates the full analysis pipeline in
# execution order. Each SECTION corresponds to one originally
# standalone script and is independently runnable (scripts read
# intermediate tables produced by earlier sections).
#
# Notes:
#  - PPMI individual-level data are controlled-access and are NOT
#    redistributed. Update hard-coded paths to your own PPMI download.
#  - Sections 05, 06, 10, 11, 12 and 23 are the revision_v2 final
#    analyses (covariate-adjusted DESeq2, adjusted GSEA/GO, GENEPARK
#    independent validation, deduplicated LMM, and the updated
#    Figure 2 with the three cell-composition-robust risk genes).
#  - Runtime environment: R 4.5.1; package dependencies are installed
#    in Section 01.
# =====================================================================

# =====================================================================
# SECTION 01 | Package dependencies
# =====================================================================

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

# =====================================================================
# SECTION 02 | Data preparation and CI-group assignment
# =====================================================================

# ============================================================
# Figure Data Prep: PC1 score, DE, VST, immune, LMM-EMM, ML
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(nnls); library(preprocessCore)
  library(lme4); library(lmerTest); library(emmeans); library(glmnet)
})
setwd("E:/PPMI帕金森数据库专用")
OUT_DIR <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2"
FD <- file.path(OUT_DIR, "figures_pc1/data"); dir.create(FD, FALSE, TRUE)

ci_df <- read.csv("complex/CI_genes_converted.csv", stringsAsFactors=FALSE)
pd_data <- read.csv(file.path(OUT_DIR, "data/PD_all_clustering_methods.csv"))
dds_bl <- readRDS("dds_pd_control_BL_object.rds")
dds_ensg <- sub("\\..*","",rownames(dds_bl))
ci_matched <- intersect(ci_df$original_id, dds_ensg)
pd_samples <- intersect(colnames(dds_bl), pd_data$SAMPLE_ID)
dds_idx <- match(ci_matched, dds_ensg)

# PC1
norm_ci <- counts(dds_bl, normalized=TRUE)[dds_idx, pd_samples, drop=FALSE]
log2_ci <- log2(norm_ci + 1)
rownames(log2_ci) <- ci_df$gene_symbol[match(ci_matched, ci_df$original_id)]
pca <- prcomp(t(log2_ci), center=TRUE, scale.=TRUE)
pc1 <- pca$x[,"PC1"]
if(cor(pc1, colSums(log2_ci)) < 0) { pc1 <- -pc1; pca$rotation[,"PC1"] <- -pca$rotation[,"PC1"]; pca$x[,"PC1"] <- -pca$x[,"PC1"] }
med <- median(pc1)
grp <- ifelse(pc1 < med, "Low", "High"); names(grp) <- pd_samples

# Save PC1 score + loadings + variance
saveRDS(pca, file.path(FD,"pca.rds"))
write.csv(data.frame(SAMPLE_ID=pd_samples, PATNO=pd_data$PATNO[match(pd_samples,pd_data$SAMPLE_ID)],
          PC1=pc1, group=grp), file.path(FD,"pc1_scores.csv"), row.names=FALSE)
write.csv(data.frame(Gene=rownames(pca$rotation), PC1_loading=pca$rotation[,"PC1"]),
          file.path(FD,"pc1_loadings.csv"), row.names=FALSE)
var_exp <- summary(pca)$importance[2,1:min(10,ncol(pca$x))]
write.csv(data.frame(PC=names(var_exp), VarExplained=var_exp), file.path(FD,"pca_variance.csv"), row.names=FALSE)
write.csv(t(log2_ci), file.path(FD,"ci_expr.csv"))  # for heatmap
cat(sprintf("PC1 done: Low=%d High=%d, var=%.1f%%\n", sum(grp=="Low"), sum(grp=="High"), var_exp[1]*100))

# ===== DESeq2 full DE (PC1) for volcano =====
cat("DESeq2...\n")
dds_pd <- dds_bl[, pd_samples]
keep <- rowSums(counts(dds_pd)>=10) >= ncol(dds_pd)*0.2
dds_pd <- dds_pd[keep,]
colData(dds_pd)$PC1_group <- factor(grp[colnames(dds_pd)], levels=c("High","Low"))
design(dds_pd) <- ~ PC1_group
dds_pd <- DESeq(dds_pd)
res <- results(dds_pd, contrast=c("PC1_group","High","Low"), alpha=0.05)
de <- as.data.frame(res); de$ENSG <- sub("\\..*","",rownames(de))
de$Symbol <- mapIds(org.Hs.eg.db, keys=de$ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
write.csv(de, file.path(FD,"de_full.csv"), row.names=FALSE)
cat(sprintf("DE saved: %d genes\n", nrow(de)))

# VST expression of panel genes for heatmap
panel <- c("MAPT","GPNMB","SCARB2","TMEM175","MMP16","SNCA","LRRK2","GBA1",
  "DMAC1","NDUFAF8","TIMMDC1","TMEM186","NDUFS6","NDUFB11",
  "GFAP","IL6","AIF1","SNAP25","SYT1","SYP","DLG4","BACE1","APP")
vsd <- vst(dds_pd, blind=TRUE)
vmat <- assay(vsd); rownames(vmat) <- sub("\\..*","",rownames(vmat))
panel_ensg <- de$ENSG[match(panel, de$Symbol)]
pv <- vmat[panel_ensg[!is.na(panel_ensg)], ]
rownames(pv) <- panel[!is.na(panel_ensg)]
write.csv(pv, file.path(FD,"panel_vst.csv"))
write.csv(data.frame(SAMPLE_ID=colnames(pv), group=grp[colnames(pv)]), file.path(FD,"heatmap_anno.csv"), row.names=FALSE)
cat("VST panel saved\n")

# ===== NNLS immune (per-sample fractions) =====
cat("Immune NNLS...\n")
ref <- as.matrix(read.table("结果/免疫浸润/refer.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE))
all_nv <- sub("\\..*","",rownames(counts(dds_bl)))
syms <- mapIds(org.Hs.eg.db, keys=all_nv, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
nc_all <- counts(dds_bl, normalized=TRUE)[, pd_samples]
valid <- !is.na(syms) & syms!="" & !duplicated(syms)
norm_sym <- nc_all[valid,]; rownames(norm_sym) <- syms[valid]
common <- intersect(rownames(ref), rownames(norm_sym))
ref_s <- ref[common,]; mix <- norm_sym[common,]
mix_qn <- normalize.quantiles(mix); rownames(mix_qn)<-rownames(mix); colnames(mix_qn)<-colnames(mix)
imm <- matrix(0, ncol(mix_qn), ncol(ref_s), dimnames=list(colnames(mix_qn), colnames(ref_s)))
for(i in 1:ncol(mix_qn)) { f<-nnls(ref_s, mix_qn[,i]); imm[i,]<-f$x/sum(f$x) }
imm_df <- as.data.frame(imm); imm_df$SAMPLE_ID <- rownames(imm); imm_df$group <- grp[rownames(imm)]
write.csv(imm_df, file.path(FD,"immune_fractions.csv"), row.names=FALSE)
# Cohen's d table
ic <- data.frame()
for(ct in colnames(imm)) {
  lo<-imm[grp[rownames(imm)]=="Low",ct]; hi<-imm[grp[rownames(imm)]=="High",ct]
  if(sd(c(lo,hi))>0) ic <- rbind(ic, data.frame(Cell=ct, d=(mean(lo)-mean(hi))/sd(c(lo,hi)),
    P=wilcox.test(lo,hi)$p.value,
    d_lo=NA, d_hi=NA, stringsAsFactors=FALSE))
}
# Bootstrap CI for Cohen's d
set.seed(1)
for(i in 1:nrow(ic)) {
  ct <- ic$Cell[i]; gg <- grp[rownames(imm)]
  ds <- replicate(500, {
    idx <- sample(seq_len(nrow(imm)), replace=TRUE)
    lo<-imm[idx,ct][gg[idx]=="Low"]; hi<-imm[idx,ct][gg[idx]=="High"]
    if(length(lo)>2 && length(hi)>2 && sd(c(lo,hi))>0) (mean(lo)-mean(hi))/sd(c(lo,hi)) else NA
  })
  ic$d_lo[i] <- quantile(ds, 0.025, na.rm=TRUE); ic$d_hi[i] <- quantile(ds, 0.975, na.rm=TRUE)
}
ic$Padj <- p.adjust(ic$P,"BH"); ic <- ic[order(ic$Padj),]
write.csv(ic, file.path(FD,"immune_cohend.csv"), row.names=FALSE)
cat(sprintf("Immune: %d/%d sig\n", sum(ic$Padj<0.05), nrow(ic)))

# ===== LMM + emmeans =====
cat("LMM emmeans...\n")
moca_l <- read_csv("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", show_col_types=FALSE) %>%
  dplyr::select(PATNO, EVENT_ID, MoCA=MCATOT) %>% filter(EVENT_ID %in% c("SC","V04","V06","V08","V10","V12"))
motor_l <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("BL","SC","V04","V06","V08","V10","V12")) %>%
  dplyr::arrange(PATNO, EVENT_ID) %>% dplyr::distinct(PATNO, EVENT_ID, .keep_all=TRUE) %>%
  mutate(PIGD=NP3GAIT+NP3PSTBL+NP3FRZGT) %>% dplyr::select(PATNO, EVENT_ID, PIGD)
pc1_pat <- data.frame(PATNO=pd_data$PATNO[match(pd_samples,pd_data$SAMPLE_ID)], group=grp)
long <- moca_l %>% full_join(motor_l, by=c("PATNO","EVENT_ID")) %>%
  inner_join(pc1_pat, by="PATNO") %>%
  mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
         VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))

emm_all <- data.frame(); int_all <- data.frame()
for(oc in c("PIGD","MoCA")) {
  sub <- long[!is.na(long[[oc]]),]
  m <- lmer(as.formula(paste(oc,"~ VISIT*group + (1|PATNO)")), data=sub, REML=TRUE)
  emm <- as.data.frame(emmeans(m, ~ VISIT|group))
  emm$Outcome <- oc; emm_all <- rbind(emm_all, emm)
  # Interaction estimates from summary
  cm <- coef(summary(m))
  for(v in c("V04","V06","V08","V10","V12")) {
    rn <- grep(paste0("VISIT",v,":groupLow"), rownames(cm), fixed=TRUE)
    if(length(rn)>0) int_all <- rbind(int_all, data.frame(Outcome=oc, Visit=v,
      Est=cm[rn,"Estimate"], SE=cm[rn,"Std. Error"], P=cm[rn,"Pr(>|t|)"]))
  }
}
write.csv(emm_all, file.path(FD,"lmm_emm.csv"), row.names=FALSE)
write.csv(int_all, file.path(FD,"lmm_interactions.csv"), row.names=FALSE)
# raw means for overlay
raw <- long %>% tidyr::pivot_longer(c(PIGD,MoCA), names_to="Outcome", values_to="val") %>%
  filter(!is.na(val)) %>% group_by(Outcome, VISIT, group) %>%
  summarise(mean=mean(val), se=sd(val)/sqrt(n()), .groups="drop")
write.csv(raw, file.path(FD,"raw_means.csv"), row.names=FALSE)
cat("LMM emmeans saved\n")

# ===== ML predictions for ROC =====
cat("ML ROC...\n")
demo <- read_csv("Demographics_03Feb2026.csv", show_col_types=FALSE) %>% dplyr::select(PATNO,SEX) %>% distinct(PATNO,.keep_all=TRUE)
age <- read_csv("Age_at_visit_17Mar2025.csv", show_col_types=FALSE) %>% filter(EVENT_ID=="BL") %>% dplyr::select(PATNO,AGE_AT_VISIT) %>% distinct(PATNO,.keep_all=TRUE)
motor3 <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>% filter(EVENT_ID %in% c("BL","SC")) %>% arrange(PATNO,EVENT_ID) %>% distinct(PATNO,.keep_all=TRUE) %>% dplyr::select(PATNO,NP3TOT)
moca3 <- read_csv("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", show_col_types=FALSE) %>% filter(EVENT_ID=="SC") %>% dplyr::select(PATNO,MCATOT) %>% distinct(PATNO,.keep_all=TRUE)
ci <- pc1_pat
ci$PC1 <- pc1  # pc1_pat and pc1 are in pd_samples order
ci <- ci %>% left_join(demo,by="PATNO") %>% left_join(age,by="PATNO") %>% left_join(motor3,by="PATNO") %>% left_join(moca3,by="PATNO") %>%
  left_join(pd_data %>% dplyr::select(PATNO,MCATOT_V12) %>% distinct(PATNO,.keep_all=TRUE), by="PATNO") %>%
  filter(!is.na(MCATOT), MCATOT>=26, !is.na(MCATOT_V12), !is.na(AGE_AT_VISIT), !is.na(NP3TOT))
ci$impair <- ci$MCATOT_V12 < 26
set.seed(42); tr <- sample(1:nrow(ci), round(nrow(ci)*0.7)); te <- setdiff(1:nrow(ci),tr)
Xc <- model.matrix(~AGE_AT_VISIT+SEX+NP3TOT+MCATOT, ci)[,-1]
Xb <- cbind(Xc, PC1=ci$PC1)
roc_df <- data.frame(impair=ci$impair[te])
cv1 <- cv.glmnet(Xc[tr,], as.factor(ci$impair[tr]), family="binomial", alpha=0)
roc_df$Clinical <- predict(cv1, Xc[te,], s="lambda.min", type="response")[,1]
cv2 <- cv.glmnet(cbind(ci$PC1,rnorm(nrow(ci),0,1e-4))[tr,], as.factor(ci$impair[tr]), family="binomial", alpha=0)
roc_df$CIscore <- predict(cv2, cbind(ci$PC1,rnorm(nrow(ci),0,1e-4))[te,], s="lambda.min", type="response")[,1]
cv3 <- cv.glmnet(Xb[tr,], as.factor(ci$impair[tr]), family="binomial", alpha=0)
roc_df$Combined <- predict(cv3, Xb[te,], s="lambda.min", type="response")[,1]
write.csv(roc_df, file.path(FD,"roc_pred.csv"), row.names=FALSE)
cat("ML ROC saved\n")

cat("\n=== ALL DATA PREP DONE ===\n")

# =====================================================================
# SECTION 03 | Distribution structure assessment (dip test, GMM, consensus clustering)
# =====================================================================

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

# =====================================================================
# SECTION 04 | PCA on 66 CI genes and eigen-correlations
# =====================================================================

# ============================================================
# Figure 1 via PCAtools: screeplot + biplot(ellipse) + loadings + eigencorplot
# ============================================================
suppressPackageStartupMessages({
  library(PCAtools); library(ggplot2); library(patchwork); library(dplyr); library(readr); library(ggplotify)
})
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
FG <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

# ---- expression matrix: genes x samples ----
expr <- read.csv(file.path(FD,"ci_expr.csv"), row.names=1, check.names=FALSE)  # samples x genes
mat <- t(as.matrix(expr))                                                       # genes x samples
pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))                                # SAMPLE_ID,PATNO,PC1,group

# ---- assemble clinical metadata (samples x variables) ----
demo <- read_csv("Demographics_03Feb2026.csv", show_col_types=FALSE) %>% dplyr::select(PATNO,SEX) %>% distinct(PATNO,.keep_all=TRUE)
age  <- read_csv("Age_at_visit_17Mar2025.csv", show_col_types=FALSE) %>% filter(EVENT_ID=="BL") %>% dplyr::select(PATNO,Age=AGE_AT_VISIT) %>% distinct(PATNO,.keep_all=TRUE)
moca <- read_csv("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("SC","V12")) %>% tidyr::pivot_wider(id_cols=PATNO, names_from=EVENT_ID, values_from=MCATOT) %>%
  dplyr::rename(MoCA_BL=SC, MoCA_Y5=V12)
m3 <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("BL","SC","V12")) %>% arrange(PATNO,EVENT_ID) %>%
  mutate(PIGD=NP3GAIT+NP3PSTBL+NP3FRZGT,
         RestTremor=NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL)
m3_bl <- m3 %>% filter(EVENT_ID %in% c("BL","SC")) %>% distinct(PATNO,.keep_all=TRUE) %>%
  dplyr::select(PATNO, UPDRSIII=NP3TOT, PIGD_BL=PIGD, RestTremor_BL=RestTremor)
m3_v12 <- m3 %>% filter(EVENT_ID=="V12") %>% distinct(PATNO,.keep_all=TRUE) %>% dplyr::select(PATNO, PIGD_Y5=PIGD)

meta <- pc1s %>% left_join(age,by="PATNO") %>% left_join(moca,by="PATNO") %>%
  left_join(m3_bl,by="PATNO") %>% left_join(m3_v12,by="PATNO")
meta$CI_group <- factor(meta$group, levels=c("Low","High"))
rownames(meta) <- meta$SAMPLE_ID
meta <- meta[colnames(mat), ]   # align to expression columns

# ---- PCAtools PCA (scale+center, same as primary analysis) ----
p <- pca(mat, metadata=meta, scale=TRUE, center=TRUE)
# enforce sign: low CI expression -> negative PC1
if(cor(p$rotated$PC1, colMeans(mat)) < 0){ p$rotated$PC1 <- -p$rotated$PC1; p$loadings$PC1 <- -p$loadings$PC1 }
v1 <- round(p$variance["PC1"],1); v2 <- round(p$variance["PC2"],1)
LOW <- "#D55E00"; HIGH <- "#0072B2"

base_thm <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.line=element_line(linewidth=0.35), axis.ticks=element_line(linewidth=0.35),
        axis.title=element_text(size=8,face="bold"), axis.text=element_text(size=7),
        plot.title=element_text(size=9,face="bold"), legend.text=element_text(size=7.5),
        legend.title=element_blank(), panel.grid=element_blank(), plot.tag=element_text(size=11,face="bold"))

# ---- A: biplot (samples, color by CI group, 95% ellipse) ----
pA <- biplot(p, colby="CI_group", colkey=c(Low=LOW, High=HIGH),
             ellipse=TRUE, ellipseLevel=0.95, ellipseFill=TRUE, ellipseAlpha=0.10,
             lab=NULL, pointSize=1.1, legendPosition="top",
             xlab=sprintf("PC1 (%.1f%%)",v1), ylab=sprintf("PC2 (%.1f%%)",v2),
             title="CI transcriptomic axis (66 complex I genes)") +
  base_thm + theme(legend.position=c(0.86,0.92)) + labs(tag="A")

# ---- B: screeplot ----
pB <- screeplot(p, components=getComponents(p,1:8), axisLabSize=7, titleLabSize=9,
                colBar="grey70", drawCumulativeSumLine=FALSE) +
  base_thm + labs(title="Variance explained", tag="B") +
  theme(axis.text.x=element_text(angle=0, size=6))

# ---- C: plotloadings (top genes on PC1) ----
pC <- plotloadings(p, components="PC1", rangeRetain=0.18, absolute=FALSE,
                   labSize=2.2, shape=21, col=c("#D95F0E","grey90","#2C7FB8"),
                   drawConnectors=TRUE, title="PC1 gene loadings", subtitle=NULL, caption=NULL,
                   legendPosition="none") +
  base_thm + labs(tag="C")

# ---- D: eigencorplot (PC1-5 vs clinical phenotypes) ----
metavars <- c("Age","UPDRSIII","PIGD_BL","RestTremor_BL","MoCA_BL","PIGD_Y5","MoCA_Y5")
pD <- eigencorplot(p, components=getComponents(p,1:5), metavars=metavars,
                   col=c("#2166AC","#92C5DE","white","#F4A582","#B2182B"),
                   cexCorval=0.7, fontCorval=1, rotLabX=45, scale=TRUE,
                   main="PC vs clinical phenotype", cexMain=1,
                   corFUN="pearson", corUSE="pairwise.complete.obs",
                   signifSymbols=c("***","**","*",""), signifCutpoints=c(0,0.001,0.01,0.05,1),
                   returnPlot=TRUE) 

# eigencorplot returns a lattice/trellis object -> convert to ggplot for patchwork
pD_gg <- ggplotify::as.ggplot(pD) + labs(tag="D")
fig_top <- pA + (pB / pC) + plot_layout(widths=c(1.5,1))
fig <- fig_top / pD_gg + plot_layout(heights=c(2,1))
w <- 183/25.4; h <- 230/25.4

ragg::agg_tiff(file.path(FG,"Figure1_PCA.tiff"), width=w, height=h, units="in", res=600)
print(fig); dev.off()
grDevices::cairo_pdf(file.path(FG,"Figure1_PCA.pdf"), width=w, height=h)
print(fig); dev.off()
ragg::agg_png(file.path(FG,"Figure1_PCA_preview.png"), width=w, height=h, units="in", res=150)
print(fig); dev.off()
cat("Saved Figure1_PCA (.pdf/.tiff/.png). PC1 var =", v1, "%\n")
cat("eigencorplot vars:", paste(metavars,collapse=", "), "\n")

# =====================================================================
# SECTION 05 | DESeq2 covariate-adjusted differential expression (final model: ~ RIN + Plate + Age + Sex + PC1_group)
# (revision_v2 final DE model)
# =====================================================================

# Strictest sensitivity: ~ RIN + Plate + Age + Sex + PC1_group
suppressPackageStartupMessages({ library(DESeq2); library(dplyr) })
BASE <- "E:/PPMI帕金森数据库专用"
OUT  <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2"

dds_bl <- readRDS(file.path(BASE, "dds_pd_control_BL_object.rds"))
pd <- read.csv(file.path(OUT, "data/PD_all_clustering_methods.csv"))
pc1s <- read.csv(file.path(OUT, "figures_pc1/data/pc1_scores.csv"))
meta <- read.csv(file.path(BASE, "metaDataIR3.csv"), check.names = FALSE)

pd_samples <- intersect(colnames(dds_bl), pd$SAMPLE_ID)
dds_pd <- dds_bl[, pd_samples]

med <- median(pc1s$PC1)
grp <- ifelse(pc1s$PC1[match(pd_samples, pc1s$SAMPLE_ID)] < med, "Low", "High")
colData(dds_pd)$PC1_group <- factor(grp, levels = c("High","Low"))

# technical + clinical covariates
bc <- as.character(meta[["Specimen Bar Code"]])
idx <- match(pd_samples, bc)
colData(dds_pd)$RIN   <- as.numeric(meta[["RIN Value"]])[idx]
colData(dds_pd)$Plate <- factor(as.character(meta[["Plate"]])[idx])
# age/sex from clinical data
age <- as.numeric(pd$AGE_AT_VISIT[match(pd_samples, pd$SAMPLE_ID)])
sex <- as.character(pd$SEX[match(pd_samples, pd$SAMPLE_ID)])
colData(dds_pd)$Age <- age
colData(dds_pd)$Sex <- factor(sex)
cat("samples:", ncol(dds_pd), " age NA:", sum(is.na(age)), " sex NA:", sum(is.na(sex)), "\n")

keep <- rowSums(counts(dds_pd) >= 10) >= ncol(dds_pd) * 0.2
d3 <- dds_pd[keep, ]
d3 <- d3[, complete.cases(colData(d3)[, c("RIN","Plate","Age","Sex")])]
cat("Model3 samples:", ncol(d3), " plates:", nlevels(colData(d3)$Plate), "\n")

design(d3) <- ~ RIN + Plate + Age + Sex + PC1_group
d3 <- tryCatch(DESeq(d3, quiet = TRUE), error = function(e) { cat("Model3 error:", conditionMessage(e), "\n"); NULL })
if (!is.null(d3)) {
  res3 <- results(d3, contrast = c("PC1_group","High","Low"), alpha = 0.05)
  res3 <- res3[!is.na(res3$padj), ]
  n3_low <- sum(res3$padj < 0.05 & res3$log2FoldChange < -1)
  n3_high <- sum(res3$padj < 0.05 & res3$log2FoldChange > 1)
  cat(sprintf("[RIN+Plate+Age+Sex] low-CI up: %d, high-CI up: %d, total: %d\n", n3_low, n3_high, n3_low + n3_high))
  write.csv(as.data.frame(res3), file.path(OUT, "FINAL_PACKAGE_2026-06-21/04_tables/DE_covariate_adjusted_full.csv"), row.names = TRUE)
  cat("Saved DE_covariate_adjusted_full.csv\n")
}
cat("DONE\n")

# =====================================================================
# SECTION 06 | KEGG-GSEA and GO enrichment on covariate-adjusted DE statistics
# (revision_v2 adjusted enrichment)
# =====================================================================

# Recompute GSEA + GO enrichment from covariate-adjusted DE (Scheme B)
suppressPackageStartupMessages({ library(dplyr); library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi) })
TAB <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"
DATA <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"

adj <- read.csv(file.path(TAB, "DE_covariate_adjusted_full.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(adj)[1] <- "ENSG"
adj$ENSG <- sub("\\..*", "", adj$ENSG)
adj$Symbol <- mapIds(org.Hs.eg.db, keys = adj$ENSG, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
adj <- adj[!is.na(adj$Symbol) & !is.na(adj$stat), ]

# ---- GSEA ranked list (same convention: rank = -Wald stat, positive NES = low-CI direction) ----
de_rank <- adj[order(-(-adj$stat)), c("Symbol", "stat")]
de_rank$rank_score <- -de_rank$stat
rnk <- setNames(de_rank$rank_score, de_rank$Symbol)
rnk <- rnk[!duplicated(names(rnk)) & !is.na(rnk)]

# KEGG GSEA
suppressPackageStartupMessages(library(msigdbr))
run_gsea_adj <- function(collection, subcollection, out_name) {
  term2gene <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection) %>%
    dplyr::select(gs_name, gene_symbol) %>% distinct()
  set.seed(20260824)
  gsea <- tryCatch(clusterProfiler::GSEA(geneList = rnk, TERM2GENE = term2gene,
                   minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0, verbose = FALSE),
                   error = function(e) { cat(out_name, "err:", conditionMessage(e), "\n"); NULL })
  if (!is.null(gsea)) {
    out <- as.data.frame(gsea)
    out <- out[order(out$p.adjust, -abs(out$NES)), ]
    write.csv(out, file.path(TAB, out_name), row.names = FALSE)
    cat(out_name, ":", nrow(out), "terms\n")
  }
  gsea
}
g_go <- run_gsea_adj("C5", "GO:BP", "GSEA_GO_BP_lowCI_ranked_ADJUSTED.csv")
g_kg <- run_gsea_adj("C2", "CP:KEGG_LEGACY", "GSEA_KEGG_lowCI_ranked_ADJUSTED.csv")
if (!is.null(g_go)) {
  gd <- as.data.frame(g_go)
  cat("GO BP FDR<0.05 low-CI:", sum(gd$p.adjust < 0.05 & gd$NES > 0, na.rm = TRUE), "\n")
  print(head(gd[gd$p.adjust < 0.05 & gd$NES > 0, c("Description","NES","p.adjust")], 8))
}
if (!is.null(g_kg)) {
  kd <- as.data.frame(g_kg)
  cat("KEGG FDR<0.05 low-CI:", sum(kd$p.adjust < 0.05 & kd$NES > 0, na.rm = TRUE), "\n")
  print(head(kd[kd$p.adjust < 0.05 & kd$NES > 0, c("Description","NES","p.adjust")], 8))
}

# ---- GO enrichment (low-CI-upregulated: padj<0.05 & log2FC < -1) ----
low_up <- adj[adj$padj < 0.05 & adj$log2FoldChange < -1, ]
cat("\nAdjusted low-CI-up DEGs for GO enrich:", nrow(low_up), "\n")
if (nrow(low_up) >= 10) {
  go_up <- enrichGO(gene = unique(low_up$ENSG), OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (!is.null(go_up)) {
    gd <- as.data.frame(go_up)
    gd$geneID <- NULL
    write.csv(gd, file.path(TAB, "GO_lowCI_up_ADJUSTED.csv"), row.names = FALSE)
    cat("GO enrich terms:", nrow(gd), "\n")
    print(head(gd[order(gd$p.adjust), c("Description","Count","p.adjust")], 10))
  }
}
cat("DONE\n")

# =====================================================================
# SECTION 07 | PD risk-gene screen (Nalls 2019)
# =====================================================================

# ============================================================
# Full PD risk gene screen in PPMI DESeq2 + GENEPARK limma
# ============================================================

# PPMI DESeq2 (High vs Low CI: negative = higher in Low CI)
ppmi <- read.csv("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/de_full.csv")
ppmi <- ppmi[!is.na(ppmi$padj) & !is.na(ppmi$Symbol), c("Symbol","log2FoldChange","pvalue","padj")]
names(ppmi) <- c("Symbol","PPMI_logFC","PPMI_pval","PPMI_padj")

# GENEPARK limma (Low vs High CI: positive = higher in Low CI)
gse <- read.csv("E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_limma_DE_PC1.csv")
gse <- gse[, c("Gene","logFC","adj.P.Val")]
names(gse) <- c("Symbol","GSE_logFC","GSE_padj")

# PD risk genes from Nalls 2019 GWAS + Mendelian genes
pd_genes <- c(
  # Mendelian
  "SNCA","LRRK2","PRKN","PINK1","PARK7","DJ1","ATP13A2","FBXO7",
  "VPS35","CHCHD2","VPS13C","PLA2G6","DNAJC6","SYNJ1","EIF4G1",
  # Strong risk
  "GBA1","GBA","TMEM175","GCH1",
  # GWAS top hits
  "MAPT","GPNMB","MMP16","SCARB2","BST1","STK39","MCCC1",
  "SYT11","ACMSD","FGF20","ITGA8","CTSB","HLA-DRB5","RAB29",
  "PM20D1","GAK","DDRGK1","UBAP2L","SIPA1L2","INPP5F",
  "BAG3","CNTNAP2","ELOVL7","SCN2A","RICTOR","MED13","BCKDK",
  "CCDC62","MEX3C","ASXL3","WNT3","STX1B","UBQLN1",
  "RIT2","SH3GL2","GLT8D1","DMPK","SATB1"
)

# Filter to genes actually in PPMI
df <- merge(ppmi, gse, by="Symbol", all=FALSE)
df <- df[df$Symbol %in% pd_genes, ]

# Label whether meets threshold
df$PPMI_sig <- (df$PPMI_padj < 0.05 & abs(df$PPMI_logFC) > 1)
df$PPMI_nominal <- (df$PPMI_padj < 0.05 & abs(df$PPMI_logFC) <= 1)
df$GSE_sig <- df$GSE_padj < 0.05

# Direction: PPMI negative=higher in Low CI, GSE positive=higher in Low CI
# Concordant if both point same way
df$concordant <- sign(df$PPMI_logFC) * sign(df$GSE_logFC) < 0

# Sort by PPMI abs logFC
df <- df[order(-abs(df$PPMI_logFC)), ]

cat("=== Full PD Risk Gene Screen (PPMI + GENEPARK) ===\n\n")
cat(sprintf("Total PD genes queried: %d\n", length(pd_genes)))
cat(sprintf("Found in both PPMI & GSE: %d\n", nrow(df)))
cat(sprintf("PPMI |logFC|>1 & padj<0.05: %d\n", sum(df$PPMI_sig)))
cat(sprintf("PPMI nominal only: %d\n", sum(df$PPMI_nominal)))
cat(sprintf("GSE padj<0.05: %d\n", sum(df$GSE_sig)))
cat(sprintf("Direction concordant: %d/%d\n\n", sum(df$concordant), nrow(df)))

cat(sprintf("%-14s %8s %10s %6s %8s %9s %6s\n",
            "Gene","PPMI_FC","PPMI_padj","Sig?","GSE_FC","GSE_padj","Dir?"))
cat(strrep("-",80),"\n")
for(i in 1:nrow(df)){
  r <- df[i,]
  sig_lab <- ifelse(r$PPMI_sig, "****",
             ifelse(r$PPMI_nominal, "nom*", "ns"))
  dir_lab <- ifelse(r$concordant, "SAME", "OPP")
  cat(sprintf("%-14s %8.2f %10.2e %6s %8.2f %9.2e %6s\n",
      r$Symbol, r$PPMI_logFC, r$PPMI_padj, sig_lab, r$GSE_logFC, r$GSE_padj, dir_lab))
}

# Also list genes NOT found
miss <- setdiff(pd_genes, df$Symbol)
cat("\nGenes NOT in expression matrix:\n")
cat(paste(miss, collapse=", "),"\n")

write.csv(df, "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data/PD_risk_gene_screen.csv", row.names=FALSE)
cat("\nSaved: PD_risk_gene_screen.csv\n")

# =====================================================================
# SECTION 08 | ssGSEA whole-transcriptome validation
# =====================================================================

# Real ssGSEA: CI 66-gene set enrichment on full transcriptome (393 PD samples)
suppressPackageStartupMessages({ library(DESeq2); library(GSVA) })
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"

dds <- readRDS("dds_pd_control_BL_object.rds")
log2_expr <- log2(counts(dds, normalized=TRUE) + 1)            # full genes x 582 samples
rownames(log2_expr) <- sub("\\..*$", "", rownames(log2_expr))  # strip ENSG version

pd <- read.csv("E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/data/PD_all_clustering_methods.csv")   # 393 PD SAMPLE_IDs
pd_ids <- intersect(pd$SAMPLE_ID, colnames(log2_expr))
mat <- log2_expr[, pd_ids]                                      # full transcriptome x 393 PD
cat(sprintf("Full PD matrix: %d genes x %d samples\n", nrow(mat), ncol(mat)))

ci <- read.csv("complex/CI_genes_converted.csv", stringsAsFactors=FALSE)
ci_set <- intersect(ci$original_id, rownames(mat))
cat(sprintf("CI gene set matched in transcriptome: %d / %d\n", length(ci_set), nrow(ci)))
gs <- list(CI = ci_set)

# ssGSEA via GSVA (new ssgseaParam API, fallback to legacy)
ss <- tryCatch({
  par <- GSVA::ssgseaParam(as.matrix(mat), gs, normalize=TRUE)
  GSVA::gsva(par)
}, error=function(e){
  cat("new API failed (", conditionMessage(e), "), trying legacy...\n")
  GSVA::gsva(as.matrix(mat), gs, method="ssgsea", ssgsea.norm=TRUE, verbose=FALSE)
})
ssgsea <- as.numeric(ss["CI", ]); names(ssgsea) <- colnames(ss)
cat("ssGSEA computed for", length(ssgsea), "samples\n")

# merge with existing scores
sc <- read.csv(file.path(FD,"three_scores.csv"))               # SAMPLE_ID, PC1, SumZ, RankComposite
sc$ssGSEA <- ssgsea[match(sc$SAMPLE_ID, names(ssgsea))]
# sign align (increase with CI expression)
if(cor(sc$ssGSEA, sc$SumZ, use="complete.obs") < 0) sc$ssGSEA <- -sc$ssGSEA

cat("\n==== 3-SCORE CONCORDANCE (now with REAL ssGSEA) ====\n")
cc <- cor(sc[,c("PC1","SumZ","ssGSEA")], use="complete.obs", method="pearson")
print(round(cc,3))
cat(sprintf("PC1 vs ssGSEA: Pearson r=%.3f, Spearman=%.3f\n",
            cor(sc$PC1, sc$ssGSEA, use="complete.obs"),
            cor(sc$PC1, sc$ssGSEA, use="complete.obs", method="spearman")))
sp <- function(x) x > median(x, na.rm=TRUE)
agr <- mean(sp(sc$PC1)==sp(sc$ssGSEA), na.rm=TRUE)
cat(sprintf("Median-split agreement PC1 vs ssGSEA: %.1f%%\n", 100*agr))

write.csv(sc, file.path(FD,"three_scores.csv"), row.names=FALSE)
cat("\nUpdated three_scores.csv with REAL ssGSEA column\n")

# =====================================================================
# SECTION 09 | k-means partition immune robustness check
# =====================================================================

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

# =====================================================================
# SECTION 10 | GENEPARK independent PC1 re-estimation and CI stratification (GSE99039)
# (revision_v2 GENEPARK validation)
# =====================================================================

# GSE99039 — PC1 stratification + transcriptome + immune validation
suppressPackageStartupMessages({library(hgu133plus2.db); library(limma)})

meta <- read.csv("E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/sample_metadata.csv", stringsAsFactors=FALSE)
keep <- meta$disease_label %in% c("IPD","CONTROL"); meta_f <- meta[keep, ]
meta_pd <- meta_f[meta_f$disease_label=="IPD", ]  # PD only for PC1 split

expr <- read.csv("E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_expr.csv",
                 row.names=1, check.names=FALSE)
idx <- match(meta_pd$geo_accession, gsub('"',"",colnames(expr)))
expr_pd <- as.matrix(expr[, idx])
colnames(expr_pd) <- meta_pd$geo_accession
cat(sprintf("PD samples: %d, probes: %d\n", ncol(expr_pd), nrow(expr_pd)))

# --- CI gene expression ---
ci_df <- read.csv("E:/PPMI帕金森数据库专用/complex/CI_genes_converted.csv", stringsAsFactors=FALSE)
gm <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(expr_pd), columns="SYMBOL", keytype="PROBEID")
gm <- gm[!is.na(gm$SYMBOL), ]
ci_gm <- gm[gm$SYMBOL %in% ci_df$gene_symbol, ]
ci_expr <- expr_pd[ci_gm$PROBEID, , drop=FALSE]

ci_genes <- unique(ci_gm$SYMBOL)
ci_gene_expr <- t(sapply(ci_genes, function(g) {
  prows <- which(ci_gm$SYMBOL==g)
  if(length(prows)==1) ci_expr[prows,] else colMeans(ci_expr[prows,,drop=FALSE])
}))
rownames(ci_gene_expr) <- ci_genes
cat(sprintf("CI genes: %d\n", nrow(ci_gene_expr)))

# --- PC1 score ---
pca <- prcomp(t(ci_gene_expr), center=TRUE, scale.=TRUE)
pc1 <- pca$x[, "PC1"]
if(cor(pc1, colMeans(ci_gene_expr)) < 0) pc1 <- -pc1
cat(sprintf("PC1 var: %.1f%%\n", summary(pca)$importance[2,1]*100))

med <- median(pc1)
pc1_group <- ifelse(pc1 < med, "Low", "High")
cat(sprintf("Low=%d, High=%d\n", sum(pc1_group=="Low"), sum(pc1_group=="High")))

# --- DESeq2-like: genome-wide DE (PC1-High vs PC1-Low) ---
cat("\n=== Genome-wide DE: PC1-High vs PC1-Low ===\n")
grp <- factor(pc1_group, levels=c("High","Low"))
design <- model.matrix(~ grp)
fit <- lmFit(expr_pd, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef="grpLow", number=Inf, sort.by="none")
tt$SYMBOL <- gm$SYMBOL[match(rownames(tt), gm$PROBEID)]

# Collapse to gene level (keep most significant probe per gene)
tt_gene <- tt[!is.na(tt$SYMBOL) & tt$SYMBOL != "", ]
tt_gene <- tt_gene[order(tt_gene$P.Value), ]
tt_gene <- tt_gene[!duplicated(tt_gene$SYMBOL), ]
cat(sprintf("DE genes (P<0.01): %d\n", sum(tt_gene$P.Value < 0.01)))
cat(sprintf("DE genes (P<0.05): %d\n", sum(tt_gene$P.Value < 0.05)))

# --- Key gene panels ---
cat("\n=== PD Risk Genes (PC1-Low vs PC1-High) ===\n")
cat(sprintf("%-12s %8s %8s %8s %8s\n","Gene","logFC","P","dir",""))
risk <- c("MAPT","GPNMB","TMEM175","SCARB2","MMP16","SNCA","LRRK2","GBA1",
          "DMAC1","NDUFAF8","TIMMDC1","TMEM186","NDUFS6","NDUFB11",
          "GFAP","IL6","TNF","SNAP25","SYP","DLG4","BACE1","APP")
for(g in risk) {
  r <- tt_gene[tt_gene$SYMBOL==g, ]
  if(nrow(r)>0) {
    cat(sprintf("%-12s %+8.3f %8.4f %8s %s\n",
                g, r$logFC, r$P.Value,
                ifelse(r$logFC>0,"Low>High","High>Low"),
                ifelse(r$P.Value<0.05,"*","")))
  }
}

# --- CI assembly factors ---
cat("\n=== CI Assembly Factors ===\n")
cat(sprintf("%-12s %8s %8s %8s\n","Gene","logFC","P",""))
ci_assembly <- c("NDUFAF8","DMAC1","TIMMDC1","TMEM186","NDUFS6","NDUFB11",
                 "NDUFAF1","NDUFAF2","NDUFAF3","NDUFAF4","NDUFAF5","NDUFAF6",
                 "FOXRED1","NUBPL","ECSIT","TMEM126B","ACAD9","COA1")
for(g in ci_assembly) {
  r <- tt_gene[tt_gene$SYMBOL==g, ]
  if(nrow(r)>0) {
    cat(sprintf("%-12s %+8.3f %8.4f %s\n",
                g, r$logFC, r$P.Value, ifelse(r$P.Value<0.05,"*","")))
  }
}

# --- CIBERSORT-like immune: simple cell-type marker enrichment ---
cat("\n=== Immune Cell Markers ===\n")
immune_markers <- list(
  Neutrophils = c("ELANE","MPO","CEACAM8","OLFM4","CXCR2","FCGR3B","CXCL8"),
  CD4_Tcells = c("CD4","CD3E","CD3D","CD2","IL7R","CCR7","CD28","LCK"),
  CD8_Tcells = c("CD8A","CD8B","GZMB","PRF1","GNLY","NKG7"),
  B_cells = c("CD19","CD79A","CD79B","MS4A1","BLK","PAX5","CR2"),
  Monocytes = c("CD14","CD68","CSF1R","ITGAM","FCGR1A","CD163"),
  NK_cells = c("KLRD1","KLRF1","NCR1","NCR3","KLRC1","KIR2DL1"),
  Dendritic = c("CD1C","CLEC10A","FCER1A","CLEC4C","IL3RA","NRP1")
)
cat(sprintf("%-15s %10s %10s %10s\n","CellType","Direction","MeanLogFC","P"))
for(ct in names(immune_markers)) {
  markers <- immune_markers[[ct]]
  found <- tt_gene[tt_gene$SYMBOL %in% markers, ]
  if(nrow(found) >= 3) {
    mean_lfc <- mean(found$logFC)
    mean_p <- mean(found$P.Value)
    dir <- ifelse(mean_lfc > 0, "Low>High", "High>Low")
    sig_markers <- found[found$P.Value < 0.05, ]
    cat(sprintf("%-15s %10s %+10.3f %10.4f (%d/%d markers P<0.05)\n",
                ct, dir, mean_lfc, mean_p,
                nrow(sig_markers), nrow(found)))
  }
}

# --- Save ---
write.csv(tt_gene, "E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_PC1_DE.csv", row.names=FALSE)
cat("\nDone.\n")

# =====================================================================
# SECTION 11 | GENEPARK risk-gene DE and immune deconvolution validation
# (revision_v2 GENEPARK validation)
# =====================================================================

# ============================================================
# GSE99039 — NNLS Immune + limma DE following PPMI methods
# Mirrors revision_v2/pipeline/03_immune_nnls.R exactly
# ============================================================
suppressPackageStartupMessages({
  library(hgu133plus2.db); library(limma)
  library(nnls); library(preprocessCore)
})

OUT_DIR <- "E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results"

# ---- Load PD-only expression (gene-level) ----
meta <- read.csv("E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/sample_metadata.csv", stringsAsFactors=FALSE)
meta_pd <- meta[meta$disease_label=="IPD", ]
expr <- read.csv(file.path(OUT_DIR, "GSE99039_expr.csv"), row.names=1, check.names=FALSE)
idx <- match(meta_pd$geo_accession, gsub('"',"",colnames(expr)))
expr_pd <- as.matrix(expr[, idx]); colnames(expr_pd) <- meta_pd$geo_accession
cat(sprintf("PD samples: %d, probes: %d\n", ncol(expr_pd), nrow(expr_pd)))

# Probe -> Symbol mapping
gm <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(expr_pd), columns="SYMBOL", keytype="PROBEID")
gm <- gm[!is.na(gm$SYMBOL) & gm$SYMBOL != "", ]
gm <- gm[!duplicated(gm$PROBEID), ]

# Collapse probes to gene-level (mean) — full transcriptome
sym_expr <- rowsum(expr_pd[gm$PROBEID, ], gm$SYMBOL) /
            as.numeric(table(gm$SYMBOL)[unique(gm$SYMBOL)])
# Safer: aggregate by symbol mean
sym_list <- split(seq_len(nrow(gm)), gm$SYMBOL)
gene_expr <- t(sapply(names(sym_list), function(s) {
  rows <- gm$PROBEID[sym_list[[s]]]
  if(length(rows)==1) expr_pd[rows, ] else colMeans(expr_pd[rows, , drop=FALSE])
}))
cat(sprintf("Gene-level: %d genes\n", nrow(gene_expr)))

# ---- PC1 stratification (66 CI genes) ----
ci_df <- read.csv("E:/PPMI帕金森数据库专用/complex/CI_genes_converted.csv", stringsAsFactors=FALSE)
ci_in <- intersect(ci_df$gene_symbol, rownames(gene_expr))
ci_mat <- gene_expr[ci_in, ]
cat(sprintf("CI genes for PC1: %d\n", length(ci_in)))

pca <- prcomp(t(ci_mat), center=TRUE, scale.=TRUE)
pc1 <- pca$x[,"PC1"]
if(cor(pc1, colMeans(ci_mat)) < 0) pc1 <- -pc1
med <- median(pc1)
grp <- factor(ifelse(pc1 < med, "Low", "High"), levels=c("High","Low"))
cat(sprintf("PC1 var=%.1f%%, Low=%d, High=%d\n",
            summary(pca)$importance[2,1]*100, sum(grp=="Low"), sum(grp=="High")))

# ============================================================
# PART 1: NNLS Immune Deconvolution (LM22) — PPMI method
# ============================================================
cat("\n========== NNLS Immune Deconvolution (LM22) ==========\n")
ref <- as.matrix(read.table("E:/PPMI帕金森数据库专用/结果/免疫浸润/refer.txt",
                            header=TRUE, row.names=1, sep="\t", check.names=FALSE))
cat(sprintf("LM22 reference: %d genes x %d cells\n", nrow(ref), ncol(ref)))

# Common genes
common <- intersect(rownames(ref), rownames(gene_expr))
cat(sprintf("LM22 genes matched: %d/%d\n", length(common), nrow(ref)))
ref_sub <- ref[common, ]
mix_sub <- gene_expr[common, ]

# Quantile normalize (PPMI method)
mix_qn <- normalize.quantiles(mix_sub)
rownames(mix_qn) <- rownames(mix_sub); colnames(mix_qn) <- colnames(mix_sub)

# NNLS per sample
ns <- ncol(mix_qn); nc <- ncol(ref_sub)
imm <- matrix(0, ns, nc, dimnames=list(colnames(mix_qn), colnames(ref_sub)))
for(i in 1:ns) {
  fit <- nnls(ref_sub, mix_qn[,i])
  imm[i,] <- fit$x / sum(fit$x)
}
cat(sprintf("NNLS: %d samples x %d cells\n", ns, nc))

# Group comparison: Wilcoxon + BH + Cohen's d (PPMI method)
imm_comp <- data.frame()
for(ct in colnames(imm)) {
  lo <- imm[grp=="Low", ct]; hi <- imm[grp=="High", ct]
  p <- wilcox.test(lo, hi)$p.value
  d <- (mean(lo) - mean(hi)) / sd(c(lo, hi))
  imm_comp <- rbind(imm_comp, data.frame(Cell=ct, LowM=mean(lo), HighM=mean(hi),
                P=p, D=d, stringsAsFactors=FALSE))
}
imm_comp <- imm_comp[!is.nan(imm_comp$D), ]  # drop zero-fraction cells
imm_comp$Padj <- p.adjust(imm_comp$P, "BH")
imm_comp <- imm_comp[order(imm_comp$Padj), ]

cat("\n=== Immune cells (sorted by Padj) ===\n")
cat(sprintf("%-28s %8s %8s %8s %8s\n","Cell","LowMean","HighMean","Cohen_d","Padj"))
for(i in 1:nrow(imm_comp)) {
  r <- imm_comp[i, ]
  cat(sprintf("%-28s %8.4f %8.4f %+8.3f %8.4f %s\n",
              r$Cell, r$LowM, r$HighM, r$D, r$Padj, ifelse(r$Padj<0.05,"*","")))
}
cat(sprintf("\nSignificant cells (Padj<0.05): %d/%d\n", sum(imm_comp$Padj<0.05), nrow(imm_comp)))

write.csv(imm_comp, file.path(OUT_DIR, "GSE99039_immune_NNLS.csv"), row.names=FALSE)

# ============================================================
# PART 2: limma DE (PC1-Low vs PC1-High) — microarray equiv of DESeq2
# ============================================================
cat("\n========== limma DE: PC1-Low vs PC1-High ==========\n")
design <- model.matrix(~ grp)
fit <- eBayes(lmFit(gene_expr, design))
tt <- topTable(fit, coef="grpLow", number=Inf, sort.by="P")
cat(sprintf("DE genes padj<0.05: %d\n", sum(tt$adj.P.Val<0.05)))
cat(sprintf("DE genes padj<0.05 & |logFC|>0.5: %d\n",
            sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>0.5)))

# Key gene panel
cat("\n=== PD Risk Genes & Biomarkers (Low vs High) ===\n")
cat(sprintf("%-12s %8s %10s %10s %s\n","Gene","logFC","P","Padj","Dir"))
panel <- c("MAPT","GPNMB","MMP16","TMEM175","SCARB2","SNCA","LRRK2","GBA1",
           "DMAC1","NDUFAF8","TIMMDC1","TMEM186","NDUFS6","NDUFB11",
           "GFAP","IL6","SNAP25","SYP","DLG4","BACE1","APP")
for(g in panel) {
  if(g %in% rownames(tt)) {
    r <- tt[g, ]
    cat(sprintf("%-12s %+8.3f %10.2e %10.2e %s %s\n",
                g, r$logFC, r$P.Value, r$adj.P.Val,
                ifelse(r$logFC>0,"Low>High","High>Low"),
                ifelse(r$adj.P.Val<0.05,"*","")))
  }
}

write.csv(data.frame(Gene=rownames(tt), tt),
          file.path(OUT_DIR, "GSE99039_limma_DE_PC1.csv"), row.names=FALSE)
cat(sprintf("\nDone: %s\n", OUT_DIR))

# =====================================================================
# SECTION 12 | Unified PIGD LMM with deduplicated longitudinal data (primary specification)
# (revision_v2 final LMM (audit 2026-08-24))
# =====================================================================

# ============================================================
# Unified LMM rerun with deduplicated longitudinal data
# Fixes: SC->BL remapping duplicate-baseline issue (audit 2026-08-24)
# Produces single source of truth for PIGD LMM results
# ============================================================
options(warn = 1)
suppressPackageStartupMessages({
  library(dplyr); library(lme4); library(lmerTest)
})

BASE   <- "E:/PPMI帕金森数据库专用"
FD     <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
TABDIR <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"

# ---------- 1. Build longitudinal data (deduplicated) ----------
pc1s <- read.csv(file.path(FD, "pc1_scores.csv"), stringsAsFactors = FALSE)
pc1s$PC1z <- as.numeric(scale(pc1s$PC1))
grp_map <- setNames(pc1s$group, pc1s$PATNO)

motor <- read.csv(file.path(BASE, "运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv"),
                  fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
motor <- motor[motor$EVENT_ID %in% c("BL","SC","V04","V06","V08","V10","V12"), ]
motor <- motor[!is.na(motor$NP3GAIT) | !is.na(motor$NP3PSTBL) | !is.na(motor$NP3FRZGT), ]
motor$PIGD <- rowSums(motor[, c("NP3GAIT","NP3PSTBL","NP3FRZGT")], na.rm = FALSE)
motor <- motor[order(motor$PATNO, motor$EVENT_ID), ]
motor <- motor[!duplicated(motor[, c("PATNO","EVENT_ID")]), ]
motor <- motor[, c("PATNO","EVENT_ID","PIGD")]

moca <- read.csv(file.path(BASE, "运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv"),
                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
moca <- moca[moca$EVENT_ID %in% c("SC","V04","V06","V08","V10","V12"), ]
moca <- moca[order(moca$PATNO, moca$EVENT_ID), ]
moca <- moca[!duplicated(moca[, c("PATNO","EVENT_ID")]), ]
moca$MoCA <- moca$MCATOT
moca <- moca[, c("PATNO","EVENT_ID","MoCA")]

long <- merge(moca, motor, by = c("PATNO","EVENT_ID"), all = TRUE)
long <- long[long$PATNO %in% pc1s$PATNO, ]
long$group <- factor(grp_map[as.character(long$PATNO)], levels = c("High","Low"))
long$PC1z  <- pc1s$PC1z[match(long$PATNO, pc1s$PATNO)]
long$VISIT <- factor(ifelse(long$EVENT_ID == "SC", "BL", long$EVENT_ID),
                     levels = c("BL","V04","V06","V08","V10","V12"))

cat("=== BEFORE dedup (per VISIT rows) ===\n")
print(table(long$VISIT, useNA = "ifany"))
cat("total rows before dedup:", nrow(long), "\n")

# Dedup: one row per PATNO+VISIT, prefer complete (SC-derived) row
long <- long[order(long$PATNO, long$VISIT,
                   !is.na(long$MoCA), !is.na(long$PIGD)), ]
long <- long[!duplicated(long[, c("PATNO","VISIT")]), ]

cat("\n=== AFTER dedup (per VISIT rows) ===\n")
print(table(long$VISIT, useNA = "ifany"))
cat("total rows after dedup:", nrow(long), "\n")
cat("subjects with any PIGD:", length(unique(long$PATNO[!is.na(long$PIGD)])), "\n")

write.csv(long, file.path(TABDIR, "longitudinal_PIGD_dedup.csv"), row.names = FALSE)

# ---------- 2. PIGD factor-time LMM (primary) ----------
sub <- long[!is.na(long$PIGD), ]
cat("\nPIGD factor-time model: n_obs =", nrow(sub), ", n_subj =", length(unique(sub$PATNO)), "\n")
m_fact <- lmer(PIGD ~ VISIT * group + (1 | PATNO), data = sub, REML = TRUE)
cm <- as.data.frame(coef(summary(m_fact)))
rns <- rownames(cm)
v12_idx <- grep("VISITV12:groupLow", rns, fixed = TRUE)
cat("\n=== PIGD factor-time interaction (V12 x group) ===\n")
cat("Estimate:", cm[v12_idx, "Estimate"], "\n")
cat("SE:", cm[v12_idx, "Std. Error"], "\n")
cat("P:", cm[v12_idx, "Pr(>|t|)"], "\n")
cat("N_obs:", nrow(sub), "N_subj:", length(unique(sub$PATNO)), "\n")

# per-visit group means for Figure 3b
emm_fact <- emmeans_simple <- NULL
if (requireNamespace("emmeans", quietly = TRUE)) {
  library(emmeans)
  emm <- as.data.frame(emmeans(m_fact, ~ VISIT | group))
  write.csv(emm, file.path(TABDIR, "LMM_PIGD_emmeans_dedup.csv"), row.names = FALSE)
}

# ---------- 3. Continuous PC1 factor-time (sensitivity) ----------
m_cont <- lmer(PIGD ~ VISIT * PC1z + (1 | PATNO), data = sub, REML = TRUE)
cmc <- as.data.frame(coef(summary(m_cont)))
rns_c <- rownames(cmc)
v12c <- grep("VISITV12:PC1z", rns_c, fixed = TRUE)
cat("\n=== Continuous PC1 factor-time interaction (V12 x PC1z) ===\n")
cat("Estimate:", cmc[v12c, "Estimate"], "\n")
cat("SE:", cmc[v12c, "Std. Error"], "\n")
cat("P:", cmc[v12c, "Pr(>|t|)"], "\n")

# ---------- 4. Random slope vs intercept (S17) ----------
cat("\n=== Random slope comparison (AIC) ===\n")
# Continuous-time random slope (2 random effects) is identifiable with sparse visits;
# full factor-time random slopes (6 per subject) exceed the deduplicated n_obs.
sub$Time <- as.numeric(sub$VISIT) - 1  # 0,1,2,3,4,5
m_ri_time <- lmer(PIGD ~ Time * group + (1 | PATNO), data = sub, REML = TRUE,
                  control = lmerControl(optCtrl = list(maxfun = 20000)))
m_rs_time <- lmer(PIGD ~ Time * group + (1 + Time | PATNO), data = sub, REML = TRUE,
                  control = lmerControl(optCtrl = list(maxfun = 20000)))
cat("RI AIC:", AIC(m_ri_time), " RS AIC:", AIC(m_rs_time),
    " Delta:", AIC(m_rs_time) - AIC(m_ri_time), "\n")
cat("RS converged:", as.character(m_rs_time@optinfo$conv$opt == 0), "\n")

# continuous-time interaction (for reference)
cmt <- as.data.frame(coef(summary(m_ri_time)))
cat("Continuous Time x group interaction p:", cmt["Time:groupLow", "Pr(>|t|)"], "\n")

# ---------- 5. Save unified interaction table ----------
interactions <- data.frame(
  Outcome = "PIGD", Visit = "V12",
  Model = c("Factor-time binary group", "Continuous PC1_z"),
  Estimate = c(cm[v12_idx, "Estimate"], cmc[v12c, "Estimate"]),
  SE = c(cm[v12_idx, "Std. Error"], cmc[v12c, "Std. Error"]),
  P = c(cm[v12_idx, "Pr(>|t|)"], cmc[v12c, "Pr(>|t|)"]),
  N_obs = nrow(sub), N_subj = length(unique(sub$PATNO))
)
write.csv(interactions, file.path(TABDIR, "LMM_PIGD_interactions_dedup.csv"), row.names = FALSE)
cat("\nSaved LMM_PIGD_interactions_dedup.csv and LMM_PIGD_emmeans_dedup.csv\n")
cat("DONE\n")

# =====================================================================
# SECTION 13 | Pre-specified longitudinal outcome sweep (43 outcomes)
# =====================================================================

# ============================================================
# LMM sweep: PC1 group x time interaction across ALL longitudinal scales
# Factor-time (V12 interaction) + continuous-time, FDR corrected
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"

# PC1 groups (PATNO -> group)
pc1 <- read.csv(file.path(FD,"pc1_scores.csv"))[,c("PATNO","group")]
VIS <- c("SC","BL","V04","V06","V08","V10","V12")
vmap <- c(BL=0,SC=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# ---- helper to load a longitudinal total ----
get_long <- function(file, scorefun, name) {
  d <- suppressWarnings(read_csv(file, show_col_types=FALSE))
  d <- d %>% filter(EVENT_ID %in% VIS)
  d$Score <- scorefun(d)
  d %>% dplyr::select(PATNO, EVENT_ID, Score) %>% filter(!is.na(Score)) %>%
    mutate(Outcome=name)
}

scales <- list()
scales[["UPDRS_I_clin"]]  <- get_long("运动症状数据/MDS-UPDRS_Part_I_31Jan2026.csv", function(d) d$NP1RTOT, "UPDRS_I_clin")
scales[["UPDRS_I_pat"]]   <- get_long("运动症状数据/MDS-UPDRS_Part_I_Patient_Questionnaire_31Jan2026.csv", function(d) d$NP1PTOT, "UPDRS_I_pat")
scales[["UPDRS_II"]]      <- get_long("运动症状数据/MDS_UPDRS_Part_II__Patient_Questionnaire_31Jan2026.csv", function(d) d$NP2PTOT, "UPDRS_II")
scales[["UPDRS_IV"]]      <- get_long("运动症状数据/MDS-UPDRS_Part_IV__Motor_Complications_29Jan2026.csv", function(d) d$NP4TOT, "UPDRS_IV")
scales[["MoCA"]]          <- get_long("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv", function(d) d$MCATOT, "MoCA")

# Part III: total, NHY, PIGD, Tremor, RestTremor
p3 <- suppressWarnings(read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% VIS) %>% arrange(PATNO, EVENT_ID) %>% distinct(PATNO, EVENT_ID, .keep_all=TRUE)
mkp3 <- function(col) p3 %>% mutate(Score=.data[[col]]) %>% dplyr::select(PATNO,EVENT_ID,Score) %>% filter(!is.na(Score)) %>% mutate(Outcome=col)
scales[["UPDRS_III"]] <- mkp3("NP3TOT"); scales[["UPDRS_III"]]$Outcome <- "UPDRS_III"
scales[["NHY"]]       <- mkp3("NHY");    scales[["NHY"]]$Outcome <- "NHY"
p3 <- p3 %>% mutate(PIGD=NP3GAIT+NP3PSTBL+NP3FRZGT,
                    Tremor=NP3PTRMR+NP3PTRML+NP3KTRMR+NP3KTRML+NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL,
                    RestTremor=NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL)
scales[["PIGD"]]       <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=PIGD) %>% filter(!is.na(Score)) %>% mutate(Outcome="PIGD")
scales[["Tremor"]]     <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=Tremor) %>% filter(!is.na(Score)) %>% mutate(Outcome="Tremor")
scales[["RestTremor"]] <- p3 %>% dplyr::select(PATNO,EVENT_ID,Score=RestTremor) %>% filter(!is.na(Score)) %>% mutate(Outcome="RestTremor")

# RBDSQ total = sum(12 items) + (any neuro disorder ? 1 : 0)
rb <- suppressWarnings(read_csv("REM_Sleep_Behavior_Disorder_Questionnaire_03Feb2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% VIS)
it12 <- c("DRMVIVID","DRMAGRAC","DRMNOCTB","SLPLMBMV","SLPINJUR","DRMVERBL","DRMFIGHT","DRMUMV","DRMOBJFL","MVAWAKEN","DRMREMEM","SLPDSTRB")
neuro <- c("STROKE","HETRA","PARKISM","RLS","NARCLPSY","DEPRS","EPILEPSY","BRNINFM","CNSOTH")
rb$Score <- rowSums(rb[,it12], na.rm=TRUE) + as.integer(rowSums(rb[,neuro], na.rm=TRUE) > 0)
scales[["RBDSQ"]] <- rb %>% dplyr::select(PATNO,EVENT_ID,Score) %>% mutate(Outcome="RBDSQ")

# SCOPA-AUT total = sum SCAU1..SCAU21 (proxy)
sc <- suppressWarnings(read_csv("SCOPA-AUT_07Feb2026.csv", show_col_types=FALSE)) %>% filter(EVENT_ID %in% VIS)
scau_cols <- paste0("SCAU", 1:21); scau_cols <- intersect(scau_cols, colnames(sc))
sc_num <- sc[,scau_cols]; sc_num[sc_num==9] <- NA   # 9 = N/A
sc$Score <- rowSums(sapply(sc_num, as.numeric), na.rm=TRUE)
scales[["SCOPA_AUT"]] <- sc %>% dplyr::select(PATNO,EVENT_ID,Score) %>% mutate(Outcome="SCOPA_AUT")

# ---- Run LMM for each scale ----
res <- data.frame()
for(nm in names(scales)) {
  d <- scales[[nm]] %>% inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & !is.na(d$Score), ]
  if(nlevels(droplevels(d$VISIT)) < 3 || nrow(d) < 50) { cat(sprintf("SKIP %s (insufficient)\n", nm)); next }

  # Factor-time RI-LMM: V12 interaction + min interaction p across visits
  v12p <- NA; v12est <- NA; minp <- NA
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)) {
    cm <- coef(summary(mf))
    int_rows <- grep(":groupLow$", rownames(cm))
    if(length(int_rows)>0) {
      ips <- cm[int_rows,"Pr(>|t|)"]; minp <- min(ips, na.rm=TRUE)
      v12r <- grep("VISITV12:groupLow", rownames(cm), fixed=TRUE)
      if(length(v12r)>0) { v12p <- cm[v12r[1],"Pr(>|t|)"]; v12est <- cm[v12r[1],"Estimate"] }
    }
  }
  # Continuous-time RS-LMM: group main + time:group interaction
  grp_p <- NA; int_p <- NA; int_est <- NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)) {
    cm <- coef(summary(mc))
    gr <- grep("^groupLow$", rownames(cm)); ir <- grep("Time:groupLow", rownames(cm), fixed=TRUE)
    if(length(gr)>0) grp_p <- cm[gr[1],"Pr(>|t|)"]
    if(length(ir)>0) { int_p <- cm[ir[1],"Pr(>|t|)"]; int_est <- cm[ir[1],"Estimate"] }
  }
  res <- rbind(res, data.frame(Outcome=nm, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)),
    V12_int_est=v12est, V12_int_p=v12p, MinVisit_int_p=minp,
    Cont_grp_p=grp_p, Cont_int_est=int_est, Cont_int_p=int_p, stringsAsFactors=FALSE))
}

# FDR correction across outcomes (on V12 interaction p)
res$V12_int_padj <- p.adjust(res$V12_int_p, "BH")
res$Cont_int_padj <- p.adjust(res$Cont_int_p, "BH")
res <- res[order(res$V12_int_p),]

cat("\n========== LMM SWEEP RESULTS (sorted by V12 interaction p) ==========\n")
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)

cat("\n=== V12 interaction: nominal p<0.05 ===\n")
print(res[!is.na(res$V12_int_p) & res$V12_int_p<0.05, c("Outcome","V12_int_est","V12_int_p","V12_int_padj")], row.names=FALSE)
cat(sprintf("\nSurvive FDR (V12 padj<0.05): %s\n",
            paste(res$Outcome[!is.na(res$V12_int_padj) & res$V12_int_padj<0.05], collapse=", ")))
cat("\n=== Continuous group main effect: nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_grp_p) & res$Cont_grp_p<0.05, c("Outcome","Cont_grp_p")], row.names=FALSE)

write.csv(res, file.path(OUT,"LMM_sweep_all_outcomes.csv"), row.names=FALSE)
cat(sprintf("\nSaved: %s/LMM_sweep_all_outcomes.csv\n", OUT))

# =====================================================================
# SECTION 14 | Full longitudinal outcome sweep
# =====================================================================

# ============================================================
# COMPREHENSIVE LMM sweep: ~43 longitudinal indicators x PC1 group
# Totals + DAT-SBR imaging + MoCA subdomains + UPDRS subdomains/items + SCOPA domains
# Factor-time (last-visit interaction) + continuous-time, FDR across ALL
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4); library(lmerTest)
})
setwd("E:/PPMI帕金森数据库专用")
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
VIS <- c("SC","BL","V04","V06","V08","V10","V12")
vmap <- c(BL=0,SC=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# Loaders (cached reads)
rd <- function(f) suppressWarnings(read_csv(f, show_col_types=FALSE)) %>% filter(EVENT_ID %in% VIS)
p1c <- rd("运动症状数据/MDS-UPDRS_Part_I_31Jan2026.csv")
p1p <- rd("运动症状数据/MDS-UPDRS_Part_I_Patient_Questionnaire_31Jan2026.csv")
p2  <- rd("运动症状数据/MDS_UPDRS_Part_II__Patient_Questionnaire_31Jan2026.csv")
p3  <- rd("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv") %>% arrange(PATNO,EVENT_ID) %>% distinct(PATNO,EVENT_ID,.keep_all=TRUE)
p4  <- rd("运动症状数据/MDS-UPDRS_Part_IV__Motor_Complications_29Jan2026.csv")
moca<- rd("运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv")
rb  <- rd("REM_Sleep_Behavior_Disorder_Questionnaire_03Feb2026.csv")
sca <- rd("SCOPA-AUT_07Feb2026.csv")
dat <- rd("运动症状数据/Xing_Core_Lab_-_Quant_SBR_23Feb2026.csv")

SS <- function(d, cols) rowSums(sapply(intersect(cols,colnames(d)), function(c) as.numeric(d[[c]])), na.rm=TRUE)
mk <- function(d, score, nm) data.frame(PATNO=d$PATNO, EVENT_ID=d$EVENT_ID, Score=score, Outcome=nm)

L <- list()
# A. Totals
L[["UPDRS_I_clin"]] <- mk(p1c, p1c$NP1RTOT, "UPDRS_I_clin")
L[["UPDRS_I_pat"]]  <- mk(p1p, p1p$NP1PTOT, "UPDRS_I_pat")
L[["UPDRS_II"]]     <- mk(p2,  p2$NP2PTOT,  "UPDRS_II")
L[["UPDRS_III"]]    <- mk(p3,  p3$NP3TOT,   "UPDRS_III")
L[["UPDRS_IV"]]     <- mk(p4,  p4$NP4TOT,   "UPDRS_IV")
L[["MoCA"]]         <- mk(moca,moca$MCATOT, "MoCA")
L[["NHY"]]          <- mk(p3,  p3$NHY,      "NHY")
L[["PIGD"]]         <- mk(p3, SS(p3,c("NP3GAIT","NP3PSTBL","NP3FRZGT")), "PIGD")
L[["Tremor"]]       <- mk(p3, SS(p3,c("NP3PTRMR","NP3PTRML","NP3KTRMR","NP3KTRML","NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL")), "Tremor")
L[["RestTremor"]]   <- mk(p3, SS(p3,c("NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL")), "RestTremor")
it12 <- c("DRMVIVID","DRMAGRAC","DRMNOCTB","SLPLMBMV","SLPINJUR","DRMVERBL","DRMFIGHT","DRMUMV","DRMOBJFL","MVAWAKEN","DRMREMEM","SLPDSTRB")
neuro <- c("STROKE","HETRA","PARKISM","RLS","NARCLPSY","DEPRS","EPILEPSY","BRNINFM","CNSOTH")
L[["RBDSQ"]] <- mk(rb, SS(rb,it12) + as.integer(SS(rb,neuro)>0), "RBDSQ")
sca_na <- sca; for(c in paste0("SCAU",1:21)) if(c %in% colnames(sca_na)) sca_na[[c]][sca_na[[c]]==9] <- NA
L[["SCOPA_total"]] <- mk(sca_na, SS(sca_na, paste0("SCAU",1:21)), "SCOPA_total")

# B. DAT-SBR imaging (lower=worse; SC/V04/V06/V10, no V12)
L[["DAT_striatum"]] <- mk(dat, as.numeric(dat$STRIATUM_REF_CWM), "DAT_striatum")
L[["DAT_caudate"]]  <- mk(dat, as.numeric(dat$CAUDATE_REF_CWM),  "DAT_caudate")
L[["DAT_putamen"]]  <- mk(dat, as.numeric(dat$PUTAMEN_REF_CWM),  "DAT_putamen")

# C. MoCA subdomains
L[["MoCA_visuospatial"]] <- mk(moca, SS(moca,c("MCAALTTM","MCACUBE","MCACLCKC","MCACLCKN","MCACLCKH")), "MoCA_visuospatial")
L[["MoCA_naming"]]       <- mk(moca, SS(moca,c("MCALION","MCARHINO","MCACAMEL")), "MoCA_naming")
L[["MoCA_attention"]]    <- mk(moca, SS(moca,c("MCAFDS","MCABDS","MCAVIGIL","MCASER7")), "MoCA_attention")
L[["MoCA_language"]]     <- mk(moca, SS(moca,c("MCASNTNC","MCAVF")), "MoCA_language")
L[["MoCA_abstraction"]]  <- mk(moca, SS(moca,c("MCAABSTR")), "MoCA_abstraction")
L[["MoCA_recall"]]       <- mk(moca, SS(moca,c("MCAREC1","MCAREC2","MCAREC3","MCAREC4","MCAREC5")), "MoCA_recall")
L[["MoCA_orientation"]]  <- mk(moca, SS(moca,c("MCADATE","MCAMONTH","MCAYR","MCADAY","MCAPLACE","MCACITY")), "MoCA_orientation")

# D. UPDRS III motor subdomains
L[["UPDRS3_rigidity"]]    <- mk(p3, SS(p3,c("NP3RIGN","NP3RIGRU","NP3RIGLU","NP3RIGRL","NP3RIGLL")), "UPDRS3_rigidity")
L[["UPDRS3_bradykinesia"]]<- mk(p3, SS(p3,c("NP3FTAPR","NP3FTAPL","NP3HMOVR","NP3HMOVL","NP3PRSPR","NP3PRSPL","NP3TTAPR","NP3TTAPL","NP3LGAGR","NP3LGAGL")), "UPDRS3_bradykinesia")
L[["UPDRS3_axial"]]       <- mk(p3, SS(p3,c("NP3SPCH","NP3GAIT","NP3FRZGT","NP3PSTBL")), "UPDRS3_axial")

# E. UPDRS I non-motor individual items
for(it in c("NP1COG","NP1HALL","NP1DPRS","NP1ANXS","NP1APAT","NP1DDS")) L[[it]] <- mk(p1c, as.numeric(p1c[[it]]), it)
for(it in c("NP1SLPN","NP1SLPD","NP1PAIN","NP1URIN","NP1CNST","NP1LTHD","NP1FATG")) L[[it]] <- mk(p1p, as.numeric(p1p[[it]]), it)

# F. SCOPA subdomains
L[["SCOPA_GI"]]       <- mk(sca_na, SS(sca_na, paste0("SCAU",1:7)),   "SCOPA_GI")
L[["SCOPA_urinary"]]  <- mk(sca_na, SS(sca_na, paste0("SCAU",8:13)),  "SCOPA_urinary")
L[["SCOPA_cardio"]]   <- mk(sca_na, SS(sca_na, paste0("SCAU",14:16)), "SCOPA_cardio")
L[["SCOPA_thermo"]]   <- mk(sca_na, SS(sca_na, c("SCAU17","SCAU18","SCAU20","SCAU21")), "SCOPA_thermo")
L[["SCOPA_pupillo"]]  <- mk(sca_na, SS(sca_na, c("SCAU19")), "SCOPA_pupillo")

cat(sprintf("Total indicators to test: %d\n", length(L)))

# ---- Run LMM ----
res <- data.frame()
for(nm in names(L)) {
  d <- L[[nm]] %>% inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & !is.na(d$Score), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT) < 3 || nrow(d) < 50) { cat(sprintf("SKIP %s\n", nm)); next }
  lastv <- tail(levels(d$VISIT), 1)

  lastp<-NA; lastest<-NA; minp<-NA
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)) { cm<-coef(summary(mf)); ir<-grep(":groupLow$", rownames(cm))
    if(length(ir)>0){ minp<-min(cm[ir,"Pr(>|t|)"],na.rm=TRUE)
      lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
      if(length(lr)>0){lastp<-cm[lr[1],"Pr(>|t|)"];lastest<-cm[lr[1],"Estimate"]} } }
  grpp<-NA; intp<-NA; intest<-NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){ cm<-coef(summary(mc)); gr<-grep("^groupLow$",rownames(cm)); ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0) grpp<-cm[gr[1],"Pr(>|t|)"]
    if(length(ir)>0){intp<-cm[ir[1],"Pr(>|t|)"];intest<-cm[ir[1],"Estimate"]} }
  res <- rbind(res, data.frame(Outcome=nm, N_obs=nrow(d), LastVisit=lastv,
    LastV_int_est=lastest, LastV_int_p=lastp, MinVisit_int_p=minp,
    Cont_grp_p=grpp, Cont_int_est=intest, Cont_int_p=intp, stringsAsFactors=FALSE))
}
res$LastV_int_padj <- p.adjust(res$LastV_int_p, "BH")
res$Cont_int_padj  <- p.adjust(res$Cont_int_p, "BH")
res$Cont_grp_padj  <- p.adjust(res$Cont_grp_p, "BH")
res <- res[order(res$LastV_int_p),]

cat(sprintf("\n========== COMPREHENSIVE LMM SWEEP (%d indicators) ==========\n", nrow(res)))
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)
cat("\n=== Last-visit interaction nominal p<0.05 ===\n")
print(res[!is.na(res$LastV_int_p)&res$LastV_int_p<0.05, c("Outcome","LastVisit","LastV_int_est","LastV_int_p","LastV_int_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (last-visit padj<0.05): %s\n", paste(res$Outcome[!is.na(res$LastV_int_padj)&res$LastV_int_padj<0.05],collapse=", ")))
cat("\n=== Continuous interaction nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_int_p)&res$Cont_int_p<0.05, c("Outcome","Cont_int_est","Cont_int_p","Cont_int_padj")], row.names=FALSE)
cat("\n=== Continuous group main effect nominal p<0.05 ===\n")
print(res[!is.na(res$Cont_grp_p)&res$Cont_grp_p<0.05, c("Outcome","Cont_grp_p","Cont_grp_padj")], row.names=FALSE)
write.csv(res, file.path(OUT,"LMM_sweep_COMPREHENSIVE.csv"), row.names=FALSE)
cat(sprintf("\nSaved LMM_sweep_COMPREHENSIVE.csv (%d indicators)\n", nrow(res)))

# =====================================================================
# SECTION 15 | Forest plot of outcome sweep
# =====================================================================

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

# =====================================================================
# SECTION 16 | Full-sweep forest plot
# =====================================================================

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

# =====================================================================
# SECTION 17 | DAT-SBR subregion sweep
# =====================================================================

# DAT-SBR sub-region longitudinal LMM sweep x PC1 group
suppressPackageStartupMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(emmeans)})
setwd("E:/PPMI帕金森数据库专用")
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
PACKAGE_TABLES <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
vmap <- c(SC=0,BL=0,V04=1,V06=2,V08=3,V10=4,V12=5)

dat <- suppressWarnings(read_csv("运动症状数据/Xing_Core_Lab_-_Quant_SBR_23Feb2026.csv", show_col_types=FALSE)) %>%
  filter(EVENT_ID %in% c("SC","BL","V04","V06","V08","V10","V12"))

# all summary sub-regions (bilateral mean)
regions <- c("STRIATUM_REF_CWM","CAUDATE_REF_CWM","PUTAMEN_REF_CWM","PRECAUDATE_REF_CWM",
  "POSCAUDATE_REF_CWM","PRECOMMISSURAL_PUTAMEN_REF_CWM","POSCOMMISSURAL_PUTAMEN_REF_CWM",
  "PREDORSALPUTAMEN_REF_CWM","PREVENTRALPUTAMEN_REF_CWM","POSDORSALPUTAMEN_REF_CWM","POSVENTRALPUTAMEN_REF_CWM")
# asymmetry indices (|L-R|/mean) for caudate & putamen
mkasym <- function(L,R) ifelse((L+R)>0, abs(L-R)/((L+R)/2), NA)
dat$PUTAMEN_asym  <- mkasym(as.numeric(dat$PUTAMEN_L_REF_CWM),  as.numeric(dat$PUTAMEN_R_REF_CWM))
dat$CAUDATE_asym  <- mkasym(as.numeric(dat$CAUDATE_L_REF_CWM),  as.numeric(dat$CAUDATE_R_REF_CWM))
dat$STRIATUM_asym <- mkasym(as.numeric(dat$STRIATUM_L_REF_CWM), as.numeric(dat$STRIATUM_R_REF_CWM))
# lowest-side putamen (most affected hemisphere) — clinically key
dat$PUTAMEN_min <- pmin(as.numeric(dat$PUTAMEN_L_REF_CWM), as.numeric(dat$PUTAMEN_R_REF_CWM), na.rm=TRUE)

outcomes <- c(regions, "PUTAMEN_asym","CAUDATE_asym","STRIATUM_asym","PUTAMEN_min")

res <- data.frame()
for(rg in outcomes) {
  d <- dat %>% mutate(Score=as.numeric(.data[[rg]])) %>%
    dplyr::select(PATNO, EVENT_ID, Score) %>% filter(!is.na(Score)) %>%
    inner_join(pc1, by="PATNO") %>%
    mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
           Time=vmap[EVENT_ID],
           VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID), levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT)<3 || nrow(d)<80) { cat(sprintf("SKIP %s\n", rg)); next }
  lastv <- tail(levels(d$VISIT),1)
  grpp<-NA;grpe<-NA;intp<-NA;inte<-NA;lastp<-NA
  mc <- tryCatch(lmer(Score ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(Score ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){cm<-coef(summary(mc));gr<-grep("^groupLow$",rownames(cm));ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0){grpe<-cm[gr[1],"Estimate"];grpp<-cm[gr[1],"Pr(>|t|)"]}
    if(length(ir)>0){inte<-cm[ir[1],"Estimate"];intp<-cm[ir[1],"Pr(>|t|)"]}}
  mf <- tryCatch(lmer(Score ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)){cm<-coef(summary(mf));lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
    if(length(lr)>0) lastp<-cm[lr[1],"Pr(>|t|)"]}
  res <- rbind(res, data.frame(Region=rg, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)), LastVisit=lastv,
    Grp_est=grpe, Grp_p=grpp, Int_est=inte, Int_p=intp, LastV_int_p=lastp, stringsAsFactors=FALSE))
}
res$Grp_padj <- p.adjust(res$Grp_p,"BH"); res$Int_padj <- p.adjust(res$Int_p,"BH")
res <- res[order(res$Grp_p),]
cat(sprintf("\n===== DAT-SBR SUB-REGION SWEEP (%d regions, follow-up to %s) =====\n", nrow(res), paste(unique(res$LastVisit),collapse="/")))
print(res %>% mutate(across(where(is.numeric), ~round(.,4))), row.names=FALSE)
cat("\n=== Group main effect nominal p<0.05 ===\n")
print(res[!is.na(res$Grp_p)&res$Grp_p<0.05,c("Region","Grp_est","Grp_p","Grp_padj")], row.names=FALSE)
cat(sprintf("SURVIVE FDR (group): %s\n", paste(res$Region[!is.na(res$Grp_padj)&res$Grp_padj<0.05],collapse=", ")))
cat("\n=== Interaction nominal p<0.05 ===\n")
print(res[!is.na(res$Int_p)&res$Int_p<0.05,c("Region","Int_est","Int_p","Int_padj")], row.names=FALSE)
cat(sprintf("SURVIVE FDR (interaction): %s\n", paste(res$Region[!is.na(res$Int_padj)&res$Int_padj<0.05],collapse=", ")))
write.csv(res, file.path(OUT,"DAT_subregion_sweep.csv"), row.names=FALSE)
write.csv(res, file.path(PACKAGE_TABLES,"DAT_subregion_sweep.csv"), row.names=FALSE)
cat("\nSaved DAT_subregion_sweep.csv\n")

# Model-estimated trajectory for the clinically representative posterior dorsal putamen.
# The full 15-region screen remains in DAT_subregion_sweep.csv.
representative_region <- "POSDORSALPUTAMEN_REF_CWM"
d_traj <- dat %>%
  mutate(Score=as.numeric(.data[[representative_region]])) %>%
  dplyr::select(PATNO, EVENT_ID, Score) %>%
  filter(!is.na(Score)) %>%
  inner_join(pc1, by="PATNO") %>%
  mutate(
    PATNO=factor(PATNO),
    group=factor(group, levels=c("High","Low")),
    VISIT=factor(ifelse(EVENT_ID=="SC","BL",EVENT_ID),
      levels=c("BL","V04","V06","V08","V10","V12"))
  ) %>%
  filter(!is.na(VISIT)) %>%
  droplevels()

m_traj <- lmer(Score ~ VISIT*group + (1|PATNO), data=d_traj, REML=FALSE)
emm_traj <- as.data.frame(emmeans(m_traj, ~ group | VISIT))
n_traj <- d_traj %>%
  count(VISIT, group, name="N")
screen_row <- res %>%
  filter(Region == representative_region)

emm_traj <- emm_traj %>%
  left_join(n_traj, by=c("VISIT","group")) %>%
  mutate(
    Region=representative_region,
    Baseline_group_p=screen_row$Grp_p[[1]],
    Linear_group_by_year_p=screen_row$Int_p[[1]],
    Last_visit_interaction_p=screen_row$LastV_int_p[[1]],
    N_subjects=nlevels(d_traj$PATNO)
  )

write.csv(emm_traj,
  file.path(PACKAGE_TABLES,"DAT_posterior_dorsal_putamen_trajectory.csv"),
  row.names=FALSE)
cat("Saved DAT_posterior_dorsal_putamen_trajectory.csv\n")

# =====================================================================
# SECTION 18 | Fluid biomarker sweep
# =====================================================================

# Biomarker LMM sweep: blood/CSF biomarkers x PC1 group (main effect + interaction)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lme4); library(lmerTest)})
OUT <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1"
pc1 <- read.csv(file.path(OUT,"data/pc1_scores.csv"))[,c("PATNO","group")]
bm <- read.csv(file.path(OUT,"data/biomarkers_long.csv"))
vmap <- c(BL=0,V04=1,V06=2,V08=3,V10=4,V12=5)

# collapse duplicate measurements (same PATNO/visit/biomarker) by mean
bm <- bm %>% group_by(PATNO, EVENT_ID, biomarker) %>%
  summarise(value=mean(value, na.rm=TRUE), .groups="drop")

# keep biomarkers with adequate longitudinal coverage
qc <- bm %>% group_by(biomarker) %>%
  summarise(n_obs=n(), n_subj=n_distinct(PATNO),
            n_visits=n_distinct(EVENT_ID), .groups="drop") %>%
  filter(n_visits>=3, n_obs>=80, n_subj>=40)
cat(sprintf("Biomarkers passing QC: %d\n", nrow(qc)))

res <- data.frame()
for(b in qc$biomarker) {
  d <- bm %>% filter(biomarker==b) %>% inner_join(pc1, by="PATNO")
  if(nrow(d)<80 || n_distinct(d$group)<2) next
  # log transform if all positive & skewed
  if(all(d$value>0, na.rm=TRUE)) d$value <- log(d$value)
  d <- d %>% mutate(PATNO=factor(PATNO), group=factor(group, levels=c("High","Low")),
                    Time=vmap[EVENT_ID],
                    VISIT=factor(EVENT_ID, levels=c("BL","V04","V06","V08","V10","V12")))
  d <- d[!is.na(d$Time) & is.finite(d$value), ]; d$VISIT <- droplevels(d$VISIT)
  if(nlevels(d$VISIT)<3 || nrow(d)<80) next
  lastv <- tail(levels(d$VISIT),1)

  # continuous-time: group main + time:group
  grpp<-NA; grpe<-NA; intp<-NA; inte<-NA
  mc <- tryCatch(lmer(value ~ Time*group + (1+Time|PATNO), data=d, REML=FALSE), error=function(e)
        tryCatch(lmer(value ~ Time*group + (1|PATNO), data=d, REML=FALSE), error=function(e2) NULL))
  if(!is.null(mc)){ cm<-coef(summary(mc)); gr<-grep("^groupLow$",rownames(cm)); ir<-grep("Time:groupLow",rownames(cm),fixed=TRUE)
    if(length(gr)>0){grpp<-cm[gr[1],"Pr(>|t|)"]; grpe<-cm[gr[1],"Estimate"]}
    if(length(ir)>0){intp<-cm[ir[1],"Pr(>|t|)"]; inte<-cm[ir[1],"Estimate"]} }
  # factor-time last visit interaction
  lastp<-NA
  mf <- tryCatch(lmer(value ~ VISIT*group + (1|PATNO), data=d, REML=FALSE), error=function(e) NULL)
  if(!is.null(mf)){ cm<-coef(summary(mf)); lr<-grep(paste0("VISIT",lastv,":groupLow"),rownames(cm),fixed=TRUE)
    if(length(lr)>0) lastp<-cm[lr[1],"Pr(>|t|)"] }
  res <- rbind(res, data.frame(Biomarker=b, N_obs=nrow(d), N_subj=nlevels(droplevels(d$PATNO)),
    LastVisit=lastv, Grp_est=grpe, Grp_p=grpp, Int_est=inte, Int_p=intp, LastV_int_p=lastp, stringsAsFactors=FALSE))
}
res$Grp_padj <- p.adjust(res$Grp_p,"BH")
res$Int_padj <- p.adjust(res$Int_p,"BH")
res <- res[order(res$Grp_p),]

cat(sprintf("\n========== BIOMARKER LMM SWEEP (%d biomarkers) ==========\n", nrow(res)))
cat("\n=== GROUP MAIN EFFECT: nominal p<0.05 (low vs high CI overall) ===\n")
print(res[!is.na(res$Grp_p)&res$Grp_p<0.05, c("Biomarker","N_subj","Grp_est","Grp_p","Grp_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (group main padj<0.05): %s\n", paste(res$Biomarker[!is.na(res$Grp_padj)&res$Grp_padj<0.05],collapse=", ")))
cat("\n=== TIME x GROUP INTERACTION: nominal p<0.05 ===\n")
print(res[!is.na(res$Int_p)&res$Int_p<0.05, c("Biomarker","Int_est","Int_p","Int_padj")], row.names=FALSE)
cat(sprintf("\nSURVIVE FDR (interaction padj<0.05): %s\n", paste(res$Biomarker[!is.na(res$Int_padj)&res$Int_padj<0.05],collapse=", ")))
write.csv(res, file.path(OUT,"biomarker_LMM_sweep.csv"), row.names=FALSE)
cat(sprintf("\nSaved biomarker_LMM_sweep.csv (%d biomarkers)\n", nrow(res)))

# =====================================================================
# SECTION 19 | Biomarker volcano figure
# =====================================================================

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

# =====================================================================
# SECTION 20 | Cell-composition sensitivity analysis
# =====================================================================

# Part 2 (standalone, no dds): PIGD year-5 LMM cell-fraction sensitivity
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(dplyr); library(readr); library(tidyr) })
setwd("E:/PPMI帕金森数据库专用")
FD <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"
TAB <- "E:/PPMI帕金森数据库专用/结果/新结果/实验/movementdisorders/new/revision_v2/FINAL_PACKAGE_2026-06-21/04_tables"

pc1s <- read.csv(file.path(FD,"pc1_scores.csv"))
imf  <- read.csv(file.path(FD,"immune_fractions.csv"), check.names=FALSE)
cf <- imf[, c("SAMPLE_ID","Neutrophils","T cells CD8","T cells CD4 memory resting",
              "T cells CD4 naive","B cells naive","NK cells resting","Monocytes")]
colnames(cf) <- c("SAMPLE_ID","Neu","CD8","CD4mr","CD4n","Bn","NKr","Mono")
cf2 <- merge(cf, pc1s[,c("SAMPLE_ID","PATNO")], by="SAMPLE_ID")

m3 <- read_csv("运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv", show_col_types=FALSE) %>%
  filter(EVENT_ID %in% c("BL","V04","V06","V08","V10","V12")) %>%
  mutate(PIGD = NP3GAIT+NP3PSTBL+NP3FRZGT) %>%
  dplyr::select(PATNO, EVENT_ID, PIGD) %>% filter(!is.na(PIGD)) %>%
  group_by(PATNO, EVENT_ID) %>% summarise(PIGD=mean(PIGD), .groups="drop")   # dedupe

ld <- m3 %>% inner_join(pc1s[,c("PATNO","PC1","group")], by="PATNO") %>%
  left_join(cf2[,c("PATNO","Neu","CD8","CD4mr","CD4n","Bn","NKr","Mono")], by="PATNO")
ld$visit <- factor(ld$EVENT_ID, levels=c("BL","V04","V06","V08","V10","V12"))
ld$group <- factor(ld$group, levels=c("Low","High"))
ld$PC1z  <- scale(ld$PC1)[,1]
ld <- ld[complete.cases(ld[,c("Neu","CD8","CD4mr")]),]
cat(sprintf("LMM long data: %d rows, %d patients\n", nrow(ld), length(unique(ld$PATNO))))

getV12 <- function(m){ s<-summary(m)$coefficients
  r <- rownames(s)[grepl("V12", rownames(s)) & grepl(":", rownames(s))]
  s[r, c("Estimate","Pr(>|t|)"), drop=FALSE] }

getV12_full <- function(m, model, adjustment, direction) {
  s <- summary(m)$coefficients
  r <- rownames(s)[grepl("V12", rownames(s)) & grepl(":", rownames(s))]
  out <- as.data.frame(s[r, , drop=FALSE], check.names=FALSE)
  out$Term <- rownames(out)
  out$Model <- model
  out$Adjustment <- adjustment
  out$Outcome <- "PIGD"
  out$Visit <- "V12"
  out$Estimate_lowCI_direction <- direction * out$Estimate
  out$N_obs <- nobs(m)
  out$N_patients <- length(unique(ld$PATNO))
  out[, c("Outcome","Visit","Model","Adjustment","Term","Estimate",
          "Estimate_lowCI_direction","Std. Error","df","t value","Pr(>|t|)",
          "N_obs","N_patients")]
}

m_bin_raw <- lmer(PIGD~group*visit+(1|PATNO),ld)
m_bin_adj <- lmer(PIGD~group*visit+Neu+CD8+CD4mr+CD4n+Bn+NKr+Mono+(1|PATNO),ld)
m_pc1_raw <- lmer(PIGD~PC1z*visit+(1|PATNO),ld)
m_pc1_adj <- lmer(PIGD~PC1z*visit+Neu+CD8+CD4mr+CD4n+Bn+NKr+Mono+(1|PATNO),ld)

cat("\n(a) binary group×visit UNADJUSTED:\n");      print(round(getV12(m_bin_raw),4))
cat("\n(b) binary group×visit CELL-ADJUSTED:\n");   print(round(getV12(m_bin_adj),4))
cat("\n(c) continuous PC1z×visit UNADJUSTED:\n");    print(round(getV12(m_pc1_raw),4))
cat("\n(d) continuous PC1z×visit CELL-ADJUSTED:\n"); print(round(getV12(m_pc1_adj),4))

cell_sens <- bind_rows(
  getV12_full(m_bin_raw, "Binary PC1 median group x visit", "Unadjusted", -1),
  getV12_full(m_bin_adj, "Binary PC1 median group x visit", "Cell-fraction adjusted", -1),
  getV12_full(m_pc1_raw, "Continuous PC1_z x visit", "Unadjusted", -1),
  getV12_full(m_pc1_adj, "Continuous PC1_z x visit", "Cell-fraction adjusted", -1)
) %>%
  mutate(
    Effect_interpretation = ifelse(
      Model == "Continuous PC1_z x visit",
      "Additional year-5 delta PIGD per 1-SD lower PC1",
      "Additional year-5 delta PIGD in low-CI vs high-CI"
    ),
    Cell_covariates = ifelse(
      Adjustment == "Cell-fraction adjusted",
      "Neutrophils, CD8 T cells, CD4 memory resting T cells, CD4 naive T cells, naive B cells, resting NK cells, monocytes",
      "None"
    )
  )
write.csv(cell_sens, file.path(TAB, "PIGD_cell_composition_sensitivity.csv"), row.names=FALSE)
cat(sprintf("Saved %s\n", file.path(TAB, "PIGD_cell_composition_sensitivity.csv")))

# =====================================================================
# SECTION 21 | PIGD latent-class trajectory (LCMM, supplementary)
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
})

if (!requireNamespace("lcmm", quietly = TRUE)) {
  stop("Package 'lcmm' is required. Install it with install.packages('lcmm') and rerun this script.")
}

library(lcmm)

PKG <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
ROOT <- "E:/PPMI帕金森数据库专用"
TAB <- file.path(PKG, "04_tables")
FIG <- file.path(PKG, "02_figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

set.seed(20260815)

event_map <- tibble(
  EVENT_ID = c("BL", "SC", "V04", "V06", "V08", "V10", "V12"),
  visit = c("BL", "BL", "V04", "V06", "V08", "V10", "V12"),
  time = c(0, 0, 1, 2, 3, 4, 5)
)

pc1 <- read_csv(file.path(TAB, "three_scores.csv"), show_col_types = FALSE) %>%
  transmute(SAMPLE_ID, PC1)
pat_map <- read_csv(file.path(TAB, "Table_clustering_methods.csv"), show_col_types = FALSE) %>%
  filter(is_PD) %>%
  distinct(SAMPLE_ID, PATNO)

ci <- pc1 %>%
  inner_join(pat_map, by = "SAMPLE_ID") %>%
  mutate(CI_group = if_else(PC1 < median(PC1, na.rm = TRUE), "Low", "High"))

motor <- read_csv(
  file.path(ROOT, "运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv"),
  show_col_types = FALSE
) %>%
  filter(EVENT_ID %in% event_map$EVENT_ID) %>%
  mutate(
    PIGD = as.numeric(NP3GAIT) + as.numeric(NP3PSTBL) + as.numeric(NP3FRZGT)
  ) %>%
  select(PATNO, EVENT_ID, PIGD) %>%
  filter(!is.na(PIGD)) %>%
  left_join(event_map, by = "EVENT_ID") %>%
  group_by(PATNO, visit, time) %>%
  summarise(PIGD = mean(PIGD, na.rm = TRUE), .groups = "drop") %>%
  inner_join(ci, by = "PATNO") %>%
  mutate(PATNO = as.integer(PATNO), time = as.numeric(time), PIGD = as.numeric(PIGD))

write_csv(motor, file.path(TAB, "PIGD_lcmm_long_input.csv"))

fit_lmm <- lmer(PIGD ~ time + (time | PATNO), data = motor, REML = TRUE)
avg_slope <- unname(fixef(fit_lmm)["time"])
progressive_threshold <- 2 * avg_slope
write_csv(
  tibble(
    model = "PIGD ~ time + (time | PATNO)",
    mean_slope_per_year = avg_slope,
    progressive_threshold_rule = "2 x cohort mean slope",
    progressive_threshold_per_year = progressive_threshold,
    n_patients = n_distinct(motor$PATNO),
    n_observations = nrow(motor)
  ),
  file.path(TAB, "PIGD_lcmm_progressive_threshold.csv")
)

fits <- vector("list", 4)
fits[[1]] <- hlme(
  fixed = PIGD ~ time,
  random = ~ time,
  subject = "PATNO",
  ng = 1,
  data = motor,
  verbose = FALSE
)
seed1 <- fits[[1]]
for (k in 2:4) {
  message("Fitting ", k, "-class LCLMM...")
  fits[[k]] <- tryCatch(
    gridsearch(
      hlme(
        fixed = PIGD ~ time,
        mixture = ~ time,
        random = ~ time,
        subject = "PATNO",
        ng = k,
        data = motor,
        B = seed1,
        verbose = FALSE
      ),
      rep = 20,
      maxiter = 30,
      minit = seed1
    ),
    error = function(e) {
      message("  failed: ", conditionMessage(e))
      NULL
    }
  )
}

mean_max_prob <- function(fit, k) {
  if (is.null(fit)) return(NA_real_)
  if (k == 1) return(1)
  prob_cols <- paste0("prob", seq_len(k))
  mean(apply(fit$pprob[, prob_cols, drop = FALSE], 1, max), na.rm = TRUE)
}

model_selection <- bind_rows(lapply(seq_along(fits), function(k) {
  fit <- fits[[k]]
  if (is.null(fit)) {
    return(tibble(n_classes = k, conv = NA_integer_, loglik = NA_real_,
                  AIC = NA_real_, BIC = NA_real_, mean_max_postprob = NA_real_))
  }
  tibble(
    n_classes = k,
    conv = fit$conv,
    loglik = fit$loglik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    mean_max_postprob = mean_max_prob(fit, k)
  )
}))

write_csv(model_selection, file.path(TAB, "PIGD_lcmm_model_selection.csv"))

eligible <- model_selection %>%
  filter(!is.na(BIC), conv == 1, mean_max_postprob >= 0.80)
if (nrow(eligible) == 0) {
  eligible <- model_selection %>% filter(!is.na(BIC), conv == 1)
}
best_k <- eligible$n_classes[which.min(eligible$BIC)]
best_fit <- fits[[best_k]]

pred_grid <- data.frame(time = seq(0, 5, by = 0.1))
pred <- predictY(best_fit, newdata = pred_grid, var.time = "time", draws = FALSE)
pred_mat <- as.data.frame(pred$pred)
names(pred_mat) <- paste0("class", seq_len(ncol(pred_mat)))
pred_plot <- bind_cols(pred_grid, pred_mat) %>%
  pivot_longer(starts_with("class"), names_to = "class_label", values_to = "PIGD_pred") %>%
  mutate(class = as.integer(sub("class", "", class_label)))

class_summary <- pred_plot %>%
  filter(time %in% c(0, 5)) %>%
  select(class, time, PIGD_pred) %>%
  pivot_wider(names_from = time, values_from = PIGD_pred, names_prefix = "time_") %>%
  mutate(
    slope_per_year = (`time_5` - `time_0`) / 5,
    progressive_binary = if_else(slope_per_year >= progressive_threshold, "Progressive", "Stable"),
    trajectory_rank = rank(slope_per_year, ties.method = "first"),
    trajectory3 = case_when(
      trajectory_rank == min(trajectory_rank) ~ "Stable",
      trajectory_rank == max(trajectory_rank) ~ "Fast progressive",
      TRUE ~ "Slow progressive"
    ),
    trajectory3_cn = recode(
      trajectory3,
      "Stable" = "稳定型",
      "Slow progressive" = "慢性进展型",
      "Fast progressive" = "快速进展型"
    )
  ) %>%
  select(
    class,
    intercept = `time_0`,
    PIGD_year5 = `time_5`,
    slope_per_year,
    progressive_binary,
    trajectory3,
    trajectory3_cn
  )

pp <- best_fit$pprob
id_col <- intersect(c("PATNO", "subject"), names(pp))[1]
if (is.na(id_col)) id_col <- names(pp)[1]

labels <- pp %>%
  as_tibble() %>%
  rename(PATNO = all_of(id_col)) %>%
  mutate(PATNO = as.integer(PATNO)) %>%
  left_join(class_summary %>% select(class, progressive_binary, trajectory3, trajectory3_cn), by = "class") %>%
  left_join(ci %>% select(PATNO, PC1, CI_group), by = "PATNO")

class_counts <- labels %>%
  count(class, progressive_binary, trajectory3, trajectory3_cn, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  left_join(class_summary, by = c("class", "progressive_binary", "trajectory3", "trajectory3_cn")) %>%
  arrange(slope_per_year)

enrichment <- table(labels$CI_group, labels$progressive_binary)
enrichment_df <- as.data.frame(enrichment) %>%
  as_tibble() %>%
  rename(CI_group = Var1, progressive_binary = Var2, n = Freq) %>%
  group_by(CI_group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()
if (nrow(enrichment) >= 2 && ncol(enrichment) >= 2) {
  enrichment_test <- fisher.test(enrichment)
  enrichment_stats <- tibble(
    test = "Fisher exact",
    p_value = enrichment_test$p.value,
    odds_ratio = unname(enrichment_test$estimate)
  )
} else {
  enrichment_stats <- tibble(
    test = "Fisher exact",
    p_value = NA_real_,
    odds_ratio = NA_real_
  )
}

labels <- labels %>%
  mutate(
    trajectory3 = factor(trajectory3, levels = c("Stable", "Slow progressive", "Fast progressive")),
    trajectory3_cn = factor(trajectory3_cn, levels = c("稳定型", "慢性进展型", "快速进展型"))
  )
enrichment3 <- table(labels$CI_group, labels$trajectory3)
enrichment3_df <- as.data.frame(enrichment3) %>%
  as_tibble() %>%
  rename(CI_group = Var1, trajectory3 = Var2, n = Freq) %>%
  group_by(CI_group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    CI_group_label = recode(CI_group, "Low" = "低CI表达组", "High" = "高CI表达组"),
    trajectory3_cn = recode(
      as.character(trajectory3),
      "Stable" = "稳定型",
      "Slow progressive" = "慢性进展型",
      "Fast progressive" = "快速进展型"
    )
  )
if (nrow(enrichment3) >= 2 && ncol(enrichment3) >= 2) {
  enrichment3_test <- fisher.test(enrichment3)
  enrichment3_stats <- tibble(
    test = "Fisher exact",
    p_value = enrichment3_test$p.value,
    table = "CI_group x 3-class LCLMM trajectory"
  )
} else {
  enrichment3_stats <- tibble(
    test = "Fisher exact",
    p_value = NA_real_,
    table = "CI_group x 3-class LCLMM trajectory"
  )
}

write_csv(class_counts, file.path(TAB, "PIGD_lcmm_class_trajectories.csv"))
write_csv(labels, file.path(TAB, "PIGD_lcmm_patient_labels.csv"))
write_csv(enrichment_df, file.path(TAB, "PIGD_lcmm_CI_enrichment.csv"))
write_csv(enrichment_stats, file.path(TAB, "PIGD_lcmm_CI_enrichment_test.csv"))
write_csv(enrichment3_df, file.path(TAB, "PIGD_lcmm_CI_trajectory3_enrichment.csv"))
write_csv(enrichment3_stats, file.path(TAB, "PIGD_lcmm_CI_trajectory3_test.csv"))

pred_plot <- pred_plot %>%
  left_join(class_summary %>% select(class, trajectory3), by = "class") %>%
  mutate(trajectory3 = factor(trajectory3, levels = c("Stable", "Slow progressive", "Fast progressive")))

p <- ggplot() +
  geom_point(
    data = motor,
    aes(time, PIGD),
    color = "grey70",
    alpha = 0.16,
    size = 0.6,
    position = position_jitter(width = 0.03, height = 0.04)
  ) +
  geom_line(
    data = pred_plot,
    aes(time, PIGD_pred, color = trajectory3, group = class),
    linewidth = 1.1
  ) +
  geom_text(
    data = class_counts,
    aes(x = 5.05, y = PIGD_year5,
        label = sprintf("Class %d: n=%d, slope=%.2f/y", class, n, slope_per_year),
        color = trajectory3),
    hjust = 0,
    size = 3.0
  ) +
  scale_color_manual(values = c(
    "Stable" = "#2F6B9A",
    "Slow progressive" = "#D99A2B",
    "Fast progressive" = "#C74343"
  )) +
  coord_cartesian(xlim = c(0, 6.15), ylim = c(0, 12), clip = "off") +
  labs(
    x = "Time (years)",
    y = "PIGD score",
    color = NULL,
    title = "Latent class linear mixed model of 5-year PIGD trajectories",
    subtitle = sprintf(
      "Selected model: %d classes; Venuto-style progressive threshold = %.2f points/year",
      best_k, progressive_threshold
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.margin = margin(8, 88, 8, 8),
    legend.position = "top"
  )

ggsave(file.path(FIG, "FigS_PIGD_lcmm_trajectories.png"), p, width = 7.2, height = 4.2, dpi = 320)
ggsave(file.path(FIG, "FigS_PIGD_lcmm_trajectories.pdf"), p, width = 7.2, height = 4.2)

message("Saved LCLMM outputs:")
message("  ", file.path(TAB, "PIGD_lcmm_progressive_threshold.csv"))
message("  ", file.path(TAB, "PIGD_lcmm_model_selection.csv"))
message("  ", file.path(TAB, "PIGD_lcmm_class_trajectories.csv"))
message("  ", file.path(TAB, "PIGD_lcmm_patient_labels.csv"))
message("  ", file.path(FIG, "FigS_PIGD_lcmm_trajectories.png"))

# =====================================================================
# SECTION 22 | Figure 1: CI axis construction and structure validation
# =====================================================================

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

# =====================================================================
# SECTION 23 | Figure 2: DE, GSEA, GO, risk-gene cross-cohort validation, immune shifts
# (revision_v2 final (3 robust genes in panel d))
# =====================================================================

# ============================================================
# Figure 2: DE, KEGG-GSEA, GO enrichment, cross-cohort, and immune panels
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(ggtext)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(sub("^--file=", "", file_arg[1]))
} else {
  getwd()
}

ROOT <- file.path(SCRIPT_DIR, "..")
FIG_DIR <- file.path(ROOT, "02_figures")
TAB_DIR <- file.path(ROOT, "04_tables")
FALLBACK_DATA_DIR <- file.path(ROOT, "github_code_release_input_files", "derived_tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

read_first <- function(paths, ...) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop("None of these input files exists: ", paste(paths, collapse = " | "))
  }
  read.csv(hit, ...)
}

LOW <- "#D55E00"
HIGH <- "#0072B2"
PPMI_COL <- "#0072B2"
GSE_COL <- "#E69F00"
SIG_COL <- "#D55E00"
NS_COL <- "grey62"

theme_pub <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.tag = element_text(face = "bold", size = 11),
    plot.margin = margin(4, 5, 4, 5)
  )

# ---- 2A: volcano (no gene labels; risk-gene story is in panels d/e) ----
de <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "de_full.csv"),
  file.path(TAB_DIR, "de_full.csv")
))
de <- de[!is.na(de$padj), ]
de$Symbol <- ifelse(is.na(de$Symbol), "", as.character(de$Symbol))
de$padj_plot <- pmax(de$padj, .Machine$double.xmin)

pd_genes <- c(
  "ASXL3", "BAG3", "BIN3", "BRIP1", "BST1", "C5orf24", "CAB39L", "CAMK2D",
  "CASC16", "CD19", "CHD9", "CHRNB1", "CLCN3", "CRHR1", "CRLS1", "CTSB",
  "CYRIB", "DLG2", "DNAH17", "DYRK1A", "ELOVL7", "FAM171A2", "FAM47E",
  "FBRSL1", "FCGR2A", "FGF20", "FYN", "GAK", "GALC", "GBA1", "GBF1", "GCH1",
  "GPNMB", "GS1-124K5.11", "HIP1R", "HLA-DRB5", "HLA-H", "IGSF9B", "INPP5F",
  "IP6K2", "ITGA8", "ITPKB", "KCNIP3", "KCNS3", "KPNA1", "KRTCAP2", "LCORL",
  "LINC00693", "LRRK2", "MAP4K4", "MBNL2", "MCCC1", "MED12L", "MEX3C",
  "MIPOL1", "NOD2", "NUCKS1", "PAM", "PMVK", "RAB29", "RETREG3", "RIMS1",
  "RIT2", "RNF141", "RPS12", "RPS6KL1", "SATB1", "SCAF11", "SCARB2",
  "SETD1A", "SH3GL2", "SIPA1L2", "SNCA", "SPPL2B", "SPTSSB", "STK39", "SYT17",
  "TMEM163", "TMEM175", "TRIM40", "UBAP2", "UBTF", "VAMP4", "VPS13C", "WNT3"
)

de$pass <- de$padj < 0.05 & abs(de$log2FoldChange) > 1
de$direction <- ifelse(de$log2FoldChange < 0, "LowUp", "HighUp")
de$group <- ifelse(!de$pass, "nonsig", de$direction)
de$group <- factor(de$group, levels = c("LowUp", "HighUp", "nonsig"))
de <- de[order(de$pass), ]

n_low_up <- sum(de$pass & de$log2FoldChange < 0, na.rm = TRUE)
n_high_up <- sum(de$pass & de$log2FoldChange > 0, na.rm = TRUE)

# gradient background: white at the centre, fading toward LOW/HIGH at the tails
mix_col <- function(col_hex, frac) {
  rgb_mix <- col2rgb(col_hex) / 255
  out <- c(1, 1, 1) * (1 - frac) + rgb_mix * frac
  rgb(out[1], out[2], out[3])
}
n_steps <- 40
max_frac <- 0.55
xs_left <- seq(-3.4, 0, length.out = n_steps)
xs_right <- seq(0, 2.4, length.out = n_steps)
grad_df <- rbind(
  data.frame(xmin = xs_left[-(n_steps)], xmax = xs_left[-1], ymin = 0, ymax = 130,
    fill = vapply(rev(seq_len(n_steps - 1)),
      function(i) mix_col(LOW, max_frac * i / n_steps), character(1))),
  data.frame(xmin = xs_right[-(n_steps)], xmax = xs_right[-1], ymin = 0, ymax = 130,
    fill = vapply(seq_len(n_steps - 1),
      function(i) mix_col(HIGH, max_frac * i / n_steps), character(1)))
)

f2a <- ggplot(de, aes(log2FoldChange, -log10(padj_plot))) +
  geom_rect(data = grad_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    inherit.aes = FALSE, alpha = 0.35) +
  scale_fill_identity(guide = "none") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.25) +
  geom_point(aes(color = group), alpha = 0.68, stroke = 0, size = 0.5) +
  scale_color_manual(values = c(
    LowUp = LOW, HighUp = HIGH, nonsig = "grey83"
  ), guide = "none") +
  coord_cartesian(xlim = c(-3.4, 2.4), ylim = c(0, 130), clip = "off") +
  annotate("text", x = -2.2, y = 128,
    label = "Low-CI expression",
    color = LOW, fontface = "bold", size = 2.8, hjust = 0.5) +
  annotate("text", x = 1.7, y = 128,
    label = "High-CI expression",
    color = HIGH, fontface = "bold", size = 2.8, hjust = 0.5) +
  labs(tag = "a", x = expression(log[2] ~ "fold change"),
    y = expression(-log[10] ~ "FDR")) +
  theme_pub +
  theme(aspect.ratio = 1)

# ---- 2B: KEGG + GO BP combined enrichment bar plot ----
kegg_gsea <- read_first(c(file.path(TAB_DIR, "GSEA_KEGG_lowCI_ranked.csv")))
kegg_keep <- c(
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_AXON_GUIDANCE",
  "KEGG_FOCAL_ADHESION",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_TIGHT_JUNCTION",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY"
)
kegg_panel <- kegg_gsea[kegg_gsea$ID %in% kegg_keep & kegg_gsea$NES > 0, ]
kegg_panel$Term <- gsub("^KEGG_", "", kegg_panel$ID)
kegg_panel$Term <- gsub("_", " ", kegg_panel$Term)
kegg_panel$Term <- tools::toTitleCase(tolower(kegg_panel$Term))
kegg_panel$Term <- gsub("Ecm", "ECM", kegg_panel$Term)
kegg_panel$Term <- gsub("Dna", "DNA", kegg_panel$Term)
kegg_panel$Term <- gsub("Parkinsons Disease", "Parkinson's Disease", kegg_panel$Term)
kegg_panel$mlog10_fdr <- -log10(pmax(kegg_panel$p.adjust, .Machine$double.xmin))
kegg_panel$Count <- vapply(strsplit(kegg_panel$core_enrichment, "/"), length, integer(1))
kegg_panel$genes <- vapply(strsplit(kegg_panel$core_enrichment, "/"),
  function(x) paste(head(x, 12), collapse = "/"), character(1))
kegg_panel$genes <- ifelse(nchar(kegg_panel$genes) > 78,
  paste0(substr(kegg_panel$genes, 1, 75), "..."), kegg_panel$genes)
kegg_panel$Category <- "KEGG"

kegg_panel$Term <- factor(kegg_panel$Term, levels = kegg_panel$Term[order(kegg_panel$mlog10_fdr)])
focus_terms <- c("Neuroactive Ligand Receptor Interaction", "Axon Guidance",
  "ECM Receptor Interaction", "Complement And Coagulation Cascades")
kegg_panel$focus <- kegg_panel$Term %in% focus_terms

x_text <- 0.5
f2b_gsea <- ggplot(kegg_panel, aes(x = mlog10_fdr, y = Term)) +
  geom_col(fill = LOW, width = 0.62, alpha = 0.25, color = NA) +
  geom_text(aes(label = Term, color = focus), x = x_text, hjust = 0,
    fontface = "bold", size = 2.05, vjust = -0.10) +
  geom_text(aes(label = genes), x = x_text, hjust = 0,
    fontface = "italic", size = 1.18, color = "grey45", vjust = 1.45) +
  scale_color_manual(values = c("TRUE" = "#CB181D", "FALSE" = "black"), guide = "none") +
  geom_point(aes(x = -0.9, size = Count),
    color = "grey25", fill = "white", shape = 21, stroke = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  scale_size_continuous(range = c(1.4, 5.5), name = "Count") +
  labs(tag = "b", x = expression(-log[10] ~ "(p.adjust)"), y = NULL) +
  coord_cartesian(xlim = c(-1.8, max(kegg_panel$mlog10_fdr) + 0.8), clip = "off") +
  theme_pub +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 5.8),
    legend.title = element_text(size = 6.0),
    legend.key.size = unit(0.16, "in"),
    axis.line.y = element_blank()
  )

# ---- 2C: GO BP enrichment horizontal bar plot ----
go_lo <- read_first(c(file.path(TAB_DIR, "GO_lowCI_up.csv")))
go_lo <- head(go_lo[order(go_lo$p.adjust), ], 10)
go_lo$Description <- ifelse(nchar(go_lo$Description) > 44,
  paste0(substr(go_lo$Description, 1, 41), "..."), go_lo$Description)
go_lo$mlog10_fdr <- -log10(pmax(go_lo$p.adjust, .Machine$double.xmin))
go_lo$Description <- factor(go_lo$Description, levels = go_lo$Description[order(go_lo$mlog10_fdr)])
focus_go <- c("axonogenesis", "axon guidance", "neuron projection guidance")
go_lo$focus <- go_lo$Description %in% focus_go

f2c_go <- ggplot(go_lo, aes(x = mlog10_fdr, y = Description)) +
  geom_col(fill = "#0072B2", width = 0.75, alpha = 0.25, color = NA) +
  geom_point(aes(x = -1.4, size = Count),
    color = "grey25", fill = "white", shape = 21, stroke = 0.5) +
  geom_text(aes(label = Description, color = focus), x = 0.5, hjust = 0,
    fontface = "bold", size = 2.5, vjust = 0.5) +
  scale_color_manual(values = c("TRUE" = "#CB181D", "FALSE" = "black"), guide = "none") +
  scale_size_continuous(range = c(1.2, 4.2), name = "Count") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  labs(tag = "c", x = expression(-log[10] ~ "(p.adjust)"), y = NULL) +
  coord_cartesian(xlim = c(-2.4, max(go_lo$mlog10_fdr) + 0.8), clip = "off") +
  theme_pub +
  theme(
    aspect.ratio = 1.35,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 5.8),
    legend.title = element_text(size = 6.0),
    legend.key.size = unit(0.16, "in"),
    axis.line.y = element_blank()
  )

# ---- 2D-E: PD-risk genes with concordant cross-cohort direction ----
risk <- read_first(c(file.path(TAB_DIR, "PD_risk_gene_screen.csv")))
risk$GSE_plot <- -risk$GSE_logFC
risk$concordant <- sign(risk$PPMI_logFC) == sign(risk$GSE_plot)
risk_conc <- risk[risk$concordant, ]
binom_p <- binom.test(nrow(risk_conc), nrow(risk), 0.5)$p.value

make_effect_panel <- function(df, tag, direction_label, xlim, breaks) {
  df <- df[order(abs(df$PPMI_logFC)), ]
  df$Gene <- factor(df$Symbol, levels = df$Symbol)
  df_long <- df %>%
    select(Gene, PPMI = PPMI_logFC, GENEPARK = GSE_plot) %>%
    pivot_longer(c(PPMI, GENEPARK), names_to = "Cohort", values_to = "logFC")
  df_long$Cohort <- factor(df_long$Cohort, levels = c("PPMI", "GENEPARK"))

  ggplot(df_long, aes(logFC, Gene, fill = Cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_col(position = position_dodge2(width = 0.62, preserve = "single"),
      width = 0.52, alpha = 0.95) +
    scale_fill_manual(values = c(PPMI = PPMI_COL, GENEPARK = GSE_COL)) +
    scale_x_continuous(breaks = breaks, expand = expansion(mult = c(0.01, 0.04))) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(tag = tag, x = direction_label, y = NULL,
      caption = "Negative values indicate higher expression in the Low-CI group") +
    theme_pub +
    theme(
      aspect.ratio = 1,
      axis.text.y = element_text(face = "italic", size = 5.4),
      axis.text.x = element_text(size = 6.2),
      axis.title.x = element_text(size = 6.6, margin = margin(t = 2)),
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(0.24, "in"),
      axis.line.y = element_blank()
    )
}

# 3 cell-composition-robust genes (survive leukocyte-fraction adjustment)
robust_genes <- c("ITGA8", "RIMS1", "DNAH17")
low_panel <- risk[risk$Symbol %in% robust_genes, ]

f2d_risk_low <- make_effect_panel(
  low_panel, "d", "log2 fold change",
  c(-3.0, 0), c(-3, -2, -1, 0)
)

cd_title <- "Cell-composition-robust PD risk genes upregulated in the low-CI expression group"
fig2de <- f2d_risk_low +
  plot_annotation(
    title = cd_title,
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  )

# ---- 2F-G: immune-cell shifts, biologically ordered across cohorts ----
ppmi_imm <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "immune_cohend.csv"),
  file.path(TAB_DIR, "immune_cohend.csv")
))
ppmi_imm <- ppmi_imm[is.finite(ppmi_imm$d), c("Cell", "d", "P", "d_lo", "d_hi", "Padj")]
names(ppmi_imm)[names(ppmi_imm) == "d"] <- "D"

gse_imm <- read_first(c(
  file.path(TAB_DIR, "GSE99039_immune_NNLS.csv"),
  "E:/PPMI帕金森数据库专用/ppmi数据表/GSE99039/results/GSE99039_immune_NNLS.csv"
))
gse_imm <- gse_imm[is.finite(gse_imm$D), c("Cell", "D", "P", "Padj")]

add_ci <- function(df) {
  if (all(c("d_lo", "d_hi") %in% names(df))) {
    return(df)
  }
  z <- abs(qnorm(df$P / 2))
  se <- abs(df$D) / z
  se[!is.finite(se)] <- NA_real_
  df$d_lo <- df$D - 1.96 * se
  df$d_hi <- df$D + 1.96 * se
  df
}
ppmi_imm <- add_ci(ppmi_imm)
gse_imm <- add_ci(gse_imm)

imm <- merge(ppmi_imm, gse_imm, by = "Cell", suffixes = c("_PPMI", "_GSE"))
imm <- imm[is.finite(imm$D_PPMI) & is.finite(imm$D_GSE), ]
imm$concordant <- sign(imm$D_PPMI) == sign(imm$D_GSE)
imm$both_sig <- imm$Padj_PPMI < 0.05 & imm$Padj_GSE < 0.05
imm$either_sig <- imm$Padj_PPMI < 0.05 | imm$Padj_GSE < 0.05
imm$ci_double_positive <- imm$d_lo_PPMI > 0 & imm$d_lo_GSE > 0
imm$ci_double_negative <- imm$d_hi_PPMI < 0 & imm$d_hi_GSE < 0
imm$ci_concordant <- imm$ci_double_positive | imm$ci_double_negative
imm$priority <- dplyr::case_when(
  imm$ci_concordant & imm$both_sig ~ 1,
  imm$ci_concordant ~ 2,
  imm$concordant & imm$either_sig ~ 3,
  imm$either_sig ~ 4,
  TRUE ~ 5
)

bio_order <- c(
  "Neutrophils",
  "Dendritic cells resting", "Dendritic cells activated",
  "Macrophages M0", "Macrophages M2", "Macrophages M1",
  "Monocytes", "Eosinophils", "NK cells resting", "Mast cells resting",
  "T cells CD8", "T cells CD4 memory resting", "T cells CD4 memory activated",
  "T cells regulatory (Tregs)", "T cells follicular helper",
  "B cells naive", "B cells memory", "Plasma cells"
)
imm$bio_rank <- match(imm$Cell, bio_order)
imm$bio_rank[is.na(imm$bio_rank)] <- length(bio_order) + seq_len(sum(is.na(match(imm$Cell, bio_order))))
imm <- imm[order(imm$priority, imm$bio_rank, -abs(imm$D_PPMI)), ]

display_cells <- imm$Cell
label_cells <- ifelse(imm$ci_concordant,
  paste0("<span style='color:", SIG_COL, "'>", display_cells, "</span>"),
  display_cells)
names(label_cells) <- display_cells

prep_immune_panel <- function(imm, cohort) {
  if (cohort == "PPMI") {
    df <- data.frame(
      Cell = imm$Cell,
      Cell_label = label_cells[imm$Cell],
      d = imm$D_PPMI,
      d_lo = imm$d_lo_PPMI,
      d_hi = imm$d_hi_PPMI,
      sig = ifelse(imm$Padj_PPMI < 0.05, "padj<0.05", "ns")
    )
  } else {
    df <- data.frame(
      Cell = imm$Cell,
      Cell_label = label_cells[imm$Cell],
      d = imm$D_GSE,
      d_lo = imm$d_lo_GSE,
      d_hi = imm$d_hi_GSE,
      sig = ifelse(imm$Padj_GSE < 0.05, "padj<0.05", "ns")
    )
  }
  df$Cell_label <- factor(df$Cell_label, levels = rev(unname(label_cells[display_cells])))
  df$sig <- factor(df$sig, levels = c("padj<0.05", "ns"))
  df
}

make_immune_forest <- function(df, tag, cohort, xlim) {
  ggplot(df, aes(d, Cell_label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_errorbar(aes(xmin = d_lo, xmax = d_hi), width = 0.3, linewidth = 0.55,
      orientation = "y", na.rm = TRUE) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("padj<0.05" = SIG_COL, "ns" = NS_COL), drop = FALSE) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(tag = tag, x = "Cohen's d", y = cohort, color = NULL) +
    theme_pub +
    theme(
      legend.position = "none",
      axis.text.y = ggtext::element_markdown(size = 5.8),
      axis.text.x = element_text(size = 6.2),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 8.5, angle = 90, margin = margin(r = 4)),
      aspect.ratio = 1,
      axis.line.y = element_blank()
    )
}

imm_ppmi <- prep_immune_panel(imm, "PPMI")
imm_gse <- prep_immune_panel(imm, "GENEPARK")
f2f_imm_ppmi <- make_immune_forest(imm_ppmi, "e", "PPMI", c(-1.15, 1.05))
f2g_imm_gse <- make_immune_forest(imm_gse, "f", "GENEPARK", c(-1.15, 1.15))

fig2fg <- (f2f_imm_ppmi | f2g_imm_gse) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Immune-cell shifts ordered by cross-cohort concordance and cell lineage",
    theme = theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  )

# ---- 2H: k-means immune robustness ----
immune_frac <- read_first(c(
  file.path(FALLBACK_DATA_DIR, "immune_fractions.csv"),
  file.path(TAB_DIR, "immune_fractions.csv"),
  file.path(TAB_DIR, "Table_immune_fractions.csv")
), check.names = FALSE)
cluster_df <- read_first(c(
  file.path(TAB_DIR, "Table_clustering_methods.csv"),
  file.path(FALLBACK_DATA_DIR, "PD_all_clustering_methods.csv")
))
cluster_df <- cluster_df[, c("SAMPLE_ID", "KM_group")]
km_dat <- merge(immune_frac, cluster_df, by = "SAMPLE_ID")
km_dat <- km_dat[km_dat$KM_group %in% c("Low", "High"), ]

cohen_ci <- function(lo, hi) {
  lo <- lo[is.finite(lo)]
  hi <- hi[is.finite(hi)]
  n_lo <- length(lo)
  n_hi <- length(hi)
  pooled_sd <- sd(c(lo, hi))
  d <- (mean(lo) - mean(hi)) / pooled_sd
  se <- sqrt((n_lo + n_hi) / (n_lo * n_hi) + d^2 / (2 * (n_lo + n_hi - 2)))
  c(d = d, d_lo = d - 1.96 * se, d_hi = d + 1.96 * se,
    p = wilcox.test(lo, hi)$p.value)
}

km_cells <- display_cells[display_cells %in% names(km_dat)]
km_imm <- do.call(rbind, lapply(km_cells, function(cell) {
  est <- cohen_ci(km_dat[km_dat$KM_group == "Low", cell],
    km_dat[km_dat$KM_group == "High", cell])
  data.frame(Cell = cell, t(est), check.names = FALSE)
}))
km_imm$padj <- p.adjust(km_imm$p, "BH")
km_imm$sig <- factor(ifelse(km_imm$padj < 0.05, "padj<0.05", "ns"),
  levels = c("padj<0.05", "ns"))
km_labels <- ifelse(km_imm$padj < 0.05,
  paste0("<span style='color:", SIG_COL, "'>", km_imm$Cell, "</span>"),
  km_imm$Cell)
km_imm$Cell_label <- factor(km_labels, levels = rev(km_labels))
write.csv(km_imm[, c("Cell", "d", "d_lo", "d_hi", "p", "padj")],
  file.path(TAB_DIR, "kmeans_immune_full.csv"), row.names = FALSE)

f2h_kmeans <- ggplot(km_imm, aes(d, Cell_label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.25) +
  geom_errorbar(aes(xmin = d_lo, xmax = d_hi), width = 0.3, linewidth = 0.55,
    orientation = "y", na.rm = TRUE) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("padj<0.05" = SIG_COL, "ns" = NS_COL), drop = FALSE) +
  coord_cartesian(xlim = c(-0.65, 0.65), clip = "off") +
  labs(tag = "g", x = "Cohen's d\n(k-means low-CI cluster)", y = "k-means", color = NULL) +
  theme_pub +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.text.y = ggtext::element_markdown(size = 5.8),
    axis.text.x = element_text(size = 6.2),
    axis.title.x = element_text(size = 6.7),
    axis.title.y = element_text(size = 8.5, angle = 90, margin = margin(r = 4))
  )

# (Panel e removed: legacy 9-gene cell-composition panel was dropped in favour of the
# 10-gene PD-risk framework; see Supplementary Table S13 for cell-composition sensitivity.)

fig2 <- patchwork::wrap_plots(
  f2a, f2b_gsea, f2c_go, f2d_risk_low,
  f2f_imm_ppmi, f2g_imm_gse, f2h_kmeans,
  ncol = 4
)

ggsave(file.path(FIG_DIR, "Figure2B_KEGG_GSEA.pdf"), f2b_gsea, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2B_KEGG_GSEA.png"), f2b_gsea, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2DE_cross_cohort.pdf"), fig2de, width = 4, height = 4.5)
ggsave(file.path(FIG_DIR, "Figure2DE_cross_cohort.png"), fig2de, width = 4, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2CD_cross_cohort.pdf"), fig2de, width = 4, height = 4.5)
ggsave(file.path(FIG_DIR, "Figure2CD_cross_cohort.png"), fig2de, width = 4, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2FG_immune_cross_cohort.pdf"), fig2fg, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2FG_immune_cross_cohort.png"), fig2fg, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2EF_immune_cross_cohort.pdf"), fig2fg, width = 7.4, height = 4.1)
ggsave(file.path(FIG_DIR, "Figure2EF_immune_cross_cohort.png"), fig2fg, width = 7.4, height = 4.1, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2H_kmeans_immune.pdf"), f2h_kmeans, width = 3.6, height = 3.6)
ggsave(file.path(FIG_DIR, "Figure2H_kmeans_immune.png"), f2h_kmeans, width = 3.6, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure2.pdf"), fig2, width = 12, height = 11)
ggsave(file.path(FIG_DIR, "Figure2.png"), fig2, width = 12, height = 11, dpi = 300)

message(sprintf(
  "Saved Figure 2 outputs to %s; PD-risk concordance: %d/%d; immune top cells: %s",
  FIG_DIR, nrow(risk_conc), nrow(risk), paste(head(display_cells, 6), collapse = ", ")
))

# =====================================================================
# SECTION 24 | Figure 3: year-5 PIGD divergence and biomarker context
# =====================================================================

# ============================================================
# Figure 3: PIGD trajectory, specificity, and biomarker context
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
FIG_DIR <- file.path(ROOT, "02_figures")
TAB_DIR <- file.path(ROOT, "04_tables")
DATA_DIR <- normalizePath(file.path(ROOT, "..", "figures_pc1", "data"), winslash = "/", mustWork = TRUE)

theme_fig3 <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#333333"),
      axis.text = element_text(colour = "#333333"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, colour = "#555555"),
      legend.title = element_blank(),
      legend.key.size = unit(0.28, "cm"),
      legend.position = "top",
      strip.background = element_rect(fill = "#f2f2f2", colour = NA),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

pal_group <- c("High-CI" = "#0072B2", "Low-CI" = "#d55e00")

fmt_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
}

fmt_q <- function(q) {
  ifelse(is.na(q), "NA", ifelse(q < 0.001, "<0.001", sprintf("%.2f", q)))
}

vmap <- c(BL = 0, V04 = 1, V06 = 2, V08 = 3, V10 = 4, V12 = 5)

classic <- tibble::tribble(
  ~Biomarker,                  ~Label,                  ~Category,
  "CSF_Alpha_synuclein_CSF",   "alpha-synuclein CSF",   "Synucleinopathy",
  "AlphaSyn_SAA_Fmax_CSF",     "CSF SAA",               "Synucleinopathy",
  "NfL_Serum",                 "NfL serum",             "Neuroaxonal / glial",
  "NFL_Plasma",                "NfL plasma",            "Neuroaxonal / glial",
  "NFL_CSF",                   "NfL CSF",               "Neuroaxonal / glial",
  "GFAP_Plasma",               "GFAP plasma",           "Neuroaxonal / glial",
  "GFAP_CSF",                  "GFAP CSF",              "Neuroaxonal / glial",
  "tTau_CSF",                  "total tau CSF",         "AD-type co-pathology",
  "pTau_CSF",                  "p-tau CSF",             "AD-type co-pathology",
  "pTau181_CSF",               "p-tau181 CSF",          "AD-type co-pathology",
  "Ptau217p_Plasma",           "p-tau217 plasma",       "AD-type co-pathology",
  "ABeta42_CSF",               "Abeta42 CSF",           "AD-type co-pathology"
)

pc1 <- read.csv(file.path(DATA_DIR, "pc1_scores.csv"), check.names = FALSE) %>%
  dplyr::select(PATNO, group)

bm_long <- read.csv(file.path(DATA_DIR, "biomarkers_long.csv"), check.names = FALSE) %>%
  filter(biomarker %in% setdiff(classic$Biomarker, "AlphaSyn_SAA_Fmax_CSF")) %>%
  group_by(PATNO, EVENT_ID, biomarker) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

saa_file <- normalizePath(file.path(ROOT, "..", "..", "..", "..", "..", "..", "..",
  "运动症状数据", "SAA_Biospecimen_Analysis_Results_23Feb2026.csv"),
  winslash = "/", mustWork = FALSE)

if (file.exists(saa_file)) {
  saa_raw <- read.csv(saa_file, check.names = FALSE)
  fmax_cols <- c("Fmax_24h_Rep1", "Fmax_24h_Rep2", "Fmax_24h_Rep3")
  saa_bm <- saa_raw %>%
    filter(CLINICAL_EVENT %in% names(vmap)) %>%
    mutate(
      across(all_of(fmax_cols), ~ suppressWarnings(as.numeric(.x))),
      value = rowMeans(as.data.frame(across(all_of(fmax_cols))), na.rm = TRUE),
      EVENT_ID = CLINICAL_EVENT,
      biomarker = "AlphaSyn_SAA_Fmax_CSF"
    ) %>%
    filter(is.finite(value)) %>%
    group_by(PATNO, EVENT_ID, biomarker) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  saa_status <- saa_raw %>%
    filter(CLINICAL_EVENT == "BL") %>%
    distinct(PATNO, .keep_all = TRUE) %>%
    inner_join(pc1, by = "PATNO") %>%
    mutate(SAA_positive = grepl("Positive", SAA_Status, ignore.case = TRUE))

  saa_tab <- table(saa_status$group, saa_status$SAA_positive)
  saa_fisher <- fisher.test(saa_tab)
  saa_status_res <- data.frame(
    Metric = "Baseline alpha-syn SAA positivity",
    N_high = sum(saa_status$group == "High"),
    N_low = sum(saa_status$group == "Low"),
    High_positive_pct = 100 * mean(saa_status$SAA_positive[saa_status$group == "High"]),
    Low_positive_pct = 100 * mean(saa_status$SAA_positive[saa_status$group == "Low"]),
    Odds_ratio_Low_vs_High = unname(saa_fisher$estimate),
    P = saa_fisher$p.value
  )
  write.csv(saa_status_res, file.path(TAB_DIR, "SAA_baseline_status.csv"), row.names = FALSE)
} else {
  saa_bm <- bm_long[0, ]
  saa_status_res <- NULL
}

bm <- bind_rows(bm_long, saa_bm) %>%
  inner_join(pc1, by = "PATNO") %>%
  mutate(
    Time = unname(vmap[EVENT_ID]),
    PATNO = factor(PATNO),
    group = factor(group, levels = c("High", "Low"))
  ) %>%
  filter(!is.na(Time), is.finite(value))

sweep_res <- read.csv(file.path(TAB_DIR, "biomarker_LMM_sweep.csv"), check.names = FALSE) %>%
  dplyr::select(Biomarker, LastVisit, Grp_padj)

fit_biomarker <- function(b) {
  d <- bm %>% filter(biomarker == b)
  if (nrow(d) < 40 || n_distinct(d$group) < 2) return(NULL)
  if (all(d$value > 0, na.rm = TRUE)) d$value <- log(d$value)
  d$zvalue <- as.numeric(scale(d$value))
  d <- d[is.finite(d$zvalue), ]

  m <- tryCatch(
    lmer(zvalue ~ Time * group + (1 + Time | PATNO), data = d, REML = FALSE),
    error = function(e) tryCatch(
      lmer(zvalue ~ Time * group + (1 | PATNO), data = d, REML = FALSE),
      error = function(e2) NULL
    )
  )
  if (is.null(m)) return(NULL)
  cm <- coef(summary(m))
  rn <- "groupLow"
  if (!rn %in% rownames(cm)) return(NULL)

  data.frame(
    Biomarker = b,
    N_obs = nrow(d),
    N_subj = n_distinct(d$PATNO),
    Estimate = cm[rn, "Estimate"],
    SE = cm[rn, "Std. Error"],
    P = cm[rn, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

# A. PIGD estimated marginal means
emm <- read.csv(file.path(DATA_DIR, "lmm_emm.csv"), check.names = FALSE) %>%
  filter(Outcome == "PIGD") %>%
  mutate(
    VISIT = factor(VISIT, levels = c("BL", "V04", "V06", "V08", "V10", "V12"),
      labels = c("BL", "Y1", "Y2", "Y3", "Y4", "Y5")),
    group = recode(group, "High" = "High-CI", "Low" = "Low-CI")
  )

p_a <- ggplot(emm, aes(VISIT, emmean, colour = group, group = group)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = group), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = pal_group) +
  scale_fill_manual(values = pal_group) +
  labs(title = "PIGD trajectory", y = "Estimated PIGD score", x = NULL) +
  theme_fig3()

# B. Visit-specific group x time interactions
ints <- read.csv(file.path(DATA_DIR, "lmm_interactions.csv"), check.names = FALSE) %>%
  filter(Outcome == "PIGD") %>%
  mutate(
    Visit = factor(Visit, levels = c("V04", "V06", "V08", "V10", "V12"),
      labels = c("Y1", "Y2", "Y3", "Y4", "Y5")),
    lcl = Est - 1.96 * SE,
    ucl = Est + 1.96 * SE,
    sig = P < 0.05
  )

p_b <- ggplot(ints, aes(Visit, Est)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#777777") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.12, linewidth = 0.45, colour = "#555555") +
  geom_point(aes(fill = sig), shape = 21, size = 2.4, colour = "#333333") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  labs(title = "Visit-specific divergence", y = "Additional \u0394PIGD\nLow-CI vs High-CI", x = NULL) +
  theme_fig3() +
  theme(legend.position = "none")

# D. Continuous PC1 sensitivity
pc1_cont <- read.csv(file.path(TAB_DIR, "PIGD_continuous_PC1_visit_sensitivity.csv"),
  check.names = FALSE) %>%
  mutate(
    Year = factor(Year, levels = c("Y1", "Y2", "Y3", "Y4", "Y5")),
    sig = P < 0.05
  )

p_c <- ggplot(pc1_cont, aes(Year, Additional_dPIGD_lowCI_per_1SD)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#777777") +
  geom_errorbar(aes(ymin = LCL_lowCI, ymax = UCL_lowCI), width = 0.12,
    linewidth = 0.45, colour = "#555555") +
  geom_point(aes(fill = sig), shape = 21, size = 2.4, colour = "#333333") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  labs(title = "Continuous PC1 sensitivity", subtitle = "Y5 p = 0.008",
    y = "Additional \u0394PIGD\nper 1-SD lower PC1", x = NULL) +
  theme_fig3() +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "Figure3D_continuous_PC1_sensitivity.pdf"), p_c,
  width = 4.1, height = 3.4, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3D_continuous_PC1_sensitivity.png"), p_c,
  width = 4.1, height = 3.4, dpi = 600)

# I. Partition sensitivity for the year-5 PIGD signal
primary_y5 <- ints %>%
  filter(Visit == "Y5") %>%
  transmute(
    Method = "PC1 median",
    Effect = Est,
    SE = SE,
    LCL = Est - 1.96 * SE,
    UCL = Est + 1.96 * SE,
    P = P,
    Source = "figures_pc1/data/lmm_interactions.csv"
  )

partition_y5 <- read.csv(file.path(TAB_DIR, "Table_LMM_factor.csv"), check.names = FALSE) %>%
  filter(Outcome == "PIGD", Year == "Y5", Track %in% c("KMeans", "Hierarchical")) %>%
  transmute(
    Method = recode(Track, "KMeans" = "k-means", "Hierarchical" = "Hierarchical"),
    Effect = as.numeric(Effect),
    SE = NA_real_,
    LCL = NA_real_,
    UCL = NA_real_,
    P = as.numeric(P),
    Source = "04_tables/Table_LMM_factor.csv"
  )

partition_sens <- bind_rows(primary_y5, partition_y5) %>%
  mutate(
    Method = factor(Method, levels = c("Hierarchical", "k-means", "PC1 median")),
    p_label = ifelse(P < 0.001, "p<0.001", sprintf("p=%.3f", P)),
    label = sprintf("\u03b2=%.2f, %s", Effect, p_label)
  )

write.csv(partition_sens, file.path(TAB_DIR, "PIGD_partition_sensitivity.csv"), row.names = FALSE)

p_part <- ggplot(partition_sens, aes(Effect, Method)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#777777") +
  geom_point(aes(fill = P < 0.05), shape = 21, size = 2.6, colour = "#333333") +
  geom_text(aes(label = label), nudge_x = 0.04, size = 2.2, hjust = 0) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  coord_cartesian(xlim = c(-0.05, max(partition_sens$Effect, na.rm = TRUE) + 0.34)) +
  labs(title = "Partition sensitivity", subtitle = "Y5 PIGD interaction",
    y = NULL, x = "Y5 additional \u0394PIGD") +
  theme_fig3() +
  theme(legend.position = "none")

# D. Clinical specificity screen
clin <- read.csv(file.path(TAB_DIR, "LMM_sweep_COMPREHENSIVE.csv"), check.names = FALSE) %>%
  mutate(
    neglog = -log10(LastV_int_padj),
    pass = LastV_int_padj < 0.05,
    Outcome2 = recode(Outcome,
      "UPDRS3_axial" = "UPDRS III axial",
      "MoCA_attention" = "MoCA attention",
      "MoCA_orientation" = "MoCA orientation",
      "SCOPA_total" = "SCOPA total")
  ) %>%
  arrange(LastV_int_padj) %>%
  slice_head(n = 12) %>%
  mutate(Outcome2 = factor(Outcome2, levels = rev(Outcome2)))

p_d <- ggplot(clin, aes(neglog, Outcome2)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = 0, xend = neglog, yend = Outcome2), linewidth = 0.45, colour = "#b8b8b8") +
  geom_point(aes(fill = pass), shape = 21, size = 2.5, colour = "#333333") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  labs(title = "Clinical specificity screen", y = NULL, x = expression(-log[10]("BH padj"))) +
  theme_fig3() +
  theme(legend.position = "none")

# E. Representative DAT-SPECT trajectory
dat_traj <- read.csv(file.path(TAB_DIR, "DAT_posterior_dorsal_putamen_trajectory.csv"),
  check.names = FALSE) %>%
  mutate(
    VISIT = factor(VISIT, levels = c("BL", "V04", "V06", "V10"),
      labels = c("BL", "Y1", "Y2", "Y4")),
    group = recode(group, "High" = "High-CI", "Low" = "Low-CI")
  )

dat_y4_p <- unique(dat_traj$Last_visit_interaction_p)[1]

p_e <- ggplot(dat_traj, aes(VISIT, emmean, colour = group, group = group)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = group),
    alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = pal_group) +
  scale_fill_manual(values = pal_group) +
  labs(
    title = "Posterior dorsal putamen DAT-SBR",
    subtitle = sprintf("Last available visit: Y4; Y4 x group p = %.3f", dat_y4_p),
    y = "Estimated DAT-SBR",
    x = NULL
  ) +
  theme_fig3() +
  theme(legend.position = "none")

# G. Established fluid biomarkers
effects <- bind_rows(lapply(classic$Biomarker, fit_biomarker)) %>%
  left_join(classic, by = "Biomarker") %>%
  left_join(sweep_res, by = "Biomarker") %>%
  mutate(
    LCL = Estimate - 1.96 * SE,
    UCL = Estimate + 1.96 * SE,
    Focused_padj = p.adjust(P, method = "BH"),
    Nominal = P < 0.05,
    Label = factor(Label, levels = rev(classic$Label)),
    Category = factor(Category, levels = c("Synucleinopathy", "Neuroaxonal / glial", "AD-type co-pathology")),
    Stat_label = sprintf("n=%d, p=%s, q=%s", N_subj, fmt_p(P), fmt_q(Focused_padj))
  )

write.csv(effects, file.path(TAB_DIR, "classic_biomarker_effects.csv"), row.names = FALSE)
write.csv(effects, file.path(TAB_DIR, "TableS10_classic_biomarkers.csv"), row.names = FALSE)

bio_x_min <- min(effects$LCL, na.rm = TRUE) - 0.08
bio_x_max <- max(effects$UCL, na.rm = TRUE) + 0.08
bio_x_text <- bio_x_max + 0.04
saa_subtitle <- if (!is.null(saa_status_res)) {
  sprintf("Focused panel; baseline SAA positivity Low vs High %.1f%% vs %.1f%%, p=%.3f",
    saa_status_res$Low_positive_pct, saa_status_res$High_positive_pct, saa_status_res$P)
} else {
  "Focused panel; SAA source file not found"
}

p_f <- ggplot(effects, aes(Estimate, Label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = LCL, xmax = UCL), width = 0.22,
    linewidth = 0.42, colour = "#555555", orientation = "y") +
  geom_point(aes(fill = Nominal), shape = 21, size = 2.3, colour = "#333333") +
  geom_text(aes(x = bio_x_text, label = Stat_label), hjust = 0, size = 2.15, colour = "#333333") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  coord_cartesian(xlim = c(bio_x_min, bio_x_text + 0.62), clip = "off") +
  labs(title = "Established fluid biomarkers", subtitle = saa_subtitle,
    x = "Low-CI vs High-CI standardized difference", y = NULL) +
  theme_fig3(base_size = 7) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 6.2),
    panel.spacing.y = unit(0.38, "lines"),
    plot.margin = margin(5.5, 82, 5.5, 5.5)
  )

ggsave(file.path(FIG_DIR, "Figure3G_classic_biomarkers.pdf"), p_f,
  width = 9.2, height = 4.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3G_classic_biomarkers.png"), p_f,
  width = 9.2, height = 4.2, dpi = 600)

# G. Compact clinical biomarker context
p_g <- p_f +
  theme(
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "none"
  )

fig3_top <- wrap_plots(p_a, p_b, p_part, p_c, p_d, p_e, ncol = 3)
fig3_mid <- wrap_plots(plot_spacer(), p_g, plot_spacer(), ncol = 3) +
  plot_layout(widths = c(0.45, 0.85, 0.45))

fig3 <- fig3_top / fig3_mid +
  plot_layout(heights = c(2.0, 0.78)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(file.path(FIG_DIR, "Figure3_PIGD_clinical_summary.pdf"), fig3,
  width = 12.5, height = 10.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3_PIGD_clinical_summary.png"), fig3,
  width = 12.5, height = 10.2, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure3.pdf"), fig3,
  width = 12.5, height = 10.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3.png"), fig3,
  width = 12.5, height = 10.2, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure3_continuousPC1_updated.pdf"), fig3,
  width = 12.5, height = 10.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3_continuousPC1_updated.png"), fig3,
  width = 12.5, height = 10.2, dpi = 600)

cat("Saved Figure 3 to:\n")
cat(file.path(FIG_DIR, "Figure3.pdf"), "\n")
cat(file.path(FIG_DIR, "Figure3.png"), "\n")
cat(file.path(FIG_DIR, "Figure3_PIGD_clinical_summary.pdf"), "\n")
cat(file.path(FIG_DIR, "Figure3_PIGD_clinical_summary.png"), "\n")
cat(file.path(FIG_DIR, "Figure3D_continuous_PC1_sensitivity.png"), "\n")
cat(file.path(FIG_DIR, "Figure3_continuousPC1_updated.png"), "\n")

# =====================================================================
# SECTION 25 | Classic biomarker figure
# =====================================================================

# ============================================================
# Established fluid biomarker check
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
FIG_DIR <- file.path(ROOT, "02_figures")
TAB_DIR <- file.path(ROOT, "04_tables")
DATA_DIR <- normalizePath(file.path(ROOT, "..", "figures_pc1", "data"), winslash = "/", mustWork = TRUE)

theme_s4 <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#333333"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#333333"),
      axis.text = element_text(colour = "#333333"),
      axis.title = element_text(colour = "#222222"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = "#555555"),
      strip.background = element_rect(fill = "#f2f2f2", colour = NA),
      strip.text.y = element_text(face = "bold", size = base_size - 1, angle = 0),
      plot.margin = margin(5.5, 80, 5.5, 5.5)
    )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "p=NA", ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)))
}

fmt_q <- function(q) {
  ifelse(is.na(q), "NA", ifelse(q < 0.001, "<0.001", sprintf("%.2f", q)))
}

vmap <- c(BL = 0, V04 = 1, V06 = 2, V08 = 3, V10 = 4, V12 = 5)

classic <- tibble::tribble(
  ~Biomarker,                  ~Label,                  ~Category,
  "CSF_Alpha_synuclein_CSF",   "alpha-synuclein CSF",   "Synucleinopathy",
  "AlphaSyn_SAA_Fmax_CSF",     "CSF SAA",               "Synucleinopathy",
  "NfL_Serum",                 "NfL serum",             "Neuroaxonal / glial",
  "NFL_Plasma",                "NfL plasma",            "Neuroaxonal / glial",
  "NFL_CSF",                   "NfL CSF",               "Neuroaxonal / glial",
  "GFAP_Plasma",               "GFAP plasma",           "Neuroaxonal / glial",
  "GFAP_CSF",                  "GFAP CSF",              "Neuroaxonal / glial",
  "tTau_CSF",                  "total tau CSF",         "AD-type co-pathology",
  "pTau_CSF",                  "p-tau CSF",             "AD-type co-pathology",
  "pTau181_CSF",               "p-tau181 CSF",          "AD-type co-pathology",
  "Ptau217p_Plasma",           "p-tau217 plasma",       "AD-type co-pathology",
  "ABeta42_CSF",               "Abeta42 CSF",           "AD-type co-pathology"
)

pc1 <- read.csv(file.path(DATA_DIR, "pc1_scores.csv"), check.names = FALSE) %>%
  dplyr::select(PATNO, group)

bm_long <- read.csv(file.path(DATA_DIR, "biomarkers_long.csv"), check.names = FALSE) %>%
  filter(biomarker %in% setdiff(classic$Biomarker, "AlphaSyn_SAA_Fmax_CSF")) %>%
  group_by(PATNO, EVENT_ID, biomarker) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

saa_file <- normalizePath(file.path(ROOT, "..", "..", "..", "..", "..", "..", "..",
  "运动症状数据", "SAA_Biospecimen_Analysis_Results_23Feb2026.csv"),
  winslash = "/", mustWork = FALSE)

if (file.exists(saa_file)) {
  saa_raw <- read.csv(saa_file, check.names = FALSE)
  fmax_cols <- c("Fmax_24h_Rep1", "Fmax_24h_Rep2", "Fmax_24h_Rep3")
  saa_bm <- saa_raw %>%
    filter(CLINICAL_EVENT %in% names(vmap)) %>%
    mutate(
      across(all_of(fmax_cols), ~ suppressWarnings(as.numeric(.x))),
      value = rowMeans(as.data.frame(across(all_of(fmax_cols))), na.rm = TRUE),
      EVENT_ID = CLINICAL_EVENT,
      biomarker = "AlphaSyn_SAA_Fmax_CSF"
    ) %>%
    filter(is.finite(value)) %>%
    group_by(PATNO, EVENT_ID, biomarker) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  saa_status <- saa_raw %>%
    filter(CLINICAL_EVENT == "BL") %>%
    distinct(PATNO, .keep_all = TRUE) %>%
    inner_join(pc1, by = "PATNO") %>%
    mutate(SAA_positive = grepl("Positive", SAA_Status, ignore.case = TRUE))

  saa_tab <- table(saa_status$group, saa_status$SAA_positive)
  saa_fisher <- fisher.test(saa_tab)
  saa_status_res <- data.frame(
    Metric = "Baseline alpha-syn SAA positivity",
    N_high = sum(saa_status$group == "High"),
    N_low = sum(saa_status$group == "Low"),
    High_positive_pct = 100 * mean(saa_status$SAA_positive[saa_status$group == "High"]),
    Low_positive_pct = 100 * mean(saa_status$SAA_positive[saa_status$group == "Low"]),
    Odds_ratio_Low_vs_High = unname(saa_fisher$estimate),
    P = saa_fisher$p.value
  )
  write.csv(saa_status_res, file.path(TAB_DIR, "SAA_baseline_status.csv"), row.names = FALSE)
} else {
  saa_bm <- bm_long[0, ]
}

bm <- bind_rows(bm_long, saa_bm) %>%
  inner_join(pc1, by = "PATNO") %>%
  mutate(
    Time = unname(vmap[EVENT_ID]),
    PATNO = factor(PATNO),
    group = factor(group, levels = c("High", "Low"))
  ) %>%
  filter(!is.na(Time), is.finite(value))

sweep_res <- read.csv(file.path(TAB_DIR, "biomarker_LMM_sweep.csv"), check.names = FALSE) %>%
  dplyr::select(Biomarker, LastVisit, Grp_padj)

fit_one <- function(b) {
  d <- bm %>% filter(biomarker == b)
  if (nrow(d) < 40 || n_distinct(d$group) < 2) return(NULL)
  if (all(d$value > 0, na.rm = TRUE)) d$value <- log(d$value)
  d$zvalue <- as.numeric(scale(d$value))
  d <- d[is.finite(d$zvalue), ]

  m <- tryCatch(
    lmer(zvalue ~ Time * group + (1 + Time | PATNO), data = d, REML = FALSE),
    error = function(e) tryCatch(
      lmer(zvalue ~ Time * group + (1 | PATNO), data = d, REML = FALSE),
      error = function(e2) NULL
    )
  )
  if (is.null(m)) return(NULL)
  cm <- coef(summary(m))
  rn <- "groupLow"
  if (!rn %in% rownames(cm)) return(NULL)

  data.frame(
    Biomarker = b,
    N_obs = nrow(d),
    N_subj = n_distinct(d$PATNO),
    Estimate = cm[rn, "Estimate"],
    SE = cm[rn, "Std. Error"],
    P = cm[rn, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

effects <- bind_rows(lapply(classic$Biomarker, fit_one)) %>%
  left_join(classic, by = "Biomarker") %>%
  left_join(sweep_res, by = "Biomarker") %>%
  mutate(
    LCL = Estimate - 1.96 * SE,
    UCL = Estimate + 1.96 * SE,
    Focused_padj = p.adjust(P, method = "BH"),
    Nominal = P < 0.05,
    Label = factor(Label, levels = rev(classic$Label)),
    Category = factor(Category, levels = c("Synucleinopathy", "Neuroaxonal / glial", "AD-type co-pathology")),
    Stat_label = sprintf("n=%d, %s, q=%s", N_subj, fmt_p(P), fmt_q(Focused_padj))
  )

write.csv(effects, file.path(TAB_DIR, "classic_biomarker_effects.csv"), row.names = FALSE)
write.csv(effects, file.path(TAB_DIR, "TableS10_classic_biomarkers.csv"), row.names = FALSE)

x_min <- min(effects$LCL, na.rm = TRUE) - 0.08
x_max <- max(effects$UCL, na.rm = TRUE) + 0.08
x_text <- x_max + 0.04

p_s4 <- ggplot(effects, aes(Estimate, Label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = LCL, xmax = UCL), width = 0.24,
    linewidth = 0.45, colour = "#555555", orientation = "y") +
  geom_point(aes(fill = Nominal), shape = 21, size = 2.5, colour = "#333333") +
  geom_text(aes(x = x_text, label = Stat_label), hjust = 0, size = 2.25, colour = "#333333") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#d55e00")) +
  coord_cartesian(xlim = c(x_min, x_text + 0.55), clip = "off") +
  labs(
    title = "Established fluid biomarkers",
    subtitle = "Visit-adjusted Low-CI vs High-CI differences; q values from focused-panel BH correction",
    x = "Standardized difference (SD units)",
    y = NULL
  ) +
  theme_s4() +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.spacing.y = unit(0.45, "lines")
  )

ggsave(file.path(FIG_DIR, "FigS4_classic_biomarkers.pdf"), p_s4,
  width = 7.8, height = 5.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "FigS4_classic_biomarkers.png"), p_s4,
  width = 7.8, height = 5.2, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure3G_classic_biomarkers.pdf"), p_s4,
  width = 7.8, height = 5.2, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3G_classic_biomarkers.png"), p_s4,
  width = 7.8, height = 5.2, dpi = 600)

cat("Saved classic biomarker figure to:\n")
cat(file.path(FIG_DIR, "FigS4_classic_biomarkers.pdf"), "\n")
cat(file.path(FIG_DIR, "FigS4_classic_biomarkers.png"), "\n")
cat(file.path(TAB_DIR, "classic_biomarker_effects.csv"), "\n")

# =====================================================================
# SECTION 26 | Supplementary Figure S1: structure assessment
# =====================================================================

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

