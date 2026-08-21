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
