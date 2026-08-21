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
