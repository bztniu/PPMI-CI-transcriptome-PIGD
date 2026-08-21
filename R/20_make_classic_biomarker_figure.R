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
