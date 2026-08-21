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
