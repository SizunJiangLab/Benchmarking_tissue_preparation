#!/usr/bin/env Rscript
# =============================================================================
# Factorial CV reanalysis - Stage 1: models + tidy tables
#
# Per-site factorial models of marker-level coefficient of variation (CV):
#   BIDMC    : CV ~ buffer + hier_duration + staining + Marker   (24 conditions)
#   Stanford : CV ~ hier_duration + staining + Marker            (12, Tris/EDTA only)
#   Roche    : CV ~ buffer + hier_duration + Marker              (6,  25 C 1 h only)
# Type II ANOVA -> partial eta-squared (effect size); emmeans -> adjusted mean CV.
#
# Consolidated R-only port of:
#   rebuttal/scripts/R1_5_factorial_model_all_sites.R   (the statistical core)
#   rebuttal/scripts/R3_build_viz_input.py              (pure relabel, no modeling)
#   rebuttal/scripts/R3_fig3_master_table.py            (ranking + traffic-light colors)
#
# Run from the repo root:  Rscript workflows/factorial_cv_model.R
# Inputs  (gitignored, produced by the main Mesmer workflow with *_all configs):
#   results/out_{BIDMC,Stanford,Roche}_all/cv_values_long.csv
#   results/out_{BIDMC,Stanford,Roche}_all/condition_summary.csv
# Outputs -> results/factorial_cv/{effect_size,adjusted_mean_cv,ranking_master_table}.csv
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(car)          # Type II ANOVA
  library(effectsize)   # partial eta-squared
  library(emmeans)      # adjusted means (EMMs)
})

out_dir <- "results/factorial_cv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- site metadata (platform labels printed in the figures) -----------------
PLATFORM <- c(
  BIDMC    = "Akoya PhenoCycler-Fusion",
  Stanford = "RareCyte Orion",
  Roche    = "Akoya PhenoCycler-Fusion"
)

# Staining display labels: EMM panel (R1_5) vs ranking table (R3_fig3_master).
STAIN_EMM <- c("4C_ON" = "4 °C overnight", "25C_1h" = "25 °C 1 h",
               "25C_ON" = "25 °C overnight", "37C_1h" = "37 °C 1 h")
STAIN_TBL <- c("4C_ON" = "4 °C, O/N", "25C_1h" = "25 °C, 1 h",
               "25C_ON" = "25 °C, O/N", "37C_1h" = "37 °C, 1 h")

# Traffic-light example-row colors (rank 1/2/3, and last = worst).
IMG_COLOR <- c("1" = "#8dc63f", "2" = "#eee827", "3" = "#f99d1c")
IMG_LAST  <- "#f26858"

require_input <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path,
         "\n  Run the main Mesmer workflow for the corresponding *_all config first.",
         call. = FALSE)
  }
  path
}

read_site_cv <- function(site) {
  f <- require_input(file.path("results", paste0("out_", site, "_all"),
                               "cv_values_long.csv"))
  read_csv(f, show_col_types = FALSE) %>%
    rename(condition_id = Staining_condition) %>%
    filter(!Is_Excluded, !is.na(CV)) %>%
    mutate(slide = as.integer(gsub(paste0(site, "_"), "", condition_id)))
}

# ---- per-site experimental-design lookups (encode the design; do not change) -
# Model lookups use the full buffer names used for the statistics.
bidmc_lookup <- data.frame(
  slide = 1:24,
  buffer = rep(c("Citraconic Anhydride", "Tris/EDTA"), each = 4, times = 3),
  hier_duration = rep(c("10", "20", "40"), each = 8),
  staining = rep(c("4C_ON", "25C_1h", "25C_ON", "37C_1h"), times = 6),
  stringsAsFactors = FALSE
)
stanford_lookup <- data.frame(
  slide = c(5:8, 13:16, 21:24),
  hier_duration = rep(c("10", "20", "40"), each = 4),
  staining = rep(c("4C_ON", "25C_1h", "25C_ON", "37C_1h"), times = 3),
  stringsAsFactors = FALSE
)
roche_lookup <- data.frame(
  slide = c(2, 6, 10, 14, 18, 22),
  buffer = rep(c("Citrate-based pH6", "EDTA-based pH9"), times = 3),
  hier_duration = rep(c("10", "20", "40"), each = 2),
  stringsAsFactors = FALSE
)

# Table lookups use the figure's display abbreviations for buffer.
tbl_lookup <- list(
  BIDMC = data.frame(
    slide = 1:24,
    hier_buffer = rep(c("CA", "Tris/EDTA"), each = 4, times = 3),
    hier_duration_min = rep(c("10", "20", "40"), each = 8),
    stain_code = rep(c("4C_ON", "25C_1h", "25C_ON", "37C_1h"), times = 6),
    stringsAsFactors = FALSE
  ),
  Stanford = data.frame(
    slide = c(5:8, 13:16, 21:24),
    hier_buffer = "Tris/EDTA",
    hier_duration_min = rep(c("10", "20", "40"), each = 4),
    stain_code = rep(c("4C_ON", "25C_1h", "25C_ON", "37C_1h"), times = 3),
    stringsAsFactors = FALSE
  ),
  Roche = data.frame(
    slide = c(2, 6, 10, 14, 18, 22),
    hier_buffer = rep(c("Bond ER 1 (Citrate-based pH6)",
                        "Bond ER 2 (EDTA-based pH9)"), times = 3),
    hier_duration_min = rep(c("10", "20", "40"), each = 2),
    stain_code = "25C_1h",
    stringsAsFactors = FALSE
  )
)

# ---- helpers to tidy model outputs ------------------------------------------
tidy_eta <- function(anova_obj, recode_map) {
  eta_squared(anova_obj, partial = TRUE) %>%
    as.data.frame() %>%
    filter(!grepl("Marker", Parameter)) %>%
    mutate(factor = dplyr::recode(Parameter, !!!recode_map)) %>%
    transmute(factor, eta2_partial = Eta2_partial,
              ci_low = CI_low, ci_high = CI_high)
}

emm_factor <- function(model, spec, factor_label, level_fun) {
  summary(emmeans(model, spec), infer = TRUE) %>%
    as.data.frame() %>%
    transmute(factor = factor_label,
              level = level_fun(.),
              adjusted_cv = emmean, lower_cl = lower.CL, upper_cl = upper.CL)
}

# =============================================================================
# Fit the three site models
# =============================================================================
effect_frames <- list()
emm_frames    <- list()

## ---- BIDMC ------------------------------------------------------------------
bidmc <- read_site_cv("BIDMC") %>%
  left_join(bidmc_lookup, by = "slide") %>%
  mutate(
    buffer = factor(buffer),
    hier_duration = factor(hier_duration, levels = c("10", "20", "40")),
    staining = factor(staining, levels = c("4C_ON", "25C_1h", "25C_ON", "37C_1h")),
    Marker = factor(Marker)
  )
contrasts(bidmc$buffer)        <- contr.sum(2)
contrasts(bidmc$hier_duration) <- contr.sum(3)
contrasts(bidmc$staining)      <- contr.sum(4)

m_bidmc <- lm(CV ~ buffer + hier_duration + staining + Marker, data = bidmc)
a_bidmc <- Anova(m_bidmc, type = "II")

effect_frames$BIDMC <- tidy_eta(a_bidmc, c(
  "buffer" = "Buffer", "hier_duration" = "HIER duration",
  "staining" = "Staining protocol")) %>%
  mutate(site = "BIDMC", platform = PLATFORM[["BIDMC"]], n_obs = nrow(bidmc))

emm_frames$BIDMC <- bind_rows(
  emm_factor(m_bidmc, ~ buffer, "Buffer", function(d) as.character(d$buffer)),
  emm_factor(m_bidmc, ~ hier_duration, "HIER duration",
             function(d) paste0(d$hier_duration, " min")),
  emm_factor(m_bidmc, ~ staining, "Staining",
             function(d) dplyr::recode(as.character(d$staining), !!!STAIN_EMM))
) %>% mutate(site = "BIDMC", platform = PLATFORM[["BIDMC"]])

## ---- Stanford ---------------------------------------------------------------
stanford <- read_site_cv("Stanford") %>%
  left_join(stanford_lookup, by = "slide") %>%
  mutate(
    hier_duration = factor(hier_duration, levels = c("10", "20", "40")),
    staining = factor(staining, levels = c("4C_ON", "25C_1h", "25C_ON", "37C_1h")),
    Marker = factor(Marker)
  )
contrasts(stanford$hier_duration) <- contr.sum(3)
contrasts(stanford$staining)      <- contr.sum(4)

m_stanford <- lm(CV ~ hier_duration + staining + Marker, data = stanford)
a_stanford <- Anova(m_stanford, type = "II")

effect_frames$Stanford <- tidy_eta(a_stanford, c(
  "hier_duration" = "HIER duration", "staining" = "Staining protocol")) %>%
  mutate(site = "Stanford", platform = PLATFORM[["Stanford"]], n_obs = nrow(stanford))

emm_frames$Stanford <- bind_rows(
  emm_factor(m_stanford, ~ hier_duration, "HIER duration",
             function(d) paste0(d$hier_duration, " min")),
  emm_factor(m_stanford, ~ staining, "Staining",
             function(d) dplyr::recode(as.character(d$staining), !!!STAIN_EMM))
) %>% mutate(site = "Stanford", platform = PLATFORM[["Stanford"]])

## ---- Roche ------------------------------------------------------------------
roche <- read_site_cv("Roche") %>%
  left_join(roche_lookup, by = "slide") %>%
  mutate(
    buffer = factor(buffer),
    hier_duration = factor(hier_duration, levels = c("10", "20", "40")),
    Marker = factor(Marker)
  )
contrasts(roche$buffer)        <- contr.sum(2)
contrasts(roche$hier_duration) <- contr.sum(3)

m_roche <- lm(CV ~ buffer + hier_duration + Marker, data = roche)
a_roche <- Anova(m_roche, type = "II")

effect_frames$Roche <- tidy_eta(a_roche, c(
  "buffer" = "Buffer", "hier_duration" = "HIER duration")) %>%
  mutate(site = "Roche", platform = PLATFORM[["Roche"]], n_obs = nrow(roche))

emm_frames$Roche <- bind_rows(
  emm_factor(m_roche, ~ buffer, "Buffer", function(d) as.character(d$buffer)),
  emm_factor(m_roche, ~ hier_duration, "HIER duration",
             function(d) paste0(d$hier_duration, " min"))
) %>% mutate(site = "Roche", platform = PLATFORM[["Roche"]])

# =============================================================================
# Write effect_size.csv + adjusted_mean_cv.csv  (folds R3_build_viz_input.py)
# =============================================================================
site_order <- c("BIDMC", "Stanford", "Roche")

effect_size <- bind_rows(effect_frames[site_order]) %>%
  select(site, platform, factor, eta2_partial, ci_low, ci_high, n_obs) %>%
  arrange(match(site, site_order), desc(eta2_partial))
write_csv(effect_size, file.path(out_dir, "effect_size.csv"))

adjusted_mean_cv <- bind_rows(emm_frames[site_order]) %>%
  select(site, platform, factor, level, adjusted_cv, lower_cl, upper_cl)
write_csv(adjusted_mean_cv, file.path(out_dir, "adjusted_mean_cv.csv"))

# =============================================================================
# Write ranking_master_table.csv  (folds R3_fig3_master_table.py)
# =============================================================================
fmt_score <- function(v) ifelse(v >= 10, sprintf("%.1f", v), sprintf("%.2f", v))

rank_one_site <- function(site) {
  plat <- PLATFORM[[site]]
  f <- require_input(file.path("results", paste0("out_", site, "_all"),
                               "condition_summary.csv"))
  cs <- read_csv(f, show_col_types = FALSE) %>%
    mutate(slide = as.integer(gsub(paste0(site, "_"), "", Staining_condition))) %>%
    inner_join(tbl_lookup[[site]], by = "slide") %>%
    arrange(Average_Score) %>%                         # ascending: rank 1 = best
    mutate(rank = row_number())

  last <- max(cs$rank)
  cs %>%
    mutate(
      tie = ifelse(duplicated(round(Average_Score, 2)) |
                     duplicated(round(Average_Score, 2), fromLast = TRUE), "*", ""),
      is_example = rank %in% c(1, 2, 3, last),
      image_color = dplyr::case_when(
        rank == last ~ IMG_LAST,
        as.character(rank) %in% names(IMG_COLOR) ~ IMG_COLOR[as.character(rank)],
        TRUE ~ ""
      ),
      ab_staining = dplyr::recode(stain_code, !!!STAIN_TBL),
      hier_duration = paste0(hier_duration_min, " min"),
      avg_score = fmt_score(Average_Score),
      avg_score_exact = Average_Score,
      site = .env$site, platform = plat
    ) %>%
    select(site, platform, rank, tie, hier_duration, hier_buffer, ab_staining,
           avg_score, avg_score_exact, is_example, image_color)
}

ranking_master <- bind_rows(lapply(site_order, rank_one_site))
write_csv(ranking_master, file.path(out_dir, "ranking_master_table.csv"))

# =============================================================================
# Report
# =============================================================================
cat("\nWrote to", out_dir, ":\n")
for (f in c("effect_size.csv", "adjusted_mean_cv.csv", "ranking_master_table.csv")) {
  cat(sprintf("  %-26s %d rows\n", f,
              nrow(read_csv(file.path(out_dir, f), show_col_types = FALSE))))
}
cat("\n--- effect_size.csv (partial eta^2) ---\n")
print(as.data.frame(effect_size), row.names = FALSE)
cat("\nDone.\n")
