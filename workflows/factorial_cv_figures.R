#!/usr/bin/env Rscript
# =============================================================================
# Factorial CV reanalysis - Stage 2: figures
#
# Reads the tidy tables from Stage 1 (results/factorial_cv/*.csv) and renders the
# three promoted Figure-3 deliverables, one per site where applicable:
#   - adjusted_mean_cv_{site}   Panel C: EMM dot + 95% CI per factor level
#                               (left factor strips, alternating bands,
#                                grand-mean dashed line, best level haloed)
#   - effect_size_bubble        Panel D: cross-site partial-eta^2 bubble matrix
#   - ranking_table_{site}      funkyheatmap "Conditions Ranked Top->Bottom" table
#
# Consolidated R-only port of A_importance_matrix.R, the build_emm panel of
# R3_stats_v3_funky_consistent.R, and R3_concept1_funky.R.
#
# Run from the repo root (after factorial_cv_model.R):
#   Rscript workflows/factorial_cv_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(forcats); library(ggplot2); library(ggtext); library(scales)
  library(svglite); library(cowplot); library(tibble); library(funkyheatmap)
})

io_dir <- "results/factorial_cv"
stopifnot(dir.exists(io_dir))

effect <- read_csv(file.path(io_dir, "effect_size.csv"), show_col_types = FALSE)
emm    <- read_csv(file.path(io_dir, "adjusted_mean_cv.csv"), show_col_types = FALSE)
master <- read_csv(file.path(io_dir, "ranking_master_table.csv"), show_col_types = FALSE)

SITE_COL   <- c(BIDMC = "#B52427", Stanford = "#FFC130", Roche = "#1E94CF")
SITE_TITLE <- c(BIDMC = "BIDMC-Harvard", Stanford = "Stanford", Roche = "Roche")
sites      <- c("BIDMC", "Stanford", "Roche")

save_fig <- function(plot, stem, w, h) {
  ggsave(paste0(stem, ".svg"), plot, width = w, height = h,
         device = svglite::svglite, fix_text_size = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = w, height = h, device = grDevices::cairo_pdf)
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
}

# =============================================================================
# PANEL C - Adjusted mean CV (EMM dot + 95% CI) per factor, per site
# Factor rows: HIER duration / Ab staining / HIER buffer; best level haloed;
# grand-mean dashed reference; alternating grey factor bands; left strip labels.
# =============================================================================
FACTOR_LABEL <- c("HIER duration" = "HIER duration",
                  "Staining" = "Ab staining", "Buffer" = "HIER buffer")
FACTOR_ORDER <- c("HIER buffer", "Ab staining", "HIER duration")  # top -> bottom

build_panel_c <- function(site) {
  col <- SITE_COL[[site]]
  d <- emm %>%
    filter(site == !!site) %>%
    mutate(factor = dplyr::recode(factor, !!!FACTOR_LABEL),
           factor = factor(factor, levels = intersect(FACTOR_ORDER, unique(factor))),
           level = str_replace(level, "overnight", "O/N"))  # match figure abbreviation

  grand <- mean(d$adjusted_cv)

  # within each factor: lowest CV (best) at the TOP -> y levels = descending CV
  # circle (halo) = any level whose adjusted CV is below the overall (grand) mean
  d <- d %>%
    mutate(is_below_mean = adjusted_cv < grand) %>%
    arrange(factor, desc(adjusted_cv)) %>%
    mutate(level = fct_inorder(level))

  # alternating background bands (1st & 3rd factor grey, 2nd white)
  bands <- tibble(factor = levels(d$factor)) %>%
    mutate(factor = factor(factor, levels = levels(d$factor)),
           fill = ifelse(seq_len(n()) %% 2 == 1, "#F2F2F2", "#FFFFFF"))

  ggplot(d, aes(x = adjusted_cv, y = level)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    scale_fill_identity() +
    geom_vline(xintercept = grand, linetype = "dashed",
               colour = "#7A7A7A", linewidth = 0.35) +
    geom_linerange(aes(xmin = lower_cl, xmax = upper_cl),
                   colour = col, linewidth = 0.7) +
    # below-mean levels: open halo ring behind the dot
    geom_point(data = dplyr::filter(d, is_below_mean),
               shape = 21, size = 4.4, fill = NA, colour = col, stroke = 0.9) +
    geom_point(size = 2.4, shape = 21, fill = col, colour = "white", stroke = 0.4) +
    facet_grid(factor ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(labels = label_number(accuracy = 0.01)) +
    labs(title = str_wrap(paste0(
           "Adjusted Mean CV with Lower and Upper Confidence Interval for Each Factor at ",
           SITE_TITLE[[site]]), 50),
         subtitle = paste0("<span style='color:", col,
                           "'>&#9673;</span> indicates lower than overall average CV (dotted line)"),
         x = "Adjusted mean CV (lower = better)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      strip.placement  = "outside",
      strip.background  = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
      panel.spacing.y   = unit(2, "pt"),
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      plot.title        = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle     = ggtext::element_markdown(hjust = 0.5, size = 9,
                                                   margin = margin(b = 4)),
      plot.title.position = "plot"
    )
}

for (s in sites) {
  message("Panel C: ", s)
  p <- build_panel_c(s)
  n <- nrow(filter(emm, site == s))
  save_fig(p, file.path(io_dir, paste0("adjusted_mean_cv_", s)),
           w = 6.0, h = 1.2 + 0.42 * n)
}

# =============================================================================
# PANEL D - "Which prep knob matters" importance matrix (port of A_importance_matrix.R)
# rows = HIER buffer / Ab staining / HIER duration; cols = 3 sites;
# dot AREA = partial eta^2, colored by site; absent cells render "n/a".
# =============================================================================
build_panel_d <- function() {
  factor_levels <- c("Buffer", "Staining", "HIER duration")
  factor_labs   <- c(Buffer = "HIER buffer", Staining = "Ab staining",
                     `HIER duration` = "HIER duration")
  site_labs <- c(
    BIDMC = "BIDMC<br><span style='font-size:6.5pt;color:#666666'>(Harvard)</span>",
    Stanford = "Stanford", Roche = "Roche")

  eff <- effect %>%
    mutate(factor = ifelse(factor == "Staining protocol", "Staining", factor)) %>%
    filter(factor %in% factor_levels, site %in% sites)

  grid <- expand.grid(factor = factor_levels, site = sites, stringsAsFactors = FALSE)
  dat <- grid %>%
    left_join(eff %>% select(site, factor, eta2_partial), by = c("factor", "site")) %>%
    mutate(site = factor(site, levels = sites),
           factor = factor(factor, levels = rev(factor_levels)),
           present = !is.na(eta2_partial))
  present_dat <- dat %>% filter(present)
  absent_dat  <- dat %>% filter(!present)

  max_eta <- max(present_dat$eta2_partial)
  size_breaks <- c(0.02, 0.05, 0.10, 0.14); size_breaks <- size_breaks[size_breaks <= max_eta * 1.05]
  present_dat <- present_dat %>%
    mutate(lab_size = scales::rescale(sqrt(eta2_partial), to = c(2.1, 3.0),
                                      from = c(0, sqrt(max_eta))))

  th <- theme_minimal(base_size = 8, base_family = "sans") +
    theme(
      panel.grid = element_blank(), axis.title = element_blank(),
      axis.text.x = element_markdown(size = 9, face = "bold", margin = margin(b = 2)),
      axis.text.y = element_text(size = 9, face = "bold", hjust = 1),
      axis.text.x.top = element_markdown(size = 9, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_markdown(size = 8, colour = "#444444", hjust = 0,
                                       margin = margin(b = 10), lineheight = 1.25),
      plot.caption = element_textbox_simple(size = 6.5, colour = "#666666", hjust = 0,
                                            halign = 0, margin = margin(t = 12), lineheight = 1.35),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = element_text(size = 7.5), legend.text = element_text(size = 7),
      legend.margin = margin(2, 6, 2, 6), plot.margin = margin(14, 18, 10, 14))

  ggplot() +
    geom_tile(data = dat, aes(x = site, y = factor), fill = NA, colour = "grey92",
              linewidth = 0.3, width = 0.96, height = 0.96) +
    geom_text(data = absent_dat, aes(x = site, y = factor), label = "n/a",
              colour = "grey75", size = 2.7, fontface = "italic") +
    geom_point(data = present_dat, aes(x = site, y = factor, size = eta2_partial, fill = site),
               shape = 21, colour = "white", stroke = 0.6, alpha = 0.95) +
    geom_text(data = present_dat, aes(x = site, y = factor, label = sprintf("%.3f", eta2_partial)),
              size = present_dat$lab_size, vjust = 0.5, fontface = "bold",
              colour = ifelse(present_dat$site == "Stanford", "#5a4500", "white")) +
    scale_x_discrete(position = "top", labels = site_labs, expand = expansion(add = 0.55)) +
    scale_y_discrete(labels = factor_labs, expand = expansion(add = 0.55)) +
    scale_fill_manual(values = SITE_COL, guide = "none") +
    scale_size_area(name = expression(paste("Partial ", eta^2, " (effect size)")),
                    max_size = 23, breaks = size_breaks, labels = sprintf("%.2f", size_breaks),
                    limits = c(0, max_eta),
                    guide = guide_legend(override.aes = list(fill = "grey55", colour = "white",
                                                             stroke = 0.5),
                                         title.position = "top", nrow = 1, order = 1,
                                         label.position = "bottom")) +
    labs(
      title = "Which prep knob matters, and where?",
      subtitle = paste0("Dot <b>area</b> = partial &eta;<sup>2</sup> (share of CV variance ",
                        "explained by each factor), colored by site.<br>",
                        "<span style='color:#B52427'>&#9679; BIDMC</span> &nbsp; ",
                        "<span style='color:#FFC130'>&#9679; Stanford</span> &nbsp; ",
                        "<span style='color:#1E94CF'>&#9679; Roche</span>"),
      caption = paste0("Buffer drives BIDMC &amp; Roche; staining dominates Stanford; HIER ",
                       "duration is minor except at Roche. Blank cells (<i>n/a</i>) = factor ",
                       "not tested at that site.<br>Partial &eta;<sup>2</sup> from per-site ",
                       "factorial models of coefficient of variation (CV).")) +
    th + coord_cartesian(clip = "off")
}

message("Panel D: importance matrix")
save_fig(build_panel_d(), file.path(io_dir, "effect_size_bubble"), w = 6.7, h = 5.5)

# =============================================================================
# FUNKY RANKING TABLE  (port of R3_concept1_funky.R)
# =============================================================================
SITE_TBL_TITLE <- c(
  BIDMC    = "Conditions Ranked from Top to Bottom for BIDMC-Harvard/Akoya PhenoCycler-Fusion",
  Stanford = "Conditions Ranked from Top to Bottom for Stanford/RareCyte Orion",
  Roche    = "Conditions Ranked from Top to Bottom for Roche/Akoya PhenoCycler-Fusion")

W_RANK <- 1.4; W_DURATION <- 3.2; W_STAINING <- 4.0; W_BAR <- 4.0; W_LABEL <- 1.4
POS_ARGS <- position_arguments(expand_xmax = 4, col_annot_offset = 4.2)

# Shared "Performance" gradient (matches the funkyheatmap bar fill) + a standalone
# horizontal colorbar legend keyed to the ACTUAL Avg Score range for each site.
PERF_COLORS <- c("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C")

make_perf_legend <- function(site) {
  rng <- range(master$avg_score_exact[master$site == site])
  p <- ggplot(tibble(v = rng), aes(x = v, y = 1, fill = v)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = PERF_COLORS, limits = rng,
      breaks = rng, labels = sprintf("%.1f", rng),
      name = "Avg Score (lower = better)",
      guide = guide_colourbar(title.position = "top", title.hjust = 0, ticks = FALSE,
                              barwidth = grid::unit(36, "pt"),
                              barheight = grid::unit(6, "pt"))) +
    theme_void(base_size = 8) +
    theme(legend.position = "bottom",
          legend.justification = "left",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text  = element_text(size = 7))
  cowplot::get_legend(p)
}

build_funky <- function(site) {
  is_bidmc <- site == "BIDMC"
  d <- master %>% filter(site == !!site) %>% arrange(rank)
  bmax <- max(d$avg_score_exact)
  d <- d %>% mutate(
    tie = ifelse(is.na(tie), "", tie),
    performance = avg_score_exact / bmax,
    avg_label = paste0("  ", avg_score),
    rank_lab = paste0(rank, tie),
    example_val = dplyr::case_when(rank == 1 ~ 1L, rank == 2 ~ 2L, rank == 3 ~ 3L,
                                   rank == max(rank) ~ 4L, TRUE ~ NA_integer_))
  buf_w <- max(nchar(d$hier_buffer)) * 0.34 + 0.6

  data <- d %>% transmute(id = rank_lab, hier_duration, hier_buffer, ab_staining,
                          performance, avg_label)

  column_info <- tribble(
    ~id,            ~name,           ~geom, ~group,      ~palette, ~width,     ~hjust,
    "id",           "Rank",          "text","condition", NA,       W_RANK,     0,
    "hier_duration","HIER duration", "text","condition", NA,       W_DURATION, 0,
    "hier_buffer",  "HIER buffer",   "text","condition", NA,       buf_w,      0,
    "ab_staining",  "Ab staining",   "text","condition", NA,       W_STAINING, 0,
    "performance",  "Avg Score",     "bar", "score",     "perf",   W_BAR,      NA,
    "avg_label",    "",              "text","score",     NA,       W_LABEL,    0)
  column_info <- column_info %>% mutate(options = lapply(seq_len(n()), function(i) {
    o <- list(width = width[i]); if (!is.na(hjust[i])) o$hjust <- hjust[i]; o
  })) %>% select(id, name, geom, group, palette, options)

  column_groups <- tribble(
    ~group,      ~palette, ~level1,
    "condition", "perf",   "Condition (HIER + antibody staining)",
    "score",     "perf",   "Performance")
  palettes <- list(perf = PERF_COLORS)

  g <- funky_heatmap(data = data, column_info = column_info, column_groups = column_groups,
                     palettes = palettes,
                     legends = list(list(palette = "perf", enabled = FALSE)),
                     position_args = POS_ARGS, add_abc = FALSE, scale_column = FALSE)
  if (!is_bidmc) return(g)

  # BIDMC: overlay the example-row color fills at the correct cells
  vdata <- funkyheatmap:::verify_data(data)
  vci   <- funkyheatmap:::verify_column_info(column_info, vdata)
  vri   <- funkyheatmap:::verify_row_info(NULL, vdata)
  vcg   <- funkyheatmap:::verify_column_groups(column_groups, vci)
  vrg   <- funkyheatmap:::verify_row_groups(NULL, vri)
  vpal  <- funkyheatmap:::verify_palettes(palettes, vci, data)
  gp <- funkyheatmap:::calculate_geom_positions(vdata, vci, vri, vcg, vrg, vpal,
                                                POS_ARGS, FALSE, FALSE)
  row_pos <- gp$row_pos %>% mutate(rank = d$rank, example_val = d$example_val,
                                   image_color = d$image_color)
  col_pos <- gp$column_pos
  zebra <- tryCatch(ggplot2::layer_data(g, 1L), error = function(e) NULL)
  if (!is.null(zebra) && all(c("xmin", "xmax") %in% names(zebra))) {
    data_xmin <- min(zebra$xmin, na.rm = TRUE); data_xmax <- max(zebra$xmax, na.rm = TRUE)
  } else {
    data_xmin <- min(col_pos$xmin); data_xmax <- max(col_pos$xmax)
  }
  ex_rows <- row_pos %>% filter(!is.na(example_val))
  row_fill <- ex_rows %>% transmute(xmin = data_xmin, xmax = data_xmax,
                                    ymin = ymin, ymax = ymax, fill = image_color)
  row_fill_layer <- ggplot2::geom_rect(
    data = row_fill, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = row_fill$fill, alpha = 1, inherit.aes = FALSE)
  g$layers <- append(g$layers, list(row_fill_layer), after = 1L)
  g
}

for (s in sites) {
  message("Funky table: ", s)
  g <- build_funky(s) + ggplot2::theme(plot.margin = ggplot2::margin(2, 16, 2, 16))
  hm_h <- g$height; w <- g$width
  title <- cowplot::ggdraw() +
    cowplot::draw_label(SITE_TBL_TITLE[[s]], fontface = "bold", size = 9.5, x = 0.018, hjust = 0)
  lg <- cowplot::ggdraw() +
    cowplot::draw_grob(make_perf_legend(s), x = 0.035, width = 0.93)
  th <- 0.45; lh <- 0.4; h <- hm_h + th + lh
  final <- cowplot::plot_grid(title, g, lg, ncol = 1, rel_heights = c(th, hm_h, lh))
  stem <- file.path(io_dir, paste0("ranking_table_", s))
  ggsave(paste0(stem, ".pdf"), final, width = w, height = h, device = grDevices::cairo_pdf)
  ggsave(paste0(stem, ".svg"), final, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(stem, ".png"), final, width = w, height = h, dpi = 600, bg = "white")
}

message("DONE")
