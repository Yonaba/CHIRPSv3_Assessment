root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 6b - Comparison of ETCCDI indices: Observations vs CHIRPS v2, CHIRPS v3 (RNL)
# and CHIRPS v3 (SAT).   (fast part: plots the output of 6a_etccdi_compute.R)
# -----------------------------------------------------------------------------
# Input : tables/4_etccdi_indices.csv
# Method:
#   1. Station-level value = mean over the years valid for all 4 series (at least
#      min_years years); p90 / p95 are already one value per station.
#   2. ggstatsplot::ggwithinstats (type = "robust": trimmed-means repeated-measures
#      ANOVA), as in 2b_boxplot_evaluation_plot.R. The station is the repeated
#      measure, so only stations with a valid value for all 4 series are kept.
#      Pairwise p-values are not displayed.
# Output: figures/6_etccdi_comparison.png
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "ggstatsplot", "ggpubr", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Settings -------------------------------------------------------------
in_file   <- "tables/4_etccdi_indices.csv"
fig_file  <- "figures/6_etccdi_comparison.png"

src_levels <- c("Observations", "CHIRPS v2", "CHIRPS v3 (RNL)", "CHIRPS v3 (SAT)")
min_years         <- 10  # min. number of valid years for a station-level mean
min_stations_plot <- 5   # min. number of stations to draw a panel

ind_order <- c("nP", "CDD", "CWD", "SDII", "R10mm", "R20mm",
               "RX1DAY", "RX5DAY", "p90", "p95", "R95pTOT", "R99pTOT")
ind_unit  <- c(nP = "days", CDD = "days", CWD = "days", SDII = "mm/day",
               R10mm = "days", R20mm = "days", RX1DAY = "mm", RX5DAY = "mm",
               p90 = "mm", p95 = "mm", R95pTOT = "mm", R99pTOT = "mm")

# Colour convention (observations in black)
source_cols <- c("Observations"    = "black",
                 "CHIRPS v2"       = "darkred",
                 "CHIRPS v3 (RNL)" = "darkblue",
                 "CHIRPS v3 (SAT)" = "darkgreen")

dir.create("figures", showWarnings = FALSE)

# ---- 2. Station-level values -----------------------------------------------------
res <- read_csv(in_file, show_col_types = FALSE)
if (!all(src_levels %in% res$Source))
  stop("Missing series in ", in_file, ": ",
       paste(setdiff(src_levels, unique(res$Source)), collapse = ", "))

# Annual indices: mean over the years valid for ALL 4 series (at least min_years)
st_ann <- res |>
  filter(!is.na(year)) |>
  pivot_wider(names_from = Source, values_from = value) |>
  filter(if_all(all_of(src_levels), is.finite)) |>
  group_by(Station, ind) |>
  summarise(n_years = n(), across(all_of(src_levels), mean), .groups = "drop") |>
  filter(n_years >= min_years) |>
  dplyr::select(-n_years)

# p90 / p95: already one value per station
st_static <- res |>
  filter(is.na(year)) |>
  dplyr::select(-year) |>
  pivot_wider(names_from = Source, values_from = value) |>
  filter(if_all(all_of(src_levels), is.finite))

idx_station <- bind_rows(st_ann, st_static) |>
  pivot_longer(all_of(src_levels), names_to = "Source", values_to = "value") |>
  mutate(Source = factor(Source, levels = src_levels))

n_st <- table(idx_station$ind) / length(src_levels)
message("Stations per index: ", paste(ind_order, n_st[ind_order], sep = "=", collapse = ", "))

# ---- 3. Figure -------------------------------------------------------------------
empty_panel <- function(ttl, msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 5) +
    labs(title = ttl) + theme_void() +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
}

within_panel <- function(d, ylab, ttl, show_x = TRUE) {
  # Repeated-measures design: only stations with a valid value for ALL groups,
  # sorted by group then station so that rows are aligned across groups
  d <- d |>
    filter(is.finite(value)) |>
    group_by(Station) |>
    filter(n() == length(src_levels)) |>
    ungroup() |>
    arrange(Source, Station) |>
    droplevels()
  
  if (nrow(d) / length(src_levels) < min_stations_plot)
    return(empty_panel(ttl, "Not enough stations"))
  
  # x-axis (product names) only on the bottom-most panels
  hide_x <- if (show_x) theme() else
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  tryCatch(
    ggwithinstats(
      data = d, x = Source, y = value,
      type = "robust",                  # trimmed means (robust repeated measures)
      pairwise.display = "none",        # no pairwise p-values
      bf.message = FALSE,
      xlab = NULL,
      ylab = ylab,
      title = ttl
    ) +
      scale_color_manual(values = source_cols) +
      scale_fill_manual(values = source_cols) +
      theme(plot.title    = element_text(colour = "black", size = 18, face = "bold", hjust = 0),
            plot.subtitle = element_text(colour = "black", size = 12, face = "italic", hjust = 0),
            axis.text     = element_text(colour = "black", size = 13),
            axis.title    = element_text(colour = "black", size = 13),
            axis.title.x  = element_blank()) +
      hide_x,
    error = function(e) {
      message("Panel ", ttl, " skipped: ", conditionMessage(e))
      empty_panel(ttl, "Test could not be computed")
    }
  )
}

n_col <- 3
n_row <- 4

plots <- lapply(seq_along(ind_order), function(k) {
  ii <- ind_order[k]
  within_panel(idx_station[idx_station$ind == ii, ],
               paste0(ii, " [", ind_unit[[ii]], "]"),
               paste0("(", letters[k], ") ", ii),
               show_x = (k > (n_row - 1) * n_col))  # last row only
})

fig <- ggarrange(plotlist = plots, ncol = n_col, nrow = n_row)

ggsave(fig_file, fig, width = 24, height = 26, units = "in", dpi = 300, bg = "white")
message("Done: ", fig_file)