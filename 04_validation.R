# 04_validation.R
# Out-of-sample check of the movement model.
#
# Fit the density model on 2023-2024 only, then predict each pitcher's 2025 Coors
# movement for each pitch type from his 2025 movement everywhere else. Compare
# against what he actually did at Coors in 2025.
#
# Three predictions per pitcher x pitch type:
#   baseline       no change: Coors movement = his non-Coors 2025 average
#   density        non-Coors average scaled by the trained rho elasticity
#   density_coors  density, plus the trained Coors-only leftover term
#
# Each prediction uses only his non-Coors 2025 pitches plus the air density of
# the Coors games he pitched in (known at game time from temperature).
#
# Noise floor: the observed Coors average is itself an estimate from 15+ pitches,
# so even a perfect model would show some error. The floor is the MAE you'd
# expect from that sampling noise alone.
#
# Inputs:  data/clean/pitches_air.parquet  (from 02)
# Outputs: output/validation_by_cell.csv
#          output/validation_summary.csv
#          output/validation_error_by_pitch_type.png

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(fixest)
  library(readr)
  library(ggplot2)
})

dir.create("output", showWarnings = FALSE)

train_years       <- 2023:2024
test_year         <- 2025
min_cell_pitches  <- 50   # same cell filter as 03
min_else_pitches  <- 30   # non-Coors pitches needed to set a 2025 starting point
min_coors_pitches <- 15   # Coors pitches needed for an observed 2025 Coors average
min_cells_chart   <- 10   # pitch types with fewer test cells stay out of the chart

pitch_groups <- tribble(
  ~pitch_type, ~pitch_group,
  "FF", "fastball", "SI", "fastball", "FC", "fastball",
  "SL", "breaking", "ST", "breaking", "SV", "breaking",
  "CU", "breaking", "KC", "breaking", "CS", "breaking",
  "CH", "offspeed", "FS", "offspeed", "FO", "offspeed", "SC", "offspeed"
)

# ---- Build frame (same definitions as 03) ------------------------------------

pitches <- read_parquet("data/clean/pitches_air.parquet")

# 01 stores the year as `season`; fall back to game_date if that ever changes
if (!"game_year" %in% names(pitches)) {
  if ("season" %in% names(pitches)) {
    pitches <- pitches |> mutate(game_year = as.integer(season))
  } else if ("game_date" %in% names(pitches)) {
    pitches <- pitches |> mutate(game_year = as.integer(format(as.Date(game_date), "%Y")))
  } else {
    stop("pitches_air.parquet has no season, game_year, or game_date column")
  }
}

# Regular season plus postseason, same as 03
if ("game_type" %in% names(pitches)) {
  pitches <- pitches |> filter(game_type %in% c("R", "F", "D", "L", "W"))
}

d <- pitches |>
  inner_join(pitch_groups, by = "pitch_type") |>
  filter(!is.na(pfx_x), !is.na(pfx_z), !is.na(release_speed),
         !is.na(release_spin_rate), !is.na(spin_axis), !is.na(rho_ratio)) |>
  mutate(
    ivb_in     = pfx_z * 12,
    hb_arm_in  = if_else(p_throws == "R", -pfx_x, pfx_x) * 12,
    brk_in     = sqrt(ivb_in^2 + hb_arm_in^2),
    spin_k     = release_spin_rate / 1000,
    sin_ax     = sin(spin_axis * pi / 180),
    cos_ax     = cos(spin_axis * pi / 180),
    at_coors_i = as.integer(at_coors),
    cell       = paste(pitcher, pitch_type, game_year, sep = "_")
  ) |>
  group_by(cell) |>
  mutate(
    n_cell  = n(),
    ivb_dir = ivb_in * sign(mean(ivb_in)),
    hb_dir  = hb_arm_in * sign(mean(hb_arm_in))
  ) |>
  ungroup() |>
  filter(n_cell >= min_cell_pitches)

train <- d |> filter(game_year %in% train_years)
test  <- d |> filter(game_year == test_year)

stopifnot(nrow(train) > 0, nrow(test) > 0, any(test$at_coors_i == 1))

cat(sprintf("Train %s: %s pitches, %s cells\n",
            paste(range(train_years), collapse = "-"),
            format(nrow(train), big.mark = ","),
            format(n_distinct(train$cell), big.mark = ",")))
cat(sprintf("Test %d: %s pitches, %s cells\n",
            test_year,
            format(nrow(test), big.mark = ","),
            format(n_distinct(test$cell), big.mark = ",")))

# ---- Fit on 2023-2024 ---------------------------------------------------------

controls <- "release_speed + spin_k + sin_ax + cos_ax"

# Point estimates only, so no clustered SEs here (03 has those).
fit_group <- function(dat, y, with_coors) {
  rhs <- paste("rho_ratio +", controls)
  if (with_coors) rhs <- paste(rhs, "+ at_coors_i")
  m <- feols(as.formula(paste(y, "~", rhs, "| cell")),
             data = dat, vcov = "iid", notes = FALSE, warn = FALSE)
  b <- coef(m)
  tibble(
    el          = b[["rho_ratio"]] * mean(dat$rho_ratio) / mean(dat[[y]]),
    coors_resid = if (with_coors) b[["at_coors_i"]] else 0
  )
}

coefs <- expand_grid(
  pitch_group = sort(unique(train$pitch_group)),
  outcome     = c("brk_in", "ivb_dir", "hb_dir"),
  model       = c("density", "density_coors")
) |>
  mutate(fit = pmap(list(pitch_group, outcome, model), function(g, y, m) {
    fit_group(train[train$pitch_group == g, ], y, m == "density_coors")
  })) |>
  unnest(fit)

coef_wide <- coefs |>
  pivot_wider(id_cols = c(pitch_group, outcome),
              names_from = model, values_from = c(el, coors_resid))

cat("\nTrained total-break terms (compare to 03's full-sample fit):\n")
coef_wide |>
  filter(outcome == "brk_in") |>
  transmute(pitch_group,
            el_density       = round(el_density, 2),
            el_with_coors    = round(el_density_coors, 2),
            coors_resid_in   = round(coors_resid_density_coors, 2)) |>
  print()

# ---- Predict 2025 Coors movement per pitcher x pitch type --------------------

cells <- test |>
  group_by(pitcher, pitch_type, pitch_group, p_throws) |>
  summarise(
    n_coors   = sum(at_coors_i == 1),
    n_else    = sum(at_coors_i == 0),
    rho_coors = mean(rho_ratio[at_coors_i == 1]),
    rho_else  = mean(rho_ratio[at_coors_i == 0]),
    across(c(brk_in, ivb_in, hb_arm_in),
           list(other     = ~ mean(.x[at_coors_i == 0]),
                coors     = ~ mean(.x[at_coors_i == 1]),
                coors_var = ~ var(.x[at_coors_i == 1]))),
    .groups = "drop"
  ) |>
  filter(n_coors >= min_coors_pitches, n_else >= min_else_pitches)

val <- cells |>
  pivot_longer(matches("^(brk_in|ivb_in|hb_arm_in)_"),
               names_to = c("measure", ".value"),
               names_pattern = "^(brk_in|ivb_in|hb_arm_in)_(other|coors|coors_var)$") |>
  mutate(coef_outcome = recode(measure, brk_in = "brk_in",
                               ivb_in = "ivb_dir", hb_arm_in = "hb_dir")) |>
  left_join(coef_wide, by = c("pitch_group", "coef_outcome" = "outcome")) |>
  mutate(
    d_rho = rho_coors / rho_else - 1,
    # Proportional scaling works on signed components too: it shrinks them toward
    # zero. The Coors leftover is in "normal break direction" units, so it gets
    # the sign of the pitch's movement (always +1 for total break).
    pred_baseline      = other,
    pred_density       = other * (1 + el_density * d_rho),
    pred_density_coors = other * (1 + el_density_coors * d_rho) +
                         coors_resid_density_coors * sign(other),
    err_baseline       = pred_baseline - coors,
    err_density        = pred_density - coors,
    err_density_coors  = pred_density_coors - coors,
    se_obs             = sqrt(coors_var / n_coors)
  )

write_csv(val, "output/validation_by_cell.csv")

# ---- Summarise errors ---------------------------------------------------------

summ <- function(x) {
  x |>
    summarise(
      cells        = n(),
      mae_base     = mean(abs(err_baseline)),
      mae_density  = mean(abs(err_density)),
      mae_full     = mean(abs(err_density_coors)),
      rmse_base    = sqrt(mean(err_baseline^2)),
      rmse_density = sqrt(mean(err_density^2)),
      rmse_full    = sqrt(mean(err_density_coors^2)),
      bias_base    = mean(err_baseline),
      bias_full    = mean(err_density_coors),
      noise_mae    = mean(sqrt(2 / pi) * se_obs),
      full_wins    = mean(abs(err_density_coors) < abs(err_baseline)),
      .groups = "drop"
    ) |>
    mutate(pct_better = 100 * (1 - mae_full / mae_base))
}

type_levels <- c(pitch_groups$pitch_type, "All")

summary_tbl <- bind_rows(
  val |> group_by(measure, pitch_group, pitch_type) |> summ(),
  val |> group_by(measure) |> summ() |> mutate(pitch_group = "all", pitch_type = "All")
) |>
  mutate(pitch_type = factor(pitch_type, levels = type_levels)) |>
  arrange(measure, pitch_type)

write_csv(summary_tbl, "output/validation_summary.csv")

measure_labels <- c(brk_in    = "Total break",
                    ivb_in    = "Induced vertical break",
                    hb_arm_in = "Horizontal break (arm side +)")

cat(sprintf(
  "\n%d pitcher x pitch type cells in %d with >= %d Coors and >= %d other pitches.",
  n_distinct(paste(val$pitcher, val$pitch_type)), test_year,
  min_coors_pitches, min_else_pitches))
cat("\nMAE in inches. bias = predicted minus observed. floor = MAE from sampling noise alone.")
cat("\npct_better and wins compare density_coors to baseline.\n")

for (m in names(measure_labels)) {
  cat("\n", measure_labels[[m]], "\n", sep = "")
  summary_tbl |>
    filter(measure == m) |>
    transmute(pitch_type, cells,
              base       = round(mae_base, 2),
              density    = round(mae_density, 2),
              dens_coors = round(mae_full, 2),
              floor      = round(noise_mae, 2),
              bias_base  = round(bias_base, 2),
              bias_full  = round(bias_full, 2),
              pct_better = round(pct_better),
              wins       = sprintf("%.0f%%", 100 * full_wins)) |>
    print(n = Inf)
}

# ---- Chart --------------------------------------------------------------------

plot_df <- summary_tbl |>
  filter(cells >= min_cells_chart | pitch_type == "All") |>
  mutate(measure = factor(measure_labels[measure], levels = measure_labels))

bars <- plot_df |>
  select(measure, pitch_type, mae_base, mae_density, mae_full) |>
  pivot_longer(starts_with("mae_"), names_to = "model", values_to = "mae") |>
  mutate(model = factor(model,
                        levels = c("mae_base", "mae_density", "mae_full"),
                        labels = c("No change (baseline)", "Density",
                                   "Density + Coors term")))

p <- ggplot(bars, aes(pitch_type, mae, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(data = plot_df,
                aes(x = pitch_type, ymin = noise_mae, ymax = noise_mae),
                inherit.aes = FALSE, width = 0.8, linewidth = 0.6) +
  facet_wrap(~ measure, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("grey70", "#8fb3d9", "#33608c")) +
  labs(
    title    = sprintf("Predicting %d Coors movement by pitch type", test_year),
    subtitle = sprintf(paste0("Model fit on %d-%d. Baseline assumes Coors movement ",
                              "equals movement elsewhere.\nBlack line: error expected ",
                              "from sampling noise alone."),
                       min(train_years), max(train_years)),
    x = NULL, y = "Mean absolute error (in)", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

ggsave("output/validation_error_by_pitch_type.png", p,
       width = 9, height = 8, dpi = 200, bg = "white")

cat("\nSaved output/validation_by_cell.csv, output/validation_summary.csv,",
    "output/validation_error_by_pitch_type.png\n")
