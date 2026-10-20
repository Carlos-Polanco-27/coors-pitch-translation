# 03_movement_model.R
# How much does air density change pitch movement, holding the pitcher and pitch fixed?
#
# Physics prior: Magnus and seam-shifted-wake forces both scale with air density,
# and so does drag. So break and velocity loss should scale roughly 1:1 with
# rho_ratio (elasticity near 1). Spin rate at release should not depend on rho.
#
# Identification: pitcher x pitch type x season fixed effects, so the rho effect
# comes from the same pitch thrown in different air (Coors vs elsewhere, plus
# temperature and elevation variation at every other park, plus Mexico City).
# Spec "rho_plus_coors" adds a Coors indicator. If density explains the Coors
# effect, that indicator should be near zero. A sizable leftover points at
# something Coors-specific (Hawk-Eye calibration, humidor, pitcher behavior).
#
# Inputs:  data/clean/pitches_air.parquet  (from 02)
# Outputs: output/movement_coefs.csv
#          output/movement_by_pitch_cell.csv

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(fixest)
  library(readr)
})

dir.create("output", showWarnings = FALSE)

min_cell_pitches  <- 50   # pitches per pitcher x pitch type x season to keep the cell
min_coors_pitches <- 15   # Coors pitches needed for an observed per-cell Coors delta

pitch_groups <- tribble(
  ~pitch_type, ~pitch_group,
  "FF", "fastball", "SI", "fastball", "FC", "fastball",
  "SL", "breaking", "ST", "breaking", "SV", "breaking",
  "CU", "breaking", "KC", "breaking", "CS", "breaking",
  "CH", "offspeed", "FS", "offspeed", "FO", "offspeed", "SC", "offspeed"
)

# ---- Load -------------------------------------------------------------------

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

# Regular season plus postseason (wild card, division, LCS, World Series)
if ("game_type" %in% names(pitches)) {
  pitches <- pitches |> filter(game_type %in% c("R", "F", "D", "L", "W"))
}

need <- c("pitcher", "game_pk", "game_year", "pitch_type", "p_throws",
          "pfx_x", "pfx_z", "release_speed", "release_spin_rate", "spin_axis",
          "rho_ratio", "at_coors")
missing_cols <- setdiff(need, names(pitches))
if (length(missing_cols) > 0) {
  stop("pitches_air.parquet is missing: ", paste(missing_cols, collapse = ", "))
}

has_traj <- all(c("vx0", "vy0", "vz0", "ax", "ay", "az") %in% names(pitches))
if (!has_traj) message("No trajectory columns (vx0..az); skipping velocity-loss outcome.")

# Speed at the front of the plate from the 9-parameter fit (y measured from home plate, ft).
plate_speed_mph <- function(vx0, vy0, vz0, ax, ay, az, y0 = 50, yf = 17 / 12) {
  # Solve y0 + vy0*t + 0.5*ay*t^2 = yf for the first positive root (vy0 < 0, ay > 0)
  disc <- vy0^2 - 2 * ay * (y0 - yf)
  t <- (-vy0 - sqrt(disc)) / ay
  vx <- vx0 + ax * t
  vy <- vy0 + ay * t
  vz <- vz0 + az * t
  sqrt(vx^2 + vy^2 + vz^2) * 3600 / 5280
}

# ---- Build analysis frame ---------------------------------------------------

d <- pitches |>
  inner_join(pitch_groups, by = "pitch_type") |>
  filter(!is.na(pfx_x), !is.na(pfx_z), !is.na(release_speed),
         !is.na(release_spin_rate), !is.na(spin_axis), !is.na(rho_ratio)) |>
  mutate(
    ivb_in    = pfx_z * 12,
    # Arm side positive for both hands (RHP arm-side run has negative pfx_x)
    hb_arm_in = if_else(p_throws == "R", -pfx_x, pfx_x) * 12,
    brk_in    = sqrt(ivb_in^2 + hb_arm_in^2),
    spin_k    = release_spin_rate / 1000,
    sin_ax    = sin(spin_axis * pi / 180),
    cos_ax    = cos(spin_axis * pi / 180),
    at_coors_i = as.integer(at_coors),
    cell      = paste(pitcher, pitch_type, game_year, sep = "_")
  )

if (has_traj) {
  d <- d |>
    mutate(
      speed50_mph = sqrt(vx0^2 + vy0^2 + vz0^2) * 3600 / 5280,
      plate_mph   = plate_speed_mph(vx0, vy0, vz0, ax, ay, az),
      velo_loss   = speed50_mph - plate_mph
    )
}

# Direction-normalized components: positive = more movement in the direction the
# pitch normally breaks, so a positive rho coefficient always means "denser air,
# more break" regardless of pitch shape or handedness.
d <- d |>
  group_by(cell) |>
  mutate(
    n_cell  = n(),
    ivb_dir = ivb_in * sign(mean(ivb_in)),
    hb_dir  = hb_arm_in * sign(mean(hb_arm_in))
  ) |>
  ungroup() |>
  filter(n_cell >= min_cell_pitches)

rho_coors <- mean(d$rho_ratio[d$at_coors_i == 1])
rho_else  <- mean(d$rho_ratio[d$at_coors_i == 0])

cat(sprintf("Model frame: %s pitches, %s cells, %s pitchers\n",
            format(nrow(d), big.mark = ","),
            format(n_distinct(d$cell), big.mark = ","),
            format(n_distinct(d$pitcher), big.mark = ",")))
cat(sprintf("Mean rho_ratio: Coors %.3f, elsewhere %.3f\n", rho_coors, rho_else))
cat(sprintf("Cells with any Coors pitches: %s\n",
            format(n_distinct(d$cell[d$at_coors_i == 1]), big.mark = ",")))

# ---- Models -----------------------------------------------------------------

outcomes <- c("brk_in", "ivb_dir", "hb_dir", if (has_traj) "velo_loss", "spin_k")
controls <- "release_speed + spin_k + sin_ax + cos_ax"

fit_one <- function(dat, y, spec) {
  # spin_k is a diagnostic outcome (should be ~flat in rho), so it gets velo only
  rhs <- if (y == "spin_k") "rho_ratio + release_speed" else paste("rho_ratio +", controls)
  if (spec == "rho_plus_coors") rhs <- paste(rhs, "+ at_coors_i")

  m <- feols(as.formula(paste(y, "~", rhs, "| cell")),
             data = dat, cluster = ~ pitcher + game_pk,
             notes = FALSE, warn = FALSE)

  ct <- coeftable(m)
  keep <- intersect(c("rho_ratio", "at_coors_i"), rownames(ct))

  tibble(
    term     = keep,
    estimate = ct[keep, "Estimate"],
    se       = ct[keep, "Std. Error"],
    n        = nobs(m),
    mean_y   = mean(dat[[y]]),
    mean_rho = mean(dat$rho_ratio)
  )
}

grid <- expand_grid(
  pitch_group = sort(unique(d$pitch_group)),
  outcome     = outcomes,
  spec        = c("rho_only", "rho_plus_coors")
)

coefs <- pmap_dfr(grid, function(pitch_group, outcome, spec) {
  dat <- d[d$pitch_group == pitch_group & !is.na(d[[outcome]]), ]
  fit_one(dat, outcome, spec) |>
    mutate(pitch_group = pitch_group, outcome = outcome, spec = spec, .before = 1)
})

coefs <- coefs |>
  mutate(
    t_stat = estimate / se,
    # Elasticity at the means; ~1 is the physics expectation for break and velo loss
    elasticity = if_else(term == "rho_ratio", estimate * mean_rho / mean_y, NA_real_),
    # Implied shift moving the average pitch from non-Coors air to Coors air
    coors_shift = if_else(term == "rho_ratio",
                          estimate * (rho_coors - rho_else),
                          estimate)  # at_coors_i coefficient is already a Coors shift
  )

write_csv(coefs, "output/movement_coefs.csv")

cat("\nRho effect by pitch group (rho_only spec):\n")
coefs |>
  filter(spec == "rho_only", term == "rho_ratio") |>
  transmute(pitch_group, outcome,
            est = round(estimate, 2), se = round(se, 2),
            elasticity = round(elasticity, 2),
            coors_shift = round(coors_shift, 2)) |>
  arrange(outcome, pitch_group) |>
  print(n = Inf)

cat("\nLeftover Coors effect after rho (rho_plus_coors spec):\n")
coefs |>
  filter(spec == "rho_plus_coors", term == "at_coors_i") |>
  transmute(pitch_group, outcome,
            coors_resid = round(estimate, 2), se = round(se, 2),
            t = round(t_stat, 1)) |>
  arrange(outcome, pitch_group) |>
  print(n = Inf)

# Soft sanity check on the headline number
brk_el <- coefs |>
  filter(spec == "rho_only", term == "rho_ratio", outcome == "brk_in")
off <- brk_el |> filter(elasticity < 0.4 | elasticity > 1.6)
if (nrow(off) > 0) {
  warning("Total-break elasticity far from the physics prior of ~1 for: ",
          paste(off$pitch_group, collapse = ", "),
          ". Check sign conventions, pfx units, and the Coors venue flag.")
}

# ---- Per-cell Coors penalty: predicted vs observed --------------------------

el_by_group <- brk_el |> select(pitch_group, el_brk = elasticity)

cells <- d |>
  group_by(cell, pitcher, pitch_type, pitch_group, game_year, p_throws) |>
  summarise(
    n_coors    = sum(at_coors_i == 1),
    n_else     = sum(at_coors_i == 0),
    brk_else   = mean(brk_in[at_coors_i == 0]),
    brk_coors  = if (any(at_coors_i == 1)) mean(brk_in[at_coors_i == 1]) else NA_real_,
    ivb_else   = mean(ivb_in[at_coors_i == 0]),
    hb_else    = mean(hb_arm_in[at_coors_i == 0]),
    rho_else_c = mean(rho_ratio[at_coors_i == 0]),
    .groups = "drop"
  ) |>
  filter(n_else > 0) |>
  left_join(el_by_group, by = "pitch_group") |>
  mutate(
    # Proportional model: bigger-breaking pitches lose more inches at Coors
    pred_delta_brk = el_brk * brk_else * (rho_coors / rho_else_c - 1),
    obs_delta_brk  = if_else(n_coors >= min_coors_pitches, brk_coors - brk_else, NA_real_)
  )

write_csv(cells, "output/movement_by_pitch_cell.csv")

val <- cells |> filter(!is.na(obs_delta_brk))
cat(sprintf(
  "\nPer-cell validation (%d cells with >= %d Coors pitches):\n  mean predicted %.2f in, mean observed %.2f in, r = %.2f\n",
  nrow(val), min_coors_pitches,
  mean(val$pred_delta_brk), mean(val$obs_delta_brk),
  cor(val$pred_delta_brk, val$obs_delta_brk)
))

cat("\nSaved output/movement_coefs.csv and output/movement_by_pitch_cell.csv\n")
