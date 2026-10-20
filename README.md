# Pitch Movement at Coors Field

## The question

Pitches break less at Coors Field. Can we predict how much a pitcher's movement will change at Coors from how his pitches move everywhere else, and does lower air density account for that change on its own?

## Data

Statcast pitch-level data from Baseball Savant, pulled with baseballr, covering the 2023 to 2025 regular seasons. Game temperatures come from the MLB Stats API and are combined with park elevation to compute air density for each game. The model is fit on 2023 and 2024 and validated on 2025.

Each validation observation is a pitcher x pitch type cell from 2025, kept only if the pitcher threw at least 15 of that pitch at Coors and at least 30 in other parks. That leaves 347 cells across 11 pitch types. Movement is measured in inches for three outcomes: total break, induced vertical break (IVB), and horizontal break (arm side positive).

## Method

For each cell, the model takes the pitcher's movement outside Coors and adjusts it for the difference in air density. The size of that adjustment is set by an elasticity, fit separately for fastballs, breaking balls, and offspeed pitches. Here the elasticity is the percent change in break for a 1% change in air density, measured at the sample averages (from a linear fit). A second version adds a Coors-specific term on top of density. Both are compared to a no-change baseline that assumes a pitch moves the same at Coors as it does elsewhere.

Total break terms fit on 2023 and 2024:

| Pitch group | Elasticity (density only) | Elasticity (with Coors term) | Coors residual (in) |
|---|---|---|---|
| Breaking | 1.04 | 0.92 | -0.32 |
| Fastball | 1.16 | 1.08 | -0.32 |
| Offspeed | 1.32 | 1.14 | -0.70 |

Metrics are mean absolute error (MAE, inches), bias (predicted minus observed), and a noise floor, which is the MAE we would expect from sampling noise alone given cell sizes. `pct_better` is the MAE reduction versus baseline and `wins` is the share of cells where the model lands closer than baseline.

## Results

**Total break MAE is 0.66 in, compared with 2.87 in for the no-change baseline. The model beats baseline in 91% of the 347 cells.**

| Outcome | Baseline MAE | Density MAE | Density + Coors MAE | Floor | % better | Wins |
|---|---|---|---|---|---|---|
| Total break | 2.87 | 0.66 | 0.67 | 0.23 | 77 | 91% |
| IVB | 1.95 | 0.80 | 0.81 | 0.25 | 58 | 74% |
| Horizontal break | 2.12 | 0.79 | 0.80 | 0.26 | 62 | 80% |

Air density explains the Coors effect on its own. Adding a Coors-specific term does not lower error on any outcome, and the density model's total break bias is close to zero (-0.07 in). Fastballs and sinkers show the largest gains, with total break error down 83% and 85% from baseline.

## Limitations

- IVB is weak for sliders, sweepers, sinkers, and changeups. Gains over baseline are only 19 to 33%, wins run 50 to 60%, and sweeper IVB is overpredicted by about 1.5 in.
- Cutter horizontal break is worse than baseline (MAE 1.06 vs 0.92, wins in 29% of cells).
- Sample sizes are small for slurves (2 cells), forkballs (1), knuckle curves (8), and splitters (8). Results for those pitch types should not be relied on.
- Overall error is still well above the noise floor, so there is real movement change the model is not capturing.

## Next steps

- Add seam-shifted wake terms, which should help most where density alone falls short (sinker and changeup IVB, cutter horizontal break).
- Translate predicted movement changes into run value to measure how much the Coors effect matters for pitch effectiveness.

## Repo structure

- `01_pull_data.R`: pulls Statcast and game temperatures
- `02_air_density.R`: adds park elevation and computes air density per game
- `03_movement_model.R`: fits movement as a function of air density by pitch group
- `04_validation.R`: fits on 2023-2024, predicts 2025, compares to a no-change baseline

Data is not included. Run `01_pull_data.R` to rebuild it.
