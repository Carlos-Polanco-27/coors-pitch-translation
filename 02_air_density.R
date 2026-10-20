# 02_air_density.R
# Coors pitch translation: estimate air density for every game and attach it
# to the pitch data. Drag and Magnus force both scale with air density, so
# this is the main environmental input for translating movement in/out of Coors.
#
# Inputs:  data/clean/pitches.parquet, data/clean/weather.csv (from 01)
# Outputs: data/clean/game_air.csv, data/clean/pitches_air.parquet

library(dplyr)
library(tidyr)
library(stringr)
library(arrow)

pitches <- read_parquet("data/clean/pitches.parquet")
weather <- read.csv("data/clean/weather.csv", stringsAsFactors = FALSE) %>%
  mutate(game_pk = as.numeric(game_pk)) %>%   # pitches$game_pk is dbl
  distinct(game_pk, .keep_all = TRUE)

# ---------- 1. Park elevations (feet) ----------
# Approximate field elevations. Anything near sea level barely matters
# (100 ft is roughly 0.4% density); Coors is the one that has to be right.
park_elev <- tribble(
  ~home_team, ~elev_ft,
  "AZ",   1086,  "ARI",  1086,
  "ATL",  1050,  "BAL",    30,  "BOS",    20,
  "CHC",   595,  "CWS",   595,  "CIN",   490,
  "CLE",   650,  "COL",  5190,  "DET",   600,
  "HOU",    40,  "KC",    750,  "LAA",   160,
  "LAD",   500,  "MIA",    10,  "MIL",   635,
  "MIN",   815,  "NYM",    10,  "NYY",    55,
  "OAK",    25,  "ATH",    25,  "PHI",    20,
  "PIT",   730,  "SD",     15,  "SF",     10,
  "SEA",    15,  "STL",   465,  "TB",     15,
  "TEX",   550,  "TOR",   270,  "WSH",    25
)

# Neutral-site games keep the "home" team code but were played elsewhere.
# Mexico City (2023 SD/SF) sits higher than Coors, so it must be caught.
neutral_sites <- tribble(
  ~pattern,        ~elev_ft,
  "Harp Hel",       7350,   # Estadio Alfredo Harp Helu, Mexico City
  "London",           20,
  "Gocheok",         130,   # Seoul
  "Tokyo Dome",       20,
  "Rickwood",        600,   # Birmingham
  "Bowman",          530,   # Williamsport
  "Bristol",        1700    # Speedway Classic
)

# ---------- 2. One row per game ----------
games <- pitches %>%
  distinct(game_pk, game_date, season, home_team) %>%
  left_join(weather, by = "game_pk") %>%
  left_join(park_elev, by = "home_team") %>%
  mutate(neutral_site = FALSE)

for (i in seq_len(nrow(neutral_sites))) {
  hit <- !is.na(games$venue) &
    str_detect(games$venue, regex(neutral_sites$pattern[i], ignore_case = TRUE))
  games$elev_ft[hit]      <- neutral_sites$elev_ft[i]
  games$neutral_site[hit] <- TRUE
}

# Venues that show up rarely for a team are usually neutral sites or
# temporary homes. Eyeball this list once to make sure nothing slipped by.
message("Rare venues by home team (check for unhandled neutral sites):")
games %>%
  count(home_team, venue) %>%
  group_by(home_team) %>%
  filter(n < 10) %>%
  ungroup() %>%
  arrange(home_team) %>%
  print(n = Inf)

unmatched <- games %>% filter(is.na(elev_ft)) %>% distinct(home_team)
if (nrow(unmatched) > 0) {
  stop("No elevation for home_team: ", paste(unmatched$home_team, collapse = ", "))
}

# ---------- 3. Clean weather fields ----------
# temp_f of 0 is the feed's placeholder for "not reported", not a real
# reading (no MLB game has been played anywhere near 0F). Anything under
# 20F is treated the same way; real cold games bottom out in the 30s.
message(sum(games$temp_f < 20, na.rm = TRUE), " games with temp_f < 20 set to missing (",
        sum(games$temp_f == 0, na.rm = TRUE), " of them exactly 0)")
games %>%
  filter(is.na(temp_f) | temp_f < 20) %>%
  count(home_team, condition, sort = TRUE) %>%
  print(n = 20)

games <- games %>%
  mutate(
    month       = as.integer(format(game_date, "%m")),
    roof_closed = str_detect(coalesce(condition, ""), regex("dome|roof closed", ignore_case = TRUE)),
    wind_mph    = as.numeric(str_extract(wind, "^\\d+")),
    wind_dir    = str_trim(str_remove(wind, "^\\d+\\s*mph,?")),
    temp_f      = if_else(temp_f < 20, NA_real_, temp_f),
    temp_raw    = temp_f,
    # Closed roof / dome with no reading: the building is climate controlled,
    # so use a typical indoor temp rather than the outdoor park median
    temp_f      = if_else(is.na(temp_f) & roof_closed, 72, temp_f)
  ) %>%
  # Fill remaining missing temps: same park and month, then same park, then overall
  group_by(home_team, month) %>%
  mutate(temp_f = coalesce(temp_f, median(temp_raw, na.rm = TRUE))) %>%
  group_by(home_team) %>%
  mutate(temp_f = coalesce(temp_f, median(temp_raw, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(temp_f       = coalesce(temp_f, median(temp_raw, na.rm = TRUE)),
         temp_imputed = is.na(temp_raw)) %>%
  select(-temp_raw)

failed_pulls <- sum(is.na(games$venue))
message(sum(games$temp_imputed), " games had temperature imputed; ",
        failed_pulls, " games have no weather row at all (failed API pulls)")

# ---------- 4. Air density ----------
# Alan Nathan's trajectory-calculator formula. Pressure is the sea-level
# corrected barometer reading, adjusted for elevation inside the formula.
# Humidity and pressure aren't in the public feed, so both use defaults.
air_density <- function(temp_f, elev_ft, rh = 0.5, baro_inhg = 29.92) {
  tc     <- (temp_f - 32) * 5 / 9
  elev_m <- elev_ft * 0.3048
  p_mmhg <- baro_inhg * 25.4
  svp    <- 4.5841 * exp((18.687 - tc / 234.5) * tc / (257.14 + tc))
  1.2929 * (273 / (tc + 273)) *
    (p_mmhg * exp(-0.0001217 * elev_m) - 0.3783 * rh * svp) / 760
}

rho_ref <- air_density(70, 0)   # 70F, sea level, 50% RH (about 1.194 kg/m^3)

games <- games %>%
  mutate(rho       = air_density(temp_f, elev_ft),
         rho_ratio = rho / rho_ref)

# ---------- 5. Sanity checks ----------
park_summary <- games %>%
  filter(!neutral_site) %>%
  group_by(home_team) %>%
  summarise(games     = n(),
            elev_ft   = first(elev_ft),
            avg_temp  = round(mean(temp_f), 1),
            avg_ratio = round(mean(rho_ratio), 3),
            .groups = "drop") %>%
  arrange(avg_ratio)
print(park_summary, n = Inf)

coors_ratio <- park_summary$avg_ratio[park_summary$home_team == "COL"]
stopifnot(coors_ratio > 0.78, coors_ratio < 0.87)
stopifnot(all(between(games$rho_ratio, 0.70, 1.10)))

# ---------- 6. Save ----------
games_out <- games %>%
  select(game_pk, game_date, season, home_team, venue, neutral_site,
         elev_ft, temp_f, temp_imputed, roof_closed, condition,
         wind_mph, wind_dir, rho, rho_ratio)

write.csv(games_out, "data/clean/game_air.csv", row.names = FALSE)

pitches_air <- pitches %>%
  left_join(games_out %>% select(game_pk, elev_ft, temp_f, roof_closed,
                                 neutral_site, rho, rho_ratio),
            by = "game_pk") %>%
  mutate(at_coors = home_team == "COL" & !neutral_site)

stopifnot(nrow(pitches_air) == nrow(pitches), !any(is.na(pitches_air$rho)))
write_parquet(pitches_air, "data/clean/pitches_air.parquet")

message("Saved ", nrow(games_out), " games and ", nrow(pitches_air),
        " pitches with air density. Coors avg density ratio: ", coors_ratio)
