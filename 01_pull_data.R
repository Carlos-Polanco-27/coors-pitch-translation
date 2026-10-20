# 01_pull_data.R
# Coors pitch translation: pull Statcast pitch-level data (2023-2025)
# and game temperatures from the MLB Stats API.
# Resumable: re-run at any time and it skips anything already saved.

library(baseballr)
library(dplyr)
library(purrr)
library(jsonlite)
library(arrow)

dir.create("data/raw/statcast", recursive = TRUE, showWarnings = FALSE)
dir.create("data/clean", recursive = TRUE, showWarnings = FALSE)

# ---------- settings ----------
seasons <- list(
  `2023` = c("2023-03-28", "2023-10-01"),
  `2024` = c("2024-03-18", "2024-09-30"),
  `2025` = c("2025-03-17", "2025-09-28")
)

keep_cols <- c(
  "game_date", "game_pk", "game_type", "home_team", "away_team",
  "pitcher", "player_name", "p_throws", "pitch_type",
  "release_speed", "release_spin_rate", "spin_axis",
  "pfx_x", "pfx_z", "release_extension", "release_pos_x", "release_pos_z",
  "plate_x", "plate_z", "description", "events",
  "launch_speed", "launch_angle", "estimated_woba_using_speedangle",
  "delta_run_exp"
)

chr_cols <- c("game_date", "game_type", "home_team", "away_team",
              "player_name", "p_throws", "pitch_type", "description", "events")

# ---------- 1. Statcast, one day at a time ----------
pull_day <- function(d) {
  f <- file.path("data/raw/statcast", paste0(d, ".parquet"))
  if (file.exists(f)) return(invisible(NULL))

  x <- tryCatch(
    statcast_search(start_date = d, end_date = d, player_type = "pitcher"),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) {
    message(d, ": no data")
    return(invisible(NULL))
  }

  x <- x %>%
    select(any_of(keep_cols)) %>%
    mutate(across(any_of(chr_cols), as.character),
           across(!any_of(chr_cols), as.numeric)) %>%
    filter(game_type == "R")

  write_parquet(x, f)
  message(d, ": ", nrow(x), " pitches")
  Sys.sleep(1)
}

for (s in names(seasons)) {
  days <- seq(as.Date(seasons[[s]][1]), as.Date(seasons[[s]][2]), by = "day")
  walk(as.character(days), pull_day)
}

files <- list.files("data/raw/statcast", full.names = TRUE)
pitches <- map_dfr(files, read_parquet) %>%
  mutate(game_date = as.Date(game_date),
         season = as.integer(format(game_date, "%Y")))

# Keep every pitch from any pitcher who threw at Coors in that season
coors_pitchers <- pitches %>%
  filter(home_team == "COL") %>%
  distinct(pitcher, season)

pitches <- pitches %>% semi_join(coors_pitchers, by = c("pitcher", "season"))
write_parquet(pitches, "data/clean/pitches.parquet")
message("Saved ", nrow(pitches), " pitches from ",
        n_distinct(pitches$pitcher), " pitchers")

# ---------- 2. Game temperatures ----------
weather_file <- "data/clean/weather.csv"
done <- if (file.exists(weather_file)) read.csv(weather_file)$game_pk else integer(0)
todo <- setdiff(unique(pitches$game_pk), done)
message(length(todo), " games need weather")

get_weather <- function(pk) {
  url <- paste0("https://statsapi.mlb.com/api/v1.1/game/", pk, "/feed/live")
  g <- tryCatch(fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
  w <- g$gameData$weather
  tibble(
    game_pk   = pk,
    temp_f    = suppressWarnings(as.numeric(w$temp %||% NA)),
    condition = as.character(w$condition %||% NA),
    wind      = as.character(w$wind %||% NA),
    venue     = as.character(g$gameData$venue$name %||% NA)
  )
}

# Save in batches of 100 so progress isn't lost
batches <- split(todo, ceiling(seq_along(todo) / 100))
for (b in batches) {
  out <- map_dfr(b, get_weather)
  write.table(out, weather_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(weather_file), append = file.exists(weather_file))
  message("Weather saved through game ", tail(b, 1))
}

message("Done.")
