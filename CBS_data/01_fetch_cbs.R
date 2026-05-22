# ============================================================================
# 01_fetch_cbs.R
#
# Run this script ONCE to download the six CBS Caribbean Netherlands tables
# to local CSVs. After that, work entirely from the cached CSVs and don't
# touch this file again unless you want fresh data.
# ============================================================================

library(cbsodataR)
library(readr)

dir.create("cbs_raw", showWarnings = FALSE, recursive = TRUE)

ids <- c(
  tourism_by_plane       = "83104NED",
  tourism_by_nationality = "83191NED",
  yachts_mooring         = "85015NED",
  aviation               = "82332NED",
  cruise_passengers      = "85007NED",
  value_per_sector       = "84769NED"
)

for (i in seq_along(ids)) {
  id <- ids[i]
  cat(sprintf("[%d/%d] Fetching %s (%s) ... ", i, length(ids), id, names(ids)[i]))

  df <- cbs_get_data(id) |>
    cbs_add_label_columns() |>   # adds *_label columns (e.g. RegioS_label = "Bonaire")
    cbs_add_date_column()        # adds Perioden_Date + Perioden_freq

  out_csv <- file.path("cbs_raw", paste0(id, ".csv"))
  write_csv(df, out_csv)
  cat(sprintf("%d rows -> %s\n", nrow(df), out_csv))

  # Also dump the metadata so you can confirm column units/titles
  meta <- cbs_get_meta(id)
  meta_path <- file.path("cbs_raw", paste0(id, "_meta.txt"))
  sink(meta_path)
  cat("Title:    ", meta$TableInfos$Title,       "\n")
  cat("Modified: ", as.character(meta$TableInfos$Modified), "\n")
  cat("Frequency:", meta$TableInfos$Frequency,   "\n")
  cat("Period:   ", meta$TableInfos$Period,      "\n\n")
  cat("Columns (DataProperties):\n")
  if (!is.null(meta$DataProperties)) {
    print(meta$DataProperties[, intersect(c("Key","Type","Title","Unit"),
                                          names(meta$DataProperties))])
  }
  sink()
}

cat("\nDone. Inspect cbs_raw/*_meta.txt for column descriptions before running 02_analyze.R.\n")
