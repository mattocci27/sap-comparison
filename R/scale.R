knn_impute <- function(rubber_raw_data_csv) {
  d <- read_csv(rubber_raw_data_csv) |>
    janitor::clean_names() |>
    mutate(date = mdy_hm(date))
  kNN(d)
}

