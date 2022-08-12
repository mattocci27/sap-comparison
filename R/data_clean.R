
d <- read_csv("data-raw/pres_tens_five_spp.csv")

d$species |> as.factor() |> as.numeric()
d$pressure |> as.factor()
