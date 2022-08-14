
my_theme <- function(){
  theme_bw() %+replace%
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
}

# theme_set(my_theme)

sma_scatter <- function(data, log = FALSE) {
#  targets::tar_load(five_spp_csv)
#  data <- five_spp_csv
  data <- read_csv(data)
  my_col <- RColorBrewer::brewer.pal(3, "Set2")
  p <- data |>
    mutate(pres_fac = paste(pressure, "MPa") |> as.factor()) |>
    ggplot(aes(x = pres_calib, y = tens_calib)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey60") +
    geom_sma(se = FALSE) +
    geom_point(aes(col = pres_fac)) +
    ylab(expression("Flow rate under tension "~(g~s^{-1}))) +
    xlab(expression("Flow rate by pressure "~(g~s^{-1}))) +
    scale_colour_manual(
      values = my_col,
      name = "") +
    stat_cor(
        label.x.npc = 0.45,
        label.y.npc = 0.01,
        vjust = 1,
        p.accuracy = 0.0001,
        aes(label = paste(..rr.label.., ..p.label.. ,sep = "~`,`~"))
      ) +
    my_theme() +
    theme(
      legend.position = c(0.2, 0.8)
    )
   if (log) {
     p +
      annotate(geom = "text", x = 3.5, y = 4,
        label = "1:1 line", angle = 45, col = "grey60") +
      scale_x_log10() +
      scale_y_log10() +
      coord_fixed()
   } else {
     p +
      annotate(geom = "text", x = 8, y = 9,
        label = "1:1 line", angle = 45, col = "grey60") +
      coord_cartesian(
        xlim = c(0,max(c(data$tens_calib, data$pres_calib))),
        ylim = c(0,max(c(data$tens_calib, data$pres_calib))))
   }
}
