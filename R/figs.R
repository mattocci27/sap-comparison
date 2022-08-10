
targets::tar_load(five_spp_csv)

my_theme <- theme_bw() %+replace%
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

theme_set(my_theme)

library(tidyverse)
library(ggsma)
library(ggpubr)
library(RColorBrewer)
d <- read_csv(five_spp_csv)

#hist(d$pres_calib |> log())

my_col <- RColorBrewer::brewer.pal(3, "Set2")
d |>
  mutate(pres_fac = paste(pressure, "MPa") |> as.factor()) |>
  ggplot(aes(x = pres_calib, y = tens_calib)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey60") +
  annotate(geom = "text", x = 8, y = 9,
    label = "1:1 line", angle = 45, col = "grey60") +
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
  coord_cartesian(
    xlim = c(0,max(c(d$tens_calib, d$pres_calib))),
    ylim = c(0,max(c(d$tens_calib, d$pres_calib)))) +
  theme(
    legend.position = c(0.2, 0.8)
  )
