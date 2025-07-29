#!/bin/bash

for ext in tiff png pdf;
do
  cp figs/pg_multi.$ext figs/fig1.$ext;
  cp figs/coef_density.$ext figs/fig2.$ext;
  cp figs/ab_points_model4.$ext figs/fig3.$ext;
  cp figs/traits_sp_points_main_re.$ext figs/fig4.$ext;
  cp figs/tr_bar_comb2.$ext figs/fig5.$ext;

  cp figs/pool_multi.$ext figs/fig_s4.$ext;
  cp figs/ab_points_four_models_all.$ext figs/fig_s5.$ext;
  cp figs/pg_ribbon_a.$ext figs/fig_s6.$ext;
  cp figs/pg_ribbon_b.$ext figs/fig_s7.$ext;
  cp figs/ab_pg_summary_bars.$ext figs/fig_s8.$ext;
  cp figs/traits_sp_points_si.$ext figs/fig_s9_1.$ext;
  cp figs/traits_seg_points_si.$ext figs/fig_s9_2.$ext;
done
