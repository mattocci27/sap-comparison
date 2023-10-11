#For general run
# targets::tar_make_clustermq(
#   workers = parallel::detectCores()
# )

#For single thread
#targets::tar_make()

#For HPC
arg = commandArgs(T)
targets::tar_make_clustermq(
                            # c(imp_points_plot, imp_points_plot2),
                #  segments_inclusive_ab_csv,
  workers = as.numeric(arg[1])
)

# For building mcmc draws, diagnostics, summary
# This saves RAM
# targets::tar_make_clustermq(
#   workers = 2
# )
