# sap flow comparison

```
# To install R packages for the first run
# apptainer exec radian.sif Rscript -e "renv::restore()"
apptainer exec --env RENV_PATHS_CACHE=/home/mattocci/renv \
	--env RENV_PATHS_PREFIX_AUTO=TRUE \
	radian.sif Rscript run.R
```
