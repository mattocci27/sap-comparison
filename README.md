# sap flow comparison

## Usage

### Running code on local

To run analysis:

```
# To install R packages for the first run
# Rscript -e "renv::restore()"

run Rscript
```

### Running code in Apptainer (Linux)

To build Apptainer containers:

```bash
sudo apptainer build radian.sif radian.def
```

```bash
apptainer exec --env RENV_PATHS_CACHE=/home/${USER}/renv \
	--env RENV_PATHS_PREFIX_AUTO=TRUE \
	radian.sif Rscript -e "renv::restore()"
```

To run analysis:

```bash
apptainer exec --env RENV_PATHS_CACHE=/home/${USER}/renv \
	--env RENV_PATHS_PREFIX_AUTO=TRUE \
	radian.sif Rscript run.R
```

`.Renviron` is ignored in git.
It is useful to add the following information to `.Renviron` to use `{renv}` caches on Linux.

```
#renv settings
RENV_PATHS_PREFIX_AUTO=TRUE
RENV_PATHS_CACHE=/home/${USER}/renv
```
