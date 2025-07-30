# Hydraulic conductivity-induced systematic parameter variation in a widely used thermal dissipation sap-flow technique

Chen Ya-Jun,
Maenpuen Phisamai,
Katabuchi Masatoshi,
Tor-Ngern Pantana,
Palmroth Sari,
Zhang Shubin,
Xiao Yunxue,
Liu Meng,
Oren Ram

Code repository to run the analysis and generate the manuscript for Chen et al. "Hydraulic conductivity-induced systematic parameter variation in a widely used thermal dissipation sap-flow technique".
https://doi.org/xxxx.

## Usage

### Running the analysis

Use the run.sh script to start the analysis.
It provides options for running locally, in Docker, or inside an Apptainer/Singularity container.

```bash
./run.sh
```

You will be prompted to select an option:

```bash
1) Run tar_make() locally (or inside Docker)
2) Run tar_make() in Apptainer
3) Enter the Apptainer container
4) Enter the Singularity container on HPC
Enter number:
```

Note: If you are working inside a Docker contianer, choose option 1.

### Building the Apptainer container (Linux)

To build the container:

```bash
sudo apptainer build radian.sif radian.def
```

### Optional: Set up .Renviron for shared renv caching (Linux)

The `.Renviron` file is excluded from version control.
To enable shared `{renv}` caches on Linux, add the following to your `.Renviron`:

```bash
# renv settings
RENV_PATHS_PREFIX_AUTO=TRUE
RENV_PATHS_CACHE=/home/${USER}/renv
```

