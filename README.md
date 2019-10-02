# Geographical Scale

Analysis of the costs of electricity autarky and small scale grids.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

1. Clone the repo. Because `euro-calliope` is added as a git submodule, you may want to clone using `git clone --recurse-submodules <link-to-this-repo>`.

2. Create an environment to run the analysis. You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    `conda env create -f environment.yaml`

3. Make sure you have a Gurobi license, or install and configure another solver.

4. You need an account at the Copernicus Climate Data Service and you need to create a `$HOME/.cdsapirc` file with your credentials, see their [How To](https://cds.climate.copernicus.eu/api-how-to) (you do not need to manually install the client).

5. Provide the input data for Euro-Calliope, as defined in "Getting Ready" in  `./euro-calliope/README.md`.

6. To run the sensitivity analysis (optional, involves manual steps to run see below), you need MATLAB and UQLab installed:

    1. Install [MATLAB R2019a](https://de.mathworks.com/products/matlab.html).

    2. [Register](https://www.uqlab.com/register) and [install](https://www.uqlab.com/install) UQLab.

    3. Add UQLab's `core` folder to the [MATLAB search path on startup](https://ch.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html), to be able to call `UQLab` from outside its own folder.

## Run the analysis

    snakemake --use-conda

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run on Euler cluster

To run on Euler, use the following command:

    snakemake --use-conda --profile config/euler [--config email=<you@provider.org>]

By providing an email address, you will be informed by mail when Snakemake finishes execution.

If you prefer working locally, you can sync this repository to Euler and receive build changes by running `snakemake send` and `snakemake receive`.

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Manual steps

At the moment, the uncertainty quantification is not in the loop and needs to be run manually. The following files have to be created manually:

* `./data/first-sobol-continental-national.csv`
* `./data/first-sobol-regional.csv`
* `./data/total-minus-first-sobol-continental-national.csv`
* `./data/total-minus-first-first-sobol-regional.csv`
* `./data/tota-sobol-continental-national.csv`
* `./data/total-sobol-regional.csv`
* `./data/pce-samples-continental-national-scales.csv`
* `./data/pce-samples-regional-scale.csv`

To replicate these files, do the following:

1. Perform necessary simulation runs:
    1. `snakemake --use-conda --configfile config/full.yaml build/output/national/uncertainty/continental-autarky-100-continental-grid/xy.csv build/output/national/uncertainty/national-autarky-100-national-grid/xy.csv build/output/regional/uncertainty/regional-autarky-100-regional-grid/xy.csv`
    2. `snakemake --use-conda --configfile config/multi-fidelity.yaml build/output/regional/uncertainty/continental-autarky-100-continental-grid/xy.csv build/output/regional/uncertainty/national-autarky-100-national-grid/xy.csv`
2. Manually merge results to necessary format for Matlab script:
    1. Merge national resolution results for continental and national scale into `xy1-national-resolution-national-and-continental-scales.csv` (also add absolute and relative cost diff).
    2. Merge regional resolution results for continental and national scale into `xy2-regional-resolution-national-and-continental-scales.csv` (also add absolute and relative cost diff).
    3. Rename regional resolution results for regional scale to `xy3-regional-resolution-regional-scale.csv`.
3. Add the following files into one single folder:
    * `./src/uncertainty/uq.m`
    * `./data/GeoScale_uncertain-parameters.csv`
    * `xy1-national-resolution-national-and-continental-scales.csv`
    * `xy2-regional-resolution-national-and-continental-scales.csv`
    * `xy3-regional-resolution-regional-scale.csv`
4. Change the value of `CALCULATE` to 1 in `uq.m`.
4. Run the `uq.m` using Matlab (make sure you have UQLab installed).
5. Copy the result files to their above mentioned locations within this repository.

## Run the tests

    snakemake test --use-conda

## Units and scaling

The default units for Euro-Calliope are `MW`, `MWh`, `EUR`, and `km2`, but you can scale all of these using the configuration values in `config/default.yaml`. Apart from convenience, this may be important to handle numerical issues with your solver.

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`. This excludes the KlinicSlab font family (all files in `./report/fonts/`) which is copyright Lost Type.
