import sys

import pytest
import calliope
import pandas as pd


def run_test(scenario_results, path_to_biofuel_potentials, path_to_output):
    exit_code = pytest.main(
        [
            f"--html={path_to_output}",
            f"--self-contained-html",
        ],
        plugins=[_create_config_plugin(scenario_results, path_to_biofuel_potentials)]
    )
    sys.exit(exit_code)


def _create_config_plugin(scenario_results, path_to_biofuel_potentials):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session", params=scenario_results)
        def model(self, request):
            return calliope.read_netcdf(request.param)

        @pytest.fixture(
            params=calliope.read_netcdf(scenario_results[0]).inputs.locs.values
        )
        def location(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def carrier_prod(self, model):
            return model.get_formatted_array("carrier_prod").squeeze("carriers")

        @pytest.fixture(scope="session")
        def biofuel_potentials(self):
            return (pd.read_csv(path_to_biofuel_potentials, index_col=0)
                      .rename(index=lambda loc: loc.replace(".", "-")))

        @pytest.fixture()
        def biofuel_potential(self, biofuel_potentials, location):
            return biofuel_potentials.loc[location, "biofuel_potential_mwh_per_year"]

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(
        scenario_results=snakemake.input.results,
        path_to_biofuel_potentials=snakemake.input.biofuel_potentials,
        path_to_output=snakemake.output[0]
    )
