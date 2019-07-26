import sys

import pytest
import calliope
import pandas as pd


def run_test(scenario_results, path_to_biofuel_potentials, scaling_factors, path_to_output, path_to_units):
    exit_code = pytest.main(
        [
            f"--html={path_to_output}",
            f"--self-contained-html",
        ],
        plugins=[_create_config_plugin(
            scenario_results=scenario_results,
            scaling_factors=scaling_factors,
            path_to_biofuel_potentials=path_to_biofuel_potentials,
            path_to_units=path_to_units)]
    )
    sys.exit(exit_code)


def _create_config_plugin(scenario_results, scaling_factors, path_to_biofuel_potentials, path_to_units):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session", params=scenario_results)
        def model(self, request):
            return calliope.read_netcdf(request.param)

        @pytest.fixture(scope="session")
        def scaling_factors(self):
            return scaling_factors

        @pytest.fixture(scope="session")
        def units(self):
            return pd.read_csv(path_to_units, index_col=0)

        @pytest.fixture(
            params=calliope.read_netcdf(scenario_results[0]).inputs.locs.values
        )
        def location(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def carrier_prod(self, model, scaling_factors):
            return model.get_formatted_array("carrier_prod").squeeze("carriers") / scaling_factors["power"]

        @pytest.fixture(scope="session")
        def carrier_con(self, model, scaling_factors):
            return model.get_formatted_array("carrier_con").squeeze("carriers") / scaling_factors["power"]

        @pytest.fixture(scope="session")
        def cost(self, model, scaling_factors):
            return model.get_formatted_array("cost").squeeze("costs") / scaling_factors["monetary"]

        @pytest.fixture(scope="session")
        def resource(self, model):
            resource = model.get_formatted_array("resource")
            timestep_resolution = model.inputs.timestep_resolution
            return resource / timestep_resolution

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
        path_to_units=snakemake.input.units,
        path_to_biofuel_potentials=snakemake.input.biofuel_potentials,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
