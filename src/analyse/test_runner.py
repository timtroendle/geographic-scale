import sys

import pytest
import calliope


def run_test(scenario_results, path_to_output):
    exit_code = pytest.main(
        [
            f"--html={path_to_output}",
            f"--self-contained-html",
        ],
        plugins=[_create_config_plugin(scenario_results)]
    )
    sys.exit(exit_code)


def _create_config_plugin(scenario_results):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session", params=scenario_results)
        def model(self, request):
            return calliope.read_netcdf(request.param)

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(
        scenario_results=snakemake.input.results,
        path_to_output=snakemake.output[0]
    )
