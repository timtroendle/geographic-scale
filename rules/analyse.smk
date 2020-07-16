localrules: overview_scenario_results

RESAMPLE_OVERRIDE_DICT = {
    "model.time.function": "resample",
    "model.time.function_options": {"resolution": config["resolution"]["time"]}
}
CALLIOPE_OVERRIDE_DICT = {
    **config["calliope-parameters"],
    **RESAMPLE_OVERRIDE_DICT
}


rule run_national:
    message: "Run the model for scenario {wildcards.scenario} on national resolution."
    input:
        model = "build/model/national/model.yaml"
    params: override_dict = CALLIOPE_OVERRIDE_DICT
    output: "build/output/national/runs/{scenario}.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --override_dict="{params.override_dict}"
        """


rule run_regional: # this is a copy of run_national which is necessary to have different cluster configs
    message: "Run the model for scenario {wildcards.scenario} on regional resolution."
    input:
        model = "build/model/regional/model.yaml"
    params: override_dict = CALLIOPE_OVERRIDE_DICT
    output: "build/output/regional/runs/{scenario}.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """ # because runs take a lot of computation time, don't fail when suboptimal
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --no_fail_when_infeasible\
            --override_dict="{params.override_dict}"
        """


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_constraints.py",
        "tests/test_assumptions.py",
        results = expand(
            "build/output/{{resolution}}/runs/{scenario}.nc",
            scenario=config["scenarios"]
        ),
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params:
        scaling_factors = config["scaling-factors"],
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule aggregated_results:
    message: "Create csv overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = expand("build/output/{{resolution}}/runs/{scenario}.nc", scenario=config["scenarios"]),
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/geo.yaml"
    output: "build/output/{resolution}/aggregation.csv"
    script: "../src/analyse/aggregation.py"


rule time_aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/time_aggregation.py",
        scenarios = expand("build/output/{{resolution}}/runs/{scenario}.nc", scenario=config["scenarios"]),
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/geo.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/time_aggregation.py"


rule plot_cost:
    message: "Plot overview over cost."
    input:
        src = "src/analyse/cost.py",
        matplotlibrc = "matplotlibrc",
        results = rules.time_aggregated_results.output[0]
    output:
        base = "build/output/{resolution}/cost.{plot_suffix}",
        special = "build/output/{resolution}/cost-special-cases.{plot_suffix}",
    conda: "../envs/default.yaml"
    script: "../src/analyse/cost.py"


rule plot_map_cost:
    message: "Create map of cost."
    input:
        src = "src/analyse/map_cost.py",
        matplotlibrc = "matplotlibrc",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        national_shapes = eurocalliope("build/data/national/units.geojson"),
        regional_shapes = eurocalliope("build/data/regional/units.geojson"),
        results = rules.time_aggregated_results.output[0]
    params:
        connected_regions = config["connected-regions"]
    output: "build/output/{resolution}/map-cost.{plot_suffix}"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map_cost.py"


rule plot_map_energy:
    message: "Create map of generation and transmission."
    input:
        src = "src/analyse/map_energy.py",
        matplotlibrc = "matplotlibrc",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        results = rules.time_aggregated_results.output[0]
    params:
        connected_regions = config["connected-regions"],
        scenarios = ["continental-autarky-100-continental-grid",
                     "national-autarky-100-continental-grid",
                     "regional-autarky-100-continental-grid"]
    output: "build/output/{resolution}/map-energy.{plot_suffix}"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map_energy.py"


rule plot_network_map:
    message: "Create map of electricity grid."
    input:
        src = "src/analyse/network.py",
        matplotlibrc = "matplotlibrc",
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        results = "build/output/{resolution}/runs/continental-autarky-100-continental-grid.nc"
    output: "build/output/{resolution}/network.{plot_suffix}"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/network.py"


rule plot_bubble_map:
    message: "Create bubble map of generation infrastructure."
    input:
        src = "src/analyse/bubble_map.py",
        matplotlibrc = "matplotlibrc",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        continent_shape = eurocalliope("build/data/continental/units.geojson"),
        results = rules.time_aggregated_results.output[0]
    params: resolution_km = 100
    output: "build/output/{resolution}/map-generation--{scenario}--{colour}--{markersize}.{plot_suffix}"
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        markersize = "((gen)|(\d\d))",
        colour = "((blue)|(yellow))"
    script: "../src/analyse/bubble_map.py"


rule graphical_abstract:
    message: "Create maps for graphical abstract."
    input:
        "build/output/regional/map-generation--continental-autarky-100-continental-grid--blue--gen.png",
        "build/output/regional/map-generation--regional-autarky-100-continental-grid--yellow--10.png",


rule plot_system_composition_per_scale:
    message: "Create plot of system composition on both scales."
    input:
        src = "src/analyse/composition.py",
        matplotlibrc = "matplotlibrc",
        results_csv = rules.aggregated_results.output[0],
        results_nc = rules.time_aggregated_results.output[0]
    params:
        transmission_today_twkm = config["report"]["interregional-transmission-capacity-today"],
        crossborder_today_tw = config["report"]["crossborder-transmission-capacity-today"]
    output: "build/output/{resolution}/composition.{plot_suffix}",
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition.py"


rule plot_system_composition_full:
    message: "Create plot of system composition of all scenarios."
    input:
        src = "src/analyse/composition_all.py",
        matplotlibrc = "matplotlibrc",
        results = rules.aggregated_results.output[0]
    output: "build/output/{resolution}/composition-all.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition_all.py"


rule plot_uncertainty_of_composition:
    message: "Create plot of uncertainty of system composition found in uncertainty analysis."
    input:
        src = "src/analyse/composition_uncertainty.py",
        matplotlibrc = "matplotlibrc",
        large_scales = "build/uncertainty/{}/pce-samples-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        small_scale = "build/uncertainty/{}/pce-samples-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
    output: "build/output/{resolution}/composition-uncertainty.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition_uncertainty.py"


rule plot_uncertainty_of_cost:
    message: "Create plot of uncertainty of system cost found in uncertainty analysis."
    input:
        src = "src/analyse/cost_uncertainty.py",
        matplotlibrc = "matplotlibrc",
        large_scales = "build/uncertainty/{}/pce-samples-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        small_scale = "build/uncertainty/{}/pce-samples-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        aggregate = rules.time_aggregated_results.output[0]
    output: "build/output/{resolution}/cost-uncertainty.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/cost_uncertainty.py"


rule plot_timeseries:
    message: "Create plot of timeseries on {wildcards.resolution} resolution."
    input:
        src = "src/analyse/timeseries.py",
        matplotlibrc = "matplotlibrc",
        result = "build/output/{resolution}/runs/regional-autarky-100-regional-grid.nc",
        units = eurocalliope("build/data/{resolution}/units.csv")
    params:
        connected_regions = config["connected-regions"],
        scaling_factors = config["scaling-factors"],
        resolution = config["report"]["timeseries-plot"]["resolution"]
    output:
        plot = "build/output/{resolution}/timeseries.{plot_suffix}",
        units = "build/output/{resolution}/timeseries-selection-{plot_suffix}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/timeseries.py"


rule plot_sobol_indices:
    message: "Create heatmap of {wildcards.order} Sobol indices of all input/output combinations."
    input:
        src = "src/analyse/sobol.py",
        matplotlibrc = "matplotlibrc",
        indices_cont_and_nat = "build/uncertainty/{}/{{order}}-sobol-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        indices_reg = "build/uncertainty/{}/{{order}}-sobol-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"])
    params: parameters = config["uncertainty"]["parameters"]
    wildcard_constraints:
        order = "((total)|(first)|(total-minus-first))",
        extent = "((all)|(diff))"
    output: "build/output/{resolution}/{order}-sobol-{extent}.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/sobol.py"


rule plot_generation_shares:
    message: "Create plot of generation shares of two scenarios."
    input:
        src = "src/analyse/generation_shares.py",
        matplotlibrc = "matplotlibrc",
        results = rules.time_aggregated_results.output[0]
    params:
        scenario1 = "continental-autarky-100-continental-grid",
        scenario2 = "regional-autarky-100-regional-grid"
    output: "build/output/{resolution}/generation-shares.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/generation_shares.py"


rule plot_use_of_bioenergy:
    message: "Create plot of use of bioenergy potentials."
    input:
        src = "src/analyse/bioenergy_use.py",
        matplotlibrc = "matplotlibrc",
        results = rules.time_aggregated_results.output[0],
        potentials = eurocalliope(
            "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv"
            .format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
        )
    params:
        scenario = "regional-autarky-100-regional-grid",
        efficiency = config["parameters"]["biofuel-efficiency"]
    output: "build/output/{resolution}/bioenergy-use.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/bioenergy_use.py"


rule overview_scenario_results:
    message: "Create table of key outputs of scenarios."
    input:
        src = "src/analyse/scenario_overview.py",
        results = rules.aggregated_results.output[0],
    output:
        table1 = "build/output/{resolution}/overview-scenario-results-1.csv",
        table2 = "build/output/{resolution}/overview-scenario-results-2.csv"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/scenario_overview.py"


rule overview_cost_assumptions:
    message: "Create table of key cost assumptions."
    input:
        src = "src/analyse/cost_assumptions.py",
        model = "build/output/{resolution}/runs/regional-autarky-70-continental-grid.nc"
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/overview-cost-assumptions.csv"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/cost_assumptions.py"
