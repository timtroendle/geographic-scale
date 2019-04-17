ruleorder: performance_run > run_national > run_regional


rule performance_run:
    message: "Run {wildcards.run_id} with feas {wildcards.feasibility_tol} and opt {wildcards.optimality_tol}."
    input: model = "build/model/national/model.yaml"
    params: scenario = "national-autarky-100-continental-grid"
    conda: "../envs/calliope.yaml"
    log: "build/logs/performance/feas_{feasibility_tol}-opt_{optimality_tol}/{run_id}.log"
    shell:
        """
        calliope run {input.model} --scenario={params.scenario}\
        --override_dict="{{run.solver_options.FeasibilityTol: {wildcards.feasibility_tol}, \
        run.solver_options.OptimalityTol: {wildcards.optimality_tol}}}" 2> {log}\

        """


rule performance_analysis:
    message: "Analyse performance runs."
    input:
        src = "src/analyse/performance.py",
        logs = expand(
            "build/logs/performance/feas_{feasibility_tol}-opt_{optimality_tol}/{run_id}.log",
            feasibility_tol=["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
            optimality_tol=["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
            run_id=range(2)
        )
    output: "build/logs/performance/analysis.log"
    conda: "../envs/default.yaml"
    script: "../src/analyse/performance.py"
