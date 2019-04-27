ruleorder: performance_run > run_national > run_regional


rule performance_run:
    message: "Run {wildcards.run_id} with {wildcards.threads} threads."
    input: model = "build/model/national/model.yaml"
    params: scenario = "national-autarky-100-continental-grid"
    conda: "../envs/calliope.yaml"
    log: "build/logs/performance/{threads}-threads/{run_id}.log"
    shell:
        """
        calliope run {input.model} --scenario={params.scenario}\
        --override_dict="{{run.solver_options.Threads: {wildcards.threads}}}" 2> {log}\
        """


rule performance_analysis:
    message: "Analyse performance runs."
    input:
        src = "src/analyse/performance.py",
        logs = expand(
            "build/logs/performance/{threads}-threads/{run_id}.log",
            threads=[1, 2, 4, 8, 12, 16],
            run_id=range(3)
        )
    output: "build/logs/performance/analysis.log"
    conda: "../envs/default.yaml"
    script: "../src/analyse/performance.py"
