rule all:
    message: "Run all steps to create concept paper."
    input:
        "build/concept/concept.html"


rule scenario_space:
    message: "Plot sketch of scenario space."
    input: src="src/concept/scenarios.py"
    output: "build/concept/scenario-space.png"
    conda: "../envs/default.yaml"
    script: "../src/concept/scenarios.py"


rule map:
    message: "Create map of results."
    input:
        src = "src/concept/map.py",
        continental_shape = "data/units-continental.geojson",
        national_shapes = "data/units-national.geojson",
        regional_shapes = "data/units-regional.geojson"
    output: "build/concept/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/concept/map.py"


rule sensitivity_heatmap:
    message: "Plot sketch of sensitivity heatmap."
    input: src="src/concept/sensitivity.py"
    output: "build/concept/parameter-sensitivity.png"
    conda: "../envs/default.yaml"
    script: "../src/concept/sensitivity.py"


rule concept:
    message: "Compile concept paper."
    input:
        "report/literature.bib",
        "report/concept.md",
        "report/pandoc-metadata.yml",
        "report/report.css",
        rules.scenario_space.output,
        rules.sensitivity_heatmap.output,
        rules.map.output
    output:
        "build/concept/concept.html"
    conda: "../envs/pdf.yaml"
    shell:
        """
        cd ./report
        {PANDOC} --self-contained --css=report.css concept.md pandoc-metadata.yml \
        --to html5 -o ../build/concept/concept.html
        """


rule pdf_concept:
    message: "Compile PDF concept."
    input: rules.concept.output # to avoid repeating all dependencies
    output: "build/concept/concept.pdf"
    conda: "../envs/pdf.yaml"
    shell:
        """
        cd ./report
        {PANDOC} --css=report.css concept.md pandoc-metadata.yml --pdf-engine weasyprint \
        -o ../build/concept/concept.pdf
        """
