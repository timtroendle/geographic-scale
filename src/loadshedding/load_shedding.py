from pathlib import Path

import jinja2

TEMPLATE = """overrides:
    load-shedding:
        locations:
            {% for id, allowance in locations.items() %}
            {{ id | replace(".", "-") }}.techs.load_shedding.exists: {{ allowance }}
            {% endfor %}
"""


def load_shedding(paths_to_feasibilities, path_to_result):
    load_shedding_at_location = {
        Path(path_to_feasibility).stem: infeasible(path_to_feasibility)
        for path_to_feasibility in paths_to_feasibilities
    }
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        locations=load_shedding_at_location
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


def infeasible(path_to_feasibility):
    with open(path_to_feasibility, "r") as f_feasibility:
        feasibility = f_feasibility.readlines()[0]
    return feasibility == "infeasibleOrUnbounded"


if __name__ == "__main__":
    load_shedding(
        paths_to_feasibilities=snakemake.input.feasibilities,
        path_to_result=snakemake.output[0]
    )
