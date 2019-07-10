import jinja2

TEMPLATE = """overrides:
    load-shedding:
        locations:
            {% for id in locations %}
            {{ id | replace(".", "-") }}.techs.load_shedding.exists: True
            {% endfor %}
"""


def load_shedding(loadshedding_locations, path_to_result):
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        locations=loadshedding_locations
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    load_shedding(
        loadshedding_locations=snakemake.params.locations,
        path_to_result=snakemake.output[0]
    )
