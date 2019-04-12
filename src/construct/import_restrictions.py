import yaml
import jinja2

TEMPLATE = """
overrides:
    {% for restriction in restrictions %}
    import_restriction_{{ restriction }}_percent:
        {% for location in locations %}
        group_constraints.import_restriction_{{ restriction }}_percent_{{ location | replace(".", "-") }}:
            locs: ["{{ location | replace(".", "-") }}"]
            techs: ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore"]
            demand_share_min:
                electricity: {{ 1 - (restriction / 100) }}
        {% endfor %}
    {% endfor %}
"""


def import_restriction(path_to_locations, restrictions_in_percent, path_to_result):
    with open(path_to_locations, "r") as f_locations:
        locations = yaml.safe_load(f_locations)["locations"].keys()
    restrictions = jinja2.Template(TEMPLATE).render(
        locations=locations,
        restrictions=restrictions_in_percent
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(restrictions)


if __name__ == "__main__":
    import_restriction(
        path_to_locations=snakemake.input.locations,
        restrictions_in_percent=snakemake.params.restrictions,
        path_to_result=snakemake.output[0]
    )
