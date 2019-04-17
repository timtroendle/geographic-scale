import jinja2
import geopandas as gpd

TEMPLATE = """
overrides:
    {% for restriction in restriction_levels %}
    {% for autarky_level, location_groups in autarky_levels.items() %}
    {{ autarky_level }}-{{ 100 - restriction }}-percent:
        {% for location_group, sublocations in location_groups.items() %}
        group_constraints.import_restriction_{{ restriction }}_percent_{{ location_group | replace(".", "-") }}:
            locs: [{{ sublocations | map("replace", ".", "-") | join(",") }}]
            techs: ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore"]
            demand_share_min:
                electricity: {{ 1 - (restriction / 100) }}
        {% endfor %}
    {% endfor %}
    {% endfor %}
"""
COUNTRY_CODE_COLUMN = "country_code"


def import_restriction(path_to_units, restrictions_in_percent, path_to_result):
    units = gpd.read_file(path_to_units).set_index("id", drop=True)

    autarky_levels = {
        "regional-autarky": _regional_autarky(units),
        "national-autarky": _national_autarky(units),
        "continental-autarky": _continental_autarky(units)
    }

    restrictions = jinja2.Template(TEMPLATE).render(
        autarky_levels=autarky_levels,
        restriction_levels=restrictions_in_percent
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(restrictions)


def _regional_autarky(units):
    # assume units are regions or larger
    return {unit_id: [unit_id] for unit_id, unit in units.iterrows()}


def _national_autarky(units):
    return {
        country_code: list(countries.index.values)
        for country_code, countries in units.groupby(COUNTRY_CODE_COLUMN)
    }


def _continental_autarky(units):
    return {
        "continental": list(units.index.values)
    }


if __name__ == "__main__":
    import_restriction(
        path_to_units=snakemake.input.units,
        restrictions_in_percent=snakemake.params.restrictions,
        path_to_result=snakemake.output[0]
    )
