import jinja2
import geopandas as gpd

TEMPLATE = """
overrides:
    {%- for restriction in restriction_levels %}
    {%- for autarky_level, location_groups in autarky_levels.items() %}
    {{ autarky_level }}-{{ 100 - restriction }}-percent:
        {%- for location_group, sublocations in location_groups.items() %}
        group_constraints.import_restriction_{{ restriction }}_percent_{{ location_group }}:
            locs:
            {%- for sublocation in sublocations %}
            - {{ sublocation }}
            {%- endfor %}
            net_import_share_max:
                electricity: {{ restriction / 100 }}
        {%- endfor %}
    {%- endfor %}
    {%- endfor %}
"""
COUNTRY_CODE_COLUMN = "country_code"


def import_restriction(path_to_units, restrictions_in_percent, connected_regions, path_to_result):
    units = gpd.read_file(path_to_units).set_index("id", drop=True).rename(index=lambda idx: idx.replace(".", "-"))

    autarky_levels = {
        "regional-autarky": _regional_autarky(units),
        "national-autarky": _national_autarky(units),
        "continental-autarky": _continental_autarky(units)
    }
    autarky_levels["regional-autarky"] = _connect_regions(autarky_levels["regional-autarky"], connected_regions)

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


def _connect_regions(autarky_dict, connected_regions):
    if all([region in autarky_dict.keys() for region in connected_regions]):
        # config is valid, and resolution is regional. apply config
        for insufficient, neighbour in connected_regions.items():
            autarky_dict[insufficient] = [insufficient, neighbour]
            del autarky_dict[neighbour]
        return autarky_dict
    elif all([region not in autarky_dict.keys() for region in connected_regions]):
        # config is valid, but resolution is not regional. do nothing
        return autarky_dict
    else:
        raise ValueError("Config of connected regions is invalid.")


if __name__ == "__main__":
    import_restriction(
        path_to_units=snakemake.input.units,
        restrictions_in_percent=snakemake.params.restrictions,
        connected_regions=snakemake.params.connected_regions,
        path_to_result=snakemake.output[0]
    )
