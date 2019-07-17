from dataclasses import dataclass

import yaml
import jinja2
import geopandas as gpd

TEMPLATE = """
overrides:
    {% for grid_size, links in restrictions.items() %}
    {{ grid_size }}:
        {% for link in links %}
        links.{{ link.location_A }},{{ link.location_B }}.exists: {{ link.is_allowed }}
        {% endfor %}
    {% endfor %}
"""


@dataclass
class Link:
    location_A: str
    location_B: str
    is_allowed: bool = True


def grid_size_restriction(path_to_units, path_to_links, connected_regions, path_to_result):
    units = gpd.read_file(path_to_units).set_index("id", drop=True)
    units.index = units.index.map(lambda unit_id: unit_id.replace(".", "-"))
    restrictions = {
        "continental-grid-size": _continental_grid_size(units, links=_read_links(path_to_links)),
        "national-grid-size": _national_grid_size(units, links=_read_links(path_to_links)),
        "regional-grid-size": _regional_grid_size(units, connected_regions, links=_read_links(path_to_links))
    }

    restrictions = jinja2.Template(TEMPLATE).render(
        restrictions=restrictions
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(restrictions)


def _regional_grid_size(units, connected_regions, links):
    # assuming units are not smaller than regions, then units can generally not be connected
    for link in links:
        if link.location_A in connected_regions.keys() and connected_regions[link.location_A] == link.location_B:
            link.is_allowed = True
        elif link.location_B in connected_regions.keys() and connected_regions[link.location_B] == link.location_A:
            link.is_allowed = True
        else:
            link.is_allowed = False
    return links


def _national_grid_size(units, links):
    # units can be connected iff they are in the same country
    for link in links:
        link.is_allowed = units.loc[link.location_A, "country_code"] == units.loc[link.location_B, "country_code"]
    return links


def _continental_grid_size(units, links):
    # all units can be connected
    for link in links:
        link.is_allowed = True
    return links


def _read_links(path_to_links):
    with open(path_to_links, "r") as f_links:
        try:
            links = yaml.safe_load(f_links)["links"].keys()
        except AttributeError:
            return [] # on continental scale, no links exist
    return [Link(location_A=location_pair[0], location_B=location_pair[1])
            for location_pair in (link.split(",") for link in links)]


if __name__ == "__main__":
    grid_size_restriction(
        path_to_units=snakemake.input.units,
        path_to_links=snakemake.input.links,
        connected_regions=snakemake.params.connected_regions,
        path_to_result=snakemake.output[0]
    )
