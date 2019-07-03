import functools
from pathlib import Path

import calliope
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import shapely
import seaborn as sns
import networkx as nx

GREY = "#EBEBEB"
EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf073"
VALUE_UNICODE = "\uf1e6"


def visualise_link_usage(path_to_units, path_to_countries, path_to_results, path_to_plot, scaling_factor):
    model = calliope.read_netcdf(path_to_results)
    units = gpd.read_file(path_to_units).set_index("id").to_crs(EPSG_3035_PROJ4)
    units.index = units.index.str.replace(".", "-")
    countries = gpd.read_file(path_to_countries).set_index("id").to_crs(EPSG_3035_PROJ4)

    graph = create_graph(model, countries, units, scaling_factor)

    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()

    plot_imports(graph, countries, axes[0])
    plot_balancing(graph, countries, axes[1])

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
        ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
        ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)

    axes[0].get_legend().remove()
    legend = axes[1].get_legend()
    for text in legend.get_texts():
        text.set_text(relabel(text.get_text()))
    fig.savefig(path_to_plot, dpi=300, transparent=False)


def create_graph(model, countries, units, scaling_factor):
    graph = nx.Graph()
    nodes = [(country_code, {"centroid": centroid}) for country_code, centroid in countries.centroid.to_dict().items()]
    graph.add_nodes_from(nodes)
    edges = [(line.split(":")[0], line.split(":")[-1]) for line in model.inputs.loc_techs_transmission.values]
    edges = filter( # transform to linkes between countries
        lambda x: x[0] != x[1],
        map(
            lambda x: (units.loc[x[0], "country_code"], units.loc[x[1], "country_code"]),
            edges
        )
    )
    edges = list(set(tuple(sorted((a, b))) for a, b in edges)) # filter duplicates
    edges = [
        (source, sink, {
            "timeseries":
                {source:
                    {sink: link_timeseries(model, units, source, sink, scaling_factor)},
                 sink:
                    {source: link_timeseries(model, units, sink, source, scaling_factor)}
                 }
        })
        for (source, sink) in edges
    ]
    for edge in edges:
        edge[2]["line"] = shapely.geometry.LineString([
            graph.nodes[edge[0]]["centroid"],
            graph.nodes[edge[1]]["centroid"]
        ])
    graph.add_edges_from(edges)
    return graph


def link_timeseries(model, units, country_code_a, country_code_b, scaling_factor):
    gen = transmission_carrier_prod(model, scaling_factor)
    nodes_a = units[units.country_code == country_code_a].index
    nodes_b = units[units.country_code == country_code_b].index
    a_to_b = gen.sel(source=nodes_a, sink=nodes_b).sum(["sink", "source"])
    b_to_a = gen.sel(source=nodes_b, sink=nodes_a).sum(["sink", "source"])
    return a_to_b - b_to_a


@functools.lru_cache(maxsize=1, typed=False)
def transmission_carrier_prod(model, scaling_factor):
    gen = model.get_formatted_array("carrier_prod").squeeze("carriers")
    gen = gen.sel(techs=gen.techs.str.contains("ac_transmission"))
    return (gen.assign_coords(techs=[tech.item().split(":")[-1] for tech in gen.techs])
               .rename(locs="sink", techs="source")) / scaling_factor


def plot_imports(graph, countries, ax):
    line_use = gpd.GeoDataFrame(
        geometry=[
            graph.edges[region, neighbour]["line"]
            for region, neighbour in graph.edges
        ],
        data=[
            abs(graph.edges[region, neighbour]["timeseries"][region][neighbour].sum("timesteps"))
            for region, neighbour in graph.edges
        ],
        crs=EPSG_3035_PROJ4
    )
    plot_line_use(countries, line_use, ax=ax, plot_name="annual")


def plot_balancing(graph, countries, ax):
    line_use = gpd.GeoDataFrame(
        geometry=[
            graph.edges[region, neighbour]["line"]
            for region, neighbour in graph.edges
        ],
        data=[
            abs(graph.edges[region, neighbour]["timeseries"][region][neighbour]).sum("timesteps")
            - abs(graph.edges[region, neighbour]["timeseries"][region][neighbour].sum("timesteps"))
            for region, neighbour in graph.edges
        ],
        crs=EPSG_3035_PROJ4
    )
    plot_line_use(countries, line_use, ax=ax, plot_name="inter-annual")


def plot_line_use(countries, line_use, plot_name, ax):
    countries.plot(ax=ax, facecolor=GREY, edgecolor="w", linewidth=0.1)
    line_use[0] = line_use[0] / line_use[0].sum() * 100
    line_use.plot(
        ax=ax,
        cmap="Blues",
        linewidth=2,
        column=0,
        legend=True,
        scheme="userdefined",
        classification_kwds={"bins": [1, 3, 8, 18]}
    )
    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.10, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{VALUE_UNICODE} ",
        xy=[0.10, 0.85],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(plot_name, xy=[0.17, 0.90], xycoords='axes fraction')
    total_line_use_twh = line_use.sum().item() / 1e6
    ax.annotate(f"{total_line_use_twh:.0f} TWh", xy=[0.17, 0.85], xycoords='axes fraction')


def relabel(label):
    number1, number2 = label.split("-")
    return f"{float(number1):.0f}-{float(number2):.0f}%"


if __name__ == "__main__":
    visualise_link_usage(
        path_to_units=snakemake.input.units,
        path_to_countries=snakemake.input.countries,
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0],
        scaling_factor=snakemake.params.scaling_factor
    )
