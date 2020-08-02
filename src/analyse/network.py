import calliope
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
import seaborn as sns
import networkx as nx

GREY = "#EBEBEB"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
LINEWIDTH = 1.5
LINECOLOR = BLUE
INTERNATIONAL_LINE_COLOUR = YELLOW
EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000


def visualise_links(path_to_units, path_to_results, path_to_plot):
    model = calliope.read_netcdf(path_to_results)
    units = (
        gpd
        .read_file(path_to_units)
        .set_index("id")
        .rename(index=lambda idx: idx.replace(".", "-"))
        .to_crs(EPSG_3035_PROJ4)
    )
    network_graph = read_network_graph(units, model)

    fig = plt.figure(figsize=(6.67, 6.67), constrained_layout=True)
    ax = fig.subplots(1, 1)

    plot_network(units, network_graph, ax)

    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)

    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def read_network_graph(shapes, model):
    g = nx.Graph()
    g.add_nodes_from((k, {"centroid": v}) for (k, v) in shapes.centroid.to_dict().items())

    try:
        lines = [(line.split(":")[0], line.split(":")[-1])
                 for line in model._model_data.loc_techs_transmission.values
                 if model.inputs.energy_cap_max.sel(loc_techs=line) > 0]
        lines = [(a, b, {"line": shapely.geometry.LineString([g.nodes[a]["centroid"], g.nodes[b]["centroid"]])})
                 for a, b in lines]
        g.add_edges_from(lines)
    except AttributeError:
        pass # there are no lines
    return g


def plot_network(units, network_graph, ax):
    units.plot(ax=ax, facecolor=GREY, edgecolor="w", linewidth=0.1)
    gpd.GeoSeries([
        network_graph.edges[region, neighbour]["line"]
        for region, neighbour in network_graph.edges
    ], crs=units.crs).plot(ax=ax, linewidth=LINEWIDTH, color=LINECOLOR)
    gpd.GeoSeries([
        network_graph.edges[region, neighbour]["line"]
        for region, neighbour in network_graph.edges if region[:3] != neighbour[:3] # hacky way to detect country code
    ], crs=units.crs).plot(ax=ax, linewidth=LINEWIDTH, color=INTERNATIONAL_LINE_COLOUR)


if __name__ == "__main__":
    visualise_links(
        path_to_units=snakemake.input.units,
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
