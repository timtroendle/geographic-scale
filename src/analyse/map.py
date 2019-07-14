from pathlib import Path

import calliope
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import shapely

RED = "#A01914"
BLUE = "#4F6DB8"
CMAP = sns.light_palette(sns.desaturate(RED, 0.85), reverse=False, as_cmap=True)
NORM = matplotlib.colors.Normalize(vmin=0, vmax=3)

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"

ALL_TECHS = ['wind_onshore_monopoly', 'demand_elec', 'load_shedding', 'biofuel',
             'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
             'wind_onshore_competing', 'battery', 'hydrogen', 'free_transmission',
             'ac_transmission', 'open_field_pv', 'hydro_run_of_river',
             'pumped_hydro']
SUPPLY_TECHS_WITHOUT_LOAD_SHEDDING = ['wind_onshore_monopoly', 'biofuel',
                                      'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
                                      'wind_onshore_competing', 'open_field_pv', 'hydro_run_of_river']


def plot_map(path_to_continental_shape, path_to_national_shapes, path_to_regional_shapes,
             path_to_continental_result, path_to_national_result, path_to_regional_result,
             path_to_shapes, path_to_plot):
    """Plot maps of results."""
    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    axes = fig.subplots(2, 2).flatten()
    shapes = (gpd.read_file(path_to_shapes)
                 .to_crs(EPSG_3035_PROJ4)
                 .rename(columns={"id": "locs"})
                 .set_index("locs")
                 .rename(index=lambda idx: idx.replace(".", "-")))

    continental_model = calliope.read_netcdf(path_to_continental_result)
    national_model = calliope.read_netcdf(path_to_national_result)
    regional_model = calliope.read_netcdf(path_to_regional_result)

    continental = gpd.read_file(path_to_continental_shape).to_crs(EPSG_3035_PROJ4).set_index("id")
    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("country_code")
    regional = (gpd.read_file(path_to_regional_shapes)
                   .to_crs(EPSG_3035_PROJ4)
                   .set_index("id")
                   .rename(index=lambda idx: idx.replace(".", "-")))

    continental["cost"] = _read_cost(shapes, continental_model, "continental")
    national["cost"] = _read_cost(shapes, national_model, "national")
    regional["cost"] = _read_cost(shapes, regional_model, "regional")

    base_costs = continental["cost"].iloc[0]
    continental["cost"] = continental["cost"] / base_costs
    national["cost"] = national["cost"] / base_costs
    regional["cost"] = regional["cost"] / base_costs

    total_costs_continental = _read_total_system_costs(continental_model) / base_costs
    total_costs_national = _read_total_system_costs(national_model) / base_costs
    total_costs_regional = _read_total_system_costs(regional_model) / base_costs

    network_continental = _read_network_graph(shapes, continental_model)
    network_national = _read_network_graph(shapes, national_model)
    network_regional = _read_network_graph(shapes, regional_model)

    _plot_layer(continental, network_continental, total_costs_continental, "continental", axes[0])
    _plot_layer(national, network_national, total_costs_national, "national", axes[2])
    _plot_layer(regional, network_regional, total_costs_regional, "regional", axes[3])
    sns.despine(ax=axes[1], top=True, bottom=True, left=True, right=True)
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    _plot_colorbar(fig, axes)
    fig.savefig(path_to_plot, dpi=300, transparent=True)


def _plot_layer(units, network_graph, total_cost, layer_name, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
        column="cost",
        vmin=NORM.vmin,
        vmax=NORM.vmax,
        cmap=CMAP,
        ax=ax
    )
    gpd.GeoSeries([
        network_graph.edges[region, neighbour]["line"]
        for region, neighbour in network_graph.edges
    ], crs=units.crs).plot(ax=ax, linewidth=0.25, color=BLUE)
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.10, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{MONEY_UNICODE} ",
        xy=[0.10, 0.85],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color=sns.desaturate(RED, 0.85)
    )
    ax.annotate(layer_name, xy=[0.17, 0.90], xycoords='axes fraction')
    ax.annotate(f"{total_cost:.1f}", xy=[0.17, 0.85], xycoords='axes fraction')


def _plot_colorbar(fig, axes):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(["{:.0f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)


def _read_cost(shapes, model, level):
    assert level in ["continental", "national", "regional"]
    assert set(model.inputs.techs.values) == set(ALL_TECHS)
    cost = model.get_formatted_array('cost').squeeze('costs').sum("techs")
    carrier_prod = (model
                    .get_formatted_array('carrier_prod')
                    .squeeze(['carriers']))
    carrier_prod = (carrier_prod
                    .sel(techs=_supply_techs(carrier_prod))
                    .sum(dim=['techs', 'timesteps']))
    if level == "continental":
        cost = cost.sum(dim="locs").item()
        carrier_prod = carrier_prod.sum(dim="locs").item()
    elif level == "national":
        cost = (cost.groupby(shapes.country_code.to_xarray())
                    .sum("locs")
                    .to_series())
        carrier_prod = (carrier_prod.groupby(shapes.country_code.to_xarray())
                                    .sum("locs")
                                    .to_series())
    elif level == "regional":
        cost = cost.to_series()
        carrier_prod.to_series()
    return (cost / carrier_prod)


def _supply_techs(carrier_prod):
    if "load_shedding" in carrier_prod.techs:
        return SUPPLY_TECHS_WITHOUT_LOAD_SHEDDING + ["load_shedding"]
    else:
        return SUPPLY_TECHS_WITHOUT_LOAD_SHEDDING


def _read_total_system_costs(model):
    return model.get_formatted_array("total_levelised_cost").item()


def _read_network_graph(shapes, model):
    g = nx.Graph()
    g.add_nodes_from((k, {"centroid": v}) for (k, v) in shapes.centroid.to_dict().items())

    try:
        lines = [(line.split(":")[0], line.split(":")[-1])
                 for line in model._model_data.loc_techs_transmission.values]
        lines = [(a, b, {"line": shapely.geometry.LineString([g.nodes[a]["centroid"], g.nodes[b]["centroid"]])})
                 for a, b in lines]
        g.add_edges_from(lines)
    except AttributeError:
        pass # there are no lines
    return g


if __name__ == "__main__":
    plot_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_continental_shape=snakemake.input.continental_shape,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_continental_result=snakemake.input.continental_result,
        path_to_national_result=snakemake.input.national_result,
        path_to_regional_result=snakemake.input.regional_result,
        path_to_plot=snakemake.output[0]
    )
