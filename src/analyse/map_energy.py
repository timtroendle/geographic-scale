from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path

import xarray as xr
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
LAYER_FONT_WEIGHT = "medium"
EDGE_WIDTH = 0.06
EDGE_COLOR = "white"

GREEN_PALETTE = sns.light_palette(GREEN, n_colors=6, reverse=True)
BLUE_PALETTE = sns.light_palette(BLUE, n_colors=6, reverse=True)
RED_PALETTE = sns.light_palette(RED, n_colors=6, reverse=False)
RED_TO_BLUE = [ # https://vis4.net/palettes/#/9|d|002d6e,a7bffa,dedfa0|dedfa0,fdad97,720000|1|0
    '#002d6e', '#2f5696', '#5f80be', '#92ace8', '#dedfa0',
    '#ed9985', '#c6685b', '#9f3831', '#720000'
]
RDBU_PALETTE = matplotlib.colors.LinearSegmentedColormap.from_list("signature-BlRd", RED_TO_BLUE)
GENERATION_COLORS = [RDBU_PALETTE(0.1), RDBU_PALETTE(0.4), RDBU_PALETTE(0.5), RDBU_PALETTE(0.6), RDBU_PALETTE(0.9)]
TRANSMISSION_COLORS = sns.light_palette(RED, n_colors=4)

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

SUPPLY_TECHS = [
    'wind_onshore_monopoly', 'biofuel',
    'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
    'wind_onshore_competing', 'open_field_pv', 'hydro_run_of_river'
]
DEMAND_TECH = "demand_elec"
TRANSMISSION_TECH = "ac_transmission"


@dataclass
class PlotData:
    shapes: gpd.GeoDataFrame
    column: str
    panel_id: str
    scale: str
    name: str
    legend_title: str
    bins: list = field(default_factory=list)
    labels: list = field(default_factory=list)
    colors: list = field(default_factory=list)


def create_map(path_to_shapes, connected_regions, scenarios, path_to_aggregated_results, path_to_plot):
    plot_datas = prepare_data(path_to_shapes, path_to_aggregated_results, connected_regions, scenarios)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=300, pil_kwargs={"compression": "tiff_lzw"})


def prepare_data(path_to_shapes, path_to_aggregated_results, connected_regions, scenarios):
    shapes = (
        gpd
        .read_file(path_to_shapes)
        .to_crs(EPSG_3035_PROJ4)
        .set_index("id")
        .rename(index=lambda idx: idx.replace(".", "-"))
    )
    results = xr.open_dataset(path_to_aggregated_results).sel(scenario=scenarios)
    rel_gen = (
        postprocess_combined_regions(results.carrier_prod.sel(techs=SUPPLY_TECHS).sum("techs"), connected_regions)
        / postprocess_combined_regions(-results.carrier_con.sel(techs=DEMAND_TECH), connected_regions)
    )
    rel_trans = (
        postprocess_combined_regions(results.carrier_con.sel(techs=TRANSMISSION_TECH), connected_regions)
        / postprocess_combined_regions(results.carrier_con.sel(techs=DEMAND_TECH), connected_regions)
    )
    letters = list("ABCDEF")
    generation = [
        PlotData(
            shapes=shapes.assign(rel_gen=rel_gen.sel(scenario=scenario).to_series()),
            column="rel_gen",
            panel_id=letters.pop(0),
            scale=scale(scenario),
            name="Generation",
            legend_title="Generation relative to demand",
            bins=[0.1, 0.75, 1.25, 10],
            labels=["≤ 0.1", "0.1 - 0.75", "0.75 - 1.25", "1.25 - 10", "≥ 10"],
            colors=GENERATION_COLORS
        )
        for scenario in scenarios
    ]
    transmission = [
        PlotData(
            shapes=shapes.assign(rel_trans=rel_trans.sel(scenario=scenario).to_series()),
            column="rel_trans",
            panel_id=letters.pop(0),
            scale=scale(scenario),
            name="Transmission",
            legend_title="Transmission relative to demand",
            bins=[1, 5, 15],
            labels=["≤ 1", "1 - 5", "5 - 15", "≥ 15"],
            colors=TRANSMISSION_COLORS
        )
        for scenario in scenarios
    ]
    return generation + transmission


def scale(scenario_name):
    return scenario_name.split("-")[-2].capitalize()


def plot_data(plot_datas):
    assert len(plot_datas) == 6
    fig = plt.figure(figsize=(8, 5.5))
    n_maps_per_row = 3
    axes = fig.subplots(
        nrows=2,
        ncols=n_maps_per_row + 1,
        sharex=False,
        sharey=False,
        gridspec_kw={'width_ratios': [2.5] * n_maps_per_row + [0.5], "wspace": 0, "hspace": 0}
    )

    for i in [0, 1, 2]:
        plot_map(axes[0][i], plot_datas[i])
        plot_map(axes[1][i], plot_datas[i + 3])
    plot_legend(axes[0][3], plot_datas[0])
    plot_legend(axes[1][3], plot_datas[3])

    axes[1][0].annotate(
        "Continental supply",
        xy=[0.5, -0.1],
        xycoords='axes fraction',
        ha='center',
        va='center',
        weight=LAYER_FONT_WEIGHT
    )
    axes[1][1].annotate(
        "National supply",
        xy=[0.5, -0.1],
        xycoords='axes fraction',
        ha='center',
        va='center',
        weight=LAYER_FONT_WEIGHT
    )
    axes[1][2].annotate(
        "Regional supply",
        xy=[0.5, -0.1],
        xycoords='axes fraction',
        ha='center',
        va='center',
        weight=LAYER_FONT_WEIGHT
    )
    fig.tight_layout()
    return fig


def plot_map(ax, plot_data):
    ax.set_aspect('equal')
    colors = plot_data.colors.copy()
    bins = plot_data.bins
    while len(bins) > 0 and plot_data.shapes[plot_data.column].min() > bins[0]:
        colors.pop(0)
        bins.pop(0)
    while len(bins) > 0 and plot_data.shapes[plot_data.column].max() < bins[-1]:
        colors.pop()
        bins.pop()
    if len(bins) == 0:
        bins.append(1)
        colors.append(RED)
    plot_data.shapes.plot(
        linewidth=EDGE_WIDTH,
        edgecolor=EDGE_COLOR,
        column=plot_data.column,
        scheme="userdefined",
        classification_kwds={"bins": plot_data.bins},
        cmap=matplotlib.colors.ListedColormap(sns.color_palette(colors).as_hex()),
        legend=False,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(plot_data.panel_id, xy=[0.1, 0.9], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def plot_legend(ax, plot_data):
    patches = [
        Line2D([0], [0], linestyle="none", marker="o", markersize=9, color=color, label=label)
        for color, label in zip(plot_data.colors, plot_data.labels)
    ]
    leg = ax.legend(
        handles=patches,
        loc="center",
        bbox_to_anchor=[0.45, 0.5],
        ncol=1,
        handletextpad=0.2,
        columnspacing=1,
    )
    leg.draw_frame(False)
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.annotate(
        plot_data.legend_title,
        xy=[2, 0.50],
        xycoords='axes fraction',
        ha='center',
        va='center',
        rotation="vertical"
    )


def postprocess_combined_regions(da, region_pairs):
    region_pairs = list(region_pairs.items())
    regions = list(chain(*region_pairs))
    if all([region in da.locs for region in regions]):
        # config is valid, and resolution is regional. apply config
        for region1, region2 in region_pairs:
            sum_data = da.sel(locs=region1) + da.sel(locs=region2)
            da.loc[dict(locs=region1)] = sum_data
            da.loc[dict(locs=region2)] = sum_data
        return da
    elif all([region not in da.locs for region in regions]):
        # config is valid, but resolution is not regional. do nothing
        return da
    else:
        raise ValueError("Config of connected regions is invalid.")


if __name__ == "__main__":
    create_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_aggregated_results=snakemake.input.results,
        connected_regions=snakemake.params.connected_regions,
        scenarios=snakemake.params.scenarios,
        path_to_plot=snakemake.output[0]
    )
