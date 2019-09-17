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
EDGE_WIDTH = 0.06
EDGE_COLOR = "white"

GREEN_PALETTE = sns.light_palette(GREEN, n_colors=6, reverse=True)
RED_PALETTE = sns.light_palette(RED, n_colors=6, reverse=False)
COLORS = [GREEN_PALETTE[1], GREEN_PALETTE[4], RED_PALETTE[1], RED_PALETTE[4]]

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

SCENARIO = "continental-autarky-100-continental-grid"
LAYER = "Continental"
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
    name: str
    bins: list = field(default_factory=list)
    labels: list = field(default_factory=list)


def create_map(path_to_shapes, connected_regions, path_to_aggregated_results, path_to_plot):
    plot_datas = prepare_data(path_to_shapes, path_to_aggregated_results, connected_regions)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def prepare_data(path_to_shapes, path_to_aggregated_results, connected_regions):
    shapes = (
        gpd
        .read_file(path_to_shapes)
        .to_crs(EPSG_3035_PROJ4)
        .set_index("id")
        .rename(index=lambda idx: idx.replace(".", "-"))
    )
    results = xr.open_dataset(path_to_aggregated_results).sel(scenario=SCENARIO)
    rel_gen = (
        postprocess_combined_regions(results.carrier_prod.sel(techs=SUPPLY_TECHS).sum("techs"), connected_regions)
        / postprocess_combined_regions(-results.carrier_con.sel(techs=DEMAND_TECH), connected_regions)
    )
    shapes["rel_gen"] = rel_gen.to_series()
    rel_trans = (
        postprocess_combined_regions(results.carrier_con.sel(techs=TRANSMISSION_TECH), connected_regions)
        / postprocess_combined_regions(results.carrier_con.sel(techs=DEMAND_TECH), connected_regions)
    )
    shapes["rel_trans"] = rel_trans.to_series()
    return [
        PlotData(
            shapes=shapes,
            column="rel_gen",
            panel_id="a",
            name="Generation",
            bins=[0.1, 1, 50],
            labels=["≤ 0.1", "0.1 - 3", "3 - 50", "≥ 50"]
        ),
        PlotData(
            shapes=shapes,
            column="rel_trans",
            panel_id="b",
            name="Transmission",
            bins=[1, 3, 50],
            labels=["≤ 1", "1 - 3", "3 - 50", "≥ 50"]
        ),
    ]


def plot_data(plot_datas):
    assert len(plot_datas) == 2
    fig = plt.figure(figsize=(8, 4.5))
    axes = fig.subplots(
        nrows=2,
        ncols=2,
        sharex=False,
        sharey=False,
        gridspec_kw={'height_ratios': [4, 0.5], "wspace": 0, "hspace": 0}
    )

    for i, plot_data in enumerate(plot_datas):
        plot_map(axes[0][i], plot_data)
        plot_legend(axes[1][i], plot_data)
    fig.tight_layout()
    return fig


def plot_map(ax, plot_data):
    ax.set_aspect('equal')
    plot_data.shapes.plot(
        linewidth=EDGE_WIDTH,
        edgecolor=EDGE_COLOR,
        column=plot_data.column,
        scheme="userdefined",
        classification_kwds={"bins": plot_data.bins},
        cmap=matplotlib.colors.ListedColormap(sns.color_palette(COLORS).as_hex()),
        legend=False,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.13, 0.95],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(LAYER, xy=[0.2, 0.95], xycoords='axes fraction')
    ax.annotate(plot_data.name, xy=[0.2, 0.90], xycoords='axes fraction')
    ax.annotate(plot_data.panel_id, xy=[0.05, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def plot_legend(ax, plot_data):
    patches = [
        Line2D([0], [0], linestyle="none", marker="o", markersize=9, color=color, label=label)
        for color, label in zip(COLORS, plot_data.labels)
    ]
    leg = ax.legend(
        handles=patches,
        loc="center",
        bbox_to_anchor=[0.45, 0.5],
        ncol=4,
        handletextpad=0.2,
        columnspacing=1
    )
    leg.draw_frame(False)
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    ax.set_xticks([])
    ax.set_yticks([])


def postprocess_combined_regions(da, region_pairs):
    region_pairs = list(region_pairs)
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
        path_to_plot=snakemake.output[0]
    )
