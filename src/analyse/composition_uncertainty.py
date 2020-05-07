from dataclasses import dataclass

import pandas as pd
import seaborn as sns
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from seaborn import utils
from seaborn.palettes import blend_palette

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
GREY = "#C0C0C0"
COLOR = BLUE
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
MAX_VALUE = 2
NUMBER_VARIABLES = 2

HYDRO_TECHS = ["hydro_run_of_river", "hydro_reservoir"]
PUMPED_HYDRO_TECHS = ["pumped_hydro"]
GW_TO_TW = 1e-3
MW_TO_TW = 1e-6

COLUMN_HEADER = [
    "y-continental-scale-cost-eur",
    "y-national-scale-cost-eur",
    "y-cost-diff-eur",
    "y-cost-diff-relative",
    "y-supply-diff-relative",
    "y-wind-diff-relative",
    "y-biostor-diff-relative",
    "y-continental-scale-pv-gw",
    "y-national-scale-pv-gw",
    "y-continental-scale-wind-gw",
    "y-national-scale-wind-gw",
    "y-continental-scale-biofuel-gw",
    "y-national-scale-biofuel-gw",
    "y-continental-scale-storage-gw",
    "y-national-scale-storage-gw",
    "y-continental-scale-storage-gwh",
    "y-national-scale-storage-gwh",
    "y-continental-scale-transmission-gwkm",
    "y-regional-scale-cost-eur",
    "y-regional-scale-pv-gw",
    "y-regional-scale-wind-gw",
    "y-regional-scale-biofuel-gw",
    "y-regional-scale-storage-gw",
    "y-regional-scale-storage-gwh",
    "y-regional-scale-transmission-gwkm"
]


@dataclass
class PlotData:
    panel_id: str
    title: str
    x: pd.Series
    y: pd.Series
    xlabel: str
    ylabel: str
    xlim: tuple
    ylim: tuple


def plot_composition_variability(path_to_large_scales, path_to_small_scale, path_to_scenario_results, path_to_plot):
    plot_datas = read_plot_data(path_to_large_scales, path_to_small_scale, path_to_scenario_results)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def read_plot_data(path_to_large_scales, path_to_small_scale, path_to_scenario_results):
    y = pd.concat([
        pd.read_csv(path_to_large_scales, index_col=None, header=None),
        pd.read_csv(path_to_small_scale, index_col=None, header=None),

    ], axis="columns") * GW_TO_TW
    y.columns = COLUMN_HEADER
    agg = xr.open_dataset(path_to_scenario_results)
    hydro_cap = pd.Series(
        agg.energy_cap.sel(techs=HYDRO_TECHS).sum(["techs", "locs"]).isel(scenario=0).item() * MW_TO_TW,
        index=y.index
    )
    pumped_hydro_cap = pd.Series(
        agg.energy_cap.sel(techs=PUMPED_HYDRO_TECHS).sum(["techs", "locs"]).isel(scenario=0).item() * MW_TO_TW,
        index=y.index
    )
    return [
        PlotData(
            panel_id="a",
            title="Supply",
            ylabel="National scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=(y["y-continental-scale-wind-gw"]
               + y["y-continental-scale-pv-gw"]
               + hydro_cap),
            y=(y["y-national-scale-wind-gw"]
               + y["y-national-scale-pv-gw"]
               + hydro_cap),
        ),
        PlotData(
            panel_id="b",
            title="Balancing",
            ylabel="National scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=y["y-continental-scale-biofuel-gw"] + y["y-continental-scale-storage-gw"] + pumped_hydro_cap,
            y=y["y-national-scale-biofuel-gw"] + y["y-national-scale-storage-gw"] + pumped_hydro_cap,
        ),
        PlotData(
            panel_id="c",
            title="Supply",
            ylabel="Regional scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=(y["y-continental-scale-wind-gw"]
               + y["y-continental-scale-pv-gw"]
               + y["y-continental-scale-biofuel-gw"]
               + hydro_cap),
            y=(y["y-regional-scale-wind-gw"]
               + y["y-regional-scale-pv-gw"]
               + y["y-regional-scale-biofuel-gw"]
               + hydro_cap),
        ),
        PlotData(
            panel_id="d",
            title="Balancing",
            ylabel="Regional scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=y["y-continental-scale-biofuel-gw"] + y["y-continental-scale-storage-gw"] + pumped_hydro_cap,
            y=y["y-regional-scale-biofuel-gw"] + y["y-regional-scale-storage-gw"] + pumped_hydro_cap,
        )
    ]


def plot_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 5))
    axes = fig.subplots(
        nrows=2,
        ncols=int(len(plot_datas) / 2),
        sharex=True,
        sharey=True
    )

    color_rgb = mpl.colors.colorConverter.to_rgb(COLOR)
    colors = [utils.set_hls_values(color_rgb, l=l)  # noqa
            for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    for i, plot_data in enumerate(plot_datas):
        ax = axes[i // NUMBER_VARIABLES][i % NUMBER_VARIABLES]
        ax.hexbin(
            x=plot_data.x,
            y=plot_data.y,
            gridsize=int((plot_data.x.max() - plot_data.x.min()) * 20),
            cmap=cmap
        )
        ax.set_aspect('equal')
        ax.set_ylim(*plot_data.ylim)
        ax.set_xlim(*plot_data.xlim)
        ax.plot(plot_data.xlim, plot_data.ylim, "--", color=GREY)
        if i // NUMBER_VARIABLES == 1:
            ax.set_xlabel(plot_data.xlabel)
        else:
            for tick in ax.xaxis.get_major_ticks():
                tick.set_visible(False)
            ax.annotate(plot_data.title, xy=(0.5, 1.2), xycoords="axes fraction",
                        size='large', ha='center', va='center', fontweight='bold')
        if i % NUMBER_VARIABLES == 0:
            ax.set_ylabel(plot_data.ylabel)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.set_visible(False)
        ax.annotate(plot_data.panel_id, xy=[-0.05, 1.05], xycoords='axes fraction',
                    fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    plt.subplots_adjust(
        left=0.2,
        right=0.8
    )
    return fig


if __name__ == "__main__":
    plot_composition_variability(
        path_to_large_scales=snakemake.input.large_scales,
        path_to_small_scale=snakemake.input.small_scale,
        path_to_scenario_results=snakemake.input.aggregate,
        path_to_plot=snakemake.output[0]
    )
