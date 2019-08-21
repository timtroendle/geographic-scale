from dataclasses import dataclass

import pandas as pd
import seaborn as sns
import numpy as np
import xarray as xr
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
from seaborn import utils
from seaborn.palettes import blend_palette

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
COLOR = GREEN
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
MAX_VALUE = 2

HYDRO_TECHS = ["hydro_run_of_river", "hydro_reservoir"]
GW_TO_TW = 1e-3
MW_TO_TW = 1e-6


@dataclass
class PlotData:
    title: str
    x: pd.Series
    y: pd.Series
    xlabel: str
    ylabel: str
    xlim: tuple
    ylim: tuple
    gridsize: int


def plot_composition_variability(path_to_xy, path_to_scenario_results, path_to_plot):
    plot_datas = read_plot_data(path_to_xy, path_to_scenario_results)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def read_plot_data(path_to_xy, path_to_scenario_results):
    xy = pd.read_csv(path_to_xy, index_col=0) * GW_TO_TW
    agg = xr.open_dataset(path_to_scenario_results)
    hydro_cap = pd.Series(
        agg.energy_cap.sel(techs=HYDRO_TECHS).sum(["techs", "locs"]).isel(scenario=0).item() * MW_TO_TW,
        index=xy.index
    )
    return [
        PlotData(
            title="Total supply",
            ylabel="National scale [TW]",
            xlabel="Continental scale [TW]",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=(xy["y-large-scale-wind-gw"]
               + xy["y-large-scale-pv-gw"]
               + xy["y-large-scale-biofuel-gw"]
               + hydro_cap),
            y=(xy["y-small-scale-wind-gw"]
               + xy["y-small-scale-pv-gw"]
               + xy["y-small-scale-biofuel-gw"]
               + hydro_cap),
            gridsize=10
        ),
        PlotData(
            title="Wind",
            ylabel="National scale [TW]",
            xlabel="Continental scale [TW]",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=xy["y-large-scale-wind-gw"],
            y=xy["y-small-scale-wind-gw"],
            gridsize=11
        ),
        PlotData(
            title="Bioenergy + storage",
            ylabel="National scale [TW]",
            xlabel="Continental scale [TW]",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=xy["y-large-scale-biofuel-gw"] + xy["y-large-scale-storage-gw"],
            y=xy["y-small-scale-biofuel-gw"] + xy["y-small-scale-storage-gw"],
            gridsize=6

        )
    ]


def plot_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    gs = gridspec.GridSpec(1, len(plot_datas))
    panel_ids = "abcdefghi"

    color_rgb = mpl.colors.colorConverter.to_rgb(COLOR)
    colors = [utils.set_hls_values(color_rgb, l=l)  # noqa
            for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    for i, plot_data in enumerate(plot_datas):
        ax = fig.add_subplot(gs[i])
        ax.hexbin(
            x=plot_data.x,
            y=plot_data.y,
            gridsize=plot_data.gridsize,
            cmap=cmap
        )
        ax.set_aspect('equal')
        ax.set_ylim(*plot_data.ylim)
        ax.set_xlim(*plot_data.xlim)
        ax.plot(plot_data.xlim, plot_data.ylim, "--", color=BLUE)
        ax.set_xlabel(plot_data.xlabel)
        if i == 0:
            ax.set_ylabel(plot_data.ylabel)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.set_visible(False)
        ax.set_title(plot_data.title)
        ax.annotate(panel_ids[i], xy=[-0.05, 1.05], xycoords='axes fraction',
                    fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    return fig


if __name__ == "__main__":
    plot_composition_variability(
        path_to_xy=snakemake.input.xy,
        path_to_scenario_results=snakemake.input.aggregate,
        path_to_plot=snakemake.output[0]
    )
