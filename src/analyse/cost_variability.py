from dataclasses import dataclass

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from seaborn import utils
from seaborn.palettes import blend_palette

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
GREY = "#C0C0C0"
COLOR = GREEN
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

HYDRO_TECHS = ["hydro_run_of_river", "hydro_reservoir"]
GW_TO_TW = 1e-3
MW_TO_TW = 1e-6

COLUMN_HEADER = [
    "y-continental-scale-cost-eur",
    "y-national-scale-cost-eur",
    "y-cost-diff-eur",
    "y-cost-diff-relative",
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


def plot_cost_variability(path_to_large_scales, path_to_small_scale, path_to_plot):
    plot_datas = read_plot_data(path_to_large_scales, path_to_small_scale)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def read_plot_data(path_to_large_scales, path_to_small_scale):
    y = pd.concat([
        pd.read_csv(path_to_large_scales, index_col=None, header=None),
        pd.read_csv(path_to_small_scale, index_col=None, header=None),

    ], axis="columns")
    y.columns = COLUMN_HEADER
    norm_value = y["y-continental-scale-cost-eur"].min()
    max_value = y["y-regional-scale-cost-eur"].max() / norm_value

    return [
        PlotData(
            panel_id="a",
            title="",
            ylabel="Normed cost national scale",
            xlabel="Normed cost continental scale",
            xlim=(0, max_value),
            ylim=(0, max_value),
            x=y["y-continental-scale-cost-eur"] / norm_value,
            y=y["y-national-scale-cost-eur"] / norm_value,
        ),
        PlotData(
            panel_id="b",
            title="",
            ylabel="Normed cost regional scale",
            xlabel="Normed cost continental scale",
            xlim=(0, max_value),
            ylim=(0, max_value),
            x=y["y-continental-scale-cost-eur"] / norm_value,
            y=y["y-regional-scale-cost-eur"] / norm_value,
        )
    ]


def plot_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    axes = fig.subplots(
        nrows=1,
        ncols=2,
        sharex=True,
        sharey=False
    )

    color_rgb = mpl.colors.colorConverter.to_rgb(COLOR)
    colors = [utils.set_hls_values(color_rgb, l=l)  # noqa
            for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    for i, plot_data in enumerate(plot_datas):
        ax = axes[i]
        ax.hexbin(
            x=plot_data.x,
            y=plot_data.y,
            gridsize=int((plot_data.x.max() - plot_data.x.min()) * 40),
            cmap=cmap
        )
        ax.set_aspect('equal')
        ax.set_ylim(*plot_data.ylim)
        ax.set_xlim(*plot_data.xlim)
        ax.plot(plot_data.xlim, plot_data.ylim, "--", color=GREY)
        ax.set_xlabel(plot_data.xlabel)
        ax.annotate(plot_data.title, xy=(0.5, 1.2), xycoords="axes fraction",
                    size='large', ha='center', va='center', fontweight='bold')
        ax.set_ylabel(plot_data.ylabel)
        ax.annotate(plot_data.panel_id, xy=[-0.05, 1.05], xycoords='axes fraction',
                    fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    plot_cost_variability(
        path_to_large_scales=snakemake.input.large_scales,
        path_to_small_scale=snakemake.input.small_scale,
        path_to_plot=snakemake.output[0]
    )
