from dataclasses import dataclass

import pandas as pd
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
COLOR = BLUE
MAX_VALUE = 2
NUMBER_VARIABLES = 2

GW_TO_TW = 1e-3
MW_TO_TW = 1e-6

COLUMN_HEADER = [
    "y-continental-scale-cost-eur",
    "y-national-scale-cost-eur",
    "y-cost-diff-eur",
    "y-cost-diff-relative",
    "y-supply-diff-relative",
    "y-wind-diff-relative",
    "y-balancing-diff-relative",
    "y-continental-scale-pv-gw",
    "y-national-scale-pv-gw",
    "y-continental-scale-wind-gw",
    "y-national-scale-wind-gw",
    "y-continental-scale-hydro-gw",
    "y-national-scale-hydro-gw",
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
    "y-regional-scale-hydro-gw",
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


def plot_composition_variability(path_to_large_scales, path_to_small_scale, path_to_plot):
    plot_datas = read_plot_data(path_to_large_scales, path_to_small_scale)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def read_plot_data(path_to_large_scales, path_to_small_scale):
    y = pd.concat([
        pd.read_csv(path_to_large_scales, index_col=None, header=None),
        pd.read_csv(path_to_small_scale, index_col=None, header=None),

    ], axis="columns") * GW_TO_TW
    y.columns = COLUMN_HEADER
    return [
        PlotData(
            panel_id="A",
            title="Supply",
            ylabel="National scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=(y["y-continental-scale-wind-gw"]
               + y["y-continental-scale-pv-gw"]
               + y["y-continental-scale-hydro-gw"]),
            y=(y["y-national-scale-wind-gw"]
               + y["y-national-scale-pv-gw"]
               + y["y-national-scale-hydro-gw"]),
        ),
        PlotData(
            panel_id="B",
            title="Balancing",
            ylabel="National scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=y["y-continental-scale-biofuel-gw"] + y["y-continental-scale-storage-gw"],
            y=y["y-national-scale-biofuel-gw"] + y["y-national-scale-storage-gw"],
        ),
        PlotData(
            panel_id="C",
            title="Supply",
            ylabel="Regional scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=(y["y-continental-scale-wind-gw"]
               + y["y-continental-scale-pv-gw"]
               + y["y-continental-scale-hydro-gw"]),
            y=(y["y-regional-scale-wind-gw"]
               + y["y-regional-scale-pv-gw"]
               + y["y-regional-scale-hydro-gw"]),
        ),
        PlotData(
            panel_id="D",
            title="Balancing",
            ylabel="Regional scale (TW)",
            xlabel="Continental scale (TW)",
            xlim=(0, MAX_VALUE),
            ylim=(0, MAX_VALUE),
            x=y["y-continental-scale-biofuel-gw"] + y["y-continental-scale-storage-gw"],
            y=y["y-regional-scale-biofuel-gw"] + y["y-regional-scale-storage-gw"],
        )
    ]


def plot_data(plot_datas):
    fig = plt.figure(figsize=(4.41, 4))
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
            ax.set_title(plot_data.title, loc="left")
        if i % NUMBER_VARIABLES == 0:
            ax.set_ylabel(plot_data.ylabel)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.set_visible(False)
        ax.set_title(plot_data.panel_id, loc="left")
    plt.subplots_adjust(
        left=0.13,
        right=0.93,
        top=0.93
    )
    return fig


if __name__ == "__main__":
    plot_composition_variability(
        path_to_large_scales=snakemake.input.large_scales,
        path_to_small_scale=snakemake.input.small_scale,
        path_to_plot=snakemake.output[0]
    )
