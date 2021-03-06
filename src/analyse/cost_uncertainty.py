from dataclasses import dataclass

import pandas as pd
import xarray as xr
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

HYDRO_TECHS = ["hydro_run_of_river", "hydro_reservoir"]
BASE_SCENARIO = "continental-autarky-100-continental-grid"
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
    x: pd.Series
    y: pd.Series
    xlabel: str
    ylabel: str
    xlim: tuple
    ylim: tuple


def plot_cost_variability(path_to_large_scales, path_to_small_scale, path_to_scenario_results, path_to_plot):
    plot_datas = read_plot_data(path_to_large_scales, path_to_small_scale, path_to_scenario_results)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def read_plot_data(path_to_large_scales, path_to_small_scale, path_to_scenario_results):
    y = pd.concat([
        pd.read_csv(path_to_large_scales, index_col=None, header=None),
        pd.read_csv(path_to_small_scale, index_col=None, header=None),

    ], axis="columns")
    y.columns = COLUMN_HEADER
    norm_value = (
        xr
        .open_dataset(path_to_scenario_results)
        .cost
        .sel(scenario=BASE_SCENARIO)
        .sum(["locs", "techs"])
        .item()
    )
    max_value = 3

    return [
        PlotData(
            panel_id="A",
            ylabel="National scale (-)",
            xlabel="Continental scale (-)",
            xlim=(0, max_value),
            ylim=(0, max_value),
            x=y["y-continental-scale-cost-eur"] / norm_value,
            y=y["y-national-scale-cost-eur"] / norm_value,
        ),
        PlotData(
            panel_id="B",
            ylabel="Regional scale (-)",
            xlabel="Continental scale (-)",
            xlim=(0, max_value),
            ylim=(0, max_value),
            x=y["y-continental-scale-cost-eur"] / norm_value,
            y=y["y-regional-scale-cost-eur"] / norm_value,
        ),
        PlotData(
            panel_id="C",
            ylabel="Regional scale (-)",
            xlabel="National scale (-)",
            xlim=(0, max_value),
            ylim=(0, max_value),
            x=y["y-national-scale-cost-eur"] / norm_value,
            y=y["y-regional-scale-cost-eur"] / norm_value,
        )
    ]


def plot_data(plot_datas):
    fig = plt.figure(figsize=(4.41, 1.9))
    axes = fig.subplots(
        nrows=1,
        ncols=3,
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
        ax.set_ylabel(plot_data.ylabel)
        ax.set_yticks([0, 1, 2, 3])
        ax.set_title(plot_data.panel_id, loc="left")
    fig.suptitle(t="System cost relative to continental-scale base case")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    plot_cost_variability(
        path_to_large_scales=snakemake.input.large_scales,
        path_to_small_scale=snakemake.input.small_scale,
        path_to_scenario_results=snakemake.input.aggregate,
        path_to_plot=snakemake.output[0]
    )
