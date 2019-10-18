from dataclasses import dataclass, field
from itertools import product

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
GREEN = "#679436"
BLUE = "#4F6DB8"
BLUE_CMAP = sns.light_palette(BLUE, n_colors=10, reverse=False, as_cmap=False)
HIGHLIGHT_COLOR = "#424242"
HIGHLIGHT_LINEWIDTH = 4
CMAP = BLUE_CMAP

ROWS = [
    "System cost (€)",
    "Solar (MW)",
    "Wind (MW)",
    "Bioenergy (MW)",
    "Storage (MW)",
    "Storage (MWh)"
]
DIFFERENCE_ROWS = [
    "System cost",
    "Total supply",
    "Wind capacity",
    "Bioenergy and\nstorage capacity"
]


@dataclass
class PlotData:
    name: str
    data: pd.DataFrame
    xticklabels: list
    yticklabels: list
    highlight_rows: list = field(default_factory=list)
    highlight_cols: list = field(default_factory=list)


def sobol(path_to_cont_and_nat_data, path_to_reg_data, all_data, parameters, path_to_plot):
    uncertain_parameters = [param["descriptive-name"] for unused, param in parameters.items()]
    if all_data:
        plot_datas = prepare_all_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters)
        fig = plot_all_data(plot_datas)
    else:
        plot_datas = prepare_diff_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters)
        fig = plot_diff_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def prepare_all_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters):
    data = pd.concat([
        pd.read_csv(path_to_cont_and_nat_data, header=None).T,
        pd.read_csv(path_to_reg_data, header=None).T
    ])
    return [
        PlotData(
            name="a – Continental scale",
            data=data.iloc[[0, 7, 9, 11, 13, 15]],
            yticklabels=ROWS,
            xticklabels=uncertain_parameters
        ),
        PlotData(
            name="b – National scale",
            data=data.iloc[[1, 8, 10, 12, 14, 16]],
            yticklabels=ROWS,
            xticklabels=uncertain_parameters
        ),
        PlotData(
            name="c – Regional scale",
            data=data.iloc[[18, 19, 20, 21, 22, 23]],
            yticklabels=ROWS,
            xticklabels=uncertain_parameters
        ),
        PlotData(
            name="d – Relative difference between continental and national scales",
            data=data.iloc[[3, 4, 5, 6]],
            yticklabels=DIFFERENCE_ROWS,
            xticklabels=uncertain_parameters
        )
    ]


def prepare_diff_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters):
    data = pd.concat([
        pd.read_csv(path_to_cont_and_nat_data, header=None).T,
        pd.read_csv(path_to_reg_data, header=None).T
    ])
    return [
        PlotData(
            name="Relative difference between continental and national scales",
            data=data.iloc[[3, 4, 5, 6]],
            yticklabels=DIFFERENCE_ROWS,
            xticklabels=uncertain_parameters
        )
    ]


def plot_all_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 12.5))
    axes = fig.subplots(
        nrows=len(plot_datas),
        ncols=1,
        sharex=True,
        sharey=False,
        gridspec_kw={'height_ratios': [len(plot_data.data.index) for plot_data in plot_datas]}
    )
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.9, 0.25, 0.015, 0.65])

    for i, plot_data in enumerate(plot_datas):
        cbar = True if i == 0 else False
        sns.heatmap(
            plot_data.data,
            ax=axes[i],
            vmin=0,
            vmax=1,
            square=True,
            linewidth=0.75,
            cbar=cbar,
            cbar_ax=cbar_ax,
            cmap=CMAP,
            yticklabels=plot_data.yticklabels,
            xticklabels=plot_data.xticklabels
        )
        axes[i].annotate(
            plot_data.name,
            xy=[0, 1 + 0.05 * 6 / len(plot_data.data.index)],
            xycoords='axes fraction',
            fontsize=PANEL_FONT_SIZE,
            weight=PANEL_FONT_WEIGHT
        )

        for (row, column) in product(plot_data.highlight_rows, plot_data.highlight_cols):
            axes[i].add_patch(
                Rectangle(
                    (row + 0.04, column + 0.04),
                    0.92,
                    0.92,
                    fill=False,
                    edgecolor=HIGHLIGHT_COLOR,
                    lw=HIGHLIGHT_LINEWIDTH
                )
            )
        axes[i].xaxis.set_tick_params(bottom=False)
        axes[i].yaxis.set_tick_params(left=False)
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    return fig


def plot_diff_data(plot_datas):
    assert len(plot_datas) == 1
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.subplots(
        nrows=len(plot_datas),
        ncols=1,
        sharex=True,
        sharey=False,
        gridspec_kw={'height_ratios': [len(plot_data.data.index) for plot_data in plot_datas]}
    )
    fig.subplots_adjust(right=0.6, bottom=0.26)
    cbar_ax = fig.add_axes([0.9, 0.25, 0.015, 0.65])

    for i, plot_data in enumerate(plot_datas):
        cbar = True if i == 0 else False
        sns.heatmap(
            plot_data.data,
            ax=ax,
            vmin=0,
            vmax=1,
            square=True,
            linewidth=0.75,
            cbar=cbar,
            cbar_ax=cbar_ax,
            cmap=CMAP,
            yticklabels=plot_data.yticklabels,
            xticklabels=plot_data.xticklabels
        )
        ax.annotate(
            plot_data.name,
            xy=[0, 1 + 0.05 * 6 / len(plot_data.data.index)],
            xycoords='axes fraction',
            fontsize=PANEL_FONT_SIZE,
            weight=PANEL_FONT_WEIGHT
        )

        for (row, column) in product(plot_data.highlight_rows, plot_data.highlight_cols):
            ax.add_patch(
                Rectangle(
                    (row + 0.04, column + 0.04),
                    0.92,
                    0.92,
                    fill=False,
                    edgecolor=HIGHLIGHT_COLOR,
                    lw=HIGHLIGHT_LINEWIDTH
                )
            )
        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False)
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    return fig


if __name__ == "__main__":
    sobol(
        path_to_cont_and_nat_data=snakemake.input.indices_cont_and_nat,
        path_to_reg_data=snakemake.input.indices_reg,
        all_data=snakemake.wildcards.extent == "all",
        parameters=snakemake.params.parameters,
        path_to_plot=snakemake.output[0]
    )
