from dataclasses import dataclass, field
from itertools import product

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
GREEN = "#679436"
GREEN_CMAP = sns.light_palette(GREEN, n_colors=5, reverse=False, as_cmap=True)
HIGHLIGHT_COLOR = "#424242"
HIGHLIGHT_LINEWIDTH = 4
CMAP = GREEN_CMAP

ROWS = [
    "System cost (€)",
    "Solar (MW)",
    "Wind (MW)",
    "Bioenergy (MW)",
    "Storage (MW)",
    "Storage (MWh)"
]
DIFFERENCE_ROWS = [
    "System cost (€)",
    "Relative system cost"
]
COLUMNS = [
    "Depreciation rate",
    "Cost solar",
    "Cost onshore wind",
    "Cost offshore wind",
    "Cost battery power",
    "Cost battery energy",
    "Cost hydrogen power",
    "Cost hydrogen energy",
    "Cost transmission",
    "Cost bioenergy",
    "Fuel cost bioenergy",
    "Availability bioenergy"
]


@dataclass
class PlotData:
    name: str
    data: pd.DataFrame
    yticklabels: str
    highlight_rows: list = field(default_factory=list)
    highlight_cols: list = field(default_factory=list)


def sobol(path_to_cont_and_nat_data, path_to_reg_data, path_to_plot):
    plot_datas = prepare_data(path_to_cont_and_nat_data, path_to_reg_data)
    fig = plot_data(plot_datas)
    fig.savefig(path_to_plot, dpi=600)


def prepare_data(path_to_cont_and_nat_data, path_to_reg_data):
    data = pd.concat([
        pd.read_csv(path_to_cont_and_nat_data, header=None).T,
        pd.read_csv(path_to_reg_data, header=None).T
    ])
    return [
        PlotData(
            name="a – Continental scale",
            data=data.iloc[[0, 4, 6, 8, 10, 12]],
            yticklabels=ROWS
        ),
        PlotData(
            name="b – National scale",
            data=data.iloc[[1, 5, 7, 9, 11, 13]],
            yticklabels=ROWS
        ),
        PlotData(
            name="c – Regional scale",
            data=data.iloc[[15, 16, 17, 18, 19, 20]],
            yticklabels=ROWS
        ),
        PlotData(
            name="d – Difference between continental and national scales",
            data=data.iloc[[2, 3]],
            yticklabels=DIFFERENCE_ROWS,
            highlight_rows=[0, 2, 9],
            highlight_cols=[1]
        )
    ]


def plot_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 11.5))
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
            xticklabels=COLUMNS
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


if __name__ == "__main__":
    sobol(
        path_to_cont_and_nat_data=snakemake.input.indices_cont_and_nat,
        path_to_reg_data=snakemake.input.indices_reg,
        path_to_plot=snakemake.output[0]
    )
