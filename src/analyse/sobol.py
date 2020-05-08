from dataclasses import dataclass, field

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
GREEN = "#679436"
BLUE = "#4F6DB8"
BLUE_CMAP = sns.light_palette(BLUE, n_colors=10, reverse=False, as_cmap=False)
HIGHLIGHT_COLOR = "#424242"
HIGHLIGHT_LINEWIDTH = 4
CMAP = BLUE_CMAP

OUTPUTS = [
    "System cost (€)",
    "Solar (MW)",
    "Wind (MW)",
    "Bioenergy (MW)",
    "Storage (MW)",
    "Storage (MWh)"
]
DIFF_OUTPUTS = [
    "System cost",
    "Total supply\ncapacity",
    "Total balancing\ncapacity"
]

ROW_INDEX = [
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
    name: str
    data: pd.DataFrame
    highlight_rows: list = field(default_factory=list)
    highlight_cols: list = field(default_factory=list)
    ax: plt.Axes = None


def sobol(path_to_cont_and_nat_data, path_to_reg_data, all_data, parameters, path_to_plot):
    uncertain_parameters = [param["descriptive-name"] for unused, param in parameters.items()]
    if all_data:
        plot_datas = prepare_all_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters)
        fig = plot_all_data(plot_datas)
    else:
        plot_data = prepare_diff_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters)
        fig = plot_diff_data(plot_data)
    fig.savefig(path_to_plot, dpi=600)


def prepare_all_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters):
    data = pd.concat([
        pd.read_csv(path_to_cont_and_nat_data, header=None).T,
        pd.read_csv(path_to_reg_data, header=None).T
    ], ignore_index=True)
    data.columns = uncertain_parameters
    data.index = ROW_INDEX

    return [
        PlotData(
            name="a – Continental scale",
            data=(
                data
                .loc[[
                    "y-continental-scale-cost-eur",
                    "y-continental-scale-pv-gw",
                    "y-continental-scale-wind-gw",
                    "y-continental-scale-biofuel-gw",
                    "y-continental-scale-storage-gw",
                    "y-continental-scale-storage-gwh"]]
                .set_index(pd.Index(OUTPUTS))
                .T
            )
        ),
        PlotData(
            name="b – National scale",
            data=(
                data
                .loc[[
                    "y-national-scale-cost-eur",
                    "y-national-scale-pv-gw",
                    "y-national-scale-wind-gw",
                    "y-national-scale-biofuel-gw",
                    "y-national-scale-storage-gw",
                    "y-national-scale-storage-gwh"]]
                .set_index(pd.Index(OUTPUTS))
                .T
            )
        ),
        PlotData(
            name="c – Regional scale",
            data=(
                data
                .loc[[
                    "y-regional-scale-cost-eur",
                    "y-regional-scale-pv-gw",
                    "y-regional-scale-wind-gw",
                    "y-regional-scale-biofuel-gw",
                    "y-regional-scale-storage-gw",
                    "y-regional-scale-storage-gwh"]]
                .set_index(pd.Index(OUTPUTS))
                .T
            )
        ),
        PlotData(
            name="d – Relative difference between continental and national scales",
            data=(
                data
                .loc[[
                    "y-cost-diff-relative",
                    "y-supply-diff-relative",
                    "y-balancing-diff-relative"]]
                .set_index(pd.Index(DIFF_OUTPUTS))
            )
        )
    ]


def prepare_diff_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters):
    all_plot_datas = prepare_all_data(path_to_cont_and_nat_data, path_to_reg_data, uncertain_parameters)
    plot_data = all_plot_datas[-1]
    plot_data.data = plot_data.data[plot_data.data.iloc[0].T.sort_values(ascending=False).index] # sort index by relevance
    plot_data.name = "Relative difference between continental and national scales"
    return plot_data


def plot_all_data(plot_datas):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 9))
    axes = fig.subplots(
        nrows=2,
        ncols=3,
        sharex=False,
        sharey=True,
        gridspec_kw={'height_ratios': [20, 10]}
    )
    cbar_ax = fig.add_axes([0.925, 0.25, 0.015, 0.65])
    gs = axes[1, 0].get_gridspec()
    plot_datas[0].ax = axes[0][0]
    plot_datas[1].ax = axes[0][1]
    plot_datas[2].ax = axes[0][2]
    for ax in axes[1, :]:
        ax.remove()
    plot_datas[3].ax = fig.add_subplot(gs[1, :])

    for i, plot_data in enumerate(plot_datas):
        cbar = True if i == 0 else False
        ax = plot_data.ax
        sns.heatmap(
            plot_data.data,
            ax=ax,
            vmin=0,
            vmax=1,
            square=True,
            linewidth=0.75,
            cbar=cbar,
            cbar_ax=cbar_ax,
            cmap=CMAP
        )
        ax.annotate(
            plot_data.name,
            xy=[0, 1 + 0.05 * 6 / len(plot_data.data.index)],
            xycoords='axes fraction',
            fontsize=PANEL_FONT_SIZE,
            weight=PANEL_FONT_WEIGHT
        )
        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False)
    fig.tight_layout(h_pad=0.7, rect=[0, 0, 0.9, 1])
    plot_datas[0].ax.invert_yaxis() # not sure why this is necessary; looks like a matplotlib bug
    plot_datas[1].ax.invert_yaxis()
    return fig


def plot_diff_data(plot_data):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.subplots(nrows=1, ncols=1)
    fig.subplots_adjust(right=0.6, bottom=0.26)
    cbar_ax = fig.add_axes([0.9, 0.25, 0.015, 0.65])

    sns.heatmap(
        plot_data.data,
        ax=ax,
        vmin=0,
        vmax=1,
        square=True,
        linewidth=0.75,
        cbar=True,
        cbar_ax=cbar_ax,
        cmap=CMAP,
    )
    ax.annotate(
        plot_data.name,
        xy=[0, 1 + 0.05 * 6 / len(plot_data.data.index)],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
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
