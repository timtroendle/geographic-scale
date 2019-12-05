import io

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import xarray as xr

RED = "#A01914"
ANTHRACITE = "#424242"
SYSTEM_SCALE_COLOR = "k"
AUTARKY_EXTENT_COLOR = "k"
PALETTE = sns.light_palette(RED)
HIGHLIGHT_COLOR = ANTHRACITE
HIGHLIGHT_LINEWIDTH = 4
HIGHLIGHT_LINESTYLE = "-"
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

BASE_SCENARIO = "continental-autarky-100-continental-grid"
DATA_INDEX = """autarky_layer,grid_scale,autarky_degree,cost
Regional,Regional,0%,
Regional,National,0%,
Regional,Continental,0%,
Regional,National,≤15%,
Regional,Continental,≤15%,
Regional,National,≤30%,
Regional,Continental,≤30%,
National,National,0%,
National,Continental,0%,
National,Continental,≤15%,
National,Continental,≤30%,
Continental,Continental,0%,
"""

AUTARKY_LEVEL_MAP = {
    "100": "0%",
    "85": "≤15%",
    "70": "≤30%"
}


def plot_costs(path_to_aggregated_results, path_to_costs):
    """Plot scenario space and results."""
    results = read_results(path_to_aggregated_results)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3.5))
    gs = gridspec.GridSpec(1, 4, width_ratios=[2, 0.4, 1, 0.1])
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[2])
    cbar_ax = fig.add_subplot(gs[3])
    base_case_box(results, ax_a, cbar_ax)
    special_case_box(results, ax_b)

    fig.text(
        s='a - Base cases',
        x=0.065,
        y=0.98,
        fontsize=PANEL_FONT_SIZE,
        verticalalignment="top",
        weight=PANEL_FONT_WEIGHT
    )
    fig.text(
        s='b - Cases with net imports',
        x=0.45,
        y=0.98,
        fontsize=PANEL_FONT_SIZE,
        verticalalignment="top",
        weight=PANEL_FONT_WEIGHT
    )
    for xy in ((0.04, 2.04), (0.04, 1.04), (1.04, 2.04)):
        ax_a.add_patch(Rectangle(
            xy,
            0.92,
            0.92,
            fill=False,
            edgecolor=HIGHLIGHT_COLOR,
            lw=HIGHLIGHT_LINEWIDTH,
            ls=HIGHLIGHT_LINESTYLE)
        )
    for xy in ((0.04, 2.04), (0.04, 1.04), (0.04, 0.04)):
        ax_b.add_patch(Rectangle(
            xy,
            0.92,
            0.92,
            fill=False,
            edgecolor=HIGHLIGHT_COLOR,
            lw=HIGHLIGHT_LINEWIDTH - 2,
            ls=HIGHLIGHT_LINESTYLE)
        )
    plt.subplots_adjust(
        left=0.05,
        bottom=0.2,
        right=0.88,
        top=0.88,
        wspace=0.3,
        hspace=0.2
    )
    fig.savefig(path_to_costs, dpi=600, transparent=False)


def base_case_box(results, ax, cbar_ax):
    cbar_ticks = np.linspace(
        start=results["cost"].min(),
        stop=results["cost"].max(),
        num=len(PALETTE) + 1
    )
    sns.heatmap(
        results[results.autarky_degree == "0%"].pivot(
            index="autarky_layer",
            columns="grid_scale",
            values="cost"
        ).reindex(index=["Continental", "National", "Regional"]),
        cbar=True,
        cbar_ax=cbar_ax,
        cbar_kws={
            "label": "Cost relative to lowest cost\n continental-scale system",
            "ticks": cbar_ticks,
            "format": "%.1f"
        },
        cmap=PALETTE,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=1.25,
        annot=True,
        square=True,
        ax=ax,
        fmt='.3g'
    )
    ax.set_xlabel("Balancing scale")
    ax.set_ylabel("Supply scale")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=90, va="center")


def special_case_box(results, ax):
    results = results.copy()
    results["case_name"] = pd.Series(-1, index=results.index)
    case1_mask = (results.grid_scale == "Continental") & (results.autarky_layer == "National")
    case2_mask = (results.grid_scale == "Continental") & (results.autarky_layer == "Regional")
    case3_mask = (results.grid_scale == "National") & (results.autarky_layer == "Regional")
    for i, mask in enumerate([case1_mask, case2_mask, case3_mask]):
        results.loc[mask, "case_name"] = i

    sns.heatmap(
        results[results.case_name.isin([0, 1, 2])].pivot(
            index="case_name",
            columns="autarky_degree",
            values="cost"
        ).reindex(columns=["0%", "≤15%", "≤30%"]),
        cmap=PALETTE,
        ax=ax,
        cbar=False,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        square=True,
        fmt='.3g'
    )
    ax.set_xlabel("Allowed net imports")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelrotation=0, reset=True)
    ax.yaxis.tick_left()
    ax.set_yticklabels([
        "National supply\ncontinental balancing",
        "Regional supply\ncontinental balancing",
        "Regional supply\nnational balancing"
    ])


def read_results(path_to_agregrated_results):
    cost = xr.open_dataset(path_to_agregrated_results)["cost"].sum(["locs", "techs"])
    relative_cost = cost / cost.sel(scenario=BASE_SCENARIO)
    results = (pd.read_csv(io.StringIO(DATA_INDEX))
                 .set_index(["autarky_layer", "grid_scale", "autarky_degree"]))
    for scenario in relative_cost.scenario:
        scenario = scenario.item()
        autarky_layer, autarky_level, grid_size = parse_scenario_name(scenario)
        autarky_level = AUTARKY_LEVEL_MAP[autarky_level]
        results.loc[autarky_layer, grid_size, autarky_level] = relative_cost.sel(scenario=scenario).item()
    return results.reset_index()


def parse_scenario_name(scenario_name):
    autarky_layer, _, autarky_level, grid_size, _ = scenario_name.split("-")
    assert autarky_layer in ["regional", "national", "continental"]
    assert grid_size in ["regional", "national", "continental"]
    return autarky_layer.capitalize(), autarky_level, grid_size.capitalize()


if __name__ == "__main__":
    plot_costs(
        path_to_aggregated_results=snakemake.input.results,
        path_to_costs=snakemake.output[0]
    )
