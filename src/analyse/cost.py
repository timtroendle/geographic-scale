import io

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import xarray as xr

GREEN = "#679436"
RED = "#A01914"
YELLOW = "#FABC3C"
SYSTEM_SCALE_COLOR = "k"
AUTARKY_EXTENT_COLOR = "k"
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
    fig = plt.figure(figsize=(8, 4))
    pal = sns.light_palette(RED)
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.15, 0.6, 2, 1, 0.2])
    ax4 = fig.add_subplot(gs[4])
    ax1 = fig.add_subplot(gs[0])

    sns.heatmap(
        results[results.autarky_layer == results.grid_scale].pivot(
            index="autarky_layer",
            columns="autarky_degree",
            values="cost"
        ),
        cmap=pal,
        ax=ax1,
        cbar=False,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=1.25,
        yticklabels=["Continental scale", "National scale", "Regional scale"],
        annot=True,
        square=True,
        fmt='.3g'
    )
    plt.ylabel("")
    plt.xlabel("")
    plt.xticks([])
    ax1.tick_params(axis="y", labelrotation=0)

    ax2 = fig.add_subplot(gs[2])
    sns.heatmap(
        results[(results.autarky_layer != results.grid_scale) & (results.grid_scale == "Continental")].pivot(
            columns="autarky_layer",
            index="autarky_degree",
            values="cost"
        ).reindex(["≤30%", "≤15%", "0%"]),
        cmap=pal,
        ax=ax2,
        cbar=False,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        square=True,
        fmt='.3g'
    )
    plt.ylabel("Allowed net imports")
    plt.xlabel("")
    ax2.tick_params(axis="y", labelrotation=0)
    ax2.set_title("Continental")

    ax3 = fig.add_subplot(gs[3])
    sns.heatmap(
        results[(results.autarky_layer != results.grid_scale) & (results.grid_scale == "National")].pivot(
            columns="autarky_layer",
            index="autarky_degree",
            values="cost"
        ).reindex(["≤30%", "≤15%", "0%"]),
        cmap=pal,
        ax=ax3,
        cbar_ax=ax4,
        cbar_kws={"label": "Cost relative to lowest cost\n continental scale system"},
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        square=True,
        fmt='.3g'
    )
    plt.ylabel("")
    plt.yticks([])
    plt.xlabel("")
    ax3.tick_params(axis="y", labelrotation=0)
    ax3.set_title("National")

    fig.text(
        s='a - Base cases',
        x=.12,
        y=0.98,
        fontsize=PANEL_FONT_SIZE,
        verticalalignment="top",
        weight=PANEL_FONT_WEIGHT
    )
    fig.text(
        s='b - Cases with lower level self-sufficiency',
        x=.4,
        y=0.98,
        fontsize=PANEL_FONT_SIZE,
        verticalalignment="top",
        weight=PANEL_FONT_WEIGHT
    )
    fig.text(
        s="─────────── System scale ────────────",
        y=0.89,
        x=0.63,
        color=SYSTEM_SCALE_COLOR,
        horizontalalignment="center",
        verticalalignment="top",
        weight="bold"
    )
    fig.text(
        s="─────── Self-sufficiency extent ───────",
        y=0.05,
        x=0.63,
        color=AUTARKY_EXTENT_COLOR,
        horizontalalignment="center",
        verticalalignment="bottom",
        weight="bold"
    )

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.6, left=0.15, top=0.85, bottom=0.1)
    fig.savefig(path_to_costs, dpi=600, transparent=False)


def cost_heatmap(data, grid_scale):
    if grid_scale == "Continental":
        columns = ["Continental", "National", "Regional"]
    elif grid_scale == "National":
        columns = ["National", "Regional"]
    elif grid_scale == "Regional":
        columns = ["Regional"]
    else:
        raise ValueError("Grid scale must be in {}.".format(
            ["Regional", "National", "Continental"])
        )

    pivoted_data = data[data["grid_scale"] == grid_scale].pivot(
        columns="autarky_layer",
        index="autarky_degree",
        values="cost"
    )
    pivoted_data.columns.name = "Self-sufficiency layer"
    pivoted_data.index.name = "Allowed net imports"
    return pivoted_data[columns].reindex(["≤30%", "≤15%", "0%"])


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
