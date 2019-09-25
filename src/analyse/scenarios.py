import io

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import xarray as xr

GREEN = "#679436"
RED = "#A01914"
HIGHLIGHT_COLORS = sns.light_palette(GREEN, n_colors=4, reverse=False)[1:]
HIGHLIGHT_LINEWIDTH = 4
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

BASE_SCENARIO = "continental-autarky-100-continental-grid"
DATA_INDEX = """autarky_layer,grid_scale,autarky_degree,cost
Regional,Regional,100%,
Regional,National,100%,
Regional,Continental,100%,
Regional,National,≥85%,
Regional,Continental,≥85%,
Regional,National,≥70%,
Regional,Continental,≥70%,
National,National,100%,
National,Continental,100%,
National,Continental,≥85%,
National,Continental,≥70%,
Continental,Continental,100%,
"""

AUTARKY_LEVEL_MAP = {
    "100": "100%",
    "85": "≥85%",
    "70": "≥70%"
}


def plot_costs(path_to_aggregated_results, path_to_costs):
    """Plot scenario space and results."""
    results = read_results(path_to_aggregated_results)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    pal = sns.light_palette(RED)
    gs = gridspec.GridSpec(1, 4, width_ratios=[3, 2, 1, 0.2])
    ax4 = fig.add_subplot(gs[3])
    ax1 = fig.add_subplot(gs[0])
    sns.heatmap(
        cost_heatmap(results, "Continental"),
        cmap=pal,
        ax=ax1,
        cbar=False,
        cbar_kws={"label": "total system costs"},
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        fmt='.3g'
    )
    ax1.annotate('a', xy=[-0.1, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    ax1.set_title("CONTINENTAL", loc="left")

    ax2 = fig.add_subplot(gs[1])
    sns.heatmap(
        cost_heatmap(results, "National"),
        cmap=pal,
        ax=ax2,
        cbar=False,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        fmt='.3g'
    )
    plt.ylabel("")
    plt.yticks([])
    ax2.annotate('b', xy=[-0.1 * 3 / 2, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    ax2.set_title("NATIONAL", loc="left")

    ax3 = fig.add_subplot(gs[2])
    sns.heatmap(
        cost_heatmap(results, "Regional"),
        cmap=pal,
        ax=ax3,
        cbar_ax=ax4,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=0.75,
        annot=True,
        fmt='.3g'
    )
    plt.ylabel("")
    plt.yticks([])
    ax3.annotate('c', xy=[-0.1 * 3, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    ax3.set_title("REGIONAL", loc="left")

    for color, ax in zip(HIGHLIGHT_COLORS, fig.axes[1:]): # highlight base cases
        ax.add_patch(Rectangle((0.04, 2.04), 0.92, 0.92, fill=False, edgecolor=color, lw=HIGHLIGHT_LINEWIDTH))

    fig.autofmt_xdate()
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.6)
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
    pivoted_data.index.name = "Self-sufficiency degree"
    return pivoted_data[columns].reindex(["≥70%", "≥85%", "100%"])


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
