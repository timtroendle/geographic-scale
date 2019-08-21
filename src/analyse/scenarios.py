import io
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import xarray as xr

GREEN = "#679436"
RED = "#A01914"
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

DATA_INDEX = """autarky_scale,grid_scale,net_exchange_potential,cost
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


def plot_costs(paths_to_results, path_to_costs):
    """Plot scenario space and results."""
    results = read_results(paths_to_results)

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
        cbar_kws={"label": "total system costs [-]"},
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
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
        annot=True,
        fmt='.3g'
    )
    plt.ylabel("")
    plt.yticks([])
    ax3.annotate('c', xy=[-0.1 * 3, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    ax3.set_title("REGIONAL", loc="left")

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
        columns="autarky_scale",
        index="net_exchange_potential",
        values="cost"
    )
    pivoted_data.columns.name = "Autarky scale"
    pivoted_data.index.name = "Net electricity import"
    return pivoted_data[columns].reindex(["≤30%", "≤15%", "0%"])


def read_results(paths_to_results):
    results = (pd.read_csv(io.StringIO(DATA_INDEX))
                 .set_index(["autarky_scale", "grid_scale", "net_exchange_potential"]))
    for path_to_results in paths_to_results:
        scenario_name = Path(path_to_results).parent.name
        autarky_layer, autarky_level, grid_size = parse_scenario_name(scenario_name)
        autarky_level = AUTARKY_LEVEL_MAP[autarky_level]
        total_system_cost = (xr.open_dataset(path_to_results)["total_levelised_cost"]
                               .squeeze(["costs", "carriers"])
                               .item())
        results.loc[autarky_layer, grid_size, autarky_level] = total_system_cost
    results = results / results.loc["Continental", "Continental", "0%"]
    return results.reset_index()


def parse_scenario_name(scenario_name):
    autarky_layer, _, autarky_level, grid_size, _ = scenario_name.split("-")
    assert autarky_layer in ["regional", "national", "continental"]
    assert grid_size in ["regional", "national", "continental"]
    return autarky_layer.capitalize(), autarky_level, grid_size.capitalize()


if __name__ == "__main__":
    plot_costs(
        paths_to_results=snakemake.input.results,
        path_to_costs=snakemake.output[0]
    )
