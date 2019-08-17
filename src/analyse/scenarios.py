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
regional,regional,0%,
regional,national,0%,
regional,continental,0%,
regional,national,≤15%,
regional,continental,≤15%,
regional,national,≤30%,
regional,continental,≤30%,
national,national,0%,
national,continental,0%,
national,continental,≤15%,
national,continental,≤30%,
continental,continental,0%,
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
        cost_heatmap(results, "continental"),
        cmap=pal,
        ax=ax1,
        cbar=False,
        cbar_kws={"label": "total system costs [-]"},
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        annot=True
    )
    ax1.annotate('a', xy=[-0.05, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    plt.title("CONTINENTAL")

    ax2 = fig.add_subplot(gs[1])
    sns.heatmap(
        cost_heatmap(results, "national"),
        cmap=pal,
        ax=ax2,
        cbar=False,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        annot=True
    )
    plt.ylabel("")
    plt.yticks([])
    ax2.annotate('b', xy=[-0.05, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    plt.title("NATIONAL")

    ax3 = fig.add_subplot(gs[2])
    sns.heatmap(
        cost_heatmap(results, "regional"),
        cmap=pal,
        ax=ax3,
        cbar_ax=ax4,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        annot=True
    )
    plt.ylabel("")
    plt.yticks([])
    ax3.annotate('c', xy=[-0.25, 1.05], xycoords='axes fraction',
                 fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    plt.title("REGIONAL")

    fig.autofmt_xdate()
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.6)
    fig.savefig(path_to_costs, dpi=600, transparent=False)


def cost_heatmap(data, grid_scale):
    if grid_scale == "continental":
        columns = ["continental", "national", "regional"]
    elif grid_scale == "national":
        columns = ["national", "regional"]
    elif grid_scale == "regional":
        columns = ["regional"]
    else:
        raise ValueError("Grid scale must be in {}.".format(
            ["regional", "national", "continental"])
        )

    pivoted_data = data[data["grid_scale"] == grid_scale].pivot(
        columns="autarky_scale",
        index="net_exchange_potential",
        values="cost"
    )
    pivoted_data.columns.name = "autarky scale"
    pivoted_data.index.name = "net electricity import"
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
    results = results / results.loc["continental", "continental", "0%"]
    return results.reset_index()


def parse_scenario_name(scenario_name):
    autarky_layer, _, autarky_level, grid_size, _ = scenario_name.split("-")
    assert autarky_layer in ["regional", "national", "continental"]
    assert grid_size in ["regional", "national", "continental"]
    return autarky_layer, autarky_level, grid_size


if __name__ == "__main__":
    plot_costs(
        paths_to_results=snakemake.input.results,
        path_to_costs=snakemake.output[0]
    )
