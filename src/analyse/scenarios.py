import io
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import xarray as xr

GREEN = "#679436"
RED = "#A01914"

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


def plot_costs(paths_to_results, path_to_costs, scaling_factor_cost):
    """Plot scenario space and results."""
    results = read_results(paths_to_results, scaling_factor_cost)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    pal = sns.light_palette(RED)
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 2, 3, 0.2])
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
    plt.title("CONTINENTAL AUTARKY")
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
    plt.title("NATIONAL AUTARKY")
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
    plt.title("REGIONAL AUTARKY")
    fig.autofmt_xdate()
    plt.tight_layout()
    fig.savefig(path_to_costs, dpi=300, transparent=True)


def cost_heatmap(data, autarky_scale):
    if autarky_scale == "regional":
        columns = ["continental", "national", "regional"]
    elif autarky_scale == "national":
        columns = ["continental", "national"]
    elif autarky_scale == "continental":
        columns = ["continental"]
    else:
        raise ValueError("Autarky scale must be in {}.".format(
            ["regional", "national", "continental"])
        )

    pivoted_data = data[data["autarky_scale"] == autarky_scale].pivot(
        columns="grid_scale",
        index="net_exchange_potential",
        values="cost"
    )
    pivoted_data.columns.name = "grid size"
    pivoted_data.index.name = "net electricity import"
    return pivoted_data[columns].reindex(["≤30%", "≤15%", "0%"])


def read_results(paths_to_results, scaling_factor_cost):
    results = (pd.read_csv(io.StringIO(DATA_INDEX))
                 .set_index(["autarky_scale", "grid_scale", "net_exchange_potential"]))
    for path_to_results in paths_to_results:
        scenario_name = Path(path_to_results).parent.name
        autarky_layer, autarky_level, grid_size = parse_scenario_name(scenario_name)
        autarky_level = AUTARKY_LEVEL_MAP[autarky_level]
        total_system_cost = (xr.open_dataset(path_to_results)["total_levelised_cost"]
                               .squeeze(["costs", "carriers"])
                               .item()) * scaling_factor_cost / 1e1 # scale from €/MWh to €ct/kWh
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
        path_to_costs=snakemake.output[0],
        scaling_factor_cost=snakemake.params.scaling_factor_cost
    )
