import io

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec


COST_DATA = """autarky_scale,grid_scale,net_exchange_potential,cost
regional,regional,0%,10
regional,national,0%,9
regional,continental,0%,7
regional,national,≤15%,8
regional,continental,≤15%,6
regional,national,≤30%,7
regional,continental,≤30%,5
national,national,0%,9
national,continental,0%,7
national,continental,≤15%,6
national,continental,≤30%,5
continental,continental,0%,4
"""


def plot_costs(path_to_costs):
    """Plot scenario space and hypothesised result."""
    results = pd.read_csv(io.StringIO(COST_DATA))

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    pal = sns.light_palette("green")
    gs = gridspec.GridSpec(1, 4, width_ratios=[3, 2, 1, 0.2])
    ax4 = fig.add_subplot(gs[3])
    ax1 = fig.add_subplot(gs[0])
    sns.heatmap(
        cost_heatmap(results, "regional"),
        cmap=pal, ax=ax1, cbar_ax=ax4,
        cbar_kws={"label": "total system costs"}
    )
    plt.title("REGIONAL AUTARKY")
    ax2 = fig.add_subplot(gs[1])
    sns.heatmap(cost_heatmap(results, "national"), cmap=pal[:-1], ax=ax2, cbar=False)
    plt.ylabel("")
    plt.yticks([])
    plt.title("NATIONAL AUTARKY")
    ax3 = fig.add_subplot(gs[2])
    sns.heatmap(cost_heatmap(results, "continental"), cmap=pal[:-2], ax=ax3, cbar=False)
    plt.ylabel("")
    plt.yticks([])
    plt.title("CONTINENTAL AUTARKY")
    fig.autofmt_xdate()
    plt.tight_layout()
    fig.savefig(path_to_costs, dpi=300, transparent=True)


def cost_heatmap(data, autarky_scale):
    if autarky_scale == "regional":
        columns = ["regional", "national", "continental"]
    elif autarky_scale == "national":
        columns = ["national", "continental"]
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


if __name__ == "__main__":
    plot_costs(path_to_costs=snakemake.output[0])
