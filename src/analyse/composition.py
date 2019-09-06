from collections import OrderedDict
from itertools import product

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec


GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
PALETTE = sns.light_palette(GREEN, n_colors=4, reverse=False)[1:]
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
ERROR_BAR_LINEWIDTH = 3.5
ERROR_BAR_COLOR = '#454545'

GENERATION_CAPACITIES = OrderedDict([
    ("Capacity|Generation total", "Total"),
    ("Capacity|Wind", "Wind"),
    ("Capacity|Solar PV", "Solar"),
    ("Capacity|Bioenergy", "Bio"),
    ("Capacity|Storage|Short term|Power", "Battery"),
    ("Capacity|Storage|Long term|Power", "Hydrogen")
])
STORAGE_CAPACITIES = OrderedDict([
    ("Capacity|Storage|Short term|Energy", "Battery"),
    ("Capacity|Storage|Long term|Energy", "Hydrogen")
])
TRANSMISSION_CAPACITIES = OrderedDict([
    ("Capacity|Transmission", "Transmission")
])
VRES = OrderedDict([
    ("Energy|Potential vres", "Potential"),
    ("Energy|Renewable curtailment|Absolute total", "Curtailment")
])
MAIN_SCENARIOS = [
    "continental-autarky-100-continental-grid",
    "national-autarky-100-national-grid",
    "regional-autarky-100-regional-grid"
]
SCALE_ORDER = OrderedDict([
    ("Continental", -0.25),
    ("National", 0),
    ("Regional", 0.25)
])


def composition(path_to_aggregated_results, path_to_output):
    """Plot system compositions for all scenarios."""
    data = pd.read_csv(path_to_aggregated_results)
    data["Scale"] = data["Scenario"].str.split("-").str.get(-2).str.capitalize()

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(2, 3, width_ratios=[2, 1, 2])

    ax = fig.add_subplot(gs[0:3])
    ax.plot([0, 0], [0, 0], color=ERROR_BAR_COLOR, lw=ERROR_BAR_LINEWIDTH, label='Range with relaxed autarky\nand interconnection strictness')
    plot_variables(data.copy(), GENERATION_CAPACITIES, ax, scaling_factor=1e-3)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[1:] + [handles[0]], labels[1:] + [labels[0]])
    ax.set_ylabel("TW")
    ax.get_legend().set_frame_on(False)
    ax.get_legend().set_title(None)
    ax.annotate('a – Generation capacities', xy=[-0.08, 1.05], xycoords='axes fraction', fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    ax = fig.add_subplot(gs[3])
    plot_variables(data.copy(), STORAGE_CAPACITIES, ax, scaling_factor=1e-3)
    ax.set_ylabel("TWh")
    ax.get_legend().remove()
    ax.annotate('b – Storage', xy=[-0.24, 1.05], xycoords='axes fraction', fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    ax = fig.add_subplot(gs[4])
    plot_variables(data.copy(), TRANSMISSION_CAPACITIES, ax)
    ax.set_ylabel("TWkm")
    ax.get_legend().remove()
    ax.annotate('c – Transmission', xy=[-0.48, 1.05], xycoords='axes fraction', fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    ax = fig.add_subplot(gs[5])
    plot_variables(data.copy(), VRES, ax=ax)
    ax.set_ylabel("TWh")
    ax.get_legend().remove()
    ax.annotate('d – Variable renewables', xy=[-0.26, 1.05], xycoords='axes fraction', fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    fig.tight_layout()
    fig.savefig(path_to_output, dpi=600)


def plot_variables(data, variables, ax, scaling_factor=1):
    data["Value"] = data["Value"] * scaling_factor
    variable_order = {var: i for i, var in enumerate(variables.keys())}
    sns.barplot(
        x="Variable",
        y="Value",
        hue="Scale",
        data=data[data["Variable"].isin(variables.keys()) & data["Scenario"].isin(MAIN_SCENARIOS)].replace(variables),
        orient="v",
        palette=PALETTE,
        order=variables.values(),
        hue_order=SCALE_ORDER,
        ax=ax
    )
    x_values = [patch.get_x() + patch.get_width() / 2 for patch in ax.patches]
    for i, (scale, variable) in enumerate(product(SCALE_ORDER.keys(), variable_order.keys())):
        group = data[(data["Scale"] == scale) & (data["Variable"] == variable)]
        min_value = group.Value.min()
        max_value = group.Value.max()
        ax.vlines(
            x=x_values[i],
            ymin=min_value,
            ymax=max_value,
            linewidth=ERROR_BAR_LINEWIDTH,
            color=ERROR_BAR_COLOR
        )
    sns.despine(ax=ax)
    ax.set_xlabel("")


if __name__ == "__main__":
    composition(
        path_to_aggregated_results=snakemake.input.results,
        path_to_output=snakemake.output[0]
    )
