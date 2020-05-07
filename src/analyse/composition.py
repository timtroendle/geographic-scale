from dataclasses import dataclass
import io

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import xarray as xr


BLUE = "#4F6DB8"
RED = "#A01914"
ANTHRACITE = "#424242"
SYSTEM_SCALE_COLOR = "k"
AUTARKY_EXTENT_COLOR = "k"
PALETTE = sns.light_palette(BLUE)
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

ONSHORE_WIND_TECHS = ["wind_onshore_competing", "wind_onshore_monopoly"]
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
VRES_TECHS_WITHOUT_HYDRO = ONSHORE_WIND_TECHS + PV_TECHS + ["wind_offshore"]
VRES_TECHS = VRES_TECHS_WITHOUT_HYDRO + ["hydro_run_of_river"]
SUPPLY_TECHS = VRES_TECHS + ["hydro_reservoir"]
STORAGE_TECHS = ["battery", "hydrogen", "pumped_hydro"]
BALANCING_TECHS = STORAGE_TECHS + ["biofuel"]


AUTARKY_LEVEL_MAP = {
    "100": "0%",
    "85": "≤15%",
    "70": "≤30%"
}


@dataclass
class PlotData:
    data: pd.DataFrame
    cbar_label: str
    fmt: str = '.3g'
    annotation_scale: float = None


def composition(path_to_aggregated_results_csv, path_to_aggregated_results_nc, path_to_output,
                transmission_capacity_today_twkm, crossborder_capacity_today_tw):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 7))
    axes = fig.subplots(2, 2).flatten()

    plot_datas = read_plot_datas(
        path_to_aggregated_results_nc,
        path_to_aggregated_results_csv,
        transmission_capacity_today_twkm,
        crossborder_capacity_today_tw
    )

    for ax, cbar_ax, plot_data in zip(axes, range(4), plot_datas):
        base_case_box(plot_data, ax, cbar_ax)

    plt.subplots_adjust(
        bottom=0.08,
        top=0.93,
        wspace=0.3,
        hspace=0.2
    )
    fig.savefig(path_to_output, dpi=300)


def read_plot_datas(path_to_aggregated_results_nc, path_to_aggregated_results_csv,
                    transmission_capacity_today_twkm, crossborder_capacity_today_tw):
    return [
        PlotData(
            data=read_total_supply_capacity(path_to_aggregated_results_nc),
            cbar_label="a - Supply capacity [TW]",
            fmt=".1f"
        ),
        PlotData(
            data=read_biostor_capacity(path_to_aggregated_results_nc),
            cbar_label="b - Balancing capacity [TW]",
            fmt=".2f"
        ),
        PlotData(
            data=read_transmission_capacity(path_to_aggregated_results_csv),
            cbar_label="c - Transmission capacity [TWkm]",
            fmt=".0f",
            annotation_scale=transmission_capacity_today_twkm,
        ),
        PlotData(
            data=read_international_transmission_capacity(path_to_aggregated_results_csv),
            cbar_label="d - Cross-border transmission capacity [TW]",
            fmt=".2g",
            annotation_scale=crossborder_capacity_today_tw
        ),
    ]


def base_case_box(plot_data, ax, cbar_ax):
    results = plot_data.data
    cbar_ticks = np.linspace(
        start=results["cost"].min(),
        stop=results["cost"].max(),
        num=len(PALETTE) + 1
    )
    heatmap_data = (
        results[results.autarky_degree == "0%"]
        .pivot(index="autarky_layer", columns="grid_scale", values="cost")
        .reindex(index=["Continental", "National", "Regional"])
    )
    if plot_data.annotation_scale:
        annot = heatmap_data.applymap(
            lambda x: f"{{:{plot_data.fmt}}}\n({x / plot_data.annotation_scale:.1f})".format(x)
        )
        fmt = "s"
    else:
        annot = True
        fmt = plot_data.fmt
    sns.heatmap(
        heatmap_data,
        annot=annot,
        cbar=True,
        cbar_kws={
            "ticks": cbar_ticks,
            "format": f"%{plot_data.fmt}",
            "aspect": 30,
            "shrink": 0.8
        },
        cmap=PALETTE,
        vmin=results["cost"].min(),
        vmax=results["cost"].max(),
        linewidth=1.25,
        square=True,
        ax=ax,
        fmt=fmt
    )
    ax.set_xlabel("Balancing scale")
    ax.set_ylabel("Supply scale")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=90, va="center")
    ax.annotate(plot_data.cbar_label, xy=[-0.2, 1.1], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def read_transmission_capacity(path_to_agregrated_results):
    da = (
        pd
        .read_csv(path_to_agregrated_results, index_col=[0, 1])
        .to_xarray()
        .rename({"Scenario": "scenario"})
        .sel(Variable="Capacity|Transmission")
        .Value
    )
    return bring_into_form(da)


def read_international_transmission_capacity(path_to_agregrated_results):
    da = (
        pd
        .read_csv(path_to_agregrated_results, index_col=[0, 1])
        .to_xarray()
        .rename({"Scenario": "scenario"})
        .sel(Variable="Capacity|Gross import national level")
        .Value
    ) / 1e3 # to TW
    return bring_into_form(da)


def read_total_supply_capacity(path_to_agregrated_results):
    da = (
        xr
        .open_dataset(path_to_agregrated_results)
        .energy_cap
        .sel(techs=SUPPLY_TECHS)
        .sum(["locs", "techs"])
    ) / 1e6 # to TW
    return bring_into_form(da)


def read_wind_capacity(path_to_agregrated_results):
    da = (
        xr
        .open_dataset(path_to_agregrated_results)
        .energy_cap
        .sel(techs=ONSHORE_WIND_TECHS + ["wind_offshore"])
        .sum(["locs", "techs"])
    ) / 1e6 # to TW
    return bring_into_form(da)


def read_biostor_capacity(path_to_agregrated_results):
    da = (
        xr
        .open_dataset(path_to_agregrated_results)
        .energy_cap
        .sel(techs=BALANCING_TECHS)
        .sum(["locs", "techs"])
    ) / 1e6 # to TW
    return bring_into_form(da)


def bring_into_form(da):
    results = (
        pd
        .read_csv(io.StringIO(DATA_INDEX))
        .set_index(["autarky_layer", "grid_scale", "autarky_degree"])
    )
    for scenario in da.scenario:
        scenario = scenario.item()
        autarky_layer, autarky_level, grid_size = parse_scenario_name(scenario)
        autarky_level = AUTARKY_LEVEL_MAP[autarky_level]
        results.loc[autarky_layer, grid_size, autarky_level] = da.sel(scenario=scenario).item()
    return results.reset_index()


def parse_scenario_name(scenario_name):
    autarky_layer, _, autarky_level, grid_size, _ = scenario_name.split("-")
    assert autarky_layer in ["regional", "national", "continental"]
    assert grid_size in ["regional", "national", "continental"]
    return autarky_layer.capitalize(), autarky_level, grid_size.capitalize()


if __name__ == "__main__":
    composition(
        path_to_aggregated_results_csv=snakemake.input.results_csv,
        path_to_aggregated_results_nc=snakemake.input.results_nc,
        transmission_capacity_today_twkm=snakemake.params.transmission_today_twkm,
        crossborder_capacity_today_tw=snakemake.params.crossborder_today_tw,
        path_to_output=snakemake.output[0]
    )
