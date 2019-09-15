from itertools import chain
from pathlib import Path

import xarray as xr
import numpy as np
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

RED = "#A01914"
BLUE = "#4F6DB8"
CMAP = 'RdBu_r'
NORM = matplotlib.colors.Normalize(vmin=-1, vmax=3)
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"

CONTINENTAL_SCENARIO = "continental-autarky-100-continental-grid"
NATIONAL_SCENARIO = "national-autarky-100-national-grid"
REGIONAL_SCENARIO = "regional-autarky-100-regional-grid"
DEMAND_TECH = "demand_elec"


def plot_map(path_to_national_shapes, path_to_regional_shapes, connected_regions,
             path_to_aggregated_results, path_to_shapes, path_to_plot):
    """Plot maps of results."""
    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()
    shapes = (gpd.read_file(path_to_shapes)
                 .to_crs(EPSG_3035_PROJ4)
                 .rename(columns={"id": "locs"})
                 .set_index("locs")
                 .rename(index=lambda idx: idx.replace(".", "-")))

    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("country_code")
    regional = (gpd.read_file(path_to_regional_shapes)
                   .to_crs(EPSG_3035_PROJ4)
                   .set_index("id")
                   .rename(index=lambda idx: idx.replace(".", "-")))

    cost_reader = CostReader(
        path_to_aggregated_results=path_to_aggregated_results,
        connected_regions=connected_regions,
        country_codes=shapes.country_code.to_xarray()
    )
    national["cost"] = cost_reader.read(scenario=NATIONAL_SCENARIO, level="national")
    regional["cost"] = cost_reader.read(scenario=REGIONAL_SCENARIO, level="regional")

    base_costs = cost_reader.read(scenario=CONTINENTAL_SCENARIO, level="continental")
    national["cost"] = national["cost"] / base_costs
    regional["cost"] = regional["cost"] / base_costs

    total_costs_national = cost_reader.read(scenario=NATIONAL_SCENARIO, level="continental") / base_costs
    total_costs_regional = cost_reader.read(scenario=REGIONAL_SCENARIO, level="continental") / base_costs

    _plot_layer(national, total_costs_national, "National", "a", axes[0])
    _plot_layer(regional, total_costs_regional, "Regional", "b", axes[1])

    _plot_colorbar(fig, axes)
    fig.savefig(path_to_plot, dpi=600, transparent=False)


def _plot_layer(units, total_cost, layer_name, panel_id, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
        column="cost",
        vmin=NORM.vmin,
        vmax=NORM.vmax,
        cmap=CMAP,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.19, 0.95],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{MONEY_UNICODE} ",
        xy=[0.19, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color=sns.desaturate(RED, 0.85)
    )
    ax.annotate(layer_name, xy=[0.26, 0.95], xycoords='axes fraction')
    ax.annotate(f"{total_cost:.1f}", xy=[0.26, 0.90], xycoords='axes fraction')
    ax.annotate(panel_id, xy=[0.10, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def _plot_colorbar(fig, axes):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cmap = s_m.get_cmap()
    colors = cmap(np.linspace(1 / 4, 1, cmap.N))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_jet', colors)
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=0, vmax=3))
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(["0", "1", "2", "≥3"])
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Relative system cost', rotation=90)


class CostReader:

    def __init__(self, path_to_aggregated_results, country_codes, connected_regions):
        self.__results = xr.open_dataset(path_to_aggregated_results)
        self.__country_codes = country_codes
        self.__connected_regions = connected_regions

    def read(self, scenario, level):
        assert level in ["continental", "national", "regional"]
        cost = self.__results.cost.sel(scenario=scenario).sum("techs")
        demand = -self.__results.energy_cap.sel(techs=DEMAND_TECH, scenario=scenario)
        if level == "continental":
            cost = cost.sum(dim="locs").item()
            demand = demand.sum(dim="locs").item()
        elif level == "national":
            cost = (
                cost
                .groupby(self.__country_codes.sortby(cost.locs))
                .sum("locs")
                .to_series()
            )
            demand = (
                demand
                .groupby(self.__country_codes.sortby(demand.locs))
                .sum("locs")
                .to_series()
            )
        elif level == "regional":
            cost = self._postprocess_combined_regions(cost).to_series()
            demand = self._postprocess_combined_regions(demand).to_series()
        return (cost / demand)

    def _postprocess_combined_regions(self, da):
        region_pairs = list(self.__connected_regions.items())
        regions = list(chain(*region_pairs))
        if all([region in da.locs for region in regions]):
            # config is valid, and resolution is regional. apply config
            for region1, region2 in region_pairs:
                sum_data = da.sel(locs=region1) + da.sel(locs=region2)
                da.loc[dict(locs=region1)] = sum_data
                da.loc[dict(locs=region2)] = sum_data
            return da
        elif all([region not in da.locs for region in regions]):
            # config is valid, but resolution is not regional. do nothing
            return da
        else:
            raise ValueError("Config of connected regions is invalid.")


if __name__ == "__main__":
    plot_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_aggregated_results=snakemake.input.results,
        connected_regions=snakemake.params.connected_regions,
        path_to_plot=snakemake.output[0]
    )
