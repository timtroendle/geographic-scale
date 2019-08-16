from dataclasses import dataclass
from itertools import chain

import calliope
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"
GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
PALETTE = sns.light_palette(GREEN, n_colors=5, reverse=True)

SUPPLY_TECHS = [
    'wind_onshore_monopoly', 'biofuel',
    'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
    'wind_onshore_competing', 'open_field_pv', 'hydro_run_of_river'
]
WIND_TECHS = ["wind_onshore_competing", "wind_onshore_monopoly", "wind_offshore"]
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_AND_PV = WIND_TECHS + PV_TECHS
VRES_TECHS = WIND_AND_PV + ["hydro_run_of_river"]
DEMAND = "demand_elec"
STORAGE_TECHS = ["hydrogen", "battery", "pumped_hydro", "hydro_reservoir"]
HYDROGEN = "hydrogen"
NS_TO_H = 1e-9 / 3600


@dataclass
class PlotData:
    name: str
    da: xr.DataArray
    mask: xr.DataArray
    ylabel: str


def timeseries(path_to_result, path_to_units, scaling_factors, connected_regions,
               unit_lcoe_threshold, biofuel_lcoe_threshold, hydrogen_lcos_threshold,
               resolution, path_to_plot, path_to_displayed_units):
    """Create plots of timeseries for regional result."""
    plot_datas = read_plot_data(path_to_result, scaling_factors, connected_regions,
                                unit_lcoe_threshold, biofuel_lcoe_threshold, hydrogen_lcos_threshold)
    fig = plot_timeseries(plot_datas, resolution)
    fig.savefig(path_to_plot, dpi=600)
    units = mask_all_units(path_to_units, plot_datas)
    units.to_netcdf(path_to_displayed_units)


def read_plot_data(path_to_result, scaling_factors, connected_regions, unit_lcoe_threshold,
                   biofuel_lcoe_threshold, hydrogen_lcos_threshold):
    model = calliope.read_netcdf(path_to_result)
    gen = postprocess_combined_regions(
        model.get_formatted_array("carrier_prod").squeeze("carriers") / scaling_factors["power"],
        connected_regions.items()
    )
    cost = postprocess_combined_regions(
        model.get_formatted_array("cost").squeeze('costs') / scaling_factors["monetary"],
        connected_regions.items()
    )
    cap = postprocess_combined_regions(
        model.get_formatted_array("energy_cap") / scaling_factors["power"],
        connected_regions.items()
    )
    stor = postprocess_combined_regions(
        model.get_formatted_array("storage_cap") / scaling_factors["power"],
        connected_regions.items()
    )
    e_stor = postprocess_combined_regions(
        model.get_formatted_array("storage") / scaling_factors["power"],
        connected_regions.items()
    )
    dem = postprocess_combined_regions(
        model.get_formatted_array("carrier_con").squeeze("carriers").sel(techs=DEMAND) / scaling_factors["power"],
        connected_regions.items()
    )
    resource = postprocess_combined_regions(
        model.get_formatted_array("resource").sel(techs=VRES_TECHS),
        connected_regions.items()
    )
    pot = (resource * cap.sel(techs=VRES_TECHS))
    unit_lcoe = (cost.sum("techs") / -dem.sum(["timesteps"]))
    tech_lcoes = (cost.sel(techs=SUPPLY_TECHS) / gen.sel(techs=SUPPLY_TECHS).sum(["timesteps"]))
    lcos = (cost.sel(techs=STORAGE_TECHS) / gen.sel(techs=STORAGE_TECHS).sum(["timesteps"]))
    resolution_in_hours = (gen.timesteps[1].item() - gen.timesteps[0].item()) * NS_TO_H
    return [
        PlotData(
            name="wind \nand \nsolar",
            ylabel="relative generation [-]",
            da=(pot.sel(techs=WIND_AND_PV).sum("techs")) / -dem,
            mask=unit_lcoe >= unit_lcoe_threshold
        ),
        PlotData(
            name="bioenergy",
            ylabel="relative generation [-]",
            da=gen.sel(techs="biofuel") / resolution_in_hours / cap.sel(techs="biofuel"),
            mask=tech_lcoes.sel(techs="biofuel") > biofuel_lcoe_threshold
        ),
        PlotData(
            name="hydrogen",
            ylabel="storage level [-]",
            da=e_stor.sel(techs=HYDROGEN) / stor.sel(techs=HYDROGEN),
            mask=lcos.sel(techs=HYDROGEN) >= hydrogen_lcos_threshold
        )
    ]


def plot_timeseries(plot_datas, resolution):
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 8))
    axes = fig.subplots(len(plot_datas), 2, sharex=True)
    plt.subplots_adjust(wspace=0)
    panel_ids = "abcdefghi"

    for i, plot_data in enumerate(plot_datas):
        axes[i][0].get_shared_y_axes().join(axes[i][0], axes[i][1])
        draw_areas(
            plot_data.da[~plot_data.mask].resample(timesteps=resolution).mean(dim="timesteps"),
            "locs",
            "timesteps",
            ax=axes[i][0],
            legend=True,
            pal=PALETTE
        )
        draw_areas(
            plot_data.da[plot_data.mask].resample(timesteps=resolution).mean(dim="timesteps"),
            "locs",
            "timesteps",
            ax=axes[i][1],
            legend=True,
            pal=PALETTE
        )
        axes[i][0].set_ylabel(plot_data.ylabel)
        axes[i][1].set_yticks([])
        if i is not len(plot_datas) - 1:
            axes[i][0].set_xlabel("")
            axes[i][1].set_xlabel("")
            axes[i][0].set_xticks([])
            axes[i][1].set_xticks([])
        sns.despine(ax=axes[i][0], top=True, right=True)
        sns.despine(ax=axes[i][1], left=True)
        axes[i][0].annotate(panel_ids[2 * i], xy=[-0.05, 1.05], xycoords='axes fraction',
                            fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
        axes[i][1].annotate(panel_ids[2 * i + 1], xy=[-0.05, 1.05], xycoords='axes fraction',
                            fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
        axes[i][0].annotate(plot_data.name, xy=(-0.2, 0.5), xycoords="axes fraction",
                            size='large', ha='right', va='center')

    leg = fig.legend(
        labels=[t.get_text() for t in axes[0][0].get_legend().get_texts()],
        loc="upper right",
        bbox_to_anchor=(0.87, 1.00),
        ncol=4
    )
    leg.draw_frame(False)
    for row in axes:
        for ax in row:
            ax.get_legend().remove()

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, top=0.90)
    return fig


def mask_all_units(path_to_units, plot_datas):
    units = (
        pd
        .read_csv(path_to_units, index_col=0)
        .rename(index=lambda x: x.replace(".", "-"))
        .rename_axis(index="locs")
        .to_xarray()
    )
    for i, plot_data in enumerate(plot_datas):
        units[f"in_mask{i}"] = plot_data.mask
    return units


def draw_areas(da, dim, time_dim, legend=True, ax=None, pal=None, lw=2.5, **kwargs):
    if not pal:
        pal = sns.color_palette("Reds_r")

    da.median(dim).plot(
        lw=lw,
        color=pal[0],
        label="50%",
        ax=ax
    )
    ax.set_title("")

    ax.fill_between(da[time_dim].values, da.min(dim), da.max(dim), color=pal[3])
    ax.fill_between(da[time_dim].values, da.quantile(0.1, dim).values, da.quantile(0.9, dim).values, color=pal[2])
    ax.fill_between(da[time_dim].values, da.quantile(0.25, dim).values, da.quantile(0.75, dim).values, color=pal[1])

    # Add legend entries
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, facecolor=pal[0], label='25% - 75%'))
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, facecolor=pal[1], label='10% - 90%'))
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, facecolor=pal[2], label='Min - Max'))

    if legend:
        leg = ax.legend(
            loc='lower center',
            bbox_to_anchor=(0, 0.85, 1, 1),
            ncol=4)
        leg.draw_frame(False)

    ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
    return ax


def postprocess_combined_regions(da, region_pairs):
    region_pairs = list(region_pairs)
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
    timeseries(
        path_to_result=snakemake.input.result,
        path_to_units=snakemake.input.units,
        path_to_plot=snakemake.output.plot,
        path_to_displayed_units=snakemake.output.units,
        scaling_factors=snakemake.params.scaling_factors,
        unit_lcoe_threshold=snakemake.params.unit_lcoe,
        biofuel_lcoe_threshold=snakemake.params.biofuel_lcoe,
        hydrogen_lcos_threshold=snakemake.params.hydrogen_lcos,
        resolution=snakemake.params.resolution,
        connected_regions=snakemake.params.connected_regions
    )
