from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr

RED = "#A01914"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"


def plot_map(path_to_continental_shape, path_to_national_shapes, path_to_regional_shapes,
             path_to_continental_result, path_to_national_result, path_to_plot):
    """Plot maps of results."""
    np.random.seed(123456789)
    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    axes = fig.subplots(2, 2).flatten()
    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.25)
    cmap = sns.light_palette(sns.desaturate(RED, 0.85), reverse=False, as_cmap=True)
    continental = gpd.read_file(path_to_continental_shape).to_crs(EPSG_3035_PROJ4).set_index("id")
    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    regional = gpd.read_file(path_to_regional_shapes).to_crs(EPSG_3035_PROJ4).set_index("id")
    regional_relaxed = regional.copy()

    continental["cost"] = _read_continental_cost(path_to_continental_result)
    national["cost"] = _read_national_cost(path_to_national_result)
    regional["cost"] = np.random.normal(loc=0.18, scale=0.04, size=len(regional))
    regional_relaxed["cost"] = np.random.normal(loc=0.12, scale=0.04, size=len(regional_relaxed))

    _plot_layer(continental, "continental", norm, cmap, axes[0])
    _plot_layer(national, "national", norm, cmap, axes[1])
    _plot_layer(regional, "regional", norm, cmap, axes[2])
    _plot_layer(regional_relaxed, "regional with trade", norm, cmap, axes[3])

    _plot_colorbar(fig, axes, norm, cmap)
    fig.savefig(path_to_plot, dpi=300, transparent=True)


def _plot_layer(units, layer_name, norm, cmap, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
        column="cost",
        vmin=norm.vmin,
        vmax=norm.vmax,
        cmap=cmap,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.10, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{MONEY_UNICODE} ",
        xy=[0.10, 0.85],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color=sns.desaturate(RED, 0.85)
    )
    ax.annotate(layer_name, xy=[0.17, 0.90], xycoords='axes fraction')
    ax.annotate(f"{units.cost.mean():.2f}€/kWh", xy=[0.17, 0.85], xycoords='axes fraction')


def _plot_colorbar(fig, axes, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(["{:.2f}€/kWh".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)


def _read_continental_cost(path_to_continental_result):
    return (xr.open_dataset(path_to_continental_result)["total_levelised_cost"]
              .squeeze(["costs", "carriers"])
              .item())


def _read_national_cost(path_to_national_result):
    results = xr.open_dataset(path_to_national_result)
    cost = results['cost']
    cost = _split_loc_techs(cost).sum(dim=['techs']).squeeze('costs')
    carrier_prod = results['carrier_prod']
    supply_only_carrier_prod = carrier_prod.sel(
        loc_tech_carriers_prod=list(results.loc_tech_carriers_supply_all.values)
    )
    carrier_prod = (
        _split_loc_techs(supply_only_carrier_prod)
        .sum(dim=['timesteps', 'techs'])
        .squeeze(['carriers'])
    )
    return (cost / carrier_prod).to_series()


def _split_loc_techs(data_var, as_='DataArray'):
    # copy from Calliope source code, copyright Calliope authors
    # I am using the copy here because Calliope does not support xarray 0.12.
    # which I need for the data handling. Thus, I cannot import calliope here.
    """
    Get a DataArray with locations technologies, and possibly carriers
    split into separate coordinates.
    Parameters
    ----------
    data_var : xarray DataArray
        Variable from Calliope model_data, to split loc_techs dimension
    as_ : string
        'DataArray' to return xarray DataArray or 'Series' to return pandas
        Series with dimensions as a MultiIndex
    Returns
    -------
    updated_data_var : xarray DataArray of pandas Series
    """

    # Separately find the loc_techs(_carriers) dimension and all other dimensions
    loc_tech_dim = [i for i in data_var.dims if 'loc_tech' in i]
    if not loc_tech_dim:
        loc_tech_dim = [i for i in data_var.dims if 'loc_carrier' in i]
    non_loc_tech_dims = list(set(data_var.dims).difference(loc_tech_dim))

    if not loc_tech_dim:
        if as_ == 'Series':
            return data_var.to_series()
        elif as_ == 'DataArray':
            return data_var
        else:
            raise ValueError('`as_` must be `DataArray` or `Series`, '
                             'but `{}` given'.format(as_))

    elif len(loc_tech_dim) > 1:
        raise ValueError("Cannot split loc_techs or loc_techs_carrier dimension "
                         "for DataArray {}".format(data_var.name))

    loc_tech_dim = loc_tech_dim[0]
    # xr.Datarray -> pd.Series allows for string operations
    data_var_df = data_var.to_series().unstack(non_loc_tech_dims)
    index_list = data_var_df.index.str.split('::').tolist()

    # carrier_prod, carrier_con, and carrier_export will return an index_list
    # of size 3, all others will be an index list of size 2
    possible_names = ['loc', 'tech', 'carrier']
    names = [i + 's' for i in possible_names if i in loc_tech_dim]

    data_var_df.index = pd.MultiIndex.from_tuples(index_list, names=names)

    # If there were no other dimensions other than loc_techs(_carriers) then
    # nothing was unstacked on creating data_var_df, so nothing is stacked now
    if isinstance(data_var_df, pd.Series):
        data_var_series = data_var_df
    else:
        data_var_series = data_var_df.stack(non_loc_tech_dims)

    if as_ == "Series":
        return data_var_series

    elif as_ == "DataArray":
        updated_data_var = xr.DataArray.from_series(data_var_series)
        updated_data_var.attrs = data_var.attrs
        updated_data_var.name = data_var.name

        return updated_data_var

    else:
        raise ValueError('`as_` must be `DataArray` or `Series`, '
                         'but `{}` given'.format(as_))


if __name__ == "__main__":
    plot_map(
        path_to_continental_shape=snakemake.input.continental_shape,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_continental_result=snakemake.input.continental_result,
        path_to_national_result=snakemake.input.national_result,
        path_to_plot=snakemake.output[0]
    )
