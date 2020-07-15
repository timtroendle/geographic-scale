import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pycountry


PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
HYDRO_TECHS = ["hydro_run_of_river", "hydro_reservoir"]
BIOFUEL_TECHS = ["biofuel"]
STORAGE_TECHS = ["hydrogen", "battery", "pumped_hydro"]
SUPPLY_TECHS = PV_TECHS + WIND_TECHS + HYDRO_TECHS + BIOFUEL_TECHS
TECHNOLOGY_MAP = {
    'open_field_pv': 'solar',
    'roof_mounted_pv': 'solar',
    'wind_offshore': 'wind',
    'wind_onshore_monopoly': 'wind',
    'wind_onshore_competing': 'wind',
    'hydro_run_of_river': 'hydro',
    'hydro_reservoir': 'hydro',
    'biofuel': 'bioenergy'
}

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
COLORS = [
    sns.desaturate(BLUE, 0.85),
    sns.desaturate(GREEN, 0.85),
    sns.desaturate(YELLOW, 0.85),
    sns.desaturate(RED, 0.85)
]
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"


def generation_shares(path_to_results, scenario1, scenario2, path_to_plot):
    shares = read_generation_shares(path_to_results)
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(121)
    (
        shares
        .sel(scenario=scenario1)
        .to_series()
        .reset_index()
        .pivot(columns="technology", index="country", values="carrier_prod")
        .sort_index(ascending=False)
        .loc[:, ["hydro", "bioenergy", "solar", "wind"]]
        .plot(kind="barh", stacked=True, ax=ax, width=0.9, color=COLORS, legend=False)
    )
    ax.annotate(
        f"A - {plot_label(scenario1)}",
        xy=[-0.1, 1.02],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )
    ax.set_ylabel("")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_xlabel("Generation share")

    ax = fig.add_subplot(122)
    (
        shares
        .sel(scenario=scenario2)
        .to_series()
        .reset_index()
        .pivot(columns="technology", index="country", values="carrier_prod")
        .sort_index(ascending=False)
        .loc[:, ["hydro", "bioenergy", "solar", "wind"]]
        .plot(kind="barh", stacked=True, ax=ax, width=0.9, color=COLORS)
    )
    ax.annotate(
        f"B - {plot_label(scenario2)}",
        xy=[-0.1, 1.02],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )
    ax.set_ylabel("")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_xlabel("Generation share")

    leg = fig.legend(
        labels=[t.get_text().capitalize() for t in ax.get_legend().get_texts()],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=4
    )
    leg.draw_frame(False)
    ax.get_legend().remove()

    sns.despine(fig)
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    fig.savefig(path_to_plot, dpi=300, pil_kwargs={"compression": "tiff_lzw"})


def read_generation_shares(path_to_results):
    ds = xr.open_dataset(path_to_results)
    supply = (
        ds
        .carrier_prod
        .groupby(ds.carrier_prod.country_code)
        .sum("locs")
        .sel(techs=SUPPLY_TECHS)
        .groupby(pd.Series(TECHNOLOGY_MAP, name="technology").rename_axis(index="techs").to_xarray())
        .sum("techs")
    )
    supply.coords["country"] = (
        supply["country_code"]
        .to_series()
        .map(lambda country_code: pycountry.countries.lookup(country_code).name)
        .replace("Bosnia and Herzegovina", value="Bosnia") # too long
        .replace("United Kingdom", value="UK") # too long
        .replace("North Macedonia", value="N Macedonia") # too long
    )
    rel_supply = (
        supply
        .to_dataframe()
        .reset_index()
    )
    rel_supply["carrier_prod"] = (
        rel_supply
        .groupby(["country_code", "scenario"])
        .transform(lambda x: x / x.sum())
    )
    return rel_supply.set_index(["country", "technology", "scenario"]).to_xarray()["carrier_prod"]


def parse_scenario_name(scenario_name):
    supply_scale, _, autarky_level, balancing_scale, _ = scenario_name.split("-")
    assert supply_scale in ["regional", "national", "continental"]
    assert balancing_scale in ["regional", "national", "continental"]
    return supply_scale, int(autarky_level), balancing_scale


def plot_label(scenario_name):
    supply_scale, autarky_level, balancing_scale = parse_scenario_name(scenario_name)

    if (autarky_level != 100) or (supply_scale != balancing_scale):
        raise NotImplementedError("Only pure cases implemented.")
    else:
        return f"{supply_scale.capitalize()} supply and balancing"


if __name__ == "__main__":
    generation_shares(
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0],
        scenario1=snakemake.params.scenario1,
        scenario2=snakemake.params.scenario2
    )
