import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import pycountry


GREEN = "#679436"
COLOR = GREEN
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"


def bioenergy_use(path_to_results, path_to_potentials, efficiency, scenario, path_to_plot):
    potentials = read_potentials(path_to_potentials, efficiency)
    use = read_use(path_to_results, scenario)

    nat_use = use.groupby("country").sum("locs")
    nat_potentials = potentials.groupby(use.country).sum("locs")

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(211)
    (
        (nat_use / nat_potentials)
        .to_series()
        .sort_values(ascending=False)
        .plot(kind="bar", width=0.9, color=COLOR, ax=ax)
    )
    ax.set_xlabel("")
    ax.set_ylabel("Share of bioenergy potential used")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_ylim(0, 1)
    ax.annotate(
        f"A - National use of bioenergy potentials",
        xy=[-0.05, 1.05],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )

    ax = fig.add_subplot(212)
    (
        (use / potentials)
        .to_dataframe(name="regional").hist(bins=20, color=COLOR, ax=ax)
    )
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_xlim(0, 1)
    ax.set_title("")
    ax.set_xlabel("Share of bioenergy potential used")
    ax.set_ylabel("Frequency")
    ax.annotate(
        f"B - Regional use of bioenergy potentials",
        xy=[-0.05, 1.05],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )

    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(path_to_plot, dpi=300, pil_kwargs={"compression": "tiff_lzw"})


def read_potentials(path_to_potentials, efficiency):
    return (
        pd
        .read_csv(path_to_potentials, index_col=0)
        .rename(index=lambda loc: loc.replace(".", "-"))
        .rename_axis(index="locs")
        .mul(efficiency)
        .to_xarray()
        .biofuel_potential_mwh_per_year
        .sortby("locs")
    )


def read_use(path_to_results, scenario):
    ds = xr.open_dataset(path_to_results)
    use = ds.carrier_prod.sel(techs="biofuel").sortby("locs").sel(scenario=scenario)
    use.coords["country"] = (
        use["country_code"]
        .to_series()
        .map(lambda country_code: pycountry.countries.lookup(country_code).name)
        .replace("Bosnia and Herzegovina", value="Bosnia") # too long
        .replace("United Kingdom", value="UK") # too long
        .replace("North Macedonia", value="N Macedonia") # too long
    )
    return use


if __name__ == "__main__":
    bioenergy_use(
        path_to_results=snakemake.input.results,
        path_to_potentials=snakemake.input.potentials,
        path_to_plot=snakemake.output[0],
        scenario=snakemake.params.scenario,
        efficiency=snakemake.params.efficiency
    )
