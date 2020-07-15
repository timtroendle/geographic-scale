import calliope
import pandas as pd


def time_diff(path_to_low_res, path_to_high_res, path_to_output):
    lowres = calliope.read_netcdf(path_to_low_res)
    highres = calliope.read_netcdf(path_to_high_res)

    energy_cap_lowres = (
        lowres
        .get_formatted_array("energy_cap")
        .sum("locs")
        .sel(techs=[tech.item() for tech in lowres.get_formatted_array("energy_cap").techs
                    if "ac_transmission" not in str(tech)])
        .to_series()
    )
    energy_cap_highres = (
        highres
        .get_formatted_array("energy_cap")
        .sum("locs")
        .sel(techs=[tech.item() for tech in highres.get_formatted_array("energy_cap").techs
                    if "ac_transmission" not in str(tech)])
        .to_series()
    )
    df = pd.DataFrame(
        data={
            "high_res": energy_cap_highres,
            "low_res": energy_cap_lowres
        }
    ).T
    df["number_timesteps"] = [model.results.timesteps.shape[0] for model in (highres, lowres)]
    df["cost"] = [model.get_formatted_array("cost").sum().item() for model in (highres, lowres)]
    df.T.to_csv(path_to_output, index=True, header=True)


if __name__ == "__main__":
    time_diff(
        path_to_low_res=snakemake.input.lowres,
        path_to_high_res=snakemake.input.highres,
        path_to_output=snakemake.output[0]
    )
