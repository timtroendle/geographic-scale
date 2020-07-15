---
title: Supplemental Information
bibliography:
    - 'literature.bib'
csl: joule.csl
link-citations: True
fignos-caption-name: Figure
tablenos-caption-name: Table
lang: en-GB
date: 2020-07-17-dev
---

# Figure S1: Transmission network

![**Possible locations of transmission capacities.** All lines visualise connections between two regions that can hold transmission capacities. International connections are coloured yellow, all others are coloured blue. The amount of capacities installed on these connections is an output of the optimisation and depends on the considered case.](report/network.png){#fig:network .class  tag="S1"}

<div class="pagebreak"> </div>

# Figure S2: System cost in cases including net imports

![**System cost of nine electricity systems in Europe with national or regional net supply, and continental or national balancing, relative to lowest cost, continental-scale system.** Cost of variations of three cases from Figure 1 in which net imports into the national or regional supply area are allowed to certain degrees (x-axis), from no net imports (0%, corresponds to cases from Figure 1) to imports covering up to 30% of national or regional electricity demand. Cost are relative to the entirely continental case in Figure 1.](report/cost-special-cases.svg){#fig:cost .class tag="S2"}

<div class="pagebreak"> </div>

# Figure S3: System composition of all cases

![**Cost-optimised technology mixes for all cases. a,b,c,d**, Europe-wide installed generation capacities (**a**), storage capacities (**b**), transmission capacities (**c**), and  potential and curtailed generation from variable renewable sources (**d**) for entirely continental-, national-, and regional-scale electricity systems with minimal system cost. The thin bars indicate the range of values including cases with supply on smaller scales. See Tables\ S2 and S3 for numerical results of all cases. **a**, The total excludes all storage capacities but includes capacities of hydropower at or below today's levels (36\ GW run of river, 103\ GW reservoirs). **b**, Europe-wide installed storage capacities for battery and hydrogen storage. Up to 97\ TWh from hydro reservoirs and 1.3\ TWh from pumped hydro storage are not shown. **c**, Europe-wide installed electricity transmission capacities, ignoring transmission capacities within regions. **d**, Europe-wide potential generation and curtailed generation from the variable renewable sources solar, wind and hydro.](report/composition-all.svg){#fig:composition .class tag="S3"}

<div class="pagebreak"> </div>

# Figure S4: Generation shares

![**Technology generation shares for all countries. a,** Generation shares when supply and balancing are continental. **b,** Generation shares when supply and balancing are regional.](report/generation-shares.svg){#fig:generation-shares .class tag="S4"}

<div class="pagebreak"> </div>

# Figure S5: Total Sobol' indices

![**Total Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/total-sobol-all.svg){#fig:total-sobol .class tag="S5"}

<div class="pagebreak"> </div>

# Figure S6: First-order Sobol' indices

![**First-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/first-sobol-all.svg){#fig:first-sobol .class tag="S6"}

<div class="pagebreak"> </div>

# Figure S7: Difference Sobol' indices

![**Total minus first-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/total-minus-first-sobol-all.svg){#fig:higher-sobol .class tag="S7"}

<div class="pagebreak"> </div>

# Figure S8: Use of bioenergy potentials

![**Use of bioenergy potentials.** All use data from the case with highest levels of deployed bioenergy capacity in which supply and balancing are regional. Potentials are based on estimations of available residual material[@RuizCastello:2015] and sum to 2400\ TWh/yr Europe-wide (1100\ TWh/yr usable electricity) for the year 2020 and reference assumptions. **a,** Use of national potentials. **b,** Use of potentials in 497 subnational regions.](report/bioenergy-use.svg){#fig:bioenergy-use .class tag="S8"}

<div class="pagebreak"> </div>

# Figure S9: Distribution of time series across regions

![**Distribution of time series for all regions of the entirely regional-scale electricity system. a,b,c,d,e,f,** Distribution of time series for all regions with low flexibility cost (combined cost of bioenergy, hydrogen, and battery storage) (**a,b,c**) and all regions whose flexibility cost lies within the highest decile (**d,e,f**). **a,b,** Combined wind and solar weekly generation potential time series relative to local demand. Seasonal fluctuations are more pronounced for the case with higher cost due to higher shares of solar electricity. Regions with high solar shares are often urban regions with low or no wind potential. **c,d,** Weekly generation time series from biomass combustion relative to demand. Generation has a more pronounced seasonality and is generally larger in regions with higher flexibility cost. **e,f,** Weekly time series of hydrogen storage levels relative to installed storage capacity. In regions with higher cost, hydrogen storage is used primarily to balance seasonal fluctuations, instead of balancing fluctuations within weeks or months as it is done for lower cost regions, leading to fewer storage cycles and higher cost.](report/timeseries.svg){#fig:timeseries .class tag="S9"}

<div class="pagebreak"> </div>

# Table S1: Input parameter uncertainty

```table
---
caption: 'Uncertain input parameters. For all parameters we assume a uniform distribution. {#tbl:overview-uncertain-parameters tag="S1"}'
alignment: LLLLLL
include: report/overview-uncertain-parameters.csv
include-encoding: UTF-8
markdown: True
---
```

<div class="pagebreak"> </div>

# Table S2: Generation and storage capacities

```table
---
caption: 'Installed capacities of photovoltaics (PV), on- and offshore wind, bioenergy, short-term (battery) and long-term (hydrogen) storage, and relative curtailment of solar, wind, and hydropower for all considered cases. Each case additionally contains fixed hydropower capacities on their current locations: 36 GW run of river, 103 GW / 97 TWh reservoirs, and 54 GW / 1.3 TWh pumped hydro storage. {#tbl:overview-scenario-results-1 tag="S2"}'
alignment: LRRRRRRRR
include: report/overview-scenario-results-1.csv
include-encoding: utf-8
markdown: True
---
```

<div class="pagebreak"> </div>

# Table S3: Transmission capacities

```table
---
caption: 'Installed transmission grid capacity, gross physical electricity flow crossing country borders, and net electricity flow imported by all countries for all cases. {#tbl:overview-scenario-results-2  tag="S3"}'
alignment: LRRR
include: report/overview-scenario-results-2.csv
include-encoding: utf-8
markdown: True
---
```

<div class="pagebreak"> </div>

# Table S4: Biomass feedstocks

```table
---
caption: 'Biomass feedstocks we consider, together with the proxy we use to derive regional from national values. {#tbl:biomass-feedstocks  tag="S4"}'
alignment: LR
include: biofuel-feedstocks.csv
include-encoding: UTF-8
markdown: True
---
```

<div class="pagebreak"> </div>

# Table S5: Technology cost assumptions

```table
---
caption: 'Assumptions on technology cost. ^AC transmission overnight cost is given in €/kW/1000km {#tbl:overview-cost-assumptions tag="S5"}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

<div class="pagebreak"> </div>

# Supplemental Experimental Procedures

## Procedure S1: Impact of net imports into supply area on cost

To confirm the inferior impact of supply options on total cost, we assess the cost of nine further cases, in which we relax the supply scale by permitting net imports to satisfy national or regional electricity demand. This relaxation has a small impact on cost (Figure\ @fig:network), with 10 percentage points or less cost reduction between net self-sufficiency (0 net imports, cases from Figure\ 1) and allowing up to 30% net imports. This reinforces the finding from Figure\ 1: geographic scale has a particularly large impact on cost because of the possibilities for balancing, not mainly because of supply options.

## Procedure S2: Today's transmission capacities

To determine today's transmission capacities, we use an approach described in ref. @Horsch:2018 together with data from the same publication. Ref. @Horsch:2018 uses an extended version of the GridKit tool[@Wiegmans:2016] to find the location, voltage level, and number of circuits of all transmission lines in Europe. They estimate the transmission capacity per line based on the voltage level. In our study area, this leads to ~190\ GW of international (cross-border) transmission capacity, ~1,600\ GW of interregional transmission capacity, and a total of ~340\ TWkm of transmission capacity.

To make the total number comparable to results from our model, we adjust it in the following way: We start by determining the transmission capacity between each pair of regions (in GW). We then ignore the actual lengths of the lines and instead map transmission capacities to transmission lines of our model, which connect centroids of all regions (see Figure\ @fig:network). Using this approach, we estimate today's interregional transmission capacity to be 215\ TWkm.

## Procedure S3: Effect of hydropower model choices on our results

The way in which we model hydropower generation in Europe leads to peculiarities in results on the regional scale. Here, we find the lowest system cost in regions with large hydropower installations. These cost peculiarities are consequences of our design choice: we assume they are amortised. In the following, we justify this model design choice and discuss its relevance for our results.

We assume hydropower capacities are amortised to avoid the need to model overnight cost. Overnight cost of hydropower capacities can vary strongly between projects and are thus difficult to model. Ignoring overnight cost can lead to low levelised cost of electricity in regions with large hydropower capacities. Cost is particularly low when dams provide local flexibility and other forms of more expensive flexibility provision can be avoided. Thus, our model choice leads to low cost in some regions and it also leads to slightly optimistic absolute system cost. But because we fix capacities to current levels for all Europe in all cases, the ignored overnight cost imply no benefit to any case and thus do not affect relative cost.

## Procedure S4: Effect of model simplifications

There are three aspects which our analysis does not consider and which may impact our findings. Some are likely to increase attractiveness of small-scale systems, others are likely to increase attractiveness of large-scale systems.

Most importantly, we do not consider flexibility from electricity demand or from electrifying the heat and transport sectors. These additional flexibilities may be especially beneficial for smaller-scale systems, whose flexibility options are more limited and expensive. However, large electricity systems can and will benefit as well and our sensitivity analysis shows that cost differences are driven largely by the cost of bioenergy; a technology mostly used to balance seasonal fluctuations of solar generation. Only if additional flexibilities can balance seasonal fluctuations, a significant impact on system cost can be expected. This is likely not the case for demand flexibility[@Aryandoust:2017] or transportation[@Brown:2018a], but it may be possible by electrifying the heat sector[@Brown:2018a].

Furthermore, we do not consider ancillary services for the distribution and transmission grids. The provision of ancillary services may be easier, i.e. less costly, for system layouts as we found them for smaller electricity systems, because not only generation but also support infrastructure is more homogeneously dispersed and thus able to provide services like frequency control or black-start everywhere. However, there is no reason to believe that this would change system cost and therefore relative system cost significantly[@Brown:2018].

We also do not model the distribution grid in any way. The cost of the distribution grid is likely to be higher for smaller systems where generation is dispersed more strongly with substantial amounts of generation from roof mounted PV embedded within the distribution grid. However, technical potentials of wind and utility-scale PV are high enough in most regions in Europe so that roof mounted PV is rarely necessary. Thus, cost of the distribution grid may be higher for smaller scales, but only if roof mounted PV is prioritised over utility-scale PV.

## Procedure S5: Temporal resolution of the model

We run the model using a temporal resolution of four hours due to high computational requirements stemming from the high spatial resolution. Our model is spatially resolved into 497 nodes, which represent all first-level administrative units in Europe. Together with the 2,190 timesteps of the four hour temporal resolution, this leads to ~1,100,000 pairs of space and time for which the optimisation algorithm must find a solution. In comparison, models with national spatial and hourly temporal resolution, which are typically used in previous studies[@Gils:2017; @Child:2019; @Zappa:2019], pose a problem consisting of 265,000 pairs of space and time. Due to the polynomial-time complexity of the optimisation algorithm, computation time does not scale linearly with the problem size. In fact, we find that solving a version of our model that is resolved on the national level runs 45 times faster than the original version of our model, although it has only 17 times fewer nodes.

Solving our model requires > 250\ GB of internal memory and > 24h of runtime for each of the 12 cases and each of the 20 high-resolution uncertainty runs we perform. This is only possible using a high performance computing cluster. Increasing the problem size by increasing the temporal resolution is challenging. In fact, we are not aware of any study with larger problem size.

The temporal resolution of our model likely impacts our results. As has been shown before, model results can change significantly with lower temporal resolutions[@Pfenninger:2017]. However, the effect can be expected to be not large for resolutions above six hours and for models that contain inter-temporal constraints like ours[@Pfenninger:2017]. We test the impact on our results by using a version of our model that is spatially resolved on the national level and has a temporal resolution of one hour. We find that total generation capacity remains nearly constant, but that total balancing capacity increases compared to a model that has a temporal resolution of four hours. This result indicates that a higher temporal resolution would likely support rather than contradict our main finding that the balancing scale drives cost more strongly than the supply scale.

## Procedure S6: Effect of disaggregation method of national electricity loads

The regional disaggregation method based on population counts and industry plants that we use yields regional electricity load time series that are likely stronger correlated than in reality. In the following, we discuss the impact this may have on our findings.

Most importantly, the stronger correlation impacts all cases with balancing scale above the regional level in the same way. As long as all regions within a country are connected through transmission lines, relative fluctuations between regions can be compensated through the grid. In all of these cases, we can expect that our model choice leads to a slight over- or underestimation of national transmission capacity but that the bias is similar in all cases.

In the single case with regional-scale balancing, relative fluctuations between regions cannot be compensated through the grid but instead must be handled in each region individually. Here, we can expect that through our method some regions suffer slightly higher cost for balancing, but that other regions in the same country enjoy lower cost of similar levels. The impact for single regions may be noticable, but the impact on total system cost can be considered low.

In summary, there is likely a small impact on the differences between the cases with balancing above the regional level and the case with balancing at the regional level. However, there is no reason to believe that this difference or any other impact on our results is large, as there is no reason to believe that correlations between regions are much lower than the ones resulting from our method.

<div class="pagebreak"> </div>

# Supplemental References
