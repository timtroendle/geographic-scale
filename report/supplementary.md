---
title: 'Fully renewable European electricity system designs to balance capacity expansion, transmission needs, and cost'
subtitle: 'Supplemental Information'
author:
    - Tim Tröndle (Lead Contact)
    - Johan Lilliestam
    - Stefano Marelli
    - Stefan Pfenninger
institute: IASS Potsdam and ETH Zürich
bibliography:
    - 'literature.bib'
csl: joule.csl
link-citations: True
fignos-caption-name: Figure
tablenos-caption-name: Table
date: 2019-12-04-dev
---

This document contains supplemental information for:

Tim Tröndle, Johan Lilliestam, Stefano Marelli and Stefan Pfenninger (<mark>YEAR</mark>). Fully renewable European electricity system designs to balance capacity expansion, transmission needs, and cost. DOI: <mark>XXXXXXXXXXXXX</mark>

<div class="pagebreak"> </div>

# Note S1: Code and data availability

The code and reproducible workflow used to perform this analysis is available online <mark>(for review process provided within submission)</mark>.

Furthermore, the resulting data is available online as well <mark>(for review process provided within submission)</mark>.

# Note S2: Effect of hydropower model choices on our results

The way in which we model hydropower generation in Europe leads to peculiarities in results on the regional scale. Here, we find the lowest and the highest system cost in regions with large hydropower installations. These cost peculiarities are consequences of two model design choices: we keep hydropower capacities fixed at today's level and we assume they are amortised. In the following, we justify these two model design choices and discuss their relevance for our results.

We assume hydropower capacities are amortised to avoid the need to model overnight cost. Overnight cost of hydropower capacities can vary strongly between projects and are thus difficult to model. Ignoring overnight cost can lead to low levelised cost of electricity in regions with large hydropower capacities. Cost is particularly low when dams provide local flexibility and other forms of more expensive flexibility provision can be avoided. Thus, our model choice leads to low cost in some regions and it also leads to slightly optimistic absolute system cost. But because we fix capacities to current levels for all Europe in all cases, the ignored overnight cost imply no benefit to any case and thus do not affect relative cost.

We fix hydropower capacities to current levels because significant capacity expansion in Europe is unlikely[@LacalArantegui:2014]. In regions in which hydropower generation exceeds local electricity demand largely, this model choice can lead to high cost in the regional-scale case. However, on average, the impact is small. While on the continental scale levelised cost of electricity of hydropower ranges from 31 to 55 EUR per MWh, it ranges from 33 to 61 EUR per MWh on the regional scale. In fact, 35% of the hydropower potential on the regional scale is curtailed. This corresponds to 2% of total electricity system cost. Thus, if we allowed for capacity reduction on the regional scale, its cost could reduce by up to 2%. This magnitude has no significant affect on our main results.

# Note S3: Effect of model simplifications

There are three aspects which our analysis does not consider and which may impact our findings. Some are likely to increase attractiveness of small-scale systems, others are likely to increase attractiveness of large-scale systems.

Most importantly, we do not consider flexibility from electricity demand or from electrifying the heat and transport sectors. These additional flexibilities may be especially beneficial for smaller-scale systems, whose flexibility options are more limited and expensive. However, large electricity systems can and will benefit as well and our sensitivity analysis shows that cost differences are driven largely by the cost of bioenergy; a technology mostly used to balance seasonal fluctuations of solar generation. Only if additional flexibilities can balance seasonal fluctuations, a significant impact on system cost can be expected. This is likely not the case for demand flexibility[@Aryandoust:2017] or transportation[@Brown:2018a], but it may be possible by electrifying the heat sector[@Brown:2018a].

Furthermore, we do not consider ancillary services for the distribution and transmission grids. The provision of ancillary services may be easier, i.e. less costly, for system layouts as we found them for smaller electricity systems, because not only generation but also support infrastructure is more homogeneously dispersed and thus able to provide services like frequency control or black-start everywhere. However, there is no reason to believe that this would change system cost and therefore relative system cost significantly[@Brown:2018].

We also do not model the distribution grid in any way. The cost of the distribution grid is likely to be higher for smaller systems where generation is dispersed more strongly with substantial amounts of generation from roof mounted PV embedded within the distribution grid. However, technical potentials of wind and utility-scale PV are high enough in most regions in Europe so that roof mounted PV is rarely necessary. Thus, cost of the distribution grid may be higher for smaller scales, but only if roof mounted PV is prioritised over utility-scale PV.

<div class="pagebreak"> </div>

# Figure S1: Total Sobol' indices

![**Total Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/total-sobol-all.svg){#fig:total-sobol .class tag="S1"}

<div class="pagebreak"> </div>

# Figure S2: First-order Sobol' indices

![**First-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/first-sobol-all.svg){#fig:first-sobol .class tag="S2"}

<div class="pagebreak"> </div>

# Figure S3: Difference Sobol' indices

![**Total minus first-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of entirely continental-, national-, and regional-scale electricity systems. a,b,c,** Sobol' indices of input parameters considering total system cost and total installed capacities (x-axis) of the continental- (**a**), national- (**b**), and regional-scale (**c**) systems. The y-axis shows the twelve input parameters included in the uncertainty analysis. **d,** Sobol' indices of input parameters considering difference in system cost between the continental- and national-scale systems. The x-axis shows the twelve input parameters included in the uncertainty analysis. The y-axis shows the model-wide result variables for which continental to national scale differences are compared.](report/total-minus-first-sobol-all.svg){#fig:higher-sobol .class tag="S3"}

<div class="pagebreak"> </div>

# Figure S4: Transmission network

![**Possible locations of transmission capacities.** All lines visualise connections between two regions that can hold transmission capacities. International connections are coloured yellow, all others are coloured blue. The amount of capacities installed on these connections is an output of the optimisation and depends on the considered case.](report/network.png){#fig:network .class  tag="S4"}

<div class="pagebreak"> </div>

# Figure S5: Distribution of time series across regions

![**Distribution of time series for all regions of the entirely regional-scale electricity system. a,b,c,d,e,f,** Distribution of time series for all regions with low flexibility cost (combined cost of bioenergy, hydrogen, and battery storage) (**a,b,c**) and all regions whose flexibility cost lies within the highest decile (**d,e,f**). **a,b,** Combined wind and solar weekly generation potential time series relative to local demand. Seasonal fluctuations are more pronounced for the case with higher cost due to higher shares of solar electricity. Regions with high solar shares are often urban regions with low or no wind potential. **c,d,** Weekly generation time series from biomass combustion relative to demand. Generation has a more pronounced seasonality and is generally larger in regions with higher flexibility cost. **e,f,** Weekly time series of hydrogen storage levels relative to installed storage capacity. In regions with higher cost, hydrogen storage is used primarily to balance seasonal fluctuations, instead of balancing fluctuations within weeks or months as it is done for lower cost regions, leading to fewer storage cycles and higher cost.](report/timeseries.svg){#fig:timeseries .class tag="S5"}

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

# Supplemental References
