---
title: 'Cost and design of fully renewable electricity supply on continental, national, and regional scales in Europe'
subtitle: 'Supplementary Information'
author:
    - Tim Tröndle
    - Johan Lilliestam
    - Stefano Marelli
    - Stefan Pfenninger
institute: IASS Potsdam and ETH Zürich
bibliography:
    - 'literature.bib'
csl: nature.csl
link-citations: True
fignos-caption-name: Supplementary Figure
tablenos-caption-name: Supplementary Table
date: draft-07-dev
---

This document contains supplementary information for:

Tim Tröndle, Johan Lilliestam, Stefano Marelli and Stefan Pfenninger (YEAR). Cost and design of fully renewable electricity supply on continental, national, and regional scales in Europe. DOI: XXXXXXXXXXXXX

# Supplementary Note 1: Code and data availability

The code and reproducible workflow used to perform this analysis is available online TODO ADD.

Furthermore, the resulting data is available online as well TODO ADD.

# Supplementary Note 2: Effect of model simplifications

There are four aspects which our analysis does not consider but which may impact our findings. Some are likely to be advantageous for smaller systems, others are likely to be advantageous for larger systems.

Most importantly, we do not consider flexibility from electricity demand or from coupling the electricity sector with the heat and transport sectors. These additional flexibilities may be especially beneficial for smaller-scale systems, whose flexibility options are more limited. However, large electricity systems can and will benefit as well and our sensitivity analysis shows that cost differences are driven largely by the cost of bioenergy; a technology mostly used to balance seasonal fluctuations of solar generation. So only if these flexibilities can balance seasonal fluctuations, a significant impact on system cost can be expected. This is likely not the case for demand flexibility [@Aryandoust:2017] or transportation [@Brown:2018a], but it may be possible by coupling the electricity and the heat sectors [@Brown:2018a].

Furthermore, we do not consider ancillary services for the distribution and transmission grids. The provision of ancillary services may be easier, i.e. less costly, for system layouts as we found them for smaller electricity systems, because not only generation but also support infrastructure is more homogeneously dispersed and thus able to provide services like frequency control or black-start everywhere. However, there is no reason to believe that this would change system cost --- especially relative system cost --- significantly [@Brown:2018].

We also do not model the distribution grid in any way. The cost of the distribution grid is likely to be higher for smaller systems where generation is dispersed more strongly, and in particular, where substantial amounts of generation is embedded within the distribution grid [@ImperialCollegeLondon:2014]. Thus, this simplification means that our model may underestimate the relative cost disadvantage of smaller-scale systems.

Lastly, we use upper-bound estimates for the potentials for solar and wind generation in our analysis, which is beneficial for systems with smaller scales. We use rather vast technical potentials which are likely not fully achievable [@Trondle:2019]. Considering social restrictions would likely lead to smaller potentials for open field photovoltaics and thus the need to install capacities on rooftops. This would increase cost especially for small-scale systems which rely more strongly on solar electricity. Furthermore, applying any restrictions on our upper-bound potential estimates would result in less regions able to cover demand locally with renewable electricity. It would thus create the need for these regions to enlarge their electricity grid and connect with neighbouring regions. Of course, this effect would be even more pronounced should the heat and transports sectors be coupled to the electricity sector.

# Supplementary Figure 1: Transmission network

![**Possible locations of transmission capacities.** All lines visualise connections between two regions that can hold transmission capacities. International connections are coloured yellow, all others are coloured blue. The amount of capacities installed on these connections is an output of the optimisation and depends on the considered case.](report/network.png){#fig:network .class}

# Supplementary Figure 2: Total Sobol' indices

![**Total Sobol' indices for combinations of all considered input uncertainties and selected model outputs of continental, national, and regional scale electricity systems.** Total Sobol' indices  determine the magnitude with which the variability of one model input explains the variability of one model output given assumptions on input variability. The x-axis comprises all twelve input parameters whose uncertainty we consider. **a,** Sobol' indices of input parameters considering total system cost and total installed capacities of the continental scale system. **b,** Sobol' indices of input parameters considering total system cost and total installed capacities of the national scale system. **c,** Sobol' indices of input parameters considering total system cost and total installed capacities of the regional scale system. **d,** Sobol' indices of input parameters considering difference in system cost between the continental and national scale systems. Framed Sobol' indices are referenced in the text.](report/total-sobol-all.png){#fig:total-sobol .class}

# Supplementary Figure 3: First-order Sobol' indices

![**First-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of continental, national, and regional scale electricity systems.** First-order Sobol' indices determine the magnitude with which the variability of one model input explains the variability of one model output given assumptions on input variability. The x-axis comprises all twelve input parameters whose uncertainty we consider. **a,** Sobol' indices of input parameters considering total system cost and total installed capacities of the continental scale system. **b,** Sobol' indices of input parameters considering total system cost and total installed capacities of the national scale system. **c,** Sobol' indices of input parameters considering total system cost and total installed capacities of the regional scale system. **d,** Sobol' indices of input parameters considering difference in system cost between the continental and national scale systems. Framed Sobol' indices are referenced in the text.](report/first-sobol-all.png){#fig:first-sobol .class}

# Supplementary Figure 4: Total minus first-order Sobol' indices

![**Total minus first-order Sobol' indices for combinations of all considered input uncertainties and selected model outputs of continental, national, and regional scale electricity systems.** Total minus first-order Sobol' indices determine the magnitude with which the variability of one model input explains the variability of one model output given assumptions on input variability. The x-axis comprises all twelve input parameters whose uncertainty we consider. **a,** Sobol' indices of input parameters considering total system cost and total installed capacities of the continental scale system. **b,** Sobol' indices of input parameters considering total system cost and total installed capacities of the national scale system. **b,** Sobol' indices of input parameters considering total system cost and total installed capacities of the regional scale system. **d,** Sobol' indices of input parameters considering difference in system cost between the continental and national scale systems. Framed Sobol' indices are referenced in the text.](report/total-minus-first-sobol-all.png){#fig:higher-sobol .class}

# Supplementary Figure 5: Distribution of time series across regions

![**Distribution of time series for all regions of the regional scale electricity system. a,b,c,d,e,f,** Distribution of time series for all regions with low flexibility cost (combined cost of bioenergy, hydrogen, and battery storage) (**a,b,c**) and all regions whose flexibility cost lies within the highest decile (**d,e,f**). **a,b,** Combined wind and solar weekly generation potential time series relative to local demand. Seasonal fluctuations are more pronounced for the case with higher cost due to higher shares of solar electricity. Regions with high solar shares are often urban regions with low or no wind potential. **c,d,** Weekly generation time series from biomass combustion relative to demand. Generation has a more pronounced seasonality and is generally larger in regions with higher flexibility cost. **e,f,** Weekly time series of hydrogen storage levels relative to installed storage capacity. In regions with higher cost, hydrogen storage is used primarily to balance seasonal fluctuations, instead of balancing fluctuations within weeks or months as it is done for lower cost regions, leading to fewer storage cycles and higher cost.](report/timeseries.png){#fig:timeseries .class}

# Supplementary Table 1: Biomass feedstocks

```table
---
caption: 'Biomass feedstocks we consider, together with the proxy we use to derive regional from national values. {#tbl:biomass-feedstocks}'
alignment: LR
include: biofuel-feedstocks.csv
include-encoding: UTF-8
markdown: True
---
```

# Supplementary Table 2: Technology cost assumptions

```table
---
caption: 'Assumptions on technology cost. ^AC transmission overnight cost is given in €/kW/1000km {#tbl:overview-cost-assumptions}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

# Supplementary Table 3: Generation and storage capacities

```table
---
caption: 'Installed capacities of photovoltaics (PV), on- and offshore wind, bioenergy, short-term (battery) and long-term (hydrogen) storage, and relative curtailment of solar, wind, and hydropower for all considered cases. Each case additionally contains fixed hydropower capacities on their current locations: 36 GW run of river, 103 GW / 97 TWh reservoirs, and 54 GW / 1.3 TWh pumped hydro storage. {#tbl:overview-scenario-results-1}'
alignment: LRRRRRRRR
include: report/overview-scenario-results-1.csv
include-encoding: utf-8
markdown: True
---
```

# Supplementary Table 4: Transmission capacities

```table
---
caption: 'Installed transmission grid capacity, gross physical electricity flow crossing country borders, and net electricity flow imported by all countries for all cases. {#tbl:overview-scenario-results-2}'
alignment: LRRR
include: report/overview-scenario-results-2.csv
include-encoding: utf-8
markdown: True
---
```

# Supplementary Table 5: Input parameter uncertainty

```table
---
caption: 'Uncertain input parameters. For all parameters we assume a uniform distribution. {#tbl:overview-uncertain-parameters}'
alignment: LLLLLL
include: report/overview-uncertain-parameters.csv
include-encoding: UTF-8
markdown: True
---
```

# Supplementary References
