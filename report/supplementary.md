---
title: 'Appropriate technology: The relationship between geographic scale, cost, and technology mix of fully renewable electricity systems in Europe'
subtitle: 'Supplementary Material'
author:
    - Tim Tröndle
    - Stefan Pfenninger
    - Stefano Marelli
    - Johan Lilliestam
institute: IASS Potsdam and ETH Zürich
bibliography:
    - 'literature.bib'
csl: nature.csl
link-citations: True
header-includes: |
    <script src="https://kit.fontawesome.com/40d87b9e43.js"></script>
date: draft-04-dev
---

# Code and data

The code and reproducible workflow used to perform this analysis is available online TODO ADD. Furthermore, the resulting data is available online as well TODO ADD.

# Model

## Biomass feedstocks

```table
---
caption: 'Biomass feedstocks we consider, together with the proxy we use to derive regional from national values. {#tbl:biomass-feedstocks tag="S1"}'
alignment: LR
include: biofuel-feedstocks.csv
include-encoding: UTF-8
markdown: True
---
```

## Technology costs

```table
---
caption: 'Assumptions on technology costs. ^AC transmission installation costs are given in [€/kW/1000km] {#tbl:overview-cost-assumptions tag="S2"}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

# System composition

```table
---
caption: 'Installed capacities of solar (<i class="fas fa-sun"></i>), wind (<i class="fas fa-wind"></i>), bioenergy (<i class="fas fa-leaf"></i>), short-term (<i class="fas fa-battery-three-quarters"></i> short) and long-term storage (<i class="fas fa-battery-three-quarters"></i> long) and relative curtailment of solar, wind, and hydro power (<i class="fas fa-traffic-light"></i>) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the geographic scale of the entire system (<i class="fab fa-connectdevelop"></i>). Each layout additionally contains fixed hydroelectricity capacities on their current locations: 36 GW run of river, 103 GW / 97 TWh reservoirs, and 54 GW / 1.3 TWh pumped hydro storage. {#tbl:overview-scenario-results-1 tag="S3"}'
alignment: LRRRRRRRR
include: report/overview-scenario-results-1.csv
include-encoding: utf-8
markdown: True
---
```

```table
---
caption: 'Installed transmission grid capacity (<i class="fab fa-connectdevelop"></i>), gross physical electricity flow crossing country borders (<i class="fas fa-shopping-cart"></i> gross), and net electricity flow imported by all countries (<i class="fas fa-shopping-cart"></i> net) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the geographic scale of the entire system (<i class="fab fa-connectdevelop"></i>). {#tbl:overview-scenario-results-2 tag="S4"}'
alignment: LRRR
include: report/overview-scenario-results-2.csv
include-encoding: utf-8
markdown: True
---
```

# Uncertainties

```table
---
caption: 'Uncertain input parameters. For all parameters we assume a uniform distribution. {#tbl:overview-uncertain-parameters tag="S5"}'
alignment: LLLLLL
include: report/overview-uncertain-parameters.csv
include-encoding: UTF-8
markdown: True
---
```

**TODO replace total with first order sobol indices**

![**First-order Sobol' indices for all combinations of considered input uncertainties and model outputs of continental and national scale electricity systems.** First-order Sobol' indices as result of a global sensitivity analysis determine the magnitude with which the variability of one model input directly explains the variability of one model output given assumptions on input variability. The x-axis comprises fifteen model outputs: absolute and relative system cost and installed generation, storage, and transmission capacity for the continental and the national scale system. The y-axis comprises all twelve input parameters whose uncertainty we consider: depreciation rate (d~r~), costs of solar (c~pv~), onshore wind (c~wind~), and offshore wind (c~offshore~) capacities, short-term storage power (c~sts,p~) and energy (c~sts,e~) capacities, long-term storage power (c~lts,p~) and energy (c~lts,e~) capacities, transmission (c~ntc~) and bioenergy (c~bio~) capacities, biomass fuel costs (c~fuel~) and biomass availability (a~bio~).](../data/total-sobol.png){#fig:first-sobol .class tag="S1"}

**TODO replace total with total-first order sobol indices**

![**Higher-order Sobol' indices for all combinations of considered input uncertainties and model outputs of continental and national scale electricity systems.** Higher-order Sobol' indices as result of a global sensitivity analysis determine the magnitude with which the variability of combinations of model inputs explain the variability of one model output given assumptions on input variability. The x-axis comprises fifteen model outputs: absolute and relative system cost and installed generation, storage, and transmission capacity for the continental and the national scale system. The y-axis comprises all twelve input parameters whose uncertainty we consider: depreciation rate (d~r~), costs of solar (c~pv~), onshore wind (c~wind~), and offshore wind (c~offshore~) capacities, short-term storage power (c~sts,p~) and energy (c~sts,e~) capacities, long-term storage power (c~lts,p~) and energy (c~lts,e~) capacities, transmission (c~ntc~) and bioenergy (c~bio~) capacities, biomass fuel costs (c~fuel~) and biomass availability (a~bio~).](../data/total-sobol.png){#fig:higher-sobol .class tag="S2"}

There are aspects which our analysis does not consider but which may impact our findings. Some are likely to be advantageous for smaller systems, others are likely to be advantageous for larger systems. Maybe most importantly, we do not consider flexibility from the electricity demand or from coupling the electricity sector with the heat and transport sectors. These additional flexibilities may be especially beneficial for smaller scale systems, whose flexibility options are more limited. However, large electricity systems can and will benefit as well and our sensitivity analysis shows that cost differences are driven largely by the cost of bioenergy; a technology mostly used to balance seasonal fluctuations of solar generation. So only if these flexibilities can balance seasonal fluctuations, a significant impact to system cost can be expected. This is likely not the case for demand flexibility [@Aryandoust:2017] or transportation [@Brown:2018a], but it may be possible by coupling the electricity and the heat sectors [@Brown:2018a].

Furthermore, we do not consider the distribution grid and ancillary services; the later neither for the distribution nor the transmission grid. The provision of ancillary services may be easier, i.e. less costly, for system layouts as we found them for smaller electricity systems, because not only generation but also support infrastructure is more homogeneously dispersed and thus able to provide services like frequency control or black-start everywhere. However, there is no reason to believe that this would change system cost --- especially relative system cost --- significantly [@Brown:2018]. The cost of the distribution grid, on the other side, is likely to be higher for smaller systems with generation being dispersed more strongly, and in particular dispersed within the distribution grid [@ImperialCollegeLondon:2014].

Lastly, we are using upper-bound estimations for the potentials for solar and wind generation in our analysis, which is beneficial for systems with smaller scales. We use rather vast technical potentials which are likely not fully achievable [@Trondle:2019]. Considering social restrictions would in our opinion lead especially to smaller potentials for open field photovoltaics and thus the need to install capacities on rooftops. This would increase cost especially for small scale systems which rely more strongly on solar electricity. Furthermore, applying any restrictions on our upper-bound potential estimation would lead to more regions being unable to supply themselves with renewable electricity locally and thus creating the need for these regions to increase their electricity grid and connect with neighbours. This effect would of course be even more pronounced should the heat and transports sectors be coupled to the electricity sector.

# Bibliography
