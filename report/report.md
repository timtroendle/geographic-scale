# Introduction

There is a strive for local electricity generation.

Opponents mention additional costs. Examples of quantifications are THIS AND THAT. However, they focus on grid structures and not on administrative division and they do not consider the idea of self-sufficiency.

How much does renewable electricity autarky cost in Europe?

## Model overview

Definition of model and scenarios here (only overview necessary to understand what we did, details in "methods section")

Model of the European power system with these scenarios:

* autarky on continental/national/regional levels, using different autarky levels (70--100% autarkic)
* continental/national/regional grid scale

Model consists of

* pv, wind, hydro, biomass on the supply side,
* pumped hydro, stylised short term and stylised long term on the storage side,

on each regional node which are connected through ac lines. Hourly resolution, one year operation, greenfield, cost optimisation.

# Results

## Continental perspective

![Total system costs of renewable electricity autarky within Europe.](build/output/report/scenario-space.png){#fig:scenario-space .class}

* costs higher for national + regional autarky by factors 1.5/3
* using grid for balancing reduces costs significantly
* allowing net import reduces costs slightly, but in cases with load shedding can change quite a bit

```table
---
caption: 'Overview over scenario results. {#tbl:overview-scenario-results}'
alignment: LRRRR
include: build/output/report/overview-scenario-results.csv
markdown: True
---
```

## Local perspective

![Costs of renewable electricity autarky within Europe, for all Europe, all countries, and all regions.](build/output/report/map.png){#fig:map .class}

Costs are generally higher in those countries/regions that ...

Costs are generally lower in those countries/regions that ...

DEEP DIVE 1: Country/region X with high costs. What are the drivers?

DEEP DIVE 2: Country/region Y with low costs. What are the drivers?

## Uncertainty Quantification

...

# Discussion

## Short comings

1. no explicit modelling of flexibility and demand from heat sector
2. no explicit modelling of flexibility and demand from transport sector
3. no explicit modelling of flexibility from demand side management
4. ignoring status quo and transition paths (including cost dynamics _on_ the transition path)
5. ignoring distribution grid and grid services
6. ignoring correlation of technology deployment and its costs (Europe is a large market)

# Conclusion
Autarky can be expensive. Using the grid for balancing can bring costs down significantly.

# Methods

## Model overview again

* cost optimisation 
* green field
* bla blub

## Transmission grid

* artificially created 
* net transfer capacities
* distribution grid ignored

TODO add plot of transmission grid here -- is that even possible in the methods section? Otherwise in supplementary.

## Photovoltaics

* determination of max installable capacity [@Trondle:2019]
* determination of time series [@Trondle:2019]
* costs, see Table @tbl:overview-cost-assumptions

```table
---
caption: 'Overview over cost assumptions. ^AC transmission installation costs are given in [€/kW/km] {#tbl:overview-cost-assumptions}'
alignment: LRRRRR
include: build/output/report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

## Wind on- and offshore

* determination of max installable capacity, [@Trondle:2019]
* determination of time series [@Trondle:2019]
* costs, see Table @tbl:overview-cost-assumptions

## Hydro run of river and reservoirs

* determination of installed capacity, [@MatteoDeFelice:2019]
* determination of time series [@Liu:2019],
* costs ignored

## Biomass

* determination of max installable capacity
* costs, see Table @tbl:overview-cost-assumptions TODO add

## Pumped hydro

* determination of installed capacities, [@MatteoDeFelice:2019]
* model description

## Stylised short-term storage (battery)

* model description
* upper bound, TODO TBD
* costs, see Table @tbl:overview-cost-assumptions

## Stylised long-term storage (hydrogen)

* model description
* upper bound, TODO TBD
* costs, see Table @tbl:overview-cost-assumptions

## Load shedding

* why necessary
* how does it work
* costs, see Table @tbl:overview-cost-assumptions

## Financial

* discussion of interest rate, WACC
* lifetime parameters, see Table @tbl:overview-cost-assumptions
* economic vs technical lifetime

## Uncertainty Analysis

...

# Bibliography


 

