# Introduction

100% renewable electricity provision through solar, wind, and hydro power is technically possible and financially viable in Europe: technologies are largely available and mature, land and water provide enough surfaces, costs are expected to be in a similar range as today [@Rasmussen:2012; @Bussar:2014; @Connolly:2016; @Gils:2017; @Plessmann:2017; @Jacobson:2017; @Brown:2018; @Zappa:2019; @Trondle:2019; @Child:2019]. However, generation of renewable sources fluctuates with the resource: weather regimes and seasons create fluctuations from days to seasons. These fluctuations can be balanced locally and temporally through generation mixes, demand flexibility, or storage. Additionally they can be balanced spatially, using a grid spanning large parts or all of Europe and institutions supporting international trade of electricity.

Relying on spatial balancing requires cooperation between regions and nations within Europe: the current electricity grid would need to be extended significantly [@Rodriguez:2014; @Schlachtberger:2017; @Child:2019] with some regions being more impacted by grid infrastructure than others (SOURCES), and institutions enabling electricity trade would need to be installed [@EuropeanCommission:2016; @Patt:2018]. However, many studies estimated lower costs for the entire electricity system when fluctuations are balanced spatially through the grid, rather than locally [@Rodriguez:2014; @Steinke:2013; @Schmid:2015, @Czisch:2005; @Schlachtberger:2017; @Child:2019]. In addition, a large grid may enable access to better resources, for example X and Y (SOURCES), decreasing system costs further [@Czisch:2005; @Schmid:2015].

Here we assess the trade-off between necessary grid infrastructure for larger electricity systems in Europe and higher costs of smaller electricity systems in detail. We quantify these impacts for electricity system sizes ranging from 502 regional systems to one continental system, and for shares of local electricity generation from 70--100%. We find that ...

## Model overview

Definition of model and scenarios here (only overview necessary to understand what we did, details in "methods section")

Model of the European power system with these scenarios:

* autarky on continental/national/regional levels, using different autarky levels (70--100% autarkic)
* continental/national/regional grid scale

Model consists of

* pv, wind, hydro, biomass on the supply side,
* pumped hydro, stylised short term and stylised long term on the storage side,

on each regional node which are connected through ac lines. Hourly resolution, one year operation, greenfield, cost optimisation.

We are assessing long term (steady state) cost of electricity supply. We are neither aiming at determining costs of a transition to such states nor are we aiming at considering developments on the global market for supply and storage technologies. Thus, we assume cost for technologies at a late stage within their learning curves, i.e. costs manifesting when technologies are deployed at large scales, as in our study.

# Results

## Continental perspective

![Total system costs of renewable electricity autarky within Europe.](report/scenario-space.png){#fig:scenario-space .class}

* costs higher for national + regional autarky by factors 1.5/3
* using grid for balancing reduces costs significantly
* allowing net import reduces costs slightly, but in cases with load shedding can change quite a bit

```table
---
caption: 'Overview over scenario results. {#tbl:overview-scenario-results}'
alignment: LRRRRRRRRRR
include: report/overview-scenario-results.csv
include-encoding: utf-8
markdown: True
---
```


![Share of total electricity flow between all countries for each interconnection for the unrestricted case with continental autarky, divided into net annual and inter-annual flows. Net annual flows cancel out imports and exports on the same interconnection and thus represent the net electricity imported and exported on any interconnection within one year. Inter-annual flows show the total amount of electricity flowing over the interconnection within one year, reduced by the net annual flows to pronounce the inter-annual effects. Both panels show the total electricity flow on all interconnections.](report/flows.png){#fig:flows .class}

## Local perspective

![Costs of renewable electricity autarky within Europe, for all Europe, all countries, and all regions. Net transfer capacities shown as blue lines connecting regions.](report/map.png){#fig:map .class}

Costs are generally higher in those countries/regions that ...

Costs are generally lower in those countries/regions that ...

DEEP DIVE 1: Country/region X with high costs. What are the drivers?

DEEP DIVE 2: Country/region Y with low costs. What are the drivers?

## Uncertainty Quantification

...

# Discussion

Large grid comes needs institutional support, and cooperation can fail.

## Short comings

1. no explicit modelling of flexibility and demand from heat sector
2. no explicit modelling of flexibility and demand from transport sector
3. no explicit modelling of flexibility from demand side management
4. ignoring status quo and transition paths (including cost dynamics _on_ the transition path)
5. ignoring distribution grid and grid services [@Brown:2018]
6. ignoring correlation of technology deployment and its costs (Europe is a large market)
7. optimistic renewable potentials (technical-potential?) in here, less optimistic ones make small scales systems impossible or more expensive

# Conclusion
Autarky can be expensive. Using the grid for balancing can bring costs down significantly. Access to better resources is not as important.

# Methods

We model the European electricity system as a set of network nodes and power flows between the nodes, with each node representing a regional administrative unit in Europe. We consider the deployment of certain renewable electricity supply and storage technologies at each node, and the deployment of transmission links between nodes, but disregard subordinate network nodes and power flows on the distribution system. We furthermore do not consider the current state of the electricity system, and thus employ a greenfield approach. We simulate a single year, 2016, with a temporal resolution of four hours. We build the model using the Calliope model framework [@Pfenninger:2018b]. The model code and all analysis steps are publicly available as a reproducible workflow (TODO publish and add DOI).

The study area comprises all countries with member organisations in the entso-e: EU-28, Norway, Switzerland, and Western Balkans countries. We exclude Iceland which is electricity autarkic, and Malta, for which no data is available. We divide the study area into 502 regional administrative units [@Trondle:2019], and add a network node at the centroid of each region. We model the transmission grid as direct net transfer capacities between network nodes, i.e. we consider net power flows on the shortest distances between nodes only. We allow transmission capacities between regional administrative units sharing a land border. Where this approach creates disconnected islands, we connect islands using the shortest connection possible to retrieve a fully connected electricity network graph as visible in Figure @fig:map. We ignore the distribution grid.

A set of supply capacities at each node, storage capacities at each node, and transmission capacities between all nodes defines an electricity system design. All system designs fulfilling Kirchhoff's law, technical constraints of components (see below), and political constraints  (see below) are possible. We use linear programming to find the design with the lowest total system costs. In the following, we describe the model of all components that can be or must be deployed at network nodes.

## Electrical load

We model electrical load as being inflexible, i.e. following a predefined profile. We determine profiles for each regional administrative unit following a method described in @Trondle:2019. First, we derive the location and annual demand of industrial facilities with highest electricity demand in Europe from emission data of the European Emission Trading Scheme [@EuropeanEnvironmentAgency:2018]. We assume industrial load to be nearly constant and thus derive flat industry profiles for each regional administrative unit.

Second, we use measured national load profiles of 2016 [@Muehlenpfordt:2019] (for Albania no data of 2016 is available and thus we use data of 2017) and subtract industrial demand to retrieve national profiles of residential and commercial load. We then assume residential and commercial load to be spatially distributed proportional to population counts and use the Global Human Settlement Population Grid with a resolution of 250 m [@JRC:2015] to allocate residential and commercial load to regional administrative units.

Finally, we sum the two time series to retrieve electricity load profiles at each network node.

## Photovoltaics

We allow generation capacities based on photovoltaics (PV) to be built on each network node. For each administrative unit, we first determine the maximum amount of capacities that can be deployed. Second, we determine the capacity factor time series which maps from installed capacity to electricity generation in each point of time. Third, we apply global cost assumptions to installed capacities.

Our model discriminates between photovoltaics deployed on open fields and those deployed on roof tops. For the maximum amount of installable capacities, we use geo spatial data with a 10 arcsecond resolution.  We allow open field PV to be built on areas of bare land [@EuropeanSpaceAgency:2010] or open vegetation [@EuropeanSpaceAgency:2010] that are not environmentally protected [@UNEP-WCMC:2018], not inhabited (i.e. < 1% of the grid cell are buildings or urban greens according to @Ferri:2017), and whose average slope [@Reuter:2007; @Danielson:2011] is 10° at maximum [@AlGarni:2018]. We assume a capacity density of 80 W/m^2^ and derive the maximum amount of installable open field PV capacity for all regional administrative units.

For the maximum installable capacities of roof-mounted PV, we consider only the areas that are inhabited (i.e. >= 1% of the grid cell are buildings or urban greens according to @Ferri:2017). Within those grid cells, we use building footprints from @Ferri:2017 as a proxy for the amount of available roof tops. Using the high-resolution sonnendach.ch data set for Switzerland [@SwissFederalOfficeofEnergy:2018], we find that within Switzerland, the ratio between building footprints from @Ferri:2017 and rooftops available for PV deployment is 0.56. We apply this ratio for all Europe to derive the maximum amount of roof space available for PV. We further discriminate the amount of roof space between flat roofs and tilted roofs based on the Swiss ratio from @SwissFederalOfficeofEnergy:2018 and assume capacity densities of 160 W/m^2^ for tilted roofs and 80 W/m^2^ for flat roofs.

We derive capacity factor time series for PV for around 2700 locations within the study area --- all those grid points of a mesh grid with 50 km edge length which fall on land. We assume a performance ratio of 90% and simulate the time series using Renewables.ninja [@Staffell:2016; @Pfenninger:2016]. Because tilt and orientation of tilted roofs have a significant impact on capacity factors, we once again discriminate tilted roofs between their tilt and orientation: we form 16 deployment situations based on east, south, west, and north facing roofs with tilts between 18° and 43°, and simulate each situation once. We then weigh and average the resulting 16 time series based on Swiss statistics from @SwissFederalOfficeofEnergy:2018 to derive one time series for roof mounted PV at each grid node. We assume capacity factor time series to be constant within the 50km^2^ grid cells and form one open field PV and one roof mounted PV time series for each administrative unit as a weighted spatial average.

Generation from photovoltaics can be curtailed, i.e. actual generation at a certain point in time can be lower than the maximum generation given by the installed capacities.

PV generation capacities come with fixed costs for installation and maintenance, whose values are given in Table @tbl:overview-cost-assumptions.

```table
---
caption: 'Economic assumptions of technologies. ^AC transmission installation costs are given in [€/kW/1000km] {#tbl:overview-cost-assumptions}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

## Wind on- and offshore

Onshore and offshore wind capacities can be deployed at each network node, and we apply a method similar to the one for photovoltaics to derive their maximum amount of installable capacities and their capacity factor time series. We apply monetary cost assumptions for wind capacities as stated in Table @tbl:overview-cost-assumptions.

We use geospatial data with 10 arcsecond resolution to derive the maximum amount of installable wind power capacities. We allow onshore wind farms to be built on areas with farmland, forests, open vegetation, and bare land [@EuropeanSpaceAgency:2010] that are not environmentally protected [@UNEP-WCMC:2018], not inhabited (i.e. < 1% of the grid cell are buildings or urban greens according to @Ferri:2017), and whose average slope [@Reuter:2007; @Danielson:2011] is 20° at maximum [@McKenna:2015]. We allow offshore wind farms to be built in offshore areas within Exclusive Economic Zones [@Claus:2018] with water depths [@Amante:2009] not below 50 m and which are not environmentally protected [@UNEP-WCMC:2018]. We assume capacity densities of 8 W/m^2^ and 15 W/m^2^ [@EuropeanEnvironmentAgency:2009] for onshore and offshore wind. Where land is available for onshore wind farms and open field PV, either technology or a mix of both technologies can be used. We allocate the installable offshore capacities to regions which share a coast with the Exclusive Economic Zone, and where there is more than one region, we allocate the capacities proportional to the length of the shared coast.

We derive capacity factor time series for on- and offshore wind on the same 50 km^2^ grid as we do for PV, resulting in around 2700 onshore locations and around 2800 offshore locations. We again use Renewables.ninja [@Staffell:2016; @Pfenninger:2016] to simulate wind generation at each grid cell, assume capacity factors to be constant within the cell, and generate a spatially weighted average to generate a capacity factor time series for each regional administrative unit. Similar to PV, wind generation can be curtailed.

## Hydro run of river and reservoirs

We assume hydro run of river and hydro reservoir potentials to be largely tapped today (ADD source) with almost no expansion potential. Thus, for hydro generation capacities we deviate from the greenfield approach and deploy the generation capacities of today. Similar to PV and wind generation, hydro generation can be curtailed. Hydro capacities come without costs.

We derive the location and installed power and storage capacities of hydro stations in Europe today from the JRC Hydro Power Database [@MatteoDeFelice:2019]. Where no storage capacity of hydro reservoirs is available, we use the median national ratio of power to storage capacity, and if that is not available, we use the median global ratio of power to storage capacity.

To derive water inflow time series for each station, we use an approach based on ERA5 runoff data [@Dee:2011] and hydrological basins [@Lehner:2013] described in [@Liu:2019]. We use Atlite [@Andresen:2015] to first determine all basins upstream of the hydro power station to be able to sum all upstream runoff while assuming a water flow speed of 1 m/s. For China, this method has shown to be able to replicate the temporal dynamics of water inflow accurately [@Liu:2019]. However, without correction factor, the method could not replicate the actual magnitude of the inflow.

To overcome this issue, and to generate power time series, we use annual national generation data from IRENA [@IRENA:2018a]. For hydro run of river plants we assume constant annual capacity factors within each country, which allows us to estimate the annual generation per plant. We use this estimation to derive electricity generation time series for each plant by scaling and capping the water inflow time series in such a way, that they sum to the annual generation without ever exceeding power capacities of the stations. For hydro reservoirs, we additionally assume they never need to spill water, i.e. their storage capacity is sufficient to use all inflowing water. We then simply scale the water inflow time series in such a way that they sum to the annual generation of the stations.

Finally, using the location data of each plant, we sum all time series and power and storage capacities per regional administrative unit to form hydro run of river (36 GW in all Europe) and hydro reservoir capacities (103 GW / 97 TWh in all Europe).

## Biofuels

We use estimations of biofuel potentials for the year 2020 and reference assumptions taken from @RuizCastello:2015, but we assume no dedicated farming for energy crops and thus consider residuals and wastes only. The data is given as national aggregates, and we use national shares of farmland [@EuropeanSpaceAgency:2010], national shares of forests [@EuropeanSpaceAgency:2010], and national shares of population [@JRC:2015] as proxies to derive regional potentials by allocating proportionally based on the proxy. Table @tbl:biofuel-feedstocks lists all feedstocks we consider together with the proxy we use.

We assume an efficiency of 40% for the combustion of all biofuels. Further, we use fuel costs of the same source [@RuizCastello:2015], which we weigh and aggregate to reach a single value for all Europe and all feedstocks, as given in Table @tbl:overview-cost-assumptions.

```table
---
caption: 'Biofuel feedstocks we consider, together with the proxy we use to derive regional values from national values. {#tbl:biofuel-feedstocks}'
alignment: LR
include: biofuel-feedstocks.csv
include-encoding: UTF-8
markdown: True
---
```

## Pumped hydro

Similar to hydro run of river and hydro reservoir capacities, we assume pumped hydro capacities in Europe to be largely tapped (ADD source) and do not allow for capacity expansion. Thus, we deploy today's pumped hydro power and storage capacities. We assume a round-trip electricity efficiency of 78% [@Schmidt:2019] and do not consider economic costs.

To determine locations, and power and storage capacities of each pumped hydro station in Europe today, we again use the JRC Hydro Power Database [@MatteoDeFelice:2019]. Where storage capacities are missing, we employ the same method as for hydro reservoirs: we assume national median ratios of power to storage capacity for all stations with missing storage capacity; and where this is not available, we assume global median ratios of power to storage capacity. Storage capacities within the JRC Hydro Power Database sum up to more than 10 TWh which is more than what other sources report (for example 1.3 TWh in [@Geth:2015]). We thus scale the storage capacities to match national data reported by @Geth:2015. Using the location data of each station, we then sum all power and storage capacities within regional administrative units to form a single pumped hydro capacity per unit.

## Short-term and long-term storage

Because it is today unknown exactly which storage technology will become dominant in a fully renewable electricity system and what their techno-economic parameters will be, we are not aiming to model specific technologies, but we are aiming at modeling technology classes: short-term and long-term storage capacities. These are available in all regional administrative units.

We are basing the two models on two technical parameters: the ratio between power and storage capacity, and the round-trip electricity efficiency. We limit short-term storage capacities to a maximum of 4 h of full power until depletion, and similarly, we limit long-term storage capacities to a minimum of 4 h of full power until depletion. Furthermore, we assume 86% (40%) of round-trip electricity efficiency for short-term (long-term) storage capacities.

Additionally, we assume that power and storage capacities can be expanded independently, only limited by above mentioned threshold and we assume distinct economic costs for each expansion as given in Table @tbl:overview-cost-assumptions. We generally assume higher costs for short-term storage capacity than for long-term storage capacity.

## Load shedding

In some regions, the local technical renewable electricity generation potential is not high enough to satisfy local electricity demand  [@Trondle:2019]. This is problematic in scenarios in which regions strive for electricity autarky. To be able to nonetheless reach a balance between demand and supply in these corner cases, we add the possibility to shed load. We assume high costs for each kWh shed so that load shedding is used as a last resort only, see Table @tbl:overview-cost-assumptions, and we allow load shedding only in scenarios with regional autarky and only in regions that could otherwise not be autarkic: Vienna, Brussels, Berlin, Oslo, and Basel.

## Financial

We are assessing long term (quasi steady state) cost of electricity supply. We are neither aiming at determining costs of a transition to such states nor are we aiming at considering disruptive developments on the global market for supply and storage technologies. Thus, we assume cost for technologies at a late stage within their learning curves, i.e. costs manifesting when technologies are deployed at large scales, as in our study. We do that by applying cost estimations for the year 2050 taken mainly from [@JRC:2014] for supply and transmission technologies, and taken from [@Schmidt:2019] for storage technologies, see Table @tbl:overview-cost-assumptions.

Costs of technologies are based on costs per installed capacities, annual maintenance costs per installed capacities, variable costs per generated electricity (Biofuel only), lifetime, and costs of capital. We use these cost components to form annuities for each technology, and variable costs for Biofuel.

We assume costs of capital to be 7.3% for all technologies and all locations based on historic average cost of capital for OECD countries [@Steffen:2019]. Extant literature has shown that costs of capital are likely technology [@Egli:2018; @Steffen:2019] and location [@Ondraczek:2015; @Steffen:2019] specific, but we consider the available data too sparse to justify a discrimination of costs.

## Political constraints

We use two types of political constraints forming the basis of our scenarios: constraints on the scale of the grid, and constraints on net imports. When we limit the scale of the grid to continental, national, or regional sizes, we force transmission capacities between continents, countries, or regions to be zero. For the continental scale this is given inherently by the scope of our study area.

When we limit net imports into continental, national, or regional areas leading to full or partial autarky, we force all or parts of the electricity demand on the continent, in the country, or in the region to be satisfied with local electricity generation from wind, sun, biomass, and water. In contrast to the constraint on the grid scale, the grid can still be used to balance fluctuations, as long as net annual imports stay below the required threshold.

## Uncertainty Analysis

...

# Bibliography
