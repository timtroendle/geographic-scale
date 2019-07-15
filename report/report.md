# Introduction

100% renewable electricity provision through solar, wind, and hydro power is technically possible and financially viable in Europe: technologies are largely available and mature, land and water provide enough surfaces, costs are expected to be in a similar range as today [@Rasmussen:2012; @Bussar:2014; @Connolly:2016; @Gils:2017; @Plessmann:2017; @Jacobson:2017; @Brown:2018; @Zappa:2019; @Trondle:2019; @Child:2019]. However, generation of renewable sources fluctuates with the resource: weather regimes and seasons create fluctuations from days to seasons. These fluctuations can be balanced locally and temporally through generation mixes, demand flexibility, or storage. Additionally they can be balanced spatially, using a grid spanning large parts or all of Europe and institutions supporting international trade of electricity.

Relying on spatial balancing requires cooperation between regions and nations within Europe: the current electricity grid would need to be extended significantly [@Rodriguez:2014; @Schlachtberger:2017; @Child:2019] with some regions being more impacted by grid infrastructure than others (SOURCES), and institutions enabling electricity trade would need to be installed [@EuropeanCommission:2016; @Patt:2018]. However, many studies estimated lower costs for the entire electricity system when fluctuations are balanced spatially through the grid, rather than locally [@Rodriguez:2014; @Steinke:2013; @Schmid:2015, @Czisch:2005; @Schlachtberger:2017; @Child:2019]. In addition, a large grid may enable access to better resources, for example X and Y (SOURCES), decreasing system costs further [@Czisch:2005; @Schmid:2015].

# Results

## An adequate mix of balancing options

Here we assess the trade-off between necessary grid infrastructure for larger electricity systems in Europe and higher costs of smaller electricity systems in detail. We quantify these impacts for a fully renewable electricity system based on solar, wind, and hydro power, and biofuels from residuals on the supply side, and pumped hydro power, short-term, and long-term electricity storage options. We consider different scales of autarky to form electricity systems on different geographical scales: First we consider a fully connected, continental system layout without any imports from outside as a large scale system with continental autarky. Second, we consider a European system in which all countries are strictly autarkic, i.e. self-sufficient, not importing or even exchanging electricity between countries as medium scale systems with national autarky. Lastly, we do the same for all regions in Europe to form small scale systems with regional autarky.

We find that total system costs can differ by more than factor 2 among these three options. Taking the fully connected continental system as a base, strict national autarky within Europe has 40% higher total system costs and the additional costs rise to 150% for a layout with strict autarky on the regional level. We additionally consider weaker autarky forms and allow for grid connections to balance renewable fluctuations within a year with and without net imports into autarkic units. We find that using the grid to balance renewable fluctuations has large impacts on total system costs, while net imports are less effective, see Figure @fig:scenario-space for the relative costs of all twelve system layouts we consider.

![Total system costs of twelve electricity system layouts in Europe, relative to a fully connected layout with continental autarky. Each layout is defined by three aspects related to autarky: first the autarky layer which defines whether the entire continent is autarkic as a whole, or whether countries or regions within Europe are autarkic each on their own. Second, the maximum level of net electricity imports which determines the strictness of autarky: strict when none of local electricity demand can be satisfied by net imports, or weak when up to 15% or 30% of demand can be satisfied by imports. Lastly, grid size defines up to which spatial scale autarkic units are connected through the electricity grid.](report/scenario-space.png){#fig:scenario-space .class}

Net imports can lower costs by generating electricity in regions with very good resources. In fact, this option is used extensively in the unrestricted case of continental autarky. Here, 45% of electricity demand is satisfied by electricity imported from abroad, with some countries relying almost entirely on imports. However, these imports are only sightly less costly than domestic electricity: when we restrict imports to 30% on the national autarky level -- reducing imports by a third -- costs rise by less than 3%. If we disallow net imports completely, costs rise to 6% above the unrestricted case. Thus, of the 40% total system cost difference between the cases of strict national autarky and unrestricted continental autarky, only 6% are related to net imports, although almost every second unit of electricity is generated abroad. For regional autarky, net imports show a more pronounced, yet still less important, effect. (WHICH is related to load shedding and shouldn't exist in my opinion -- but what can _I_ do about it?).

The differences in costs are thus mainly driven by the available options to balance renewable fluctuations and their respective costs. In the large scale system layout, all balancing options are available everywhere: renewable over-capacities with curtailed peaks, flexible biofuel combustion, short- and long-term electricity storage, pumped hydro storage, hydro reservoirs, and spatial balancing of renewable generation and electricity demand. Consequently, all of these options are used to form a least cost mix, typically with lowest contribution per option compared to system layouts with smaller scale grid, see Table @tbl:overview-scenario-results-1. Apart from hydro electricity which we fix to current levels (see Methods), there is almost 1 TW of generation capacity from wind and solar --- roughly two times European peak demand --- in the large scale system which is curtailed by 5%. Other layouts on smaller scales demand up to 1.5 TW of solar and wind capacities of which up to 10% are curtailed. Relatively expensive biofuel capacity completes the generation mix on the continental scale with less than 15 GW, while the smallest scale layout comprises more than 10 times this amount. To balance daily and seasonal fluctuations, the electricity system of the large scale layout includes 0.2 TWh of short-term, and 3.6 TWh of long-term storage capacity. For smaller scale systems with limited balancing options, this rises to almost 1 TWh and 16 TWh respectively. The availability of all balancing options on layouts with continental scale grid allows to choose most adequate solutions for all possible situations and thus allows to find a cost optimal mix, in contrast to layouts with smaller grids, where the range of options is limited.

```table
---
caption: 'Installed capacities of solar (<i class="fas fa-sun"></i>), wind (<i class="fas fa-wind"></i>), biofuel (<i class="fas fa-leaf"></i>), short-term (<i class="fas fa-battery-three-quarters"></i> short) and long-term storage (<i class="fas fa-battery-three-quarters"></i> long) and relative curtailment of solar, wind, and hydro power (<i class="fas fa-traffic-light"></i>) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the size of the electricity grid surrounding autarkic units (<i class="fab fa-connectdevelop"></i>). Each layout additionally contains fixed hydroelectricity capacities: 36 GW run of river, 103 GW / 97 TWh reservoirs, and 54 GW / 1.3 TWh pumped hydro storage. {#tbl:overview-scenario-results-1}'
alignment: LRRRRRRRR
include: report/overview-scenario-results-1.csv
include-encoding: utf-8
markdown: True
---
```

Contrary to all other balancing options, we of course find the highest reliance on the transmission grid in all system layouts with continental grid size, see Table @tbl:overview-scenario-results-2. This is the case for transmission capacity, as well as for the usage of that capacity. The largest scale system layout contains high voltage transmission capacity of nearly 400 TW km, which is roughly 4 times current capacity [@ENTSO-E:2019] (CURRENT NUMBER CORRECT?). Using this capacity, about 3600 TWh electricity are crossing country borders, which, again, is roughly 4 times current international electricity flows [@ENTSO-E:2019]. Of those 3600 TWh, 1500 TWh remain as net imports within countries, which is more than 7 times the amount of today [@ENTSO-E:2019]. All these numbers are smaller for national and regional autarky with continental grid size, but their transmission capacity and international electricity flows are still 2-3 times today's values. The cost benefits promised by layouts with continental grid size thus require a significant  increase of transmission grid infrastructure and, even more so, a strong increase of international trade.

```table
---
caption: 'Installed transmission grid capacity (<i class="fab fa-connectdevelop"></i>), gross physical electricity flow crossing country borders (<i class="fas fa-shopping-cart"></i> gross), net electricity flow imported by all countries (<i class="fas fa-shopping-cart"></i> net), and the relative amount of electricity demand that is shed, i.e. remains unserved (<i class="fas fa-truck-loading"></i>) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the size of the electricity grid surrounding autarkic units (<i class="fab fa-connectdevelop"></i>). {#tbl:overview-scenario-results-2}'
alignment: LRRRR
include: report/overview-scenario-results-2.csv
include-encoding: utf-8
markdown: True
---
```

Generally, there is a trend in the generation capacity composition shifting from wind to solar capacities when moving along the axis of grid size towards smaller electricity grids, see Table @tbl:overview-scenario-results-1. In the case of continental autarky, the share of solar to wind capacities is 1 to 3 despite the higher installation costs and the higher levelised costs of wind capacities. The reason is that wind generation in all Europe is less correlated than that of solar and thus wind generation can be balanced to large parts using a sufficient transmission grid only; an effect that has been shown in extant literature [@Grams:2017]. In electricity system layouts with disconnected, or insufficient transmission grids this characteristic of wind generation cannot be exploited and thus more balanced, but also more expensive generation capacities compositions are reached --- with solar capacities being even slightly higher than wind capacities in the case of strict regional autarky.

Autarkic units can have very different levelised cost of electricity, depending on the renewable resources they have available, and depending on magnitude and shape of their electricity demand, see Figure @fig:map for the relative cost compared to the baseline of continental autarky. On the regional level, there are five regions with levelised costs 10 times above the baseline. These regions are urban areas with high population density and renewable potentials insufficient for autarky: Brussels, the Swiss canton Basel-Stadt, Berlin, Vienna, and Oslo. Their very high levelised costs stem from the fact that they have to shed load, an option that we allow only in these five regions and for which we consider high costs (see Methods). Four other urban regions --- Geneva, Prague, Budapest, and Bucharest --- have very high levelised costs of 4-6 times the baseline despite having sufficient potentials. Their potential is however almost exclusive limited to roof mounted photovoltaics, whose generation has strong seasonalities. To balance the fluctuations locally within the region, the cities need to build up large stocks of long-term storage capacities which drive up their electricity costs.

Roughly 50% of the remaining regions has levelised costs in the range of 2-3 times that of the baseline. The regions in this group fall in one of two categories: on one hand urban regions, whose costs are driven by similar effects as described above, yet less pronounced. On the other hand regions with often large potentials and vast amounts of available land, but whose electricity profile does not fit well to the solar and wind generation profiles and where consequently more investments into balancing options need to be made. Their costs are often dominated by rather expensive biofuel combustion capacities. Several baltic and Scandinavian regions fall into this category. Those regions with costs below the baseline all have significant amounts of hydro electricity generation or pumped hydro storage for which we do not consider costs (see Methods).

![Levelised costs of renewable electricity for strict continental autarky, strict national autarky, and strict regional autarky in relation to the baseline of continental autarky. Strict autarky does not allow for any form of exchange of electricity between autarkic units. Each panel shows the autarky level and the relative total system costs. Blue lines visualise net transfer capacities connecting regions in the cases of continental and national autarky.](report/map.png){#fig:map .class}

And what about Figure @fig:flows? Why is it here, and what does it tell us?

No one knows.

![Share of total electricity flow between all countries for each interconnection for the unrestricted case with continental autarky, divided into net annual and inter-annual flows. Net annual flows cancel out imports and exports on the same interconnection and thus represent the net electricity imported and exported on any interconnection within one year. Inter-annual flows show the total amount of electricity flowing over the interconnection within one year, reduced by the net annual flows to pronounce the inter-annual effects. Both panels show the total electricity flow on all interconnections.](report/flows.png){#fig:flows .class}

## Uncertainty Quantification

Because small scale systems are restricted versions of large scale systems, they can only be equally or more expensive. Above, we show the mechanics driving cost differences: mainly, it is the ability of large grids to balance fluctuating generation of renewables. However, this mechanism, but also the exact cost difference of layouts of electricity systems may be sensitive to a range of parameters, for which exact values are not known. Thus, we assess the uncertainty of our results based on input uncertainty of two important groups of inputs: meteorological conditions, and assumptions on costs of technologies. Because this analysis is computationally heavy, we perform it with a reduced version of our model by running it on national instead of regional resolution on which results deviate up to 7% (MAY BE LOWER, need to rerun regional resolution) compared to the full resolution.

Any supply system based on 100% renewable electricity is largely dependent on meteorological conditions: when, where, and how strong the sun shines and the wind blows impacts the suitability of electricity system layouts. On top of that, the conditions vary between years and can thus impact annual performance (CITE Stefan). To understand whether our main result --- the cost difference between large and small scale system layouts --- depends on meteorological conditions, we use meteorological data of ten instead of only one year (see Methods). We find that the additional cost of national autarky compared to continental autarky is unaffected: the difference to the case with only one year is negligible (< 1‰). This is not to say that the electricity system layout is not sensitive to meteorological conditions. In fact, we find that total system costs are generally slightly higher, and more wind and biofuel capacities are deployed in exchange of solar capacities. However, large scale and small scale system layouts are impacted similarly and so the difference between both is unaffected by the longer time duration of considered meteorological conditions.

We furthermore assess the uncertainty of our result stemming from uncertainty of technology costs. While we do know current costs and we do know that costs are likely to fall with deployment due to learning effects, we do not know exact future costs with certainty. This is mainly because it is unknown exactly how much technologies are deployed, and exactly how much costs will fall with deployment. Because in our analysis we are assuming technologies are deployed at large scale, the first aspect does not lead to much uncertainty. Uncertainty of technology costs in our analysis thus stems largely from the uncertain relationship of deployment and cost decreases. Being a cost optimisation, the absolute total system costs of any assessed electricity system layout is sensitive to costs of technologies (SOURCE). Here, we assess the sensitivity of the cost difference between large and small scale system layouts.

We consider ten technology cost parameters and the weighted cost of capital as uncertain and describe their uncertainty by an uniform distribution over ranges taken from literature (see Methods). We perform a global sensitivity analysis of the additional costs of national autarky compared to continental autarky in this eleven dimensional space. We find that national autarky is 20%--60% (FAKE RESULT) more expensive and that this range is largely driven by parameters X, Y, and Z. Higher values of X lead to higher relative costs of national autarky, while higher values of Y lead to lower relative costs of national autarky. High values of Z lead to higher relative costs of national autarky, but only if X is high as well.

# Discussion

## Findings

1. Regional or national scale grids come with significantly higher costs than one that spans the continent.

2. (net) Autarky must not be much more expensive, as long as a sufficient transmission grid exists for balancing.

3. The scale of the grid defines the cost optimal composition of supply and support infrastructure: large grids come with wind, some storage, and some flexible biofuel. Small grids come with less wind, more solar, more storage and more flexible biofuel. (I would like to say: benefit of grid can be tapped only if the other infrastructure fits, but we do not know how much worse 'non-fitting' mixes are, see conclusions.)

## Short comings

1. no explicit modelling of flexibility and demand from heat sector
2. no explicit modelling of flexibility and demand from transport sector
3. no explicit modelling of flexibility from demand side management
4. ignoring status quo and transition paths
5. ignoring distribution grid and grid services [@Brown:2018]
6. optimistic renewable potentials (technical-potential?) in here, less optimistic ones make small scales systems impossible or more expensive
7. uncertainty analysis on national resolution
8. We did not investigate the necessary transmission grid capacities, see paragraph below.

For a trade off between costs and transmission grid, not only the scale of the grid must be taken into account, but also the "strength" / capacity (a fourth dimension in our scenario space, we do not do that). While we did not investigate that we speculate that costs (and necessary infrastructure of course) lie between those layouts we looked at. Others have shown favourable non linear relationships between transmission capacity and benefits: large parts of the benefits can be reached with fractions of the capacity [@Rodriguez:2014; @Schlachtberger:2017]. (Should we investigate this? Or is this [FRIN](https://en.wikipedia.org/wiki/Further_research_is_needed)?)

## (Policy) Conclusion

Large grids offer cost reductions but are much larger than today, need infrastructure  expansion and institutional support; with risks of cooperation failing.

A sensible compromise between costs and grid infrastructure seem to be a system layout without large net imports but with intra-annual exchange.

A long term strategy (and cooperation) for Europe can reduce costs and help achieve adequate compositions of supply and support infrastructure: large wind infrastructure unfolds its economic potential with a strong grid; while large solar infrastructure would not benefit from a strong grid as much.
(BUT we do not know how much worse a large grid with dominant solar, or a small grid with dominant wind is. Should we investigate this? Or is this [FRIN](https://en.wikipedia.org/wiki/Further_research_is_needed)?)

# Methods

We model the European electricity system as a set of network nodes and power flows between the nodes, with each node representing a regional administrative unit in Europe. We consider the deployment of certain renewable electricity supply and storage technologies at each node, and the deployment of transmission links between nodes, but disregard subordinate network nodes and power flows on the distribution system. We furthermore do not consider the current state of the electricity system, and thus employ a greenfield approach.

Using the Calliope model framework [@Pfenninger:2018b], we build a linear programming model that simultaneously optimises electricity system design and operation for a single year, 2016, with a temporal resolution of four hours. We choose a single year to reduce problem size, but we test our choice in a sensitivity analysis, see below. The objective function of the model is to find the design with the lowest total system costs. An electricity system design is defined by a set of supply capacities at each node, storage capacities at each node, and transmission capacities between all nodes. All system designs fulfilling Kirchhoff’s law, the technical constraints of all possible technology components, and political constraints (see below) are possible.

The model code and all analysis steps are publicly available as a reproducible workflow (TODO publish and add DOI). While the following text describes the modelling choices and assumptions made, the entirety of the model and all embedded assumptions are available for inspection through the repository linked to above.

## Geographic scope and transmission grid

The study area comprises all countries within member organisations ENTSO-E: EU-28, Norway, Switzerland, and Western Balkans countries. We exclude Iceland which is electricity autarkic, and Malta, for which insufficient data are available. We divide the study area into 502 regional administrative units [@Trondle:2019], and add a network node at the centroid of each region. We model the transmission grid as direct net transfer capacities between network nodes, i.e. we consider net power flows on the shortest distances between nodes only, and assume the distribution network within each node is able to handle distribution load. We allow transmission capacities between regional administrative units sharing a land border. Where this approach creates disconnected islands, we connect islands using the shortest connection possible to retrieve a fully connected electricity network graph as visible in Figure @fig:map.

## Electrical load

We determine electrical load profiles for each regional administrative unit following a method described in @Trondle:2019. First, we derive the location and annual demand of industrial facilities with highest electricity demand in Europe from emission data of the European Emission Trading Scheme [@EuropeanEnvironmentAgency:2018]. We assume industrial load to be nearly constant and thus derive flat industry profiles for each regional administrative unit.

Second, we use measured national load profiles of 2016 [@Muehlenpfordt:2019] (for Albania no data of 2016 is available and thus we use data of 2017) and subtract industrial demand to retrieve national profiles of residential and commercial load. We then assume residential and commercial load to be spatially distributed proportional to population counts. Using the Global Human Settlement Population Grid with a resolution of 250 m [@JRC:2015] we allocate residential and commercial load to regional administrative units.

Finally, we sum the two industrial and residential time series in each administrative unit to retrieve electricity load profiles at each network node.

## Photovoltaics

Photovoltaics (PV) can be built at each network node. For each administrative unit, we first determine the maximum amount of capacities that can be deployed. Then, we determine the capacity factor time series which maps from installed capacity to electricity generation in each point of time.

Our model differentiates between photovoltaics deployed on open fields and those deployed on roof tops. For the maximum amount of installable capacities, we use geo spatial data with a 10 arcsecond resolution.  We allow open field PV to be built on areas of bare land [@EuropeanSpaceAgency:2010] or open vegetation [@EuropeanSpaceAgency:2010] that are not environmentally protected [@UNEP-WCMC:2018], not inhabited (i.e. < 1% of the grid cell are buildings or urban greens according to @Ferri:2017), and whose average slope [@Reuter:2007; @Danielson:2011] is 10° at maximum [@AlGarni:2018]. We assume a capacity density of 80 W/m^2^ and derive the maximum amount of installable open field PV capacity for all regional administrative units.

For the maximum installable capacities of roof mounted PV, we consider inhabited areas only (i.e. >= 1% of the grid cell are buildings or urban greens according to @Ferri:2017). Within those grid cells, we use building footprints from @Ferri:2017 as a proxy for the amount of available roof tops. Using the high-resolution Sonnendach.ch data set for Switzerland [@SwissFederalOfficeofEnergy:2018], we find that within Switzerland, the ratio between building footprints from @Ferri:2017 and rooftops available for PV deployment is 0.56. Due to the lack of comparable data for other countries, we apply this ratio for all Europe to derive the maximum amount of roof space available for PV. We further differentiate between roof space on flat roofs and on tilted roofs based on the ratio from @SwissFederalOfficeofEnergy:2018 and assume capacity densities of 160 W/m^2^ for tilted roofs and 80 W/m^2^ for flat roofs.

We derive capacity factor time series for roof mounted and open-field PV on a mesh grid with 50 km edge length, resulting in around 2700 locations within the study area. We assume a performance ratio of 90% and simulate the time series using Renewables.ninja [@Staffell:2016; @Pfenninger:2016]. For roof mounted PV, because tilt and orientation of tilted roofs have a significant impact on capacity factors, we model 16 different deployment situations covering roofs facing east, south, west, and north, with tilts between 18° and 43°. We then weight and average the resulting 16 time series based on the distribution of roofs from @SwissFederalOfficeofEnergy:2018 to derive a single time series for roof mounted PV at each grid node. For open field PV we optimise the tilt based on location [@Jacobson:2018]. Eventually, we form one open field PV and one roof mounted PV time series for each administrative unit as a weighted spatial average from the ~2700 grid nodes.

Generation from photovoltaics can be curtailed, i.e. actual generation at a certain point in time can be lower than the maximum generation given by the installed capacities.

## Wind on- and offshore

Onshore and offshore wind capacities can be deployed at each network node, and we apply a method similar to the one for photovoltaics to derive their maximum amount of installable capacities and their capacity factor time series.

We use geospatial data with 10 arcsecond resolution to derive the maximum amount of installable wind power capacities. We allow onshore wind farms to be built on areas with farmland, forests, open vegetation, and bare land [@EuropeanSpaceAgency:2010] that are not environmentally protected [@UNEP-WCMC:2018], not inhabited (i.e. < 1% of the grid cell are buildings or urban greens according to @Ferri:2017), and whose average slope [@Reuter:2007; @Danielson:2011] is 20° at maximum [@McKenna:2015]. We allow offshore wind farms to be built in offshore areas within Exclusive Economic Zones [@Claus:2018] with water depths [@Amante:2009] not below 50 m and which are not environmentally protected [@UNEP-WCMC:2018]. We assume capacity densities of 8 W/m^2^ and 15 W/m^2^ [@EuropeanEnvironmentAgency:2009] for onshore and offshore wind. Where land is available for onshore wind farms and open field PV, either technology or a mix of both technologies can be used. We allocate the installable offshore capacities to regions which share a coast with the Exclusive Economic Zone, and where there is more than one region, we allocate the capacities proportional to the length of the shared coast.

We derive capacity factor time series for on- and offshore wind on the same 50 km^2^ grid as we do for PV, resulting in around 2700 onshore locations and around 2800 offshore locations. We again use Renewables.ninja [@Staffell:2016; @Pfenninger:2016] to simulate wind generation at each grid cell, assume capacity factors to be constant within the cell, and generate a spatially weighted average to generate a capacity factor time series for each regional administrative unit. Similar to PV, wind generation can be curtailed.

## Hydro run of river and reservoirs

We assume hydro run of river and hydro reservoir potentials to be largely tapped today (ADD source) with almost no expansion potential. Thus, for hydro generation capacities we deviate from the greenfield approach and fix today's capacities. Similar to PV and wind generation, hydro generation can be curtailed.

We derive the location and installed power and storage capacities of hydro stations in Europe today from the JRC Hydro Power Database [@MatteoDeFelice:2019]. Where no storage capacity of hydro reservoirs is available, we use the median national ratio of power to storage capacity, and if that is not available, we use the median Europe-wide ratio of power to storage capacity.

To create power generation time series for each station, we use a two-stage approach. First, we derive water inflow time series for each station using an approach based on ERA5 runoff data [@Dee:2011] and hydrological basins [@Lehner:2013] described and validated for China in [@Liu:2019]. We use Atlite [@Andresen:2015] to first determine all basins upstream of the hydro power station to be able to sum all upstream runoff while assuming a water flow speed of 1 m/s. 

Second, we are applying bias correction factors based on annual generation necessary for this method to represent the actual magnitude of the inflow and thus accurately model power generation. As we do not have data per station, we use national generation data from IRENA [@IRENA:2018a]. For hydro run of river plants we assume constant annual capacity factors within each country, which allows us to estimate the annual generation per plant. We use this estimation to derive electricity generation time series for each plant by scaling and capping the water inflow time series in such a way, that they sum to the annual generation without ever exceeding power capacities of the stations. For hydro reservoirs, we additionally assume they never need to spill water, i.e. their storage capacity is sufficient to use all inflowing water. We then scale the water inflow time series in such a way that they sum to the annual generation of the stations.

Using location data of each plant, we sum up time series as well as power and storage capacities per regional administrative unit. Our total resulting capacities are 36 GW for run of river and 103 GW / 97 TWh for reservoirs.

## Biofuels

We use estimations of biofuel potentials for the year 2020 and reference assumptions taken from @RuizCastello:2015, but we assume no dedicated farming for energy crops and thus consider residuals and wastes only. The data is given as national aggregates, and we use national shares of farmland [@EuropeanSpaceAgency:2010], national shares of forests [@EuropeanSpaceAgency:2010], and national shares of population [@JRC:2015] as proxies to derive proportionally allocated potentials per regional administrative unit. Table @tbl:biofuel-feedstocks lists all feedstocks we consider together with the allocation proxy we use.

We assume an efficiency of 45% for the combustion of all biofuels.

```table
---
caption: 'Biofuel feedstocks we consider, together with the proxy we use to derive regional from national values. {#tbl:biofuel-feedstocks}'
alignment: LR
include: biofuel-feedstocks.csv
include-encoding: UTF-8
markdown: True
---
```

## Pumped hydro

Similar to hydro run of river and hydro reservoir capacities, we assume pumped hydro capacities in Europe to be largely tapped (ADD source) and do not allow for capacity expansion. Thus, we deploy today's pumped hydro power and storage capacities. We assume a round-trip electricity efficiency of 78% [@Schmidt:2019].

To determine location, power and storage capacity of each pumped hydro station in Europe today, we also use the JRC Hydro Power Database [@MatteoDeFelice:2019]. Where storage capacities are missing, we employ the same method as for hydro reservoirs: we assume national median ratios of power to storage capacity for all stations with missing storage capacity; and where this is not available, we assume Europe-wide median ratios of power to storage capacity. Storage capacities within the JRC Hydro Power Database sum up to more than 10 TWh which is more than what other sources report (for example 1.3 TWh in [@Geth:2015]). We thus scale the storage capacities to match national data reported by @Geth:2015. Using location data of each station, we then sum all power and storage capacities within regional administrative units to form a single pumped hydro capacity per unit.

## Short-term and long-term storage

Because it is not yet known which storage technology will become dominant in a fully renewable electricity system and what its techno-economic parameters will be, we do not model specific technologies, but two distinct technology classes: short-term and long-term storage. We assume that these can be deployed in all regional administrative units.

We model both classes models with two technical parameters: the ratio between power and storage capacity, and the round-trip efficiency. Short-term storage is permitted to a maximum capacity of 4 h of full power, while long-term storage has a minimum of 4 h capacity at full power. We assume 86% of round-trip efficiency for short-term and 40% for long-term storage.

Additionally, we assume that power and storage capacities can be expanded independently, only limited by the above mentioned minimum/maximum storage capacities.

## Load shedding

In some regions, local technical renewable electricity generation potential is not high enough to satisfy local electricity demand  [@Trondle:2019]. This is problematic in scenarios in which regions strive for electricity autarky.

For the model to remain feasible in these corner cases, we permit load shedding with high cost per kWh shed (see Table @tbl:overview-cost-assumptions), and we allow load shedding only in scenarios with regional autarky and only in regions that could otherwise not be autarkic: Vienna, Brussels, Berlin, Oslo, and Basel.

## Technology costs

We assess long term (quasi steady state) cost of electricity supply. We aim neither to determine costs of a transition to such states nor to consider disruptive developments on the global market for supply and storage technologies. Thus, we assume expected learning-rate based costs once renewable supply and electricity storage are deployed at the large scales consistent with our study. Cost estimations for the year 2050 are primarily from [@JRC:2014] for supply and transmission technologies, from [@Schmidt:2019] for storage technologies, and from [@RuizCastello:2015] for fuel costs of biofuels.

Technology costs are modelled as installed capacities costs, annual maintenance costs based on installed capacity, and, for biofuels only, variable costs per unit of generated electricity. We do not model variable operation and maintenance costs for technologies other than biofuels because it is the only fuel-burning technology we model, so annual per-capacity maintenance costs are sufficient for the other technologies. Technology lifetime and costs of capital are used to derive annuities for each technology. We do not consider costs for hydroelectricity. See Table @tbl:overview-cost-assumptions for an overview of all cost assumptions.

We assume costs of capital to be 7.3% for all technologies and all locations based on historic average cost of capital for OECD countries [@Steffen:2019]. Some recent literature suggests costs of capital are likely specific to technology [@Egli:2018; @Steffen:2019] and location [@Ondraczek:2015; @Steffen:2019], but we consider the data available so far too sparse to provide a solid basis on which to model this.

```table
---
caption: 'Assumptions on technology costs. Hydroelectricity has no cost. ^AC transmission installation costs are given in [€/kW/1000km] {#tbl:overview-cost-assumptions}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

## Political constraints

We use two types of political constraints as the basis for our scenarios: constraints on the scale of the electricity grid, and constraints on net electricity imports. When we limit the scale of the grid to continental, national, or regional sizes, we force transmission capacities between continents, countries, or regions to be zero. For the continental scale this is given inherently by the scope of our study area.

When we limit net imports into continental, national, or regional areas leading to full or partial autarky, we force all or parts of the electricity demand on the continent, in the country, or in the region to be satisfied with local electricity generation from wind, sun, biomass, and water. In contrast to the constraint on the grid scale, the grid can still be used to balance fluctuations, as long as net annual imports stay below the required threshold.

## Sensitivity to meteorological conditions

We assess the sensitivity of our main result to meteorological conditions impacting generation from wind and solar power. While we use only a single year of meteorological conditions, 2016, in the scenario analysis, we use ten years, 2007--2016, in this sensitivity analysis. In this way, we are not assessing the variability between meteorological years, but we are assessing how much the result changes when considering a long duration of meteorological conditions. In this way, we find electricity system layouts that are optimal for the entire range of years, not only for single years.

We assess the sensitivity of relative total system costs of national autarky with continental autarky as a baseline. Because computational requirements to solve a model with regional spatial resolution and temporal resolution of 4h for ten years are very high, we perform this sensitivity analysis using a model with national spatial resolution while keeping temporal resolution the same. Comparing the results between national and regional resolution for the case with only one year of meteorological conditions, we find a difference of 7% for the relative costs of national autarky.

## Sensitivity to technology costs

Ignore correlation of uncertain inputs because "end of learning" type of costs are assumed, so the link between cost parameters based on political decisions can be ignored. There is still a link through developments on markets of raw material for example, but we ignore that.

* interest rate [i]: 1.9—13.5% max and min for PV and wind between 2009 and 2017 (Steffen, 2019)
* installation costs of PV [c_pv]: 280—580 €/kW (JRC, 2014) (Table 7)
* installation costs wind onshore [c_wind]: 800-1700 €/kW (Fraunhofer ISE, 2015)(JRC, 2014) (Table 4)
* installation costs wind offshore [c_offshore]: 1790–3270 €/kW (JRC, 2014) (Table 5)
* installation costs short term storage power [c_sts_p]: 110-160 €/kW (JRC, 2014) (Table 61)
* installation costs short term storage energy [c_sts_e]: 128.37–336.72 €/kWh (Schmidt et al., 2017) (JRC, 2014) (Table 61)  (could be correlated with power)
* installation costs long term storage power [c_lts_p]: 238—363 €/kW (Schmidt et al., 2017)
* installation costs long term storage energy [c_lts_e]: €/kWh::(REMOVE? NO SOURCES)::
* installation costs net transfer capacities [c_ntc]: 700—1080 € / kW / km (JRC, 2014)
* installation costs biofuel [c_bio]: 1380—3450 €/kW (JRC, 2014) (Table 48)
* fuel costs biofuel [c_biofuel]: 25.02–30.58 €/MWh (Ruiz Castello et al., 2015)

# Bibliography
