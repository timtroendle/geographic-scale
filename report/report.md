# Introduction

A large body of work has shown that 100% or near-100% renewable electricity provision through solar, wind, and hydro power is technically possible and economically viable in Europe. The necessary technologies are available and increasingly mature (refs), there is enough land for the infrastructure required (refs), costs can be expected to be in a similar range as today (refs) [@Rasmussen:2012; @Bussar:2014; @Connolly:2016; @Gils:2017; @Plessmann:2017; @Jacobson:2017; @Brown:2018; @Zappa:2019; @Trondle:2019; @Child:2019], and a range of options are available to balance the variability of renewable generation on time scales ranging from hours to years (refs). However, the different strategies to deal with this variability each come with drawbacks and with dramatically different political and social implications across areas as diverse as system ownership, governance and political control, and landscape impacts ([@Lilliestam:2016], MORE refs).  On one end of the spectrum are small-scale, decentralised solutions, often coupled with visions of local or regional autarky. HERE’S HOW THEY WORK (Cite: various decentralists). On the other end are very large-scale solutions, relying on stochastic smoothing of supply through very large and effective grids, spanning an entire continent (CITATIONS) or more (CITATIONS: EU-NA; GOBITEC; IEA-GLOBAL GRID).

Here we investigate how moving along this spectrum changes the technical feasibility and the cost of getting to a 100% renewable-powered Europe. We do this using a single, internally consistent framework to model different degrees of autarky across Europe, while taking into consideration the spatial-temporal variability of renewable generation.

There is much reason to expect that, from a techno-economic perspective, a densely meshed continental grid is cheaper than smaller, autarkic systems. Nevertheless, the idea of citizen power, small-is-beautiful, and local, democratic control over the essential commodity electricity and its infrastructure has great appeal to many, citizens and decision-makers alike (Cite: political actors stating this; [@Scheer:2012] Fridaysforfuture). As a renewable future is possible in both locally autarkic and continental grids but the techno-political paths to achieve them are very different, the question of which path to take is a real one. In this paper, we investigate the trade-off between geographical scale and cost: we know from previous research that a continental system has lower costs than systems with smaller scales, but we do not yet know how much and if local systems are techno economically feasible when looking at fluctuations. Here, we investigate the system performance of renewable electricity autarky on all political levels in Europe, from the continent-wide, via national to the regional level.

# Results

## Grid scale drives system costs

We first assess the impact of transmission grid scale on total system costs for a fully renewable electricity system in Europe based on solar, wind, and hydro power, and bioenergy from residuals, with pumped hydro power, short-term, and long-term electricity storage options. When we determine system layouts with minimal cost for each grid scale we find that system cost increases when grid scale decreases. For a continental scale transmission grid, we find system costs corresponding to levelised cost of electricity of slightly below 5 €ct/kWh. However, when we consider transmission grids span countries only, with 33 isolated national grids, the costs rise by nearly 40%. When we decrease grid scale even further and consider no exchange between 497 regions in Europe, costs rise to 70% above costs for the system with continental scale grid.

These system cost differences are based on expected technology costs for which exact values are not know. In our study we assume technologies are deployed at large scale in Europe and globally and thus uncertainty on technology costs stems mainly from uncertainty about cost decreases during technology deployment. When we consider the uncertainty of technology costs and cost of capital in a global sensitivity analysis, we find that the relationship between grid scale and system cost is robust: the continental grid always has cost lower than those of national grids. However, the extra cost of national scale grids varies and can be higher than 140%, but also below 110% of continental grid system costs in certain cases.

Lower cost of systems with larger grids can be explained by the ability to spatially balance renewable electricity supply and demand, and by the ability to trade electricity and thus generate electricity in places with beneficial resources and correspondingly lower costs of generation. To understand the impact of these two aspects on system cost we consider nine additional scenarios in which we restrict net imports and enforce local generation up to the point at which annual electricity demand in countries or regions is fulfilled entirely by local generation yet the transmission grid is still available to balance intra-annual fluctuations, see Figure @fig:scenario-space. We find that net imports have only minor impact on system cost, while the impact of spatial balancing is significant. Prohibiting net imports entirely leads to cost only 40% and 70% (SHOULD BE 6% and 15%) above the unconstrained case. When we allow to import up to 30% of annual electricity demand into countries or regions, cost decrease only slightly. Thus, generating electricity locally, in countries or regions, can be done with system costs very close to those of an electricity system in which electricity is generated in places with best wind and solar resources --- as long as a continental scale transmission grid exists that can be used sufficiently to balance fluctuations of renewable generation.

![Total system costs of twelve electricity system layouts in Europe, relative to an unconstrained layout with a transmission grid spanning the entire continent. Each layout is defined by three aspects: First, the scale of the transmission grid; whether it spans the entire continent (continental grid), does not cross country borders (national grids), or exists only within regions (regional grids). Second, the autarky scale which defines where all, or the majority of electricity is generated: anywhere on the continent (continental), or within countries (national) or regions (regional). Lastly, the maximum level of net electricity imports which determines the strictness of autarky: strict when none of local electricity demand can be satisfied by net imports, or weak when up to 15% or 30% of demand can be satisfied by imports.](report/scenario-space.png){#fig:scenario-space .class}

## Grid scale drives system composition and operation

Cost optimal electricity systems with different scales of the transmission grid also have different compositions and operations of supply and support technologies, with some technologies being more important for the large grids and some technologies being more important for smaller grids. Indeed, the differences in total system costs are driven by these system compositions.

Electricity systems with larger grids of course compose of more transmission capacities but also make use of these capacities more, see Table @tbl:overview-scenario-results-2. The continental scale system layout contains high voltage transmission capacity of nearly 400 TW km, which is roughly 4 times current capacity [@ENTSO-E:2019] (CURRENT NUMBER CORRECT?). Using this capacity, about 3000 TWh electricity are crossing country borders, which, again, is roughly 4 times current international electricity flows [@ENTSO-E:2019]. Of those 3000 TWh, 1400 TWh remain as net imports within countries, which is more than 7 times the amount of today [@ENTSO-E:2019]. All these numbers are smaller when electricity is generated locally, but transmission capacity and international electricity flows in those cases are still 2-3 times today's values. While not covered by our model, electricity systems with such high amount of cooperation require institutional support in the form of markets, grid governance, and long-term trade agreements.

```table
---
caption: 'Installed transmission grid capacity (<i class="fab fa-connectdevelop"></i>), gross physical electricity flow crossing country borders (<i class="fas fa-shopping-cart"></i> gross), and net electricity flow imported by all countries (<i class="fas fa-shopping-cart"></i> net) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the size of the electricity grid surrounding autarkic units (<i class="fab fa-connectdevelop"></i>). {#tbl:overview-scenario-results-2}'
alignment: LRRR
include: report/overview-scenario-results-2.csv
include-encoding: utf-8
markdown: True
---
```

Because of the ability to use the transmission grid to balance demand and renewable supply, electricity systems with larger grids require fewer other flexibility options. Indeed we find that the system with continental scale grid contains only moderate amounts of short-term and long-term electricity storage, bioenergy, and requires less curtailment of renewable electricity than systems with smaller grids, see Table @tbl:overview-scenario-results-1. Systems with larger grid extract the majority of their renewable electricity from onshore wind resources despite the higher installation costs and the higher levelised costs of wind capacities compared to solar electricity. Continental grids allow to exploit the weather pattern of wind all over Europe --- an effect that has been shown in extant literature [@Grams:2017] --- which simply do not exist for photovoltaics.

```table
---
caption: 'Installed capacities of solar (<i class="fas fa-sun"></i>), wind (<i class="fas fa-wind"></i>), biofuel (<i class="fas fa-leaf"></i>), short-term (<i class="fas fa-battery-three-quarters"></i> short) and long-term storage (<i class="fas fa-battery-three-quarters"></i> long) and relative curtailment of solar, wind, and hydro power (<i class="fas fa-traffic-light"></i>) for all electricity system layouts without net imports into autarkic units. Each layout is characterised by the autarky scale (<i class="fas fa-layer-group"></i>), the maximum level of annual net imports into autarkic units (<i class="fas fa-shield-alt"></i>), and the size of the electricity grid surrounding autarkic units (<i class="fab fa-connectdevelop"></i>). Each layout additionally contains fixed hydroelectricity capacities on their current locations: 36 GW run of river, 103 GW / 97 TWh reservoirs, and 54 GW / 1.3 TWh pumped hydro storage. {#tbl:overview-scenario-results-1}'
alignment: LRRRRRRRR
include: report/overview-scenario-results-1.csv
include-encoding: utf-8
markdown: True
---
```

Generally, there is a trend in generation capacity composition shifting from wind to solar capacities when moving along the axis of grid scale towards smaller transmission grids, see Table @tbl:overview-scenario-results-1. However, that is not the case because wind generation capacity is not needed, but because solar generation capacity is needed on top of that. Partly to compensate higher losses from electricity storage, but partly also because over-capacities of renewables together with curtailment are cost efficient options to balance fluctuations of renewable generation. Consequently, the system with regional grid size contains 50% more renewable generation capacities and curtails one and half times more potential renewable electricity than in the system with continental scale grid. Furthermore, small scale systems come with much higher capacities for biomass combustion, which are used to balance seasonalities of renewable electricity, mostly solar. Thus, biomass combustion competes in its use case with long-term storage and further solar capacities, albeit being more cost efficient in many circumstances.

The structural differences of systems on different grid scales explain differences in total system cost: large transmission grids are not necessarily the most cost efficient option to balance renewable fluctuations by themselves, but they add one more tool to the bouquet of options. No single option is the best, but it's an adequate mix of several options for different situations that is cost efficient.

The structural differences explain as well in parts the sensitivity of additional cost mentioned earlier. In fact, in the global sensitivity analysis we find that the uncertainty of three parameters is causing almost the entire variability of the additional cost of systems with national scale grid compared to the system with continental scale grid: cost of capital (total Sobol' index 0.4), cost of bioenergy combustion capacity (total Sobol' index 0.3), and cost of onshore wind capacity (total Sobol' index of 0.15). While the cost of capital impacts all components of all systems, it especially increases cost of transmission capacities due to their long lifetime. Large scale systems are relatively more sensitive to changes in cost of wind capacities, and small scale systems are more sensitive to changes in cost of bioenergy capacities. Thus, higher cost of capital and higher cost of wind capacities decrease relative costs of small scale systems, and relative cost of bioenergy capacity increase relative costs of small scale systems.

## Grid scale impacts local electricity systems

Locally, in countries and regions in Europe, the composition of the electricity system is very diverse but the effects we see are again dependent on the scale of the transmission grid. For the case of the continental grid, some countries and regions are equipped with much more generation capacities than others: Ireland, Lithuania, Estonia, and Albania all generate more than 4 times of what they actually consume. Large parts of their land are used to farm electricity for export. This effect is even more pronounced in single regions: several Irish counties facing the Atlantic ocean --- like Mayo, Kerry, or Cork ---, the Swedish island Gotland in the baltic sea, or Tulcea County at the Romanian shore of the black sea generate more than 50 times than what they need for themselves. Some of these exports are transmitted over long distances to their recipients, crossing many regions that merely pass through the electricity. In extreme cases, this amounts to more than 250 times local demand, for example in the central Ireland counties Laoighis, Tipperary, or Meath, where the electricity from onshore wind turbines at the Atlantic ocean is on its way to United Kingdom and eventually central Europe. These Irish regions, but also other regions in Romania, Estonia, Montenegro, Switzerland and other countries passing through large amounts of electricity are faced with large installations of transmission capacities which they do not, or only marginally benefit from directly. On the receiving end of these electricity trades are regions and countries that rely heavily, or sometimes completely, on imports: Belgium, Czech Republic, and Germany all feed less than 10% of their electricity demand with electricity generated locally.

On the other end of the scale, in a European electricity system with isolated regional transmission grids, levelised cost of electricity varies strongly between regions, depending on the renewable resources they have available, and depending on magnitude and shape of their electricity demand, see Figure @fig:map. About 10% of the regions has two times higher costs than in the continental grid case, and for half of them costs are four times those of the reference. Electricity supply in these regions is very expensive. There are three types of regions with high costs: regions with large hydro capacities, regions with low or no potential for wind capacities, and --- with the lowest cost --- some regions with often large potentials and vast amounts of available land, but whose electricity profile does not fit well to the solar and wind generation profiles and where consequently more investments into balancing options need to be made.

Because we enforce hydro capacities where they exist today, and only there (see Methods), regions with hydro generation exceeding local demand have high generation cost. This is the case for several Swiss cantons (Graubünden, Uri, Valais, Glarus), other alpine regions like Valle d'Aosta in Italy, but also in other parts of Europe, for example for Nikšic in Montenegro, or Jämtland in Sweden. Of course, this is an artefact of regional autarky, for which hydro stations in these regions have not been built for.

Several regions have low or no potential for wind capacities and thus their main source of electricity is the sun. But because solar generation has pronounced seasonal fluctuations in Europe, these regions have a high demand for balancing and consequently require more long-term electricity storage or flexible generation from bioenergy of hydroelectricity. Where the latter two are not available, generation costs are particularly high due to the high cost of long-term storage. This is the case in urban regions like Geneva, Prague, Budapest, and Bucharest. In Brussels, the Swiss canton Basel-Stadt, Berlin, Vienna, and Oslo the generation potential is too low for the region to be self-sufficient in the first place. To be able to fulfil their electricity demand, we connect them with their encompassing region or a neighbouring region, but even then some of them have high generation costs.

![Levelised costs of renewable electricity for strict continental autarky, strict national autarky, and strict regional autarky in relation to the baseline of continental autarky. Strict autarky does not allow for any form of exchange of electricity between autarkic units. Each panel shows the autarky level and the relative total system costs. Blue lines visualise net transfer capacities connecting regions in the cases of continental and national autarky.](report/map.png){#fig:map .class}

# Discussion

## Findings

1. Regional or national scale grids come with significantly higher costs than one that spans the continent.

2. (net) Autarky must not be much more expensive, as long as a sufficient transmission grid exists for balancing.

3. The scale of the grid defines the cost optimal composition of supply and support infrastructure: large grids come with wind, some storage, and some flexible biofuel. Small grids come with less wind, more solar, more storage and more flexible biofuel. (I would like to say: benefit of grid can be tapped only if the other infrastructure fits, but we do not know how much worse 'non-fitting' mixes are, see conclusions.)

ADD somehere here: "The availability of all balancing options on layouts with continental scale grid allows to choose most adequate solutions for all possible situations and thus allows to find a cost optimal mix, in contrast to layouts with smaller grids, where the range of options is limited."

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

The study area comprises all countries within member organisations ENTSO-E: EU-28, Norway, Switzerland, and Western Balkans countries. We exclude Iceland which is electricity autarkic, Cyprus, which is not directly connected to the rest of the study area, and Malta, for which insufficient data are available. We divide the study area into 497 regional administrative units [@Trondle:2019], and add a network node at the centroid of each region. We model the transmission grid as direct net transfer capacities between network nodes, i.e. we consider net power flows on the shortest distances between nodes only, and assume the distribution network within each node is able to handle distribution load. We allow transmission capacities between regional administrative units sharing a land border. We use currently existing sea connections and those that are currently under construction to connect regions that do not share a land border [@ENTSO-E:2019a]. We furthermore connect the islands Hiiu and Saare to the Estonian mainland to retrieve a fully connected electricity network graph as visible in Figure @fig:map.

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

We assume hydro run of river and hydro reservoir potentials to be largely tapped today (ADD source maybe JRC Technology Map) with almost no expansion potential. Thus, for hydro generation capacities we deviate from the greenfield approach and fix today's capacities. Similar to PV and wind generation, hydro generation can be curtailed.

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

Similar to hydro run of river and hydro reservoir capacities, we assume pumped hydro capacities in Europe to be largely tapped (ADD source, maybe JRC Technology Map) and do not allow for capacity expansion. Thus, we deploy today's pumped hydro power and storage capacities. We assume a round-trip electricity efficiency of 78% [@Schmidt:2019].

To determine location, power and storage capacity of each pumped hydro station in Europe today, we also use the JRC Hydro Power Database [@MatteoDeFelice:2019]. Where storage capacities are missing, we employ the same method as for hydro reservoirs: we assume national median ratios of power to storage capacity for all stations with missing storage capacity; and where this is not available, we assume Europe-wide median ratios of power to storage capacity. Storage capacities within the JRC Hydro Power Database sum up to more than 10 TWh which is more than what other sources report (for example 1.3 TWh in [@Geth:2015]). We thus scale the storage capacities to match national data reported by @Geth:2015. Using location data of each station, we then sum all power and storage capacities within regional administrative units to form a single pumped hydro capacity per unit.

## Short-term and long-term storage

We assume that short-term and long-term storage capacities can be deployed in all regional administrative units. We model short-term as Lithium-ion batteries and long-term as hydrogen storage as they are likely to become the dominant technology in their respective applications [@Schmidt:2019]. We models are based on two technical parameters: the ratio between power and storage capacity, and the round-trip efficiency. Short-term storage is permitted to a maximum capacity of 4 h of full power, while long-term storage has a minimum of 4 h capacity at full power. We assume 86% of round-trip efficiency for short-term and 40% for long-term storage.

Additionally, we assume that power and storage capacities can be expanded independently, only limited by the above mentioned minimum/maximum storage capacities.

## Insufficient potentials

In some regions, local technical potential for renewable electricity is not high enough to satisfy local electricity demand [@Trondle:2019]. This is problematic in system layouts in which regions strive for electricity autarky. To provide sufficient electricity supply in these regions, we connect them with a neighbouring or the encompassing region: Vienna with Lower Austria, Brussels with Flanders, Berlin with Brandenburg, Oslo with Akershus, and the canton of Basel-City with the canton of Basel-Country. So when we require regional autarky, we actually require autarky of each combined region in these five corner cases.

## Technology costs

We assess long term (quasi steady state) cost of electricity supply. We aim neither to determine costs of a transition to such states nor to consider disruptive developments on the global market for supply and storage technologies. Thus, we assume expected learning-rate based costs once renewable supply and electricity storage are deployed at the large scales consistent with our study. Cost estimations for the year 2050 are primarily from [@JRC:2014] for supply and transmission technologies, from [@Schmidt:2019] for storage technologies, and from [@RuizCastello:2015] for fuel costs of biofuels.

Technology costs are modelled as installed capacities costs, annual maintenance costs based on installed capacity, and variable costs per unit of generated electricity. For solar and wind we assume small variable costs of 0.1 €ct / kWh only to enforce curtailment whenever generation potential is higher than demand and storage capacities. We subtract these variable costs from the fixed operation and maintenance cost based on average capacity factors to not increase the overall cost of the technologies. Technology lifetime and costs of capital are used to derive annuities for each technology. For hydroelectricity, we consider annual maintenance and variable costs only as we assume capacities are already built today and as installation costs of hydroelectricity have no impact on our results. See Table @tbl:overview-cost-assumptions for an overview of all cost assumptions.

We assume costs of capital to be 7.3% for all technologies and all locations based on historic average cost of capital for OECD countries [@Steffen:2019]. Some recent literature suggests costs of capital are likely specific to technology [@Egli:2018; @Steffen:2019] and location [@Ondraczek:2015; @Steffen:2019], but we consider the data available so far too sparse to provide a solid basis on which to model this.

```table
---
caption: 'Assumptions on technology costs. ^AC transmission installation costs are given in [€/kW/1000km] {#tbl:overview-cost-assumptions}'
alignment: LRRRRRR
include: report/overview-cost-assumptions.csv
include-encoding: UTF-8
markdown: True
---
```

## Political constraints

We use two types of political constraints as the basis for our system layouts: constraints on the scale of the electricity grid, and constraints on net electricity imports. When we limit the scale of the grid to continental, national, or regional sizes, we force transmission capacities between continents, countries, or regions to be zero. For the continental scale this is given inherently by the scope of our study area.

When we limit net imports into continental, national, or regional areas leading to full or partial autarky, we force all or parts of the electricity demand on the continent, in the country, or in the region to be satisfied with local electricity generation from wind, sun, biomass, and water. In contrast to the constraint on the grid scale, the grid can still be used to balance fluctuations, as long as net annual imports stay below the required threshold.

## Sensitivity to meteorological conditions

We assess the sensitivity of our main result to meteorological conditions impacting generation from wind and solar power. While we use only a single year of meteorological conditions, 2016, in the analysis of system layouts, we use ten years, 2007--2016, in this sensitivity analysis. In this way, we are not assessing the variability between meteorological years, but we are assessing how much the result changes when considering a long duration of meteorological conditions. This will let us find electricity system layouts that are optimal for the entire range of years, not only for single years.

We assess the sensitivity of relative total system costs of national autarky with continental autarky as a baseline. Because computational requirements to solve a model with regional spatial resolution and temporal resolution of 4h for ten years are very high, we perform this sensitivity analysis using a model with national spatial resolution while keeping temporal resolution the same. Comparing the results between national and regional resolution for the case with only one year of meteorological conditions, we find a difference of 6% for the relative costs of national autarky.

The additional cost of national autarky compared to continental autarky, however, is unaffected by the longer time duration: the difference to the case with only one year is negligible (< 3‰). This is not to say that the electricity system layout is not sensitive to meteorological conditions. In fact, we find that total system costs are generally slightly higher, and more wind and biofuel capacities are deployed in exchange of solar capacities. However, large scale and small scale system layouts are impacted similarly and so the difference between both is unaffected by the longer time duration of considered meteorological conditions. These results justify the use of only a single meteorological year.

## Sensitivity to technology costs

Ignore correlation of uncertain inputs because "end of learning" type of costs are assumed, so the link between cost parameters based on political decisions can be ignored. There is still a link through developments on markets of raw material for example, but we ignore that.

* depreciation rate [i]: 1.9—13.5% max and min for PV and wind between 2009 and 2017 [@Steffen:2019]
* installation costs of PV [c_pv]: 280—580 €/kW [@JRC:2014] (Table 7)
* installation costs wind onshore [c_wind]: 800-1700 €/kW [@JRC:2014] (Table 4)
* installation costs wind offshore [c_offshore]: 1790–3270 €/kW [@JRC:2014] (Table 5)
* installation costs short term storage power [c_sts_p]: 611.324 * (0.05--0.23) = 31 -- 141 €/kW [@Schmidt:2019]
* installation costs short term storage energy [c_sts_e]: 723.130 * (0.05--0.23) = 36--166 [@Schmidt:2019]  
* installation costs long term storage power [c_lts_p]: 4884.287 * (0.23--0.43) = 1123-2100
 €/kW [@Schmidt:2019]
* installation costs long term storage energy [c_lts_e]: 27.951 * (0.23--0.43) = 6 - 12 €/kWh [@Schmidt:2019]
* installation costs net transfer capacities [c_ntc]: 700—1080 € / kW / km [@JRC:2014]
* installation costs biofuel [c_bio]: 1380—3450 €/kW [@JRC:2014](Table 48)
* fuel costs biofuel [c_biofuel]: 25.02–30.58 €/MWh [@RuizCastello:2015]
* availability biofuels [a_biofuel]: [low/reference/high availability] [@RuizCastello:2015]

We furthermore assess the uncertainty of our result stemming from uncertainty of technology costs. While we do know current costs and we do know that costs are likely to fall with deployment due to learning effects, we do not know exact future costs with certainty. This is mainly because it is unknown exactly how much technologies are deployed, and exactly how much costs will fall with deployment. Because in our analysis we are assuming technologies are deployed at large scale, the first aspect does not lead to much uncertainty. Uncertainty of technology costs in our analysis thus stems largely from the uncertain relationship of deployment and cost decreases. Being a cost optimisation, the absolute total system costs of any assessed electricity system layout is sensitive to costs of technologies (SOURCE). Here, we assess the sensitivity of the cost difference between large and small scale system layouts.

We consider ten technology cost parameters and the weighted cost of capital as uncertain and describe their uncertainty by an uniform distribution over ranges taken from literature (see Methods). We perform a global sensitivity analysis of the additional costs of national autarky compared to continental autarky in this eleven dimensional space. We find that national autarky is 20%--60% (FAKE RESULT) more expensive and that this range is largely driven by parameters X, Y, and Z. Higher values of X lead to higher relative costs of national autarky, while higher values of Y lead to lower relative costs of national autarky. High values of Z lead to higher relative costs of national autarky, but only if X is high as well.

# Bibliography
