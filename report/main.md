# Introduction

Europe is on its path towards a future of a fully renewable electricity system. While the target is less and less questioned, exactly how this future should look like and what the consequences of decisions for one of the possible futures are, is not completely known. One such unknown aspect of electricity systems of the future is their geographical scale: should Europe focus on local structures, or should it instead think big, ensure strong collaboration, and embrace trade between all countries from Portugal to Poland and from Sweden to Spain? While Europe is yet making early steps on its path, we are analysing this aspect of future electricity systems.

The electricity system comprises consumers and producers of electrical energy, but also technical infrastructure and institutions allowing for transmission, distribution, and trade of the energy. It forms a network of connected consumers and producers. From the perspective of a single stakeholder within the system, parts of the total system can be divided into local elements and global elements. For example, seen from Slovenia, all Slovenian consumers and producers are local, and all others are global. Furthermore, Slovenia may decide to trade with other countries, but not with all other countries in Europe. For Slovenia, two questions arise: Should it trade at all, or rather minimise trade and focus on local electricity generation? And if it decides to trade, with whom should it trade, i.e. how large should its trade network be?

We are investigating the cost aspect of these questions, i.e. how total energy system costs change with different scales of the trade network of Slovenia. Total costs may be lower with more excessive trading, as electricity may be generated with lower costs outside of Slovenia. We investigate the costs not only for Slovenia but for all Europe and all scales, so that in addition we can also quantify the change in system costs when not countries like Slovenia are minimising trade, but regions within Slovenia and all Europe. Geographic scale of the electricity system thus comprises of two dimensions: the size of the subsystems, and the size of the trade networks surrounding the subsystems. Using this definition, we are answering the question: how does geographical scale of renewable electricity systems impact total system costs?

# Methods

We will quantify the impact of geographical scale on system costs using a model of the European electricity system. The model optimises the design of the system by deciding on capacities for generation of renewable electricity, capacities for electricity storage, and capacities for transmission of electricity. As it minimises costs, it is a best case consideration given the constraining factors. One of the most important constraining factors is the geographical scale, which we vary to form scenarios. Because many other factors are uncertain, we perform an uncertainty quantification to verify the robustness of the results.

## Model

The model is simple representation of the European electricity system, where each node represents a European region (502 in all Europe) including its demand for electricity. Electricity demand and all other system processes are modelled with a time resolution of one hour. Neighbouring nodes can be connected with transmission lines which allow energy to flow between the nodes. Applying Kirchhoff's first law, the model makes sure that the sum of all energy flows at each node is zero at each point in time.

On each node, renewable electricity generation capacity can be build: wind on- and offshore, open field and roof-mounted photovoltaics (PV), hydro power, and biomass incineration. For all types, time series are taken from [@Trondle:2019]. We take maximum generation capacities for wind and solar power from the same study, and amend it with potential estimations for hydro power and biomass incineration.

Furthermore, storage capacities can be build on each node: Pumped hydro storage can be build where it already exists today (assuming its potential is largely exceeded today). And in addition, two idealised models of other storage forms can be deployed. First, a model for the provision of flexibility to balance daily fluctuations with an energy to power ratio of 3. Second, a model for the provision of flexibility to balance yearly seasonalities with an energy to power ratio of 500. These idealised models are to represent, respectively, battery storage or similar and power-to-gas or similar. Their technical parameter values like are derived from [@FfE:2016].

The optimiser of our model thus has three types of decision variables: renewable generation capacities, electricity storage capacities, and electricity transmission capacities. All capacities come with installation costs and some with variable costs during operation. Based on the decision variables, the optimiser is to find the system design with lowest total system costs.

## Scenarios

To assess the impact of geographical scale, we form scenarios and quantify the minimal system costs given a geographical scale. As we are using administrative divisions to form subsystems of the electricity system, we consider the following scales in this study:

* continental,
* national (34 within Europe),
* and regional (502 within Europe).

We use these scales for sub dividing the European system, but also for defining the trade network for each sub division: for the regional sub division for example, we have one scenario with a regional size trade network (i.e. no trade), a national trade network (i.e. each region can trade with all other regions in the same country), and a continental trade network.

In addition, we also consider different intensities of trade, ranging from allowing net zero trade to allowing up to 30% electricity to be imported. Together with the size of the divisions and the size of the trade network, this gives our scenario space a third dimension. The scenario space and hypothesised results for each scenario are shown in Figure @fig:scenario-space.

![Sketch of scenario space including hypothesised total system costs for each scenario.
](../build/scenario-space.png){#fig:scenario-space .class}


## Uncertainty Quantification

# Bibliography
