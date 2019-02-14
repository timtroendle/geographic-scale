# Introduction

Europe is on its path towards a future of a fully renewable electricity system. While the target is less and less questioned, exactly how this future should look like and what the consequences of decisions for one of the possible futures are, is not completely known. One such unknown aspect of electricity systems of the future is their geographical scale: should Europe focus on local structures, or should it instead think big, ensure strong collaboration, and embrace trade between all countries from Portugal to Poland and from Sweden to Spain? While Europe is yet making early steps on its path, we are analysing this aspect of future electricity systems.

The electricity system comprises consumers and producers of electrical energy, but also technical infrastructure and institutions allowing for transmission, distribution, and trade of the energy. It forms a network of connected consumers and producers. From the perspective of a single stakeholder within the system, parts of the total system can be divided into local elements and global elements. For example, seen from Slovenia, all Slovenian consumers and producers are local, and all others are global. Furthermore, Slovenia may decide to trade with other countries, but not with all other countries in Europe. For Slovenia, two questions arise: Should it trade at all, or rather minimise trade and focus on local electricity generation? And if it decides to trade, with whom should it trade, i.e. how large should its trade network be?

We are investigating the cost aspect of these questions, i.e. how total energy system costs change with different scales of the trade network of Slovenia. Total costs may be lower with more excessive trading, as electricity may be generated with lower costs outside of Slovenia. We investigate the costs not only for Slovenia but for all Europe and all scales, so that in addition we can also quantify the change in system costs when not countries like Slovenia are minimising trade, but regions within Slovenia and all Europe. Geographic scale of the electricity system thus comprises two dimensions: the size of the subsystems, and the size of the trade networks surrounding the subsystems. Using this definition, we are answering the question: how does geographical scale of renewable electricity systems impact total system costs?

# Methods

We will quantify the impact of geographical scale on system costs using a model of the European electricity system. The model optimises the design of the system by deciding on capacities for generation of renewable electricity, capacities for electricity storage, and capacities for transmission of electricity. As it minimises costs, it is a best case consideration given the constraining factors. One of the most important constraining factors is the geographical scale, which we vary to form scenarios. Because many other factors are uncertain, we perform an uncertainty quantification to verify the robustness of the results.

## Model

The model is simple representation of the European electricity system, where each node represents a European region (502 in all Europe) including its demand for electricity. Electricity demand and all other system processes are modelled with a time resolution of one hour. Neighbouring nodes can be connected with transmission lines which allow energy to flow between the nodes. Applying Kirchhoff's first law, the model makes sure that the sum of all energy flows at each node is zero at each point in time.

On each node, renewable electricity generation capacity can be build: wind on- and offshore, open field and roof-mounted photovoltaics (PV), hydro power, and biomass incineration. For all types, time series are taken from [@Trondle:2019]. We take maximum generation capacities for wind and solar power from the same study, and amend it with potential estimations for hydro power and biomass incineration.

Furthermore, storage capacities can be build on each node: Pumped hydro storage can be build where it already exists today (assuming its potential is largely exceeded today). In addition, two idealised models of other storage forms can be deployed. First, a model for the provision of flexibility to balance daily fluctuations with an energy to power ratio of 3. Second, a model for the provision of flexibility to balance yearly seasonalities with an energy to power ratio of 500. These idealised models are to represent, respectively, battery storage or similar and power-to-gas or similar. Their technical parameter values like are derived from [@FfE:2016].

The optimiser of our model thus has three types of decision variables: renewable generation capacities, electricity storage capacities, and electricity transmission capacities. All capacities come with installation costs and some with variable costs during operation. Based on the decision variables, the optimiser is to find the system design with lowest total system costs.

## Scenarios

To assess the impact of geographical scale, we form scenarios and quantify the minimal system costs given a geographical scale. As we are using administrative divisions to form subsystems of the electricity system, we consider the following scales in this study:

* continental,
* national (34 within Europe),
* and regional (502 within Europe).

We use these scales for sub dividing the European system, but also for defining the trade network for each sub division: for the regional sub division for example, we have one scenario with a regional size trade network (i.e. no trade), a national trade network (i.e. each region can trade with all other regions in the same country), and a continental trade network.

In addition, we also consider different intensities of trade, ranging from allowing net zero trade to allowing up to 30% electricity to be imported. Together with the size of the divisions and the size of the trade network, this gives our scenario space a third dimension and 12 scenarios in total. The scenario space and hypothesised results for each scenario are shown in Figure @fig:scenario-space.

![Sketch of scenario space including hypothesised total system costs for each scenario.
](../build/scenario-space.png){#fig:scenario-space .class}

## Uncertainty Quantification

Qualitatively, and because of the method described above, the answer to our research question is known even before we start our work: small scale designs of the electricity system will always be equally costly, or -- more likely -- more costly than large scale designs. That is because theoretically in our model, small scale designs add constraints to an optimisation model which can make the solution only worse. So we already know parts of the answer, but we do not know exactly how much more costly small scale will be.

But we cannot find out exactly how much more costly small scale will be. That is because we know neither how costly small scale will be, nor do we know how costly large scale will be. This is mostly due to parametric uncertainty: for many parameters, the correct value (whatever that means) is not know. Some of those unknown parameters will have a large impact: the costs of generation, storage, and transmission capacity for example. Solar photovoltaics may cost around 500 €/kW in Europe in the foreseeable future, but it may cost as well more than 1000 €/kW. Considering that a third or more of European electricity could be produced with photovoltaics, this will impact total system costs strongly. Unfortunately the value is not only unknown (as in: try harder finding out), but it is also partly in our hands (as a society or as consumers) and based on our future decisions the value may stay at above 1000 €/kW or may come down to 500 €/kW.

This means our study cannot tell how much more costly small scale will be in the future compared to large scale. But it can do three things, each adding a little more information. First, it can tell how much more costly small scale will be, depending on some assumed parameter values. Second, assuming we know reasonable ranges for the unknown parameters and assuming we can run the model sufficiently often, we can show sensitivities of the results to parameter values (something like: small scale is always much more expensive, but not if battery prices drop below 100 €/kW). Third, assuming we can furthermore attach probability densities to the parameter ranges, we can determine probabilities of cost differences between small and large scale (something like: with 90% probability, small scale will cost 0.05 €/kWh more than large scale). Some more details about these three types of findings below.

### Isolated sets of parameter values

Because we assume the uncertainty of our model results stems mainly from the parameters, not the model itself, we can always choose values for all parameters and we will receive a result which is accurate. But this won't allow us to say anything about the probability of the result being close to reality, or about how strongly results vary should parameter values vary.

### Parameter sensitivity

If we are able to run the model many times, we can determine the sensitivity of the result to changes in the parameter values --- changes of the values of single parameters, but also to combined changes of several parameters. This will allow us to rank parameters -- solely or in combination. For example, we will be able to say that the most important parameters are the cost of batteries and the cost of transmission grid, but that the costs of batteries are only important when the costs for photovoltaics are below 750 €/kW. We can also visualise the relationship of the most important parameters on the result, see Figure @fig:parameter-sensitivity. Limiting the parameter values to possible values only, will furthermore allow us to quantify the range within which the result will lie: small scale will at least be 0.02 €/kWh more costly than large scale, but not more than 0.09 €/kWh.

![Sketch of the relationship of the two most important parameters on the additional costs of small scale over large scale electricity grids in Europe [€/kWh].
](../build/parameter-sensitivity.png){#fig:parameter-sensitivity .class}

Considering that we can easily identify 10 parameters (or more) for which we would like to analyse the sensitivity, to be able to actually analyse it, we likely need 100,000s of model runs or more even when using smarter methods than brute force. Our model takes hours to run, maybe even days. So reasonably, we can run it in the order of 10s of times. A sensitivity analysis as described above is thus only possible with a surrogate model: a model that resembles ours in its input-output behaviour, but runs orders of magnitude faster.

### Probability of results

Should we know the probability density of the input parameters, we can additionally determine the probability of extra costs for small scale designs using Monte Carlo and the surrogate model demanded above. Additionally to the range 0.02 €/kWh to 0.09 €/kWh from above, we could for example also say that with 90% probability the additional cost of small scale are more than 0.05 €/kWh and with only 10% probability they are above 0.07 €/kWh.

Can we assume any reasonable probability density of the input parameters? I do not know. But indeed we do know a few things on their shapes. First, we do know ranges: all components are very unlikely to be more costly in the future than they are now, and they will always stay above 0. Second, we also know that the extreme ends of that range are less likely than values in the center of the range. Maybe we can find out more about the probability density, and maybe that is already all we need to know.

# Bibliography
