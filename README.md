# COVID workplace simulation
This code constructs and implements a SEIRS driven simulation. The simulation produces the percentage infected,  peak percentage hospitalized and fatality rate in office enviroments during a disease outbreak. The simulation is tuned to the disease parameters of the current COVID-19 pandemic but can be easily tweaked to study other (fictive)  diseases and the effect of certain interventions on the percentage infected, peak percentage hospitalized and fatality rate in office or other enviroments. The code is an extension of the wonderful work that ryansmcgee performed. 

## Overview

#### Full documentation of this package's models, code, use cases, examples, and more can be found on [the wiki](https://github.com/ryansmcgee/seirsplus/wiki/)

[**Basic SEIRS Model**](https://github.com/ryansmcgee/seirsplus/wiki/SEIRS-Model-Description) | [**Extended SEIRS Model**](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description)
:-----:|:-----:
<img src="https://github.com/ryansmcgee/seirsplus/blob/master/images/BasicSEIRS_compartments_padded.png" width="400"> | <img src="https://github.com/ryansmcgee/seirsplus/blob/master/images/ExtSEIRS_compartments.png" width="400">
[Model Description](https://github.com/ryansmcgee/seirsplus/wiki/SEIRS-Model-Description) | [Model Description](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description)
[`SEIRSNetworkModel` docs](https://github.com/ryansmcgee/seirsplus/wiki/SEIRSModel-class)<br>[`SEIRSModel` docs](https://github.com/ryansmcgee/seirsplus/wiki/SEIRSNetworkModel-class) | [`ExtSEIRSNetworkModel` docs](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-class)

### SEIRS Dynamics

The foundation of the models in this package is the classic SEIR model of infectious disease. The SEIR model is a standard compartmental model in which the population is divided into **susceptible (S)**, **exposed (E)**, **infectious (I)**, and **recovered (R)** individuals. A susceptible member of the population becomes infected (exposed) when making a transmissive contact with an infectious individual and then progresses to the infectious and finally recovered states. In the SEIRS model, recovered individuals may become resusceptible some time after recovering (although re-susceptibility can be excluded if not applicable or desired). 
<p align="center">
  <img src="https://github.com/ryansmcgee/seirsplus/blob/master/images/BasicSEIRS_compartments_resus.png" width="400"></div>
</p>

### Extended SEIRS Model

This model extends the classic SEIRS model of infectious disease to represent pre-symptomatic, asymptomatic, and severely symptomatic disease states, which are of particular **relevance to the SARS-CoV-2 pandemic**. In this extended model, the infectious subpopulation is subdivided into **pre-symptomatic (*I<sub>pre</sub>*)**, **asymptomatic (*I<sub>asym</sub>*)**, **symptomatic (*I<sub>sym</sub>*)**, and **hospitalized (severely symptomatic, *I<sub>H</sub>*)**. All of these *I* compartments represent contagious individuals, but transmissibility, rates of recovery, and other parameters may vary between these disease states.

<p align="center">
  <img src="https://github.com/ryansmcgee/seirsplus/blob/master/images/ExtSEIRS_compartments.png" width="600"></div>
</p>


### Testing, Tracing, & Isolation

The effect of isolation-based interventions (e.g., isolating individuals in response to testing or contact tracing) are modeled by introducing compartments representing quarantined individuals. An individual may be quarantined in any disease state, and every disease state has a corresponding quarantine compartment (with the exception of the hospitalized state, which is considered a quarantine state for transmission and other purposes). Quarantined individuals follow the same progression through the disease states, but the rates of transition or other parameters may be different. There are multiple methods by which individuals can be moved into or out of a quarantine state in this framework.

<p align="center">
  <img src="https://github.com/ryansmcgee/seirsplus/blob/master/images/BothSEIRS_compartments_quarantine.png" width="800"></div>
</p>

<a name="model-network"></a>


### Possible implementations of the code

The code can be modified to your own liking to simulate situations that are tailored to your research. Parameters are based on literature but can be tweaked if literature is outdated or own found literature gives reason to doubt about chosen parameters. The disease parameters are specified for the COVID-19 virus but can be varied to research different diseases in a working enviroment.
