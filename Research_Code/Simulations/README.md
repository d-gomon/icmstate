# Simulations
This folder contains the code used to perform the Simulation study.

The `Simulation_Code` folder contains the code used to run the simulation studies on a high performance computing clusted and the `ScenarioX` folders contain the code used to extract and plot statistics from the resulting simulation studies.



## Description of the folders

**Scenario1** 

In Scenario 1 we consider:

* Illness-death model
* All subjects start in state 1. Alive
* No transitions are exactly observed
* Transition intensities are Exponential with rates $0.1, 0.05$ and $0.1$ respectively ($1 \to 2, 1 \to 3, 2 \to 3$).


**Scenario2** 

In Scenario 2 we consider:

- Illness-death model
- Subjects start with equal probability in state 1. Alive or 2. Illness.
- No transitions are exactly observed
- Transition intensities are Exponential with rates $0.1, 0.05$ and $0.1$ respectively ($1 \to 2, 1 \to 3, 2 \to 3$).

**Scenario3** 

In Scenario 3 we consider:

* Illness-death model
* All subjects start in state 1. Alive
* No transitions are exactly observed
* Transition intensities are Weibull with scale/shape $0.5/\frac{1}{\sqrt{5}}, 0.5/\frac{1}{\sqrt{10}}5$ and $2/\left( \frac{\Gamma(1.5)}{10}  \right)^2$ respectively ($1 \to 2, 1 \to 3, 2 \to 4$).

**Scenario4** 

In Scenario 4 we consider:

* Extended Illness-death model
* All subjects start in state 1. Alive
* Transitions to state 3. Death and 4. Death after illness are exactly observed
* Transition intensities are Exponential with rates $0.1, 0.05$ and $0.1$ respectively ($1 \to 2, 1 \to 3, 2 \to 4$).


