# Swiss distribution network system
The Swiss distribution network system is derived from an actual distribution network in Zurich. We are sorry that it is inconvenient to disclose the data of distribution network. The data could be found upon request.

In the calculation of our work, only the loadings of low-voltage cables exceeded the nameplate rating, so we only calculated the capacities of low-voltage cables in Scenario DTR.
The parameters of cables and transformers are included in **Swiss case-equipment parameters.xlsx**.

cabletem758.m is the code to calculate the temperatures of cables with nameplate capacity of 0. 5514 MW in each time period (nameplate capacity of 0. 5514 MW is equivalent to 758 A).

In Scenario DTR, we perform power flow calculation with the maximum available output power of PV systems and calculate the temperatures of underground cables and transformers at each time period.
The temperatures of underground cables are calculated by **cabletem758.m**, while the codes to calculate the temperatures of transformers are in **f_DTR.m**.
If the temperatures of underground cables and transformers exceed the upper limits, it means that the equipment loading exceeds the DTR of the equipment and the PV output power needs to be curtailed. If the nodal voltages exceed the upper limits, the PV output power also needs to be curtailed. We gradually increase the curtailment ratio and calculate the corresponding equipment temperatures and nodal voltages until all operating constraints are satisfied.

