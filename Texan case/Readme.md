# Texan distribution network system
The Texan test distribution network system is selected by a three-phase balanced distribution feeder in the Combined Transmission and Distribution Synthetic Dataset (https://electricgrids.engr.tamu.edu/) as the Texan test system, named p5uhs27_1247--p5udt26583.

**CaseTexan.mat** is data file in MATPOWER format, including the parameters of the distribution network, such as the bus load, bus voltage limit, bus voltage level, branch impedance and branch nameplate capacity (STR).

Since the source dataset does not provide detailed hourly load data for a year, we assume that the hourly load factor (the ratio of actual load to the maximum load in one year) of each node is the same, which is included in **Pload.mat**.

The branches numbered 1 to 38 are lines, which are assumed to employ underground cables.
The lines with nameplate capacities of 0.4573 MW and 0.3076 MW are low-voltage cables, while other lines are medium-voltage cables.
The branches numbered 39 to 52 are transformers.

In the calculation of our work, only the loadings of low-voltage cables exceeded the nameplate capacity, so we only calculated the capacities of low-voltage cables in Scenario DTR.

The parameters of cables and transformers are included in **Texan case-equipment parameters .xlsx**.

**cabletem370.m** and **cabletem550.m** are the codes to calculate the temperatures of cables with nameplate capacities of 0. 3076 MW or 0. 4573 MW in each time period (nameplate capacity of 0.3076 MW is equivalent to 370A, while nameplate capacity of 0. 4573 MW is equivalent to 550A).

In Scenario DTR, we perform power flow calculation with the maximum available output power of PV systems and calculate the temperatures of underground cables and transformers at each time period. The temperatures of underground cables are calculated by **cabletem370.m** and **cabletem550.m**, while the codes to calculate the temperatures of transformers are in **f_DTR.m**.

If the temperatures of underground cables and transformers exceed the upper limits, it means that the equipment loading exceeds the DTR of the equipment and the PV output power needs to be curtailed. If the nodal voltages exceed the upper limits, the PV output power also needs to be curtailed. We gradually increase the curtailment ratio and calculate the corresponding equipment temperatures and nodal voltages until all operating constraints are satisfied.
