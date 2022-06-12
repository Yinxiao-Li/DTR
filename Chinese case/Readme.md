# Chinese distribution network system
The Chinese distribution network system is taken from https://doi.org/10.16081/j.issn.1006-6047.2019.07.007.
The original distribution network system consists of three feeders, which represent the residential load, industrial load and commercial load.
We chose the feeder that represents the residential load to conduct the case study.

**CaseChinese.mat** is the data file in MATPOWER format, including the parameters of the distribution network, such as the bus load, bus voltage limit, bus voltage level, branch impedance and branch nameplate capacity (STR).

The source dataset does not provide detailed hourly load data for a year and only provides the load data of typical days (spring rainy, spring sunny, summary rainy, summer sunny, fall rainy, fall sunny, winter rainy and winter sunny).
The type of day is determined based on the date and precipitation data.

It is assumed that there are only medium-voltage lines and no low-voltage lines in the Chinese distribution network because the source paper does not provide the data of low-voltage lines. 
The capacities of secondary transformers are set according to the peak loads of customers.

The branches numbered 1 to 48 are lines, which are assumed to employ overhead lines.
The branches numbered 49 to 96 are secondary transformers.

**DLR.m** is the code to calculate the capacity of overhead lines in Scenario DTR.

The parameters of overhead lines and transformers are included in **Chinese case-equipment parameters.xlsx**.

In Scenario DTR, we perform power flow calculation with the maximum available output power of PV systems and calculate the loading of overhead lines and the temperatures of transformers at each time period. 
The codes to calculate the temperatures of transformers are in **f_DTR.m**.
If the loading of overhead lines exceeds the upper limit calculated by DLR.m or the temperatures of transformers exceed the upper limits, it means that the equipment loading exceeds the DTR of the equipment and the PV output power needs to be curtailed. If the nodal voltages exceed the upper limits, the PV output power also needs to be curtailed. We gradually increase the curtailment ratio and calculate the corresponding equipment loading and temperatures and nodal voltages until all operating constraints are satisfied.

