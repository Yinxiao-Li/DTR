# Codes and Data of Improving Distributed PV Integration with Dynamic Thermal Rating of Power Distribution Equipment

In the folder of each case, there are 4 folders named **code**, **weather_from_2011_to_2020**, **weather_2025** and **weather_2100**.

The code folder includes the parameters of distribution systems and the codes to perform the one-year operating simulation and calculate the installed capacities and investment benefits of customers’ PV systems.

**test_main.m** is the code file to calculate the installed capacities and investment benefits of customers’ PV systems.

**f_STR.m** and **f_DTR.m** are the code files to perform the one-year operating simulation and calculated PV curtailment ratio in Scenario STR and Scenario DTR, respectively.

In the **test_main.m** file, the code **weather = xlsread('weather_2020.xlsx');** refers to reading weather data in 2020. 

**weather_2020.xlsx** can be replaced to other weather data files to calculate the installed capacities and investment benefits of customers’ PV systems using weather data of other years.

The weather data from 2011 to 2020 are included in the folder named **weather_from_2011_to_2020**.
The weather data in 2025 are included in the folder named **weather_2025**.
The weather data in 2100 are included in the folder named **weather_2100**.

**Capacity2100.m** is the code file to calculate the STR of equipment in 2100 according to the weather data in 2025 and 2100.

To calculate the installed capacities and investment benefits of customers’ PV systems using weather data in 2100 in Scenario STR (assuming the equipment capacity as STR), the STR of equipment in 2100 should be calculated first and then brought into the distribution network parameters file and then run the **test_main.m**.

**policy.m** is the code file to calculate the installed capacities and investment benefits of customers’ PV systems under the other two PV policies.

**ESS.m** is the code file to estimate the capacity of energy storage systems to achieve the same improving effects as the application of DTR.

The detailed descriptions about the parameters of distribution networks are included in the Readme file of each folder.
