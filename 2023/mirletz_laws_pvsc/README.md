This repository contains code and data used in the analysis of
Mirletz, Brian T.; Laws, Nicholas D. Impacts of Dispatch Strategies and Forecast Errors on the Economics of
Behind-the-Meter PV-Battery Systems, IEEE PVSC 2023

Running the code in this repository requres:
1. System Advisor Model (SAM) https://sam.nrel.gov/ - version 2022.11.22r2 was used in the analysis
2. Julia (1.8.5 was used)
3. PySAM (4.2 was used), package name nrel-pysam

- system_sizing.ipynb was used to size the systems and generate full-year dispatch
- run_mpc.jl and the SAM gui were used to generate other dispatch plans