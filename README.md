# roughWF_UDF
This repository contains the implementation of the adaptive roughness parameter calculation presented in INSERT CITATION, along with two tutorial cases.
The methodology presented in the above paper was implemented in Ansys Fluent, release 2023 R1, using User-Defined Functions (UDF).
The contents of this repository:
- README.md
- UDF.c
- tutorials:
  - case_A: model-scale channel flow with a power-law velocity distribution, a logarithmic-polynomial TKE distribution and vertical blending using a modified vertical coordinate
  - case_E: full-scale channel flow with a logarithmic velocity distribution and a homogeneous TKE distribution

How to run the tutorial cases:
1. Download the complete "tutorials/case_A" or "tutorials/cases_E" directory
2. By default, the case and data files are set up for using Method 4 with the parameter selection Ks+ = 1000. To run the simulation with these settings, the user can read the "run_cases.jou" journal file into Fluent, which loads the initialised case and data files, recompiles and loads the libudf directory, runs 2000 iterations in steady-state, saves the case and data files for the final iteration and exports the distribution of the streamwise velocity and the turbulent kinetic energy at the outlet. To run the simulation with a different method or a Ks+ value, the user has to follow steps 2.1 - ...
2.1. To change the fixed Ks+ value, the desired value of the "ks_plus_array" variable must be changed in line 66 of the "case_A/ABL_adaptive_roughness_case_A.c" UDF file or line 43 of the "case_E/ABL_adaptive_roughness_case_E.c" UDF file
3. change the target Ks+ in line 66 of the "case_A/ABL_adaptive_roughness_case_A.c" UDF file or line 43 of the "case_E/ABL_adaptive_roughness_case_E.c" UDF file
