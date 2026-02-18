# roughWF_UDF
This repository contains the implementation of the adaptive roughness parameter calculation presented in INSERT CITATION, along with two tutorial cases.
The methodology presented in the above paper was implemented in Ansys Fluent, release 2023 R1, using User-Defined Functions (UDF).
The "src" directory contains the UDF files, which have been set up for cases A and E, respectively.
The "tutorials" directory contains two tutorial cases:
- case_A: model-scale channel flow with a power-law velocity distribution, a logarithmic-polynomial TKE distribution and vertical blending using a modified vertical coordinate (case A in the referenced publication)
- case_E: full-scale channel flow with a logarithmic velocity distribution and a homogeneous TKE distribution (case E in the referenced publication)

Each tutorial case contains the following files and directories:
- 2D_channel.msh, Silsoe_2D_channel.msh: mesh files for cases A and E, respectively
- ABL_adaptive_roughness_case_A.c or ABL_adaptive_roughness_case_E.c: UDF file, which specifies:
  -  the inflow conditions for the k-omega and k-epsilon model families
  -  on-demand executables to patch the inflow distribution for the whole domain during initialisation
  -  the adaptive calculation of the Ks and Cs roughness parameters using all four methods presented in the publication
- libudf: compiled UDF directory
- tutorial_case_initialized.cas.h5: initialised case file
- tutorial_case_initialized.dat.h5: initialised data file
- setup_case.jou: journal file containing text commands for setting up the simulation
- run_case.jou: journal file containing text commands for running the simulation

How to run the tutorial cases with the default settings for the roughness parameter calculations:
1. Download the complete "tutorials/case_A" or "tutorials/cases_E" directory
2. By default, the case and data files are set up for using Method 4 with the parameter selection Ks+ = 1000. To run the simulation with these settings, the user can read the "run_cases.jou" journal file into Fluent, which loads the initialised case and data files, recompiles and loads the libudf directory, runs 2000 iterations in steady-state, saves the case and data files for the final iteration and exports the distribution of the streamwise velocity and the turbulent kinetic energy at the outlet.
3. The vertical distributions of the streamwise velocity and the TKE are exported to the "mean_u_velocity_profiles.txt" and "tke_profiles.txt" text files in the working directory.

To run the simulation with a different method or Ks+ value, follow these steps:
1. Download the complete "tutorials/case_A" or "tutorials/cases_E" directory
2. Open the "ABL_adaptive_roughness_case_A.c" or "ABL_adaptive_roughness_case_E.c" file in a suitable text editor and change the value of the "ks_plus_array" variable in line 66 for case A or line 43 for case E, then save the file
3. To use Methods 1, 2 or 3, open Ansys Fluent and read the "case_setup_M1.jou", "case_setup_M2.jou" or "case_setup_M3.jou" journal files, respectively, into Fluent to initialise the case
4. Read the "run_case.jou" journal file to run the initialised case
5. The vertical distributions of the streamwise velocity and the TKE are exported to the "mean_u_velocity_profiles.txt" and "tke_profiles.txt" text files in the working directory.
