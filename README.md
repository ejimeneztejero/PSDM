# PSDM: Prestack Depth Migration for Marine Data

- This software is a specialized geophysical toolset designed for the Prestack Depth Migration of marine seismic data.
- It transform seismic reflections into depth images of the subsurface with three different migration algorithms.
- It is tested for MCS data but it is also prepared for WAS data (still to be tested).

## Features

The software offers three distinct migration algorithms to suit different computational and imaging needs:

- Pure Eikonal Calculation: Primarily used for debugging and quality control of velocity models.
- Kirchhoff Migration: A classic, efficient approach for structural imaging.
- Reverse Time Migration (RTM): High-fidelity imaging for complex geological structures.

HPC is implemented and paralelized with the number of shotgathers.

> [!IMPORTANT]
> Field data pre-processing and final PSDM image post-processing must be performed by the user using external tools.

## Prerequisites

Before installing, ensure you have the following dependencies:
- Fortran Compilers (gfortran/mpif90/mpirun).
- Seismic Unix (SU): The software handles data in `.su` format.
- Make: For build automation.

## Installation
- Navigate to the source directory and compile:
- cd PSDM/src
- make

## Inputs
- The software requires MCS seismic marine field data in Seismic Unix (.su) format.
- Navigation: An ASCII file containing the geometry/navigation information.
- Parameter File: A configuration file (read during execution) where you specify which algorithm to use (eikonal, kirchoff or RTM), velocity models, grid dimensions, and algorithm-specific parameters. An input file "parfile" is included in this folder as an example.

## Run 
- Example (run with 40 MPI processes):
- mpirun -np 40 psdm_run parfile
  
## Documentation and testing
- Manual: A detailed PDF manual is currently under construction and will be added to the repository soon.
- Test Data: Sample datasets for testing the code will be hosted on the Zenodo database (link forthcoming).

## Author & Acknowledge
- Author: Clara Estela Jiménez Tejero.
- Institution: Barcelona Center for Subsurface Imaging, ICM-CSIC.
- While this software is freely available, we would be grateful if you could notify us of its use. Please send a brief email to ejimenez@icm.csic.es to acknowledge your usage, which helps us justify continued development and support.
