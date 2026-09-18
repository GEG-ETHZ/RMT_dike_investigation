# RMT_dike_investigation
This repository contains all code used in the Master's thesis "The Potential of Radiomagnetotellurics for Detecting Dike Instabilities". It includes both original code developed for this thesis and modified code based on existing software or scripts.

The Folder mtcode contains all scripts from the [mtcode](https://github.com/darcycordell/mtcode) repository by `darcycordell` which were adapted for the use of RMT:
> Cordell, D., Lee, B., Unsworth, M.J., 2022.  
> mtcode: A repository of MATLAB scripts for magnetotelluric data  
> analysis, data editing, model building, and model viewing,  
> doi:10.5281/zenodo.6784201

## Workflow
The Workflow of the Project is as follows
1) Raw Data visualization with python files
2) Topography file creation with Python/TIF_to_DEM.py
3) Model creation with mtcode
4) Model editing with python files
5) Data editing with mtcode
6) Forward/Inverse runs with ModEM
7) Data and model visualization with MATLAB files
