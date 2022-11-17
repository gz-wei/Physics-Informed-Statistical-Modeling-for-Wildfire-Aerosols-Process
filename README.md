# Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Process

## Environment
The codes have been tested on the following environment.
```
R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.5.1

(required packages)
bnstruct_1.0.12 
expm_0.999-6
ggplot2_3.3.6  
Rfast_2.0.6
SpatialVx_1.0 
tidyr_1.2.1
tidyverse_1.3.2
Matrix_1.5-1 
matrixStats_0.62.0 
MCMCpack_1.6-3  
ncdf4_1.19
```

## Data 
The complete AOD data can be downloaded in the following websites, and the sample AOD data (16.nc, 17.nc) are given at the [data](https://github.com/gz-wei/Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Process/tree/main/data) folder for testing the codes. 

- GOES-16 AOD data: [https://console.cloud.google.com/storage/browser/gcp-public-data-goes-16/ABI-L2-AODC;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false.](https://console.cloud.google.com/storage/browser/gcp-public-data-goes-16/ABI-L2-AODC;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false)
- GOES-17 AOD data: [https://console.cloud.google.com/storage/browser/gcp-public-data-goes-17/ABI-L2-AODC?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false](https://console.cloud.google.com/storage/browser/gcp-public-data-goes-17/ABI-L2-AODC?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false)


## Usage  
- Codes are completed in R project. After cloning the repository, this project can be openned with the ```Project.Rproj``` file. 
- The ```util.R``` file contains the required functions, and it does not give any results. 
- The ```illustrative_figures.R```, ```simulation_study.R```, and ```application.R``` files need to be run line by line and they can be used seperatelly. 

## Reproducibility
The figures and tables of the papaer can be reproduced in the following way.
- [illustrative_figures.R](https://github.com/gz-wei/Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Process/blob/main/illustrative_figures.R): Figures 2, 5, 6.
- [simulation_study.R](https://github.com/gz-wei/Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Process/blob/main/simulation_study.R): Figures 13, 14, 15 and Table 1. 
- [application.R](https://github.com/gz-wei/Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Process/blob/main/application.R): Figures 8, 9, 10, 11, 12.
