# Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Propagation

## Environments
The codes have been testing on the following environment
```
R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.5.1

(packages required)
ggplot2_3.3.6  
tidyr_1.2.1
expm_0.999-6
MCMCpack_1.6-3  
Rfast_2.0.6
matrixStats_0.62.0 
Matrix_1.5-1 
tidyverse_1.3.2
SpatialVx_1.0 
bnstruct_1.0.12  
ncdf4_1.19
```

## User guidence 
Codes are completed in R project. After cloning the repository, one can open the project with “Project.Rproj” file. The codes can be successfully run in the "main.R" file. One can obtain the results after a short time running. In "main.R", the highest retained frequency N is set to 10 considering the computational time. However, N = 10 gives relatively poor results as shown in the "results" folder. I also upload the case when N = 20, which shows a much better performance. 
