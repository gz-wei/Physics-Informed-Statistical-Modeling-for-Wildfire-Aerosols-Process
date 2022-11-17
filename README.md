# Physics-Informed-Statistical-Modeling-for-Wildfire-Aerosols-Propagation

## Required environment 
software: R version 3.6.2 \
packages: tidyverse, expm, Rfast, MCMCpack

R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.5.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base 

other attached packages:
 [1] bnstruct_1.0.12       igraph_1.3.4          bitops_1.0-7          SpatialVx_1.0         turboEM_2021.1       
 [6] quantreg_5.94         SparseM_1.81          numDeriv_2016.8-1.1   doParallel_1.0.17     iterators_1.0.14     
[11] foreach_1.5.2         smatr_3.4-8           smoothie_1.0-3        fields_14.1           viridis_0.6.2        
[16] viridisLite_0.4.1     spam_2.9-1            spatstat_2.3-4        spatstat.linnet_2.3-2 spatstat.core_2.4-4  
[21] rpart_4.1.16          nlme_3.1-157          spatstat.random_2.2-0 spatstat.geom_2.4-0   spatstat.data_2.2-0  
[26] matrixStats_0.62.0    Rfast_2.0.6           RcppZiggurat_0.1.6    Rcpp_1.0.9            MCMCpack_1.6-3       
[31] MASS_7.3-57           coda_0.19-4           expm_0.999-6          Matrix_1.5-1          ncdf4_1.19           
[36] forcats_0.5.2         stringr_1.4.1         dplyr_1.0.10          purrr_0.3.4           readr_2.1.2          
[41] tidyr_1.2.1           tibble_3.1.8          tidyverse_1.3.2       ggplot2_3.3.6        

loaded via a namespace (and not attached):
 [1] googledrive_2.0.0     colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2        rprojroot_2.0.3      
 [6] fs_1.5.2              rstudioapi_0.13       farver_2.1.1          MatrixModels_0.5-1    fansi_1.0.3          
[11] lubridate_1.8.0       xml2_1.3.3            codetools_0.2-18      splines_4.2.1         polyclip_1.10-0      
[16] jsonlite_1.8.0        mcmc_0.9-7            broom_1.0.1           dbplyr_2.2.1          spatstat.sparse_2.1-1
[21] compiler_4.2.1        httr_1.4.4            backports_1.4.1       assertthat_0.2.1      gargle_1.2.1         
[26] cli_3.4.1             tools_4.2.1           dotCall64_1.0-1       gtable_0.3.1          glue_1.6.2           
[31] maps_3.4.0            cellranger_1.1.0      vctrs_0.4.1           fastcluster_1.2.3     rvest_1.0.3          
[36] lifecycle_1.0.2       googlesheets4_1.0.1   goftest_1.2-3         distillery_1.2-1      scales_1.2.1         
[41] hms_1.1.2             spatstat.utils_2.3-1  gridExtra_2.3         stringi_1.7.8         boot_1.3-28          
[46] rlang_1.0.6           pkgconfig_2.0.3       lattice_0.20-45       tensor_1.5            labeling_0.4.2       
[51] tidyselect_1.1.2      here_1.0.1            magrittr_2.0.3        R6_2.5.1              generics_0.1.3       
[56] DBI_1.1.3             pillar_1.8.1          haven_2.5.1           withr_2.5.0           mgcv_1.8-40          
[61] survival_3.3-1        abind_1.4-5           modelr_0.1.9          crayon_1.5.1          utf8_1.2.2           
[66] waveslim_1.8.4        tzdb_0.3.0            grid_4.2.1            readxl_1.4.1          CircStats_0.2-6      
[71] reprex_2.0.2          digest_0.6.29         munsell_0.5.0      


## User guidence 
Codes are completed in R project. After cloning the repository, one can open the project with “Project.Rproj” file. The codes can be successfully run in the "main.R" file. One can obtain the results after a short time running. In "main.R", the highest retained frequency N is set to 10 considering the computational time. However, N = 10 gives relatively poor results as shown in the "results" folder. I also upload the case when N = 20, which shows a much better performance. 
