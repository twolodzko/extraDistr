
## Update - bug fix

Removed erf, erfc, inv_erf from shared.h that are not used at this moment
and caused problems when compiling on Fedora and Solaris
(as noted in https://cran.r-project.org/web/checks/check_results_extraDistr.html). 