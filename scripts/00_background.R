#################################################################################
#
#               Microbiome variation across bodysites in hyenas
#                      
#              Rojas et al 2020.Body-site specific variation reflect sex and 
#                       age-class among wild spotted hyenas
#
#                           By: Connie Rojas
#                           Created: 20 Jan 2019
#                           Last updated: 8 May 2021
#
################################################################################

##CODE FOR: configuring R workspace and printing R version and package versions
#for reader

################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE) ;

#load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "lme4","lmtest","multcomp","grid","phyloseq","gridExtra",
               "Biostrings","QsRutils","pheatmap");


################################################################################
#             2. Communicate the R version and package versions to reader                 
################################################################################
print("This code was developed with R version 3.6.2");

print("The packages used and their versions were: pheatmap_1.0.12 | gridExtra_2.3|
QsRutils_0.1.4| Biostrings_2.54.0| multcomp_1.4-15| lmtest_0.9-38| 
lme4_1.1-26| ggplot2_3.3.3| vegan_2.5-7| reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| 
MASS_7.3-53| car_3.0-10");

