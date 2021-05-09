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

##CODE FOR: formatting metadata factors for analyses

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1.  Load sample metadata  and format factors           
################################################################################
#load metadata
meta=read.csv("data/00_sample_metadata.csv", stringsAsFactors = F);
meta$Group=as.character(meta$Group);
meta$age_cat<-factor(meta$age_cat, levels=c("adult","juvenile"));
meta$sex<-factor(meta$sex, levels=c("female","male"));
meta$bodysite<-factor(meta$bodysite);
meta=meta[order(meta$Group),];

#save metadata file
save(meta, file="data/02_sample_metadata_formatted.Rdata");
