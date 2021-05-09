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

##CODE FOR: removing contaminant OTUs from dataset using sequence data from 
#extraction kit controls

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load OTU Table, OTU Taxonomy, and sample metadata                 
################################################################################

#load OTU table output by mothur SOP tutorial 
#https://mothur.org/wiki/miseq_sop/
#file is too large, so download from here: 
# https://drive.google.com/file/d/1p4AjWhk4ZZNe4wIrgEMS5MCEXPu4ooem/view?usp=sharing
#then add it to the data directory of this R project
otu=read.table("data/00_OTU_table.txt",sep="\t", header=T, stringsAsFactors = F);

#load sample metadata 
meta=read.csv("data/00_sample_metadata.csv",header=T, stringsAsFactors = F);

#load OTU taxonomy file
#ASVs classified as Eukarya, Mitochondria, Chloroplast, and Unknown were removed
tax=read.csv("data/00_OTU_taxonomy.txt", sep="\t",header=T); 


################################################################################
#             2.clean up OTU table                 
################################################################################
#remove unecessary columns
otu=otu[,c(2, 4:ncol(otu))]

#remove samples that did not amplify well [had <100 reads] and unwanted samples
bad=c("100","135", "165", "184", "208", "219","168","96");
meta=meta[!meta$Group %in% bad,];
otu=otu[!otu$Group %in% bad,];


################################################################################
#             2.  Remove OTUs >1% abundance in control samples 
################################################################################
#subset df to only control samples
blanks=c("BEE","BEEE2","BL","BLK1","BS","BSS","H20");
controls=otu[otu$Group %in% blanks,];
rownames(controls)=controls$Group; controls$Group=NULL;

#retain OTUs that have >1% relative abundance across control samples
controls=apply(controls, 1, function(i) (i/sum(i))*100);
controls=as.data.frame((controls))
colSums(controls);
controls$Abun=rowMeans(controls); 
controls=controls[controls$Abun>1,];

#view the taxonomic classification of these ASVs
View(tax[tax$OTU
         %in% row.names(controls),]);

#we wont remove OTU6
controls=controls[2:15,];


################################################################################
#             4.  Remove the contaminant OTUs from OTU table
#################################################################################

#remove control samples from OTU table
otu=otu[!otu$Group %in% blanks,];

#remove contaminant OTUs
asvf=otu; rownames(asvf)=asvf$Group; asvf$Group=NULL;
asvf=as.data.frame(t(asvf))
asvf=asvf[!row.names(asvf) %in% 
            row.names(controls),];

##also remove ASVs classified as chloroplast, mitochondria, and unknown
asvf=asvf[row.names(asvf) %in% tax$OTU,];

#remove singleton and doubleton ASVs
temp=asvf; 
temp=(temp>0)*1; 
temp=as.data.frame(temp);
temp$sm=rowSums(temp);
temp=temp[temp$sm>0,]; 
asvf=asvf[rownames(asvf) %in% rownames(temp),]
save(asvf, file="data/01_OTU_table_filtered_withsindoub.Rdata");

temp=temp[temp$sm>2,];
asvf=asvf[rownames(asvf) %in% rownames(temp),]
save(asvf, file="data/01_OTU_table_filtered.Rdata");

gg=asvf[ , order(names(asvf))]
