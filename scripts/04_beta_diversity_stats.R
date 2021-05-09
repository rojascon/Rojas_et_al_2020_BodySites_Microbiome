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

##CODE FOR: calculating distance matrices for PERMANOVA beta-diversity analyses
#Bray-Curtis (bray), and Jaccard (jacc)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered OTU table and sample metadata         
################################################################################
load("data/01_OTU_table_filtered.Rdata");
load("data/02_sample_metadata_formatted.Rdata");


###############################################################################
#             2. Generate the 2 types of distance matrices
################################################################################
#transpose ASV table so ASVs are columns
asvf=t(asvf);

###BRAY-CURTIS distance
bray<-apply(asvf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(asvf>0)*1;
print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

#convert both to data frames
bray=as.data.frame(as.matrix(bray.dist));
jacc=as.data.frame(as.matrix(jac.dist));


################################################################################
#         3. subset distance matrices to only adults; to only juveniles
################################################################################
#only retain juvenile samples
js=meta[meta$age_cat=="juvenile",];
jbray=bray[rownames(bray) %in% js$Group,
          colnames(bray) %in% js$Group];
jjacc=jacc[rownames(jacc) %in% js$Group,
          colnames(jacc) %in% js$Group];

#only retain adult samples
as=meta[meta$age_cat=="adult",];
abray=bray[rownames(bray) %in% as$Group,
           colnames(bray) %in% as$Group];
ajacc=jacc[rownames(jacc) %in% as$Group,
           colnames(jacc) %in% as$Group];

################################################################################
#                     4. conduct PERMANOVA -- MODEL 1
#                    microbiota similarity ~ bodysite
#                   microbiota dispersion ~ bodysite (PERMDISP only)
################################################################################
#have 4 distance matrices so 4 tests to run
mydist=list(jbray, jjacc, abray, ajacc);
names=c("Bray-Curtis", "Jaccard", "Bray-Curtis", "Jaccard");
met=c("bray","jaccard","bray","jaccard");
mysam=list(js, js, as,as);

#run for loop to conduct the 4 PERMANOVA tests
#1- juvenile bray, 2- juvenile jaccard, 3-adult bray, 4- adult jaccard
for(i in 1:4)
{
  print(paste("PERMANOVA test, across bodysites, using:", names[i]));
  print(adonis(mydist[[i]]~     
                 bodysite+
                 hyenaID,
               data=mysam[[i]],
               method = met[i],
               permutations = 999));
};

#PERMDISP (heterogeneity analysis) - JUVENILES
bdisper<-with(js, betadisper(as.dist(jbray), bodysite));  #jbray or jjacc
anova(bdisper);
TukeyHSD(bdisper);

#PERMDISP (heterogeneity analysis) - ADULTS
bdisper<-with(as, betadisper(as.dist(abray), bodysite)); #abray or ajacc
anova(bdisper);
TukeyHSD(bdisper);


################################################################################
#                     5. conduct PERMANOVA -- MODEL 2
#              microbiota similarity ~ sex [juveniles only]
#               separate test for each body site
################################################################################
#have 6 bodysites, so 6 tests to run, create for loop
bsite=levels(meta$bodysite);

#run for loop to conduct the 6 PERMANOVA tests
for(i in 1:6)
{
  A=js[js$bodysite==bsite[i],];
  B=jbray[rownames(jbray) %in% A$Group,
              colnames(jbray) %in% A$Group];
  print(paste("PERMANOVA test,sex effects, using Bray-Curtis:", bsite[i]));
  print(adonis(B~sex,   
               data=A,
               method = bray,
               permutations = 999));
};

##for jaccard, repeat the same code except changes lines 113-114: replace "jbray" 
#with "jjacc" and in line 118- also change bray to jaccard


################################################################################
#                     6. conduct PERMANOVA -- MODEL 3
#                     microbiota similarity ~ age_class
#               separate test for each body site
################################################################################
#remove hyenas from that are not from the same social group
#only want Talek hyenas adult females and juvenile females

aj=meta[meta$clan=="Talek",]; aj=aj[aj$sex=="female",];
ajbray=bray[rownames(bray) %in% aj$Group,
           colnames(bray) %in% aj$Group];
ajjacc=jacc[rownames(jacc) %in% aj$Group,
           colnames(jacc) %in% aj$Group];

#have 6 bodysites, so 6 tests to run, create for loop
bsite=levels(meta$bodysite);

#run for loop to conduct the 6 PERMANOVA tests
for(i in 1:6)
{
  A=aj[aj$bodysite==bsite[i],];
  B=ajbray[rownames(ajbray) %in% A$Group,
          colnames(ajbray) %in% A$Group];
  print(paste("PERMANOVA test,age effects, using Bray-Curtis:", bsite[i]));
  print(adonis(B~age_cat,   
               data=A,
               method = bray,
               permutations = 999));
};

##for jaccard, repeat the same code except changes lines 147-148: replace "ajbray" 
#with "ajjacc" and in line 152 change bray to jaccard


################################################################################
#                     7. conduct mantel test
#                     microbiota similarity ~ social rank [juveniles]
################################################################################
#retrieve social rank data
mr=data.frame(js[, (colnames(js) %in% c("Group","stan.soc.rank"))]);
rownames(mr)=mr$Group; mr$Group=NULL;
mr[mr==0]=0.0000000001;
dist.mat2=vegdist(mr, method="bray");
dist.mat2b=vegdist(mr, method="jaccard");

#run mantel test on bray curtis distances and jaccard distances
mantel(jbray,dist.mat2, method="spear", permutations=999);
mantel(jjacc,dist.mat2b, method="spear", permutations=999);

################################################################################
#                     8.save distance matrices for plotting later
################################################################################
#only Bray-Curtis distances
dist=save(jbray, abray, ajbray, file="data/04_distances.Rdata");


