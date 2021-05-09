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

##CODE FOR: running linear models on gut microbiota alpha-diversity
# 3 metrics of alpha-diversity: 
## Chao1 Richness, Shannon diversity, simpson's index

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load mothur file and sample metadata         
################################################################################
#load alpha-diversity output by mothur (https://mothur.org/wiki/miseq_sop/)
#step 1: sub.sample(shared=00_OTU_table.txt, size=13340, persample=true)
#step 2: summary.single(shared=00_OTU_table.opti_mcc.txt, calc=nseqs-coverage-
#sobs-invsimpson, subsample=F)

alpha=read.table("data/00_OTU_alphadiversity.txt",sep="\t", header=T, 
                 stringsAsFactors = F);

load("data/02_sample_metadata_formatted.Rdata");


################################################################################
#             2. Append metadata to alpha diversity data        
################################################################################
am=inner_join(alpha[,,], meta[,,], by="Group");

#subset data to samples from adults; to samples from juveniles
ama=am[am$age_cat=="adult",];
amj=am[am$age_cat=="juvenile",];


################################################################################
#             3. Run linear mixed models
#     microbiota alpha diversity ~ bodysite + (1| hyenaID)
################################################################################
#adults
m1=lmer(log(chao)~bodysite + (1|hyenaID), data=ama);
m2=lmer(npshannon~bodysite + (1|hyenaID), data=ama);
m3=lmer(simpson~bodysite + (1|hyenaID), data=ama);
Anova(m1);Anova(m2);Anova(m3);

#juveniles
m1=lmer(log(chao)~bodysite + (1|hyenaID), data=amj);
m2=lmer(npshannon~bodysite + (1|hyenaID), data=amj);
m3=lmer(simpson~bodysite + (1|hyenaID), data=amj);
Anova(m1);Anova(m2);Anova(m3);


################################################################################
#             4. Run linear mixed models
#     microbiota alpha diversity ~ sex + (1| hyenaID) in JUVENILES
################################################################################
m1=lmer(log(chao)~sex+ (1|hyenaID), data=amj);
m2=lmer(npshannon~sex+ (1|hyenaID), data=amj);
m3=lmer(simpson~sex + (1|hyenaID), data=amj);
Anova(m1);Anova(m2);Anova(m3);


################################################################################
#             5. Run linear mixed models
#     microbiota alpha diversity ~ ageclass + (1| hyenaID) 
################################################################################
am2=am[am$clan=="Talek",];
m1=lmer(log(chao)~age_cat+ (1|hyenaID), data=am2);
m2=lmer(npshannon~age_cat+ (1|hyenaID), data=am2);
m3=lmer(simpson~age_cat + (1|hyenaID), data=am2);
Anova(m1);Anova(m2);Anova(m3);


################################################################################
#             6. make boxplots of alpha diversity ~ bodysite
#                                 ADULTS
################################################################################
mycol=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02");

#create plot - chao 1 richness
abox=ggplot(data=ama, 
                mapping=aes(x=bodysite,y=chao, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Chao 1 Richness")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(abox);

#create plot - shannon diversity
abox2=ggplot(data=ama, 
            mapping=aes(x=bodysite,y=npshannon, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(abox2);

#create plot - simpson's index
abox3=ggplot(data=ama, 
            mapping=aes(x=bodysite,y=1-simpson, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Simpson's index (1-D)")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(abox3);


################################################################################
#             6. make boxplots of alpha diversity ~ bodysite
#                             JUVENILES
################################################################################
#create plot - chao 1 richness
jbox=ggplot(data=amj, 
            mapping=aes(x=bodysite,y=chao, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Chao 1 Richness")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(jbox);

#create plot - shannon diversity
jbox2=ggplot(data=amj, 
             mapping=aes(x=bodysite,y=npshannon, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Shannon Diversity")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(jbox2);

#create plot - simpsons index
jbox3=ggplot(data=amj, 
             mapping=aes(x=bodysite,y=1-simpson, fill=bodysite))+
  geom_boxplot()+
  theme_bw()+ 
  labs(x = "",
       y = "Simpson's index (1-D)")+
  #ylim(3.4,7)+
  scale_fill_manual(values=mycol)+
  theme(legend.position="none", 
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 0.66),
        strip.text = element_text(size =11, face="bold")); plot(jbox3);


################################################################################
#             6. save plots
################################################################################
pbox=arrangeGrob(abox,abox2,abox3, jbox, jbox2,jbox3, nrow=2);

ggsave(filename="06_boxplots.pdf",
       device="pdf",path="./figures",
       plot=pbox,
       width=10,
       height=6,
       units="in",
       dpi=500);