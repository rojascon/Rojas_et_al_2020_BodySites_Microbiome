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

##CODE FOR: generating stacked bar plots of gut microbiota composition
# at the bacterial phylum, family, and genus taxonomic level

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         formatted metadata table                 
################################################################################
tax=read.csv("data/00_OTU_taxonomy.txt", header=T, sep="\t");
load("data/01_OTU_table_filtered.Rdata");
load("data/02_sample_metadata_formatted.Rdata");

rownames(tax)=tax$OTU; tax$OTU=NULL;

#attach taxonomy to the OTU table 
otu_tax=merge(asvf, tax,by="row.names"); 

#make a vector of the samples you want for analysis
samples=as.character(meta$Group);

#remove chlorplast
otu_tax=otu_tax[otu_tax$Phylum!="Cyanobacteria_Chloroplast",]


################################################################################
#             2. Create Phylum level composition barplots                  
################################################################################
#select bacterial taxonomic rank
phylum=otu_tax[,which(names(otu_tax) 
                      %in% c(samples, "Phylum"))];
colnames(phylum)[ncol(phylum)]="taxa";

#calculate ASV relative abundances 
phylum=aggregate(.~taxa, phylum, sum);  
phylum[,-1] <- lapply(phylum[,-1], function(x) (x/sum(x))*100);
print(colSums(phylum[-1]));

#keep phyla >1% relative abundance across samples
phylum$AVG=rowMeans(phylum[,-1]);
phylum=phylum[phylum$AVG>1,];
phylum$AVG=NULL;

#denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(phylum[2:ncol(phylum)])); 
phylum=rbind(phylum, newrow); 
phylum$taxa=as.character(phylum$taxa);
phylum[nrow(phylum),1]="Other";

#melt data frame for ggplot
pbar<-reshape2::melt(phylum, id.vars="taxa",value.name = "abun");
colnames(pbar)[2]="Group";
pbar=merge(pbar, meta, by="Group");

#color-palette
phy_col=c("#8c6bb1","#6baed6","#2171b5","#74c476",
          "grey83","coral","lightgoldenrod");

#create plot- ADULTS
pbar1=pbar[pbar$age_cat=="adult",];
barphy=ggplot(data=pbar1, 
              mapping=aes(x=Group,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~bodysite, scales="free_x")+ 
  theme_bw()+ 
  labs(title="Adults",
       x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =10, face="bold"));plot(barphy);

#create plot- JUVENILES
pbar2=pbar[pbar$age_cat=="juvenile",];
barphy2=ggplot(data=pbar2, 
              mapping=aes(x=Group,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~bodysite, scales="free_x")+ 
  theme_bw()+ 
  labs(title="Juveniles",
       x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =10, face="bold"));plot(barphy2);


################################################################################
#             3. Save phylum barplots                 
################################################################################
#combine the two plots and save
pone=arrangeGrob(barphy,barphy2, nrow=2);

ggsave(filename="03_barplot_bacterial_phylum.pdf",
       device="pdf",path="./figures",
       plot=pone,
       width=14,
       height=10,
       units="in",
       dpi=500);


################################################################################
#             4. Create Family level composition barplots                 
################################################################################
#select bacterial taxonomic rank 
fam=otu_tax[,which(names(otu_tax) 
                   %in% c(samples, "Family"))];
colnames(fam)[ncol(fam)]="taxa";

#calculate ASV relative abundances 
fam=aggregate(.~taxa, fam, sum);  
fam[,-1] <- lapply(fam[,-1], function(x) (x/sum(x))*100);
print(colSums(fam[-1]));

#keep families >0.7% relative abundance across samples
fam$AVG=rowMeans(fam[,-1]);
fam=fam[fam$AVG>1.4,];
fam$AVG=NULL;

#denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(fam[2:ncol(fam)])); 
fam=rbind(fam, newrow); 
fam$taxa=as.character(fam$taxa);
fam[nrow(fam),1]="Other";

#shorten a few taxon names
fam$taxa[fam$taxa=="Bacteroidales _unclassified"]="Bacteroidales_unclass"
fam$taxa[fam$taxa=="Clostridia _unclassified"]="Clostridia_unclass"
fam$taxa[fam$taxa=="Clostridiales _unclassified"]="Clostridiales_unclass"
fam$taxa[fam$taxa=="Clostridiales_Incertae_Sedis_XI"]="Clostridiales_IncSed_XI"
fam$taxa[fam$taxa=="Firmicutes _unclassified"]="Firmicutes_unclass"
fam$taxa[fam$taxa=="unclassified _unclassified"]="unclass_unclass"

#melt data frame for ggplot
fbar<-reshape2::melt(fam, id.vars="taxa",value.name = "abun");
colnames(fbar)[2]="Group";
fbar=merge(fbar, meta, by="Group");

#color-palette
fam_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", 
          "#737373","#77CCCC","#117744","#88CCAA", "#777711",  "#DDDD77", 
          "grey","#774411", 
          "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","black")

#create plot- ADULTS
fbar1=fbar[fbar$age_cat=="adult",];
barfam=ggplot(data=fbar1, 
              mapping=aes(x=Group,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~bodysite, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial family",
       title="Adults")+
  scale_fill_manual(values=fam_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =10, face="bold"));plot(barfam);

#create plot- JUVENILES
fbar2=fbar[fbar$age_cat=="juvenile",];
barfam2=ggplot(data=fbar2, 
              mapping=aes(x=Group,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~bodysite, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial family",
       title="Juveniles")+
  scale_fill_manual(values=fam_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =10, face="bold"));plot(barfam2);

################################################################################
#             4. save Family level composition barplots                 
################################################################################
#combine the two plots and save
ptwo=arrangeGrob(barfam,barfam2, nrow=2);

ggsave(filename="03_barplot_bacterial_family.pdf",
       device="pdf",path="./figures",
       plot=ptwo,
       width=16,
       height=10,
       units="in",
       dpi=500);