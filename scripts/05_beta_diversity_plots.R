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

##CODE FOR: plotting PCoAs based on bray curtis and jaccard distances

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load beta diversity distances and sample metadata         
################################################################################
load("data/04_distances.Rdata"); 
load("data/02_sample_metadata_formatted.Rdata");

#abray= adults
#jbray= juveniles
#ajbray= adult and juvenile females


################################################################################
#             2. make a PCoA ordinations color-coded by bodysite
################################################################################
########## ADULTS ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(abray, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set color-palette for bodysites
mycol=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

##plot the PCoA
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Adults")+
  scale_fill_manual(values=mycol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa1)

########## JUVENILES ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(jbray, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa2=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Juveniles")+
  scale_fill_manual(values=mycol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa2);


################################################################################
#             3. save PCoA ordinations color-coded by bodysite
################################################################################
#combine the two plots and save
pone=arrangeGrob(pcoa1,pcoa2, nrow=1);

ggsave(filename="05_pcoa_bodysites.pdf",
       device="pdf",path="./figures",
       plot=pone,
       width=9.5,
       height=5,
       units="in",
       dpi=500);



################################################################################
#                     4. subset distance matrices 
#             (ears rectum prepuce gland) VS. (nasal, oral)
################################################################################
#subset adult and juvenile distance matrices
atemp=meta$Group[(meta$age_cat=="adult") & 
                   meta$bodysite %in% c("oral","nasal")]
jtemp=meta$Group[(meta$age_cat=="juvenile") & 
                   meta$bodysite %in% c("oral","nasal")]

a1=abray[!rownames(abray) %in% atemp,
                     !colnames(abray) %in% atemp];
a2=abray[rownames(abray) %in% atemp,
         colnames(abray) %in% atemp];

j1=jbray[!rownames(jbray) %in% jtemp,
         !colnames(jbray) %in% jtemp];
j2=jbray[rownames(jbray) %in% jtemp,
         colnames(jbray) %in% jtemp];

#set color-palettes for plots
fourcol=c("#1b9e77","#e7298a", "#66a61e", "#e6ab02");
twocol=c("#d95f02", "#7570b3");


################################################################################
#             3. make a PCoA ordinations color-coded by bodysite
#                  but nasal and oral separate from the rest
#                           ADULTS
################################################################################
########## ADULTS - ears, prepuce, rectum, gland ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(a1, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa3=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Adults")+
  scale_fill_manual(values=fourcol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa3)

########## ADULTS - nasal, oral ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(a2, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa4=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="")+
  scale_fill_manual(values=twocol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa4);


################################################################################
#             3. make a PCoA ordinations color-coded by bodysite
#                  but nasal and oral separate from the rest
#                           JUVENILES
################################################################################
########## JUVENILES - ears, prepuce, rectum, gland ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(j1, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa5=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Juveniles")+
  scale_fill_manual(values=fourcol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa5)

########## JUVENILES - nasal, oral ########
#calculate coordinates for PCoA
pcoa_dec=cmdscale(j2, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa6=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=bodysite),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="")+
  scale_fill_manual(values=twocol)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"))+
  guides(fill=guide_legend(nrow=1)); plot(pcoa6);


################################################################################
#             3. save PCoA ordinations  from above
#               (nasal and oral are separate from rest)
################################################################################
A=arrangeGrob(pcoa3,pcoa4, nrow=1);
B=arrangeGrob(pcoa5,pcoa6, nrow=1);

ggsave(filename="05_pcoa_bodysites_adults.pdf",
       device="pdf",path="./figures",
       plot=A,
       width=7.5,
       height=4,
       units="in",
       dpi=500);

ggsave(filename="05_pcoa_bodysites_juveniles.pdf",
       device="pdf",path="./figures",
       plot=B,
       width=7.5,
       height=4,
       units="in",
       dpi=500);


################################################################################
#             3. make a PCoA ordinations scent gland only
#                 juvenile females vs juvenile males
################################################################################
jt=meta$Group[(meta$age_cat=="juvenile") & 
                   meta$bodysite=="scentgland"];

jtd=jbray[rownames(jbray) %in% jt,
         colnames(jbray) %in% jt];

#get pcoa coordinates
pcoa_dec=cmdscale(jtd, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa7=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=sex),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Scent Gland")+
  scale_fill_manual(values=c("tan1","#a6d96a"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa7)

ggsave(filename="05_pcoa_scentgland.pdf",
       device="pdf",path="./figures",
       plot=pcoa7,
       width=4,
       height=4,
       units="in",
       dpi=500);


################################################################################
#             3. make a PCoA ordinations prepuce and rectum
#                 adult females vs juvenile females
################################################################################
dfp=meta$Group[(meta$clan=="Talek") & 
                meta$bodysite=="prepuce"];

dfr=meta$Group[(meta$clan=="Talek") & 
                meta$bodysite=="rectum"];


ajt=ajbray[rownames(ajbray) %in% dfp,
          colnames(ajbray) %in% dfp];

ajt2=ajbray[rownames(ajbray) %in% dfr,
          colnames(ajbray) %in% dfr];

#####  PREPUCE #######
#get pcoa coordinates
pcoa_dec=cmdscale(ajt, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa8=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=age_cat),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Prepuce")+
  scale_fill_manual(values=c("#8dd3c7","#bc80bd"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa8)

#####  PREPUCE #######
#get pcoa coordinates
pcoa_dec=cmdscale(ajt2, eig=TRUE);  
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group");
pcoa_met=merge(pcoa,meta,by="Group"); 

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

##plot the PCoA
pcoa9=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=age_cat),
             size = 3,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="Rectum")+
  scale_fill_manual(values=c("#8dd3c7","#bc80bd"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));plot(pcoa9)

##save plots
C=arrangeGrob(pcoa8,pcoa9, nrow=1);
ggsave(filename="05_pcoa_prepuce_rectum.pdf",
       device="pdf",path="./figures",
       plot=C,
       width=7.7,
       height=4,
       units="in",
       dpi=500);
