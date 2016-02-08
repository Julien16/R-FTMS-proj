#setwd("D:/Documents/R-FTMS-proj")
formula_list <- read.table("C_634_FORMULAE.dat",skip=15,header=T)
library(extremevalues)
library(plyr)
library(ggplot2)
source("classify_mol.R")
source("filter_1.R")
source("filter_2.R")
source("filter_3.R")

#here adapted to the molecular constraints: Hx-C100-O80-N2-S1_________####
mol_type <- classify_mol(formula_list)

attach(formula_list)
HC <- HIon/C
OC <- O/C
NC <- N/C
detach(formula_list)

#________________________Check high Intensity_____________####
L <- getOutliers(formula_list$Intensity, method="I", distribution="lognormal")
# outlierPlot(formula_list$Intensity, L, mode="qq")
#L$iRight
#L$limit['Right']->can be used for Intensity_max

#________________________N=0, mass even: N0_mass_even_____####
#set NtoC_max
NtoC_max <- 1
N0_mass_even <- ifelse(NtoC_max==0 , 
                      ifelse(mod(round(ExpMass+1, digits=0), 2)==0,
                             ifelse(N==0, 1, 0),
                             0),
                      1)

#________________________Nrule___________________________####
Nrule <- 1

#________________________filter 1_________________________####
filter1 <- filter_1(formula_list,
                    OtoC_min = 0, OtoC_max = 1, 
                    Intensity_Max = L$limit['Right'], 
                    Intensity_Min = 0, Intensity = formula_list$Intensity, 
                    ExpMass_Min = 140, ExpMass_Max = 800,
                    OC = formula_list$O/formula_list$C)

#________________________filter 2_________________________####
#set DBEtoC_min & DBEtoC_max & OplusNtoC_max
filter2 <- filter_2(formula_list, DBEtoC_min = 0, DBEtoC_max = 5, OplusNtoC_max = 3)

#________________________filter 3_________________________####
#set AI min & max
filter3 <- filter_3(formula_list, AI_max = 1, AI_min = -20)

attach(formula_list)

#________________________H/C charge condition:HCcc________####
#set charge... I write this way to obtain a vector with the right length. maybe not the right way to do...
charge <- ifelse(formula_list$ExpMass!=0,-1,0)
HCcc <- ifelse(Nrule*filter1*filter2*filter3==1, 
               ifelse(charge==(-1),((HIon+1)/C),((HIon-1)/C)),
               0)

#________________________H/C filtered: HCf________________####
#set HC_min & HC_max
HC_min <- 0
HC_max <- 2.5  
HCf <- ifelse(HCcc<HC_max, 
              ifelse(HCcc>HC_min,HCcc,0), 
              0)

#________________________O/C filtered: OCf________________####
OCf <- ifelse(Nrule*filter1*filter2*filter3==1,OC,0)

#________________________N/C filtered: NCf________________####
NCf <- ifelse(Nrule*filter1*filter2*filter3==1,NC,0) 

#________________________DBE______________________________####
dbe <- ifelse(HCf>0, 
              ifelse(Nrule*filter1*filter2*filter3==1,(1+.5*(2*C-HIon+N-1)),0), 
              0)

#________________________AI_______________________________####
AI <- ifelse((C-O-N-S)>0, 
             ifelse((1+C-O-0.5*(HIon+1)-S)>0, 
                    ifelse(HCf>0, 
                           ifelse(Nrule*filter1*filter2*filter3==1,(1+C-O-.5*(HIon+1)-S)/(C-O-N-S), 
                                  0), 
                           0), 
                    0), 
             0)

#________________________AI_mod_______________________________####
AI_mod <- ifelse((C-.5*O-N-S)>0, 
                 ifelse((1+C-0.5*O-.5*(HIon+1)-S)>0, 
                        ifelse(HCf>0, 
                               ifelse(Nrule*filter1*filter2*filter3==1,(1+C-.5*O-.5*(HIon+1)-S)/(C-.5*O-N-S), 
                                      0), 
                               0), 
                        0), 
                 0)

#________________________Xc_______________________________####
#set m & n values
m <- 1
n <- 1
Xc <- ifelse((3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))<=0,0,
             ifelse(dbe<=(m*O+n*S),0,(3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))))

#________________________Xc_mod_______________________________####
#set m & n values
m <- .5
n <- .5
Xc_mod <- ifelse((3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))<=0,0, 
                 ifelse(dbe<=(m*O+n*S),0,(3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))))

#________________________KMD______________________________####
nommass <- ifelse(HCcc==0,0, 
                  ifelse(Nrule*filter1*filter2*filter3==1, 
                         floor(formula_list$ExpMass), 
                         0))
kmd <- ifelse(HCcc==0,0, 
              ifelse(Nrule*filter1*filter2*filter3==1, formula_list$ExpMass-nommass, 
                     0))

#________________________KMD CH2__________________________####
km_CH2 <- (formula_list$ExpMass+1.007825-.000549)*14/14.01565
nommass_CH2 <- ifelse(Nrule*filter1*filter2*filter3==1,
                      ifelse(HCcc==0,0,ceiling(km_CH2)),0)
kmd_CH2 <- ifelse(Nrule*filter1*filter2*filter3==1,
                  ifelse(HCcc==0,0,(nommass_CH2-km_CH2)),0)

#________________________KMD COO__________________________####
km_COO <- (formula_list$ExpMass+1.007825-.000549)*44/43.989829
nommass_COO <- ifelse(Nrule*filter1*filter2*filter3==1,
                      ifelse(HCcc==0,0,floor(km_COO)),0)
kmd_COO <- ifelse(Nrule*filter1*filter2*filter3==1,
                  ifelse(HCcc==0,0,(nommass_COO-km_COO)),0)

#________________________compile in dataframe_____________####
dfVK_raw <- data.frame(formula_list, mol_type, 
                       HCf, OCf, NCf, 
                       dbe, AI, AI_mod, Xc, Xc_mod,
                       nommass_CH2, kmd_CH2,
                       nommass_COO, kmd_COO)
dfVK_filter <- subset(subset(subset(subset(dfVK_raw,mol_type!="NA"), 
                                      OCf!=0), 
                               HCf!=0), 
                        NCf<1)
detach(formula_list)

#____________________________Graph_1:Mass spectra_____________________####
 
ggplot(formula_list, aes(x=ExpMass, y=Intensity)) +
  geom_segment(aes(xend=ExpMass), yend=0, colour="royalblue1", size=.1) +
  scale_x_continuous(name="m/z", limits=c(149, 1000)) +
  scale_y_continuous(name="Intensity (counts)", limits=c(0, max(formula_list$Intensity)+1), expand=c(0, 0)) +
  theme_classic() +
  theme(axis.line=element_line(colour="black", size=.2),
        axis.ticks=element_line(size=.1),
        axis.text.x=element_text(size=6, colour="black"),
        axis.text.y=element_text(size=6, colour="black"),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        panel.background=element_blank(),
        plot.background=element_blank()) -> spectra_plot

#ggasve(filename="spectra_plot.eps", plot=spectra_plot, width=8.4, height=5.6872, units="cm", dpi=75)
#ggasve(filename="spectra_plot.pdf", plot=spectra_plot, width=8.4, height=5.6872, units="cm", dpi=75)


#____________________________Graph_1*:Mass spectra filtered_____________________####

ggplot(dfVK_filter, aes(x=ExpMass, y=Intensity)) + 
  geom_segment(aes(xend=ExpMass), yend=0, colour="royalblue1", size=.1) +
  scale_x_continuous(name="m/z", limits=c(149, 1000)) +
  scale_y_continuous(name="Intensity (counts)", limits=c(0, max(dfVK_filter$Intensity)+1), expand=c(0, 0)) +
  theme_classic() +
  theme(axis.line=element_line(colour="black", size=.2),
        axis.ticks=element_line(size=.1),
        axis.text.x=element_text(size=6, colour="black"),
        axis.text.y=element_text(size=6, colour="black"),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        panel.background=element_blank(),
        plot.background=element_blank()) -> spectra_filtered_plot

#ggasve(filename="spectra_filtered_plot.eps", plot=spectra_plot, width=8.4, height=5.6872, units="cm", dpi=75)
#ggasve(filename="spectra_filtered_plot.pdf", plot=spectra_plot, width=8.4, height=5.6872, units="cm", dpi=75)


#____________________________Graph_2:counts of formula________________####
ggplot(dfVK_filter, aes(x=mol_type)) +
geom_bar(stat="count", width=.6, fill=c("deepskyblue1", "orange", "olivedrab4", "red"), colour="black", size=.1) +
         #probleme pour gerer les couleurs lorsqu'une ou plusieurs classes de "class_mol_type" et "_mol_type" sont vides, je recois:
         #Error: Aesthetics must be either length 1 or the same as the data (3): colour, fill
         #Error: Insufficient values in manual scale. 4 needed but only 3 provided.
         scale_x_discrete(limits=c("CHO", "CHON", "CHOS", "CHONS"), name="") +
scale_y_continuous(name="counts", expand = c(0, 0)) +
                   #, limits=c(0, 2500), breaks=c(0, 500, 1000, 1500, 2000, 2500)
theme_classic() +
theme(axis.line.x=element_blank(),
      panel.background=element_blank(),
      plot.background=element_blank(),
      panel.grid.major = element_line(colour = "grey", size=.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_line(colour = "azure2", size=.05),
      panel.grid.minor.x = element_blank(),
      axis.line.y=element_line(size=.1, colour="black"),
      axis.ticks.y=element_line(size=.1),
      axis.ticks.x=element_blank(),
      axis.text.x=element_text(face="bold", colour=c("deepskyblue3", "orange3", "olivedrab4", "red")))->  count_plot

#ggasve(filename="count_plot.eps", plot=count_plot, width=8.4, height=8.4, units="cm", dpi=75)
#ggasve(filename="count_plot.pdf", plot=count_plot, width=8.4, height=8.4, units="cm", dpi=75)

#____________________________Graph_3: van Krevelen diagram: H/C vs O/C__####
ggplot(dfVK_filter, aes(x=OCf, y=HCf, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=3) +
  scale_x_continuous(name="O/C", limits=c(0, 1), breaks=c(0, .2, .4, .6, .8, 1.0), expand = c(0, 0)) +
  scale_y_continuous(name="H/C", limits=c(0, 2.5), breaks=c(0, .5, 1.0, 1.5, 2.0, 2.5), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=.25/1) -> vk_plot

#ggasve(filename="vk_plot.eps", plot=vk_plot, width=8.4, height=5.8, units="cm", dpi=75)
#ggasve(filename="vk_plot.pdf", plot=vk_plot, width=8.4, height=5.8, units="cm", dpi=75)

#____________________________Graph_4: van Krevelen diagram: H/C vs m/z__####
ggplot(dfVK_filter, aes(x=ExpMass, y=HCf, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=3) +
  scale_x_continuous(name="m/z", limits=c(140, 810), breaks=c(200, 300, 400, 500, 600, 700, 800), expand = c(0, 0)) +
  scale_y_continuous(name="H/C", limits=c(0, 2.5), breaks=c(0, .5, 1.0, 1.5, 2.0, 2.5), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=0.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=335/2) -> vk_mz_plot

#ggasve(filename="vk_mz_plot.eps", plot=vk_mz_plot, width=8.4, height=5.8, units="cm", dpi=75)
#ggasve(filename="vk_mz_plot.pdf", plot=vk_mz_plot, width=8.4, height=5.8, units="cm", dpi=75)

#____________________________Graph_5: AI vs C__________________________####
ggplot(dfVK_filter, aes(x=C, y=AI, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=5) +
  scale_x_continuous(name="C atoms", limits=c(0, 40)) +
  scale_y_continuous(name="AI", limits=c(0, 1)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=80/1) -> AI_plot

#ggasve(filename="AI_plot.eps", plot=AI_plot, width=8.4, height=12, units="cm", dpi=75)
#ggasve(filename="AI_plot.pdf", plot=AI_plot, width=8.4, height=12, units="cm", dpi=75)

#____________________________Graph_5_mod: AI_mod vs C__________________####
ggplot(dfVK_filter, aes(x=C, y=AI_mod, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=5) +
  scale_x_continuous(name="C atoms", limits=c(0, 40)) +
  scale_y_continuous(name="AI mod", limits=c(0, 1)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=80/1) -> AI_mod_plot

#ggasve(filename="AI_mod_plot.eps", plot=AI_mod_plot, width=8.4, height=12, units="cm", dpi=75)
#ggasve(filename="AI_mod_plot.pdf", plot=AI_mod_plot, width=8.4, height=12, units="cm", dpi=75)

#____________________________Graph_6: Xc vs C__________________________####
ggplot(dfVK_filter, aes(x=C, y=Xc, fill=mol_type, size=Intensity)) +
scale_size_area(max_size=5) +
scale_x_continuous(name="C atoms", limits=c(0, 40), expand = c(0, 0)) +
scale_y_continuous(name="Xc", limits=c(2.4, 3), expand = c(0, 0)) +
geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
theme_classic() +
theme(legend.position="none", 
      panel.background=element_blank(),
      plot.background=element_blank(),
      axis.line=element_line(size=.1, colour="black"), 
      axis.text.x=element_text(size=6), 
      axis.text.y=element_text(size=6), 
      axis.ticks=element_line(size=.1), 
      axis.title.x=element_text(size=8, vjust=-.5), 
      axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=80/.75) -> Xc_plot

#ggasve(filename="Xc_plot.eps", plot=Xc_plot, width=8.4, height=12, units="cm", dpi=75)
#ggasve(filename="Xc_plot.pdf", plot=Xc_plot, width=8.4, height=12, units="cm", dpi=75)
#____________________________Graph_6_mod: Xc_mod vs C__________________####
ggplot(dfVK_filter, aes(x=C, y=Xc_mod, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=5) +
  scale_x_continuous(name="C atoms", limits=c(0, 40), expand = c(0, 0)) +
  scale_y_continuous(name="Xc mod", limits=c(2.4, 3), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=8, vjust=2.5)) +
  coord_fixed(ratio=80/.75) -> Xc_mod_plot

#ggasve(filename="Xc_mod_plot.eps", plot=Xc_mod_plot, width=8.4, height=12, units="cm", dpi=75)
#ggasve(filename="Xc_mod_plot.pdf", plot=Xc_mod_plot, width=8.4, height=12, units="cm", dpi=75)

#____________________________Graph_7: KMD CH2__________________________####
ggplot(dfVK_filter, aes(x=nommass_CH2, y=kmd_CH2, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=2) +
  scale_x_continuous(name="m/z", limits=c(140, 810), breaks=c(200, 300, 400, 500, 600, 700, 800), expand = c(0, 0)) +
  scale_y_continuous(name="Kendrick mass defect (CH2)", limits=c(0, .8), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=6, vjust=2.5)) +
  coord_fixed(ratio=335/.64) -> KMD_CH2_plot

#ggasve(filename="KMD_CH2_plot.eps", plot=KMD_CH2_plot, width=8.4, height=5.6872, units="cm", dpi=75)
#ggasve(filename="KMD_CH2_plot.pdf", plot=KMD_CH2_plot, width=8.4, height=5.6872, units="cm", dpi=75)
#____________________________Graph_8: KMD COO__________________________####
ggplot(dfVK_filter, aes(x=nommass_COO, y=kmd_COO, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=2) +
  scale_x_continuous(name="m/z", limits=c(140, 810), breaks=c(200, 300, 400, 500, 600, 700, 800), expand = c(0, 0)) +
  scale_y_continuous(name="Kendrick mass defect (COO)", limits=c(-.7, 0), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=6, vjust=2.5)) +
  coord_fixed(ratio=335/.56) -> KMD_COO_plot

#ggasve(filename="KMD_COO_plot.eps", plot=KMD_COO_plot, width=8.4, height=5.6872, units="cm", dpi=75)
#ggasve(filename="KMD_COO_plot.pdf", plot=KMD_COO_plot, width=8.4, height=5.6872, units="cm", dpi=75)


#____________________________Graph_7: KMD CH2 zoom__________________________####
ggplot(dfVK_filter, aes(x=nommass_CH2, y=kmd_CH2, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=2) +
  scale_x_continuous(name="m/z", limits=c(300, 500), breaks=c(300, 400, 500), expand = c(0, 0)) +
  scale_y_continuous(name="Kendrick mass defect (CH2)", limits=c(.28, .33), breaks=c(.28, .29, .3, .31, .32, .33), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=6, vjust=2.5)) +
  coord_fixed(ratio=200/.05) -> KMD_CH2_plot

#ggasve(filename="KMD_CH2_plot_zoom.eps", plot=KMD_CH2_plot, width=8.4, height=10, units="cm", dpi=75)
#ggasve(filename="KMD_CH2_plot_zoom.pdf", plot=KMD_CH2_plot, width=8.4, height=10, units="cm", dpi=75)
#____________________________Graph_8: KMD COO zoom__________________________####
ggplot(dfVK_filter, aes(x=nommass_COO, y=kmd_COO, fill=mol_type, size=Intensity)) +
  scale_size_area(max_size=2) +
  scale_x_continuous(name="m/z", limits=c(300, 500), breaks=c(300, 400, 500), expand = c(0, 0)) +
  scale_y_continuous(name="Kendrick mass defect (COO)", limits=c(-.2, -.15), breaks=c(-.2, -0.19, -.18, -.17, -.16, -.15), expand = c(0, 0)) +
  geom_point(aes(fill=mol_type, size=Intensity), shape=21, colour='black', stroke=.1) +
  scale_fill_manual(values=c("deepskyblue1", "orange", "red", "olivedrab4")) +
  theme_classic() +
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.background=element_blank(),
        axis.line=element_line(size=.1, colour="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=.1),
        axis.title.x=element_text(size=8, vjust=-.5),
        axis.title.y=element_text(size=6, vjust=2.5)) +
  coord_fixed(ratio=200/.05) -> KMD_COO_plot

#ggasve(filename="KMD_COO_plot_zoom.eps", plot=KMD_COO_plot, width=8.4, height=10, units="cm", dpi=75)
#ggasve(filename="KMD_COO_plot_zoom.pdf", plot=KMD_COO_plot, width=8.4, height=10, units="cm", dpi=75)





