##############################
##    Hardwood content impacts the parasitoid community associated with Eastern spruce budworm (Lepidoptera: Tortricidae)
##
##  Christopher J. Greyson-Gaito, Sarah J. Dolson, Glen Forbes, Rosanna Lamb,
##  Wayne E. MacKinnon, Kevin S. McCann, M. Alex Smith, Eldon S. Eveleigh
##
##    R script for statistical analysis and figure plotting
##
##############################


# Packages ---------------------------------------------------------------
rm(list=ls())

theme_simple <- function () { 
  theme_grey() %+replace% 
    theme(
      axis.line=element_line(colour="black"),
      panel.grid.minor=element_blank(), 
      panel.grid.major=element_blank(),
      panel.background=element_blank(), 
      axis.title=element_text(size=28,face="bold"),
      axis.text.x=element_text(size=24, colour="Black"),
      axis.text.y=element_text(size=24, colour="Black"),
      axis.ticks.length=unit(0.5,"cm"),
      axis.ticks=element_line(size=0.5, colour="Black"),
      panel.border=element_rect(fill=FALSE,size=0.5),
      legend.title=element_text(size=15),
      legend.key=element_blank()
    )
  
}

library("tidyverse"); theme_set(theme_simple())
library("readxl")
library("picante")
library(viridis)
library(vegan)
library(goeveg)
library(permute)
library(nlme)

# Data Input and Cleaning -----------------------------------------------------------

## 2016 malaise caught parasitoids plus 2015 reared parasitoids (that were DNA barcoded)
ASSBW_ASBNAmetadata <- read_csv("data/ASSBW&ASBNA.csv") %>%
  select(ProcessID, SampleID, HWGrad, Plot, Method, CONTIG, BIN)

ASSBW_ASBNAphy <- read.tree("data/ASSBW&ASBNA.nwk")

## 1980s reared parasitoids (that were DNA barcoded)
ASSPP_ASSPQmetadata <- read_csv("data/ASSPP&ASSPQ.csv") %>%
  separate(CollectionNotes, c("Plot", "Composition"), sep="_", remove=FALSE)

ASSPP_ASSPQphy <- read.tree("data/ASSPP&ASSPQ.nwk")

## 1980s Malaise Samples - Stable Isotope
malcatfol <- read_csv("data/SIdataChrisGreysonGaito2017.csv")%>%
  mutate(DummyIdent=Identifier, DummyIdent=gsub(" A| B| C","",DummyIdent))%>% #remove A B C to prep for find average values of repeats of a few samples
  group_by(DummyIdent)%>% # set up averaging by group by a dummy identifier
  summarise(d13C=mean(cald13C), percentC=mean(calpercentC), d15N=mean(cald15N), percentN=mean(calpercentN), CNRatio=mean(CNratio))%>% #find the average values of the repeats of a few samples and return the same values of the other samples that did not have SI repeats
  rename(Identifier=DummyIdent) %>% #rename the DummyIdent column to Identifier
  mutate(Year = case_when(
    !grepl("GR", Identifier) & !grepl("17", Identifier) ~ 2017,
    grepl("17", Identifier) ~ 2017,
    grepl("82", Identifier) ~ 1982,
    grepl("83", Identifier) ~ 1983,
    grepl("86", Identifier) ~ 1986,
    grepl("87", Identifier) ~ 1987
  ))%>%
  mutate(Plot=case_when(
    grepl("01", Identifier) & !grepl("GR", Identifier) ~ 1,
    grepl("02", Identifier) & !grepl("GR", Identifier) ~ 2,
    grepl("03", Identifier) & !grepl("GR", Identifier) ~ 3,
    grepl("04", Identifier) & !grepl("GR", Identifier) ~ 4,
    grepl("05", Identifier) & !grepl("GR", Identifier) ~ 5,
    grepl("06", Identifier) & !grepl("GR", Identifier) ~ 6,
    grepl("07", Identifier) & !grepl("GR", Identifier) ~ 7,
    grepl("08", Identifier) & !grepl("GR", Identifier) ~ 8,
    grepl("09", Identifier) & !grepl("GR", Identifier) ~ 9
  ))%>%
  mutate(PlotType=as.factor(case_when(
    Plot %in% c(1,2,3) ~ "BFBF",
    Plot %in% c(4,5,6) ~ "BFMX",
    Plot %in% c(7,8,9) ~ "HWBF"
  )))%>%
  mutate(SamplingPeriod=as.factor(case_when(
    grepl("JUN",Identifier) ~ "SBWOUT",
    grepl("JULY",Identifier) | grepl("AUG",Identifier)  ~ "SBWGONE"
  )))%>%
  mutate(FunctionalTrophicPosition=as.factor(case_when(
    grepl("BFP",Identifier) | grepl("HWP",Identifier) | grepl("WRAUG",Identifier) | grepl("RHOAUG",Identifier) ~ "Foliage",
    grepl("ALT",Identifier) | grepl("SBW", Identifier) ~ "Caterpillar",
    grepl("GR", Identifier) ~ "Parasitoid"
  )))%>%
  mutate(TreeType=as.factor(case_when(
    grepl("BF", Identifier) ~ "BF",
    grepl("HW", Identifier) ~ "HW",
    grepl("WRAUG",Identifier) | grepl("RHOAUG",Identifier) ~ "Shrub"
  )))%>%
  mutate(SBWALT=as.factor(case_when(
    grepl("SBW", Identifier) ~ "SBW",
    grepl("ALT", Identifier) ~ "ALT"
  )))%>%
  mutate(ParasitoidGroup=as.factor(case_when( 
    grepl("GR01", Identifier) | grepl("GR02", Identifier) ~ "1_2",
    grepl("GR03", Identifier) ~ "3",
    grepl("GR04", Identifier) ~ "4",
    grepl("GR05", Identifier) ~ "5"
  )))

# Parasitoid community differences along a hardwood gradient -----------------------------

#Parasitoid community composition
nmdsdata<-ASSBW_ASBNAmetadata %>%
  group_by(HWGrad,Plot, BIN) %>%
  summarise(nind=length(BIN)) %>%
  ungroup() %>%
  mutate(HWGrad=as.factor(HWGrad))

nmdscommat<-nmdsdata%>%
  spread(BIN,nind)%>%
  ungroup()#creates a community matrix (but still including the factors) of the dataframe nmds

factors<-tibble(
  Plot=nmdscommat$Plot,
  HWGrad=nmdscommat$HWGrad,
)#creates a dataframe of the different factors that will be used in the nMDS and the PERMANOVA

nmdscommatfin <- nmdscommat %>%
  select(-c(HWGrad,Plot))#removes the factors from the dataframe to create a true community matrix

nmdscommatfin[is.na(nmdscommatfin)]<-0#replaces any NA values with 0. these NA values were created when the long format data was spread into a wide format dataframe and where taxa were not found for a certain Peak or Plot.

dimcheckMDS(nmdscommatfin, distance = "bray", k = 6, trymax = 30, autotransform=TRUE) #after exmining screeplot from 6 dimensions to 1 dimension, chosen to use 2 dimensions because additional dimensions provide small reductions in stress

set.seed(1223)
HWNMDS <- metaMDS(nmdscommatfin, distance = "bray", k = 2, trymax = 100, autotransform=TRUE)

stressplot(HWNMDS)  ##to test if distances as represented in ordination are correlated with actual distances

stress10<-c(0.1594808,0.169071,0.1457156,0.2536508,0.1457156,0.1989812,0.1605333,0.1605333,0.1750966,0.2383071)
sd(stress10)

#Create nMDS plot
HWdatascores<-as.data.frame(scores(HWNMDS))
HWdatascores$Plot<-factors$Plot
HWdatascores$HWGrad<-factors$HWGrad

# function for creating ellipses in nmds plot
#adapted from  http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipsenew <- function (x, scale = 1, npoints = 100) 
{
  cov <- cov.wt(cbind(x$NMDS1,x$NMDS2),wt=rep(1/length(x$NMDS1),length(x$NMDS1)))$cov 
  center <- c(mean(x$NMDS1),mean(x$NMDS2))
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  HWGradunique<-unique(x$HWGrad)
  tibble(
    NMDS1=t(center + scale * t(Circle %*% chol(cov)))[,1],
    NMDS2=t(center + scale * t(Circle %*% chol(cov)))[,2],
    HWGrad=HWGradunique
  )
}

df_ell.HWGrad<-HWdatascores%>%
  split(.$HWGrad)%>%
  map(veganCovEllipsenew)%>%
  bind_rows()

HWNMDSplot<-ggplot(HWdatascores)+
  geom_point(aes(NMDS1,NMDS2, colour=HWGrad),size=5)+
  geom_path(data=df_ell.HWGrad, aes(x=NMDS1, y=NMDS2, colour=HWGrad), size=1.5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.justification=c(1,0.99), 
        legend.position=c(1,0.99),
        legend.box = "horizontal")+
  coord_fixed(ratio = 1)+scale_color_manual(name="Dominant Tree\nType", breaks=c("BFBF","BFMX","HWBF"),labels=c("Balsam Fir","Mixed","Hardwood"), values=c("#39568CFF","#1F968BFF","#B8DE29FF"))

HWNMDSplot

ggsave("figs/nmdsHW.pdf",plot=HWNMDSplot,width=10,height=10) #Figure 1


# Testing the separation of the parasitoid communities between before, during, and after the budworm peak, using the adonis function for a permanova. - Guide to formulating the permutation structure and permutation test from the following webpage http://thebiobucket.blogspot.ca/2011/04/repeat-measure-adonis-lately-i-had-to.html#more
adonis(vegdist(nmdscommatfin)~HWGrad, factors, permutations=999)


## Phylogenetic clustering

### 2010s Malaise traps
ASSBWfinal<-ASSBW_ASBNAmetadata %>%
  filter(Method=="Malaise") %>%
  select(HWGrad,BIN)

ASSBWBINpresenceabsence<- ASSBWfinal %>%
  group_by(HWGrad,BIN) %>%
  summarise(abun=length(HWGrad)) %>%
  mutate(pa=ifelse(abun>0,1,0))%>%
  select(-abun)%>%
  spread(key=BIN, value=pa)%>%
  ungroup%>%
  select(-HWGrad)

ASSBWBINpresenceabsence[is.na(ASSBWBINpresenceabsence)]<-0

ASSBWBINPAmatrix<-as.matrix(ASSBWBINpresenceabsence)
row.names(ASSBWBINPAmatrix)<-c("BFBF","BFMX","HWBF")

#Picante

## Making the Phylogeny for Picante

# Making pruned phylogeny with the matrix version of our data
ASSBWprunedphy <- prune.sample(ASSBWBINPAmatrix, ASSBW_ASBNAphy)

# To view only 1 plot per output screen
par(mfrow = c(1, 1))

# This is to view the phylogeny that we will use for the picante analysis

plot(ASSBWprunedphy)


# To view 3 plots per output screen
par(mfrow = c(1, 3))

#Shows a red dot on the tips where a species is present at that site (where in community matrix it is >0)

for (i in row.names(ASSBWBINPAmatrix)) {
  plot(ASSBWprunedphy, show.tip.label = FALSE, main = i, cex.main = 3) 
  tiplabels(tip = which(ASSBWprunedphy$tip.label%in%names(which(ASSBWBINPAmatrix[i, ]  > 0))), pch = 18, cex = 1, col="red")} # Figure 3 A

## Species Richness & Phylogenetic Community Structure 

# Calculating phylogenetic diversity across our phylogeny
ASSBWpd.result <- pd(ASSBWBINPAmatrix, ASSBWprunedphy, include.root = FALSE)

# This will give you the site names and corresponding phylogenetic diversity (PD) and the species richness (SR)
ASSBWpd.result

# phydist object is measure of distance within the phylogeny
ASSBWphydist <- cophenetic(ASSBWprunedphy) 

# This is the test to determine if any sites are phylogenetically clustered or dispersed, as well as get a general species richness count
# The number of the runs is not fixed, but the standard for ses.mntd is ~999, If you need a trial then use fewer runs, ~100, so it runs faster 
# In order to randomize community matrix we run null models. Here we chose the null model "taxa.labels" because it is the default null model to randomize all tips across all sites
#full description with citations available at <http://picante.r-forge.r-project.org/picante-intro.pdf>
set.seed(42); ASSBWses.mntd.result <- ses.mntd(ASSBWBINPAmatrix, ASSBWphydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999) 

# This will give a table of results 
# ntaxa is the species richness per site
# mntd.obs.z is the standard effect size
# mntd.obs.p is the p-value; <0.05 is significantly clustered, >0.95 is significantly dispersed 
ASSBWses.mntd.result



## 1980s reared parasitoids

ASSPPPQfinal<-ASSPP_ASSPQmetadata %>%
  select(Plot,BIN)

ASPPPQBINpresenceabsence<- ASSPPPQfinal %>%
  group_by(Plot,BIN) %>%
  summarise(abun=length(Plot)) %>%
  mutate(pa=ifelse(abun>0,1,0))%>%
  select(-abun)%>%
  spread(key=BIN, value=pa)%>%
  ungroup%>%
  select(-Plot)

ASPPPQBINpresenceabsence[is.na(ASPPPQBINpresenceabsence)]<-0

ASSPPPQBINPAmatrix<-as.matrix(ASPPPQBINpresenceabsence)
row.names(ASSPPPQBINPAmatrix)<-c("Plot 1", "Plot 2", "Plot 3")

#Picante

## Making the Phylogeny for Picante

# Making pruned phylogeny with the matrix version of our data
ASSPPPQprunedphy <- prune.sample(ASSPPPQBINPAmatrix, ASSPP_ASSPQphy)

# To view only 1 plot per output screen
par(mfrow = c(1, 1))

# This is to view the phylogeny that we will use for the picante analysis

plot(ASSPPPQprunedphy)

plot(ASSPP_ASSPQphy)
# To view 3 plots per output screen
par(mfrow = c(1, 3))

#Shows a red dot on the tips where a species is present at that site (where in community matrix it is >0)
# We want no tip labels so the phylogeny is clean 

for (i in row.names(ASSPPPQBINPAmatrix)) {
  plot(ASSPPPQprunedphy, show.tip.label = FALSE, main = i, cex.main = 3) 
  tiplabels(tip = which(ASSPPPQprunedphy$tip.label%in%names(which(ASSPPPQBINPAmatrix[i, ]  > 0))), pch = 18, cex = 1, col="red")} # Figure 3 B

## Species Richness & Phylogenetic Community Structure 

# Calculating phylogenetic diversity across our phylogeny
ASSPPPQpd.result <- pd(ASSPPPQBINPAmatrix, ASSPPPQprunedphy, include.root = FALSE)

# This will give you the site names and corresponding phylogenetic diversity (PD) and the species richness (SR)
ASSPPPQpd.result

# phydist object is measure of distance within the phylogeny
ASSPPPQphydist <- cophenetic(ASSPPPQprunedphy) 

# This is the test to determine if any sites are phylogenetically clustered or dispersed, as well as get a general species richness count
# The number of the runs is not fixed, but the standard for ses.mntd is ~999, If you need a trial then use fewer runs, ~100, so it runs faster 
# In order to randomize community matrix we run null models. Here we chose the null model "taxa.labels" because it is the default null model to randomize all tips across all sites
#full description with citations available at <http://picante.r-forge.r-project.org/picante-intro.pdf>
set.seed(42); ASSPPPQses.mntd.result <- ses.mntd(ASSPPPQBINPAmatrix, ASSPPPQphydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999) 


# This will give a table of results 
# ntaxa is the species  richness per site
# mntd.obs.z is the standard effect size
# mntd.obs.p is the p-value; <0.05 is significantly clustered, >0.95 is significantly dispersed 
ASSPPPQses.mntd.result

# Parasitoid community balsam fir/hardwood trophic relationships --------

## 1980s Malaise Samples - Stable Isotope

### Tests to see difference in d15N between foliage, caterpillars and parasitoids

ggplot(data = malcatfol, aes(FunctionalTrophicPosition, d15N))+
  geom_jitter(aes(colour = TreeType),size=4)+
  theme(axis.line=element_line(colour="black"), panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), panel.background=element_blank(), 
        axis.title=element_text(size=28,face="bold"), axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        axis.text.x=element_text(size=24, colour="Black"), 
        axis.text.y=element_text(size=24, colour="Black"),axis.ticks.length=unit(0.5,"cm"),axis.ticks=element_line(size=0.5, colour="Black"),panel.border=element_rect(fill=FALSE,size=0.5)) #not a figure in the text


### Tests to see differences in d13C for foliage of balsam fir versus hardwoods
malcatfoltrees<-malcatfol%>%
  filter(FunctionalTrophicPosition=="Foliage", !TreeType=="Shrub") #shrubs were removed

ggplot(data = malcatfoltrees, aes(TreeType, d13C))+
  geom_point(aes(colour=SamplingPeriod),size=4)+
  theme(axis.line=element_line(colour="black"), panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(), panel.background=element_blank(), 
        axis.title=element_text(size=28,face="bold"), axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        axis.text.x=element_text(size=24, colour="Black"), 
        axis.text.y=element_text(size=24, colour="Black"),axis.ticks.length=unit(0.5,"cm"),axis.ticks=element_line(size=0.5, colour="Black"),panel.border=element_rect(fill=FALSE,size=0.5)) #not a figure in the text

t.test(d13C~TreeType,malcatfoltrees)

### Tests to see differences in d13C for caterpillars on balsam firs versus hardwoods

malcatfolcaterpillar <- malcatfol%>%
  filter(FunctionalTrophicPosition=="Caterpillar")

ggplot(malcatfolcaterpillar, aes(TreeType, d13C))+
  geom_point(aes(colour=SamplingPeriod,shape=SBWALT)) # not a figure in the text

t.test(d13C~TreeType,malcatfolcaterpillar)

# Comparison of d13C between years, time periods (May/June or July/August/September), and functional groups
carbonparagroup<-function(paragrp){
  para<-malcatfol %>%
    filter(FunctionalTrophicPosition=="Parasitoid")%>%
    filter(ParasitoidGroup==paragrp)%>%
    group_by(Year,SamplingPeriod)%>%
    summarise(meand13C=mean(d13C))%>%
    ungroup%>%
    mutate(ParasitoidGroup=paragrp)
  
  return(para)
} # This function takes a parasitoid group as a variable and then takes the main malcatfol dataset, filters for the parasitoid group specified, then calculates the mean d13C for each year and time period

carbonparagrp<-bind_rows(carbonparagroup("1_2"),carbonparagroup("3"),carbonparagroup("5")) %>% #NOTE group12 or 1_2 denotes previous grouping. grouping for these para in this manuscript is group 1
  mutate(ParaLabel=as.factor(ifelse(ParasitoidGroup=="1_2","one",
                                    ifelse(ParasitoidGroup=="3","two","three"))),smallyr=as.numeric(substr(Year,3,4)))

carbonparagrp$ParaLabel <- factor(carbonparagrp$ParaLabel, levels=c("one","two","three"))

group3_meand13C<-filter(carbonparagrp,ParaLabel=="three")%>%
  group_by(SamplingPeriod)%>%
  summarise(av13Cgrp3=mean(meand13C)) #-24.975 (balsam fir) #-28.025 (hardwood) for group 3 parasitoids (used for dotted lines in figure 1)

group3_meand13C[2,2]/group3_meand13C[1,2] #difference of 12.2% -  used in Results - Alternative hardwood-softwood parasitoids hypothesis

mean(malcatfoltrees$d13C)/mean(carbonparagrp$meand13C) #enriched by 16% (used in first sentence of statistical analysis for alternating hardwood-softwood parasitoids hypothesis)

carbonparagrpplot<-ggplot(carbonparagrp)+
  geom_hline(yintercept=-24.975, linetype="dashed")+
  geom_hline(yintercept=-28.025, linetype="dashed")+
  geom_point(aes(smallyr,meand13C, colour=SamplingPeriod),size=5)+
  facet_grid(.~ParaLabel)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=25),
        axis.title.y=element_text(hjust=0.5, vjust=1.5), panel.spacing = unit(1, "lines"),legend.text=element_text(size=14))+
  ylab(expression(delta*"13C"))+xlab("Year")+scale_color_viridis(name="budworm\nlarvae", breaks=c("SBWOUT","SBWGONE"),labels=c("May/June","July/August/September"), alpha = 1, begin = 0, end = 1, direction = 1, discrete = TRUE, option = "D")

ggsave("figs/carbonparagrp.pdf",carbonparagrpplot,width=10,height=4) # Figure 1


##Linear mixed effects models of year, time period, and d13C for all three parasitoid functional groups
#Data Exploration
#1.Outliers in the response variables (can't have outliers in explanatory variable of year).
dotchart(carbonparagrp$meand13C, main = "meand13C", group=carbonparagrp$SamplingPeriod)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(meand13C~Year, data=carbonparagrp, groups=SamplingPeriod)

#Creating Model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
carbonparagrpfulllinear<-lm(meand13C~Year*SamplingPeriod*ParaLabel,data=carbonparagrp)
carbonparagrp$resid<-NA
carbonparagrp$resid<-carbonparagrpfulllinear$residuals
#Test normality and homogeniety of residuals
plot(carbonparagrpfulllinear, add.smooth = FALSE, which = 1,col=carbonparagrp$SamplingPeriod)
hist(carbonparagrp$resid, xlab = "Residuals", main = "")
plot(carbonparagrp$Year,carbonparagrp$resid,xlab="Year",ylab="Residuals")
boxplot(carbonparagrp$resid~carbonparagrp$SamplingPeriod,xlab="SamplingPeriod",ylab="Residuals")
acf(residuals(carbonparagrpfulllinear), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
#Obviously violation of independence of x values (time) and maybe heterogeneity of variation (between sampling periods)

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
carbonparagrpfullgls<-gls(meand13C~Year*SamplingPeriod*ParaLabel,data=carbonparagrp, method="REML")
#3. Choose an appropriate variance structure depending on graphical analysis above.

carbonparagrpvarIdent<-gls(meand13C~Year*SamplingPeriod*ParaLabel,weights=varIdent(form=~1|SamplingPeriod),data=carbonparagrp, method="REML")


carbonparagrpvarIdentautocorr<-gls(meand13C~Year*SamplingPeriod*ParaLabel,weights=varIdent(form=~1|SamplingPeriod),correlation=corAR1(form=~Year|ParaLabel/SamplingPeriod),data=carbonparagrp, method="REML")

carbonparagrpautocorr<-gls(meand13C~Year*SamplingPeriod*ParaLabel,correlation=corAR1(form=~Year|ParaLabel/SamplingPeriod),data=carbonparagrp, method="REML")

anova(carbonparagrpfullgls,carbonparagrpvarIdent,carbonparagrpvarIdentautocorr,carbonparagrpautocorr)

#check assumptions graphically for model with autocorrelation
carbonparagrp$fmnormresid<-NA
carbonparagrp$fmnormresid<-resid(carbonparagrpvarIdent, type="normalized")
plot(carbonparagrpvarIdent, add.smooth = FALSE,col=carbonparagrp$SamplingPeriod)
hist(carbonparagrp$fmnormresid, xlab = "Normalised Residuals", main = "")
plot(carbonparagrp$Year,carbonparagrp$fmnormresid,xlab="Year",ylab="Normalised Residuals")
boxplot(carbonparagrp$fmnormresid~carbonparagrp$SamplingPeriod,xlab="SamplingPeriod",ylab="Normalised Residuals")
coplot(fmnormresid ~ Year | SamplingPeriod,
       ylab = "Normalised residuals", data = carbonparagrp)
acf(carbonparagrp$fmnormresid, na.action = na.pass,
    main = "Auto-correlation plot for normalised residuals")

#AIC lowest for model with varIdent and without autocorr (ACF shows little improvement from autocorr, boxplot little improvement from varIdent)
#decided to not include autocorr (but will include varIdent)

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#no random effects

#5. Compare the new gls model with the earlier results using AIC.


#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
carbonparagrpfullglsml<-gls(meand13C~Year*SamplingPeriod*ParaLabel,weights=varIdent(form=~1|SamplingPeriod),data=carbonparagrp, method="ML")

carbonparagrpfullglsdrop3interactionml<-update(carbonparagrpfullglsml,.~.-Year:SamplingPeriod:ParaLabel)

anova(carbonparagrpfullglsml,carbonparagrpfullglsdrop3interactionml) #drop three way interaction because not significant

carbonparagrpfullglsdropYRSPinteractionml<-update(carbonparagrpfullglsdrop3interactionml,.~.-Year:SamplingPeriod)
anova(carbonparagrpfullglsdrop3interactionml,carbonparagrpfullglsdropYRSPinteractionml) 

carbonparagrpfullglsdropYRPLinteractionml<-update(carbonparagrpfullglsdrop3interactionml,.~.-Year:ParaLabel)
anova(carbonparagrpfullglsdrop3interactionml,carbonparagrpfullglsdropYRPLinteractionml) 

carbonparagrpfullglsdropSPPLinteractionml<-update(carbonparagrpfullglsdrop3interactionml,.~.-SamplingPeriod:ParaLabel)
anova(carbonparagrpfullglsdrop3interactionml,carbonparagrpfullglsdropSPPLinteractionml) 

#drop Year:SamplingPeriod because not significant #do not drop Year:ParaLabel in this step because significant #do not drop SamplingPeriod:ParaLabel in this step because significant

carbonparagrpfullglsdropYRPL2interactionml<-update(carbonparagrpfullglsdropYRSPinteractionml,.~.-Year:ParaLabel)
anova(carbonparagrpfullglsdropYRSPinteractionml,carbonparagrpfullglsdropYRPL2interactionml) #do not drop Year:ParaLabel in this step because significant

carbonparagrpfullglsdropSPPL2interactionml<-update(carbonparagrpfullglsdropYRSPinteractionml,.~.-SamplingPeriod:ParaLabel)
anova(carbonparagrpfullglsdropYRSPinteractionml,carbonparagrpfullglsdropSPPL2interactionml) #do not drop SamplingPeriod:ParaLabel in this step because significant

#cannot drop single terms because of the significant interactions between Year:ParaLabel and SamplingPeriod:ParaLabel

carbonparagrpfinalmodelREML<-gls(meand13C~Year+SamplingPeriod+ParaLabel+Year:ParaLabel+SamplingPeriod:ParaLabel,weights=varIdent(form=~1|SamplingPeriod),data=carbonparagrp,method="REML")
summary(carbonparagrpfinalmodelREML)
anova(carbonparagrpfinalmodelREML)

gr12yeardiff<-carbonparagrp %>% #NOTE group12 or 1_2 denotes previous grouping. grouping for these para in this manuscript is group 1
  filter(ParasitoidGroup=="1_2") %>%
  group_by(Year) %>%
  summarise(Avd13C=mean(meand13C))

(((gr12yeardiff$Avd13C[4]/gr12yeardiff$Avd13C[1]) ^ (1/3)) - 1 ) * 100 # about 0.5% more negative each year for group 1 (used in Results for Alternating hardwood-softwood parasitoids hypothesis)

gr12seasondiff<-carbonparagrp %>% #NOTE group12 or 1_2 denotes previous grouping. grouping for these para in this manuscript is group 1
  filter(ParasitoidGroup=="1_2") %>%
  group_by(SamplingPeriod) %>%
  summarise(Avd13C=mean(meand13C))

(gr12seasondiff$Avd13C[1]/gr12seasondiff$Avd13C[2] - 1) *100 #2.4% more negative between Juliy/August/September and May/June (used in Results for Alternating hardwood-softwood parasitoids hypothesis)

gr3yeardiff<-carbonparagrp %>% #NOTE group3 or 3 denotes previous grouping. grouping for these para in this manuscript is group 2
  filter(ParasitoidGroup=="3") %>%
  group_by(Year) %>%
  summarise(Avd13C=mean(meand13C))

(1 - ((gr3yeardiff$Avd13C[4]/gr3yeardiff$Avd13C[1]) ^ (1/3))) * 100 # 1.6% less negative each year (used in Results for Alternating hardwood-softwood parasitoids hypothesis)



# Supporting Information --------------------------------------------------
maldipichabund <- read_csv("data/Malaise_data/maldipich_allyears_long.csv") %>%
  filter(!year %in% c(84,85,88)) %>%
  filter(!Species %in% c("Agathis_males", "Agathis_females", "Apanteles_other_sp_females", "Apanteles_other_sp_males","Charmon_extensor_males","Charmon_extensor_females","Choristoneura_fumiferana_females","Choristoneura_fumiferana_males", "Itoplectis_females","Itoplectis_males","Phaeogenes_females","Phaeogenes_males", "Ephialtes_ontario_females","Ephialtes_ontario_males"))

maldipichabund$Species <- gsub("_males","",maldipichabund$Species)
maldipichabund$Species <- gsub("_females","",maldipichabund$Species)
maldipichabund$Species <- gsub("aurifrons","auricaudata",maldipichabund$Species)
maldipichabund$Species <- gsub("tortricis","parva",maldipichabund$Species)
maldipichabund$Species <- gsub("setifacies","fumipennis",maldipichabund$Species)
maldipichabund$Species <- gsub("Pseudoperichaeta","Nilea",maldipichabund$Species)
maldipichabund$Species <- gsub("Pseudosarcophaga","Agria",maldipichabund$Species)
maldipichabund$Species <- gsub("Winthemia","Smidtia",maldipichabund$Species)
maldipichabund$Species <- gsub("psyte","pyste",maldipichabund$Species)
maldipichabund$Species <- gsub("Hemistermia","Hemisturmia",maldipichabund$Species)
maldipichabund$Species <- gsub("_"," ",maldipichabund$Species)

# 1980s Malaise Samples - Community Analysis

maldipichabundtotal <- maldipichabund %>%
  mutate(Sampling_Period=ifelse(collection_date<182,"SBWOUT","SBWGONE"), grp=case_when(
    Species %in% c("Apanteles fumiferanae", "Glypta fumiferanae", "Smidtia fumiferanae","Lypha fumipennis") ~ 1,
    Species %in% c("Actia interrupta", "Eumea caesar", "Sarcophaga aldrichi","Nilea erecta","Hemisturmia parva","Agria affinis","Compsilura concinnata", "Tachinomyia nigricans") ~ 2,
    Species %in% c("Meteorus trachynotus", "Ceromasia auricaudata", "Nemorilla pyste","Phryxe pecosensis","Madremyia saundersii","Tranosema rostrale") ~ 3
  ))%>%
  group_by(year, grp, Sampling_Period)%>%
  summarise(total=sum(Number, na.rm = TRUE))
  

maldipichabundsum <- maldipichabund %>%
  mutate(Sampling_Period=ifelse(collection_date<182,"SBWOUT","SBWGONE"))%>%
  group_by(year, Sampling_Period,Species)%>%
  summarise(abund=sum(Number, na.rm = TRUE)) %>%
  mutate(abund=ifelse(is.na(abund),0,abund),grp=case_when(
    Species %in% c("Apanteles fumiferanae", "Glypta fumiferanae", "Smidtia fumiferanae","Lypha fumipennis") ~ 1,
    Species %in% c("Actia interrupta", "Eumea caesar", "Sarcophaga aldrichi","Nilea erecta","Hemisturmia parva","Agria affinis","Compsilura concinnata", "Tachinomyia nigricans") ~ 2,
    Species %in% c("Meteorus trachynotus", "Ceromasia auricaudata", "Nemorilla pyste","Phryxe pecosensis","Madremyia saundersii","Tranosema rostrale") ~ 3
  ))%>%
  left_join(maldipichabundtotal, by = c("year", "grp", "Sampling_Period")) %>%
  mutate(prop = abund/total)

maldipichabundsum$Species<-gsub(" "," \n",maldipichabundsum$Species)

sbwout <- filter(maldipichabundsum, Sampling_Period=="SBWOUT")
sbwgone <- filter(maldipichabundsum, Sampling_Period=="SBWGONE")

maldipichabundgrp1bar<-ggplot()+
  geom_bar(data = filter(sbwout, grp==1), aes(year - (0.175+0.04),prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  geom_bar(data = filter(sbwgone, grp==1), aes(year + (0.175+0.04) ,prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=20),
        axis.title.y=element_text(hjust=0.5, vjust=1.5), panel.spacing = unit(1, "lines"),legend.text=element_text(size=14))+
  scale_x_continuous(breaks=c(82,83,84,85,86,87))+
  scale_fill_viridis(alpha = 1, begin = 0.23, end = 0.87, direction = 1, discrete = TRUE, option = "D")+
  ylab("Proportion")+xlab("Year")

ggsave("figs/maldipichpropgrp1.pdf",maldipichabundgrp1bar, width=6, height=4) #Figure S1

maldipichabundgrp2bar<-ggplot()+
  geom_bar(data = filter(sbwout, grp==2), aes(year - (0.175+0.04),prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  geom_bar(data = filter(sbwgone, grp==2), aes(year + (0.175+0.04) ,prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=20),
        axis.title.y=element_text(hjust=0.5, vjust=1.5), panel.spacing = unit(1, "lines"),legend.text=element_text(size=14))+
  scale_x_continuous(breaks=c(82,83,84,85,86,87))+
  scale_fill_viridis(alpha = 1, begin = 0.23, end = 0.87, direction = 1, discrete = TRUE, option = "D")+
  ylab("Proportion")+xlab("Year")

ggsave("figs/maldipichpropgrp2.pdf",maldipichabundgrp2bar, width=6.5, height=4.5) #Figure S2

maldipichabundgrp3bar<-ggplot()+
  geom_bar(data = filter(sbwout, grp==3), aes(year - (0.175+0.04),prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  geom_bar(data = filter(sbwgone, grp==3), aes(year + (0.175+0.04) ,prop,fill = Species), position="stack", stat="identity", width = 0.35, colour="black")+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=20),
        axis.title.y=element_text(hjust=0.5, vjust=1.5), panel.spacing = unit(1, "lines"),legend.text=element_text(size=14))+
  scale_x_continuous(breaks=c(82,83,84,85,86,87))+
  scale_fill_viridis(alpha = 1, begin = 0.23, end = 0.87, direction = 1, discrete = TRUE, option = "D")+
  ylab("Proportion")+xlab("Year")

ggsave("figs/maldipichpropgrp3.pdf",maldipichabundgrp3bar, width=6, height=4) #Figure S3

##### Extra useful information --------------

allyrsreared %>%
  group_by(Yr, HWGrad) %>%
  summarise(mnSBW = mean(NoSBWReared)) #Comparing number of budworm sampled each year and along hardwood gradient. Definitely differences but looks to be mostly random.

allyrsSBWaov <- aov(NoSBWReared~HWGrad*Yr,data=allyrsreared)

plot(allyrsSBWaov)

allyrsSBWaov
summary(allyrsSBWaov)

effectsize::cohens_f(allyrsSBWaov)
pwr.2way(a = 3, b = 3, alpha = 0.05, size.A = 9, size.B = 9, f.A = 0.33, f.B = 0.08) # Power.A = 0.202
