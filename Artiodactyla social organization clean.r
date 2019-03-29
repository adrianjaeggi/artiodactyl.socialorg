## load packages
library(ape)
library(wesanderson)
library(brms)
library(rethinking)
library(readxl)
library(paleotree)
library(phytools)
library(plyr)
library(car)
library(MCMCglmm)
library(mice)
library(phytools)

## read data file
d<- as.data.frame(read_excel("Artiodacytla_Dataset_1Dec_Final_LDHMM.xlsx", sheet=1, skip=3)) # assuming working directory is set to the folder containing this file

# remove empty columns
d<- d[,-c(55:56)]

# fix continent misspellings
d$Continent[as.factor(d$Continent)=="africa"]<- "Africa"	
d$Continent[as.factor(d$Continent)=="Europe (Introduced)"]<- "Europe"	
d$Continent[as.factor(d$Continent)=="North America (Introduced)"]<- "North America"	
d$Continent[as.factor(d$Continent)=="Oceania (Introduced)"]<- "Oceania"	
d$Continent[as.factor(d$Continent)=="Norther America"]<- "North America"	
# that covers all populations - summary(as.factor(d[is.na(d[,"Common Name"]),"Continent"]))
# no need to recode species-summaries for now
#d$Continent[as.factor(d$Continent)=="Africa; Australia; Europe; Asia"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="Australia; Europe; Asia; North America; South America"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="Europe; Asia"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="Europe; Asia; North America"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="North America; Europe; Asia"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="North America; South America"]<- "Cosmopolitan"	
#d$Continent[as.factor(d$Continent)=="South America; North America"]<- "Cosmopolitan"	
	
	
# exclude species with missing data and species summaries	
d2<- d[d$no.studies>0,]
# d2<- d2[d2$no.SO>0,]


## read phylogeny - mammal supertree
mamm<- read.tree("mammaltree.tre") # assuming working directory is set to the folder containing this file
plot(mamm)

## check species names match and trim trees to dataset
setdiff(d$Genus_species, mamm$tip.label) #    [1] "Beatragus_hunteri"   "Capricornis_crispus" "Nanger_granti"       "Oryx_beisa"          NA                    "Lama_glama_guanicoe"
# [7] "Rucervus_duvaucelii" "Rusa_unicolor"       "Moschus_cupreus"     "Moschus_leucogaster"

## see iucn.org for synonyms or subspecies names
# Beatragus_hunteri = Damaliscus_hunteri
# Capricornis_crispus = Naemorhedus_crispus
# Nanger_granti = Gazella_granti
# Oryx_beisa = missing -> replace with dammah or leucoryx (both equally related to gazella)
# Lama_glama_guanicoe = Lama_glama
# Rucervus_duvaucelii = Cervus_duvaucelii
# Rusa_unicolor = Cervus_unicolor
# Moschus_leucogaster = Moschus_cupreus = Moschus_chrysogaster
d2$Genus_species[d2$Genus_species=="Beatragus_hunteri"]<- "Damaliscus_hunteri"
d2$Genus_species[d2$Genus_species=="Capricornis_crispus"]<- "Naemorhedus_crispus"
d2$Genus_species[d2$Genus_species=="Nanger_granti"]<- "Gazella_granti"
d2$Genus_species[d2$Genus_species=="Eudorcas_thomsonii"]<- "Gazella_thomsonii"
d2$Genus_species[d2$Genus_species=="Oryx_beisa"]<- "Oryx_dammah"
d2$Genus_species[d2$Genus_species=="Lama_glama_guanicoe"]<- "Lama_glama"
d2$Genus_species[d2$Genus_species=="Rucervus_duvaucelii"]<- "Cervus_duvaucelii"
d2$Genus_species[d2$Genus_species=="Rusa_unicolor"]<- "Cervus_unicolor"
d2$Genus_species[d2$Genus_species=="Moschus_leucogaster"]<- "Moschus_chrysogaster"
d2$Genus_species[d2$Genus_species=="Moschus_cupreus"]<- "Moschus_chrysogaster"
setdiff(d2$Genus_species, mamm$tip.label) # "Eudorcas_thomsonii" NA      
mamm.prun<- drop.tip(mamm, setdiff(mamm$tip.label,d2$Genus_species))
plot(mamm.prun, type="fan", cex=0.5)
is.ultrametric(mamm.prun) # TRUE
is.rooted(mamm.prun) # TRUE
is.binary.tree(mamm.prun) # FALSE
mamm.prun.di <- multi2di(mamm.prun) # this forces the tree to be binary (as required for analysis)
is.binary.tree(mamm.prun.di) # TRUE
mamm.prun.di.min<- minBranchLength(mamm.prun.di, 1e-10) # makes the artificially introduced branches so small that it doesn't matter
mamm.prun.di.min$node.label<- NULL

inv.phylo.ma<- inverseA(mamm.prun.di.min, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo.ma$Ainv)
rownames(A) <- rownames(inv.phylo.ma$Ainv)



		### Ancestral state estimation for predictors ###
		#################################################
	

bm<- log(d.spec[,"F_BM (kg)"])
names(bm)<- d.spec$Genus_species

bm<- bm[which(!is.na(bm))]

# reprune tree
mamm.prun3<- drop.tip(mamm, setdiff(mamm$tip.label,names(bm)))
plot(mamm.prun3, type="fan", cex=0.5)
is.ultrametric(mamm.prun3)
is.rooted(mamm.prun3)
is.binary.tree(mamm.prun3) # FALSE
mamm.prun3.di <- multi2di(mamm.prun3)
is.binary.tree(mamm.prun3.di) # TRUE
mamm.prun3.di.min<- minBranchLength(mamm.prun3.di, 1e-10)
mamm.prun3.di.min$node.label<- NULL

# fastAnc
ace.bm.fastAnc = fastAnc(mamm.prun3.di.min, bm, CI = TRUE)
ace.bm.fastAnc
# --> at root (node 79): 4.187368, 95%CI = 1.778210 6.596526
# ace
ace.bm.ace = ace(bm, mamm.prun3.di.min, type="continuous")
ace.bm.ace
# --> at root (node 79): 4.187368, 95%CI = 1.755541 6.619196
# virtually the same
# --> center body size on ace.bm.ace$ace[1] 

## sexual dimorphism M:F_BM
dimo<- d.spec[,"M:F_BM"]
names(dimo)<- d.spec$Genus_species

dimo<- dimo[which(!is.na(dimo))]

# reprun tree
mamm.prun4<- drop.tip(mamm, setdiff(mamm$tip.label,names(dimo)))
plot(mamm.prun4, type="fan", cex=0.5)
is.ultrametric(mamm.prun4)
is.rooted(mamm.prun4)
is.binary.tree(mamm.prun4) # FALSE
mamm.prun4.di <- multi2di(mamm.prun4)
is.binary.tree(mamm.prun4.di) # TRUE
mamm.prun4.di.min<- minBranchLength(mamm.prun4.di, 1e-10)
mamm.prun4.di.min$node.label<- NULL

# fastAnc
ace.dimo.fastAnc = fastAnc(mamm.prun4.di.min, dimo, CI = TRUE)
ace.dimo.fastAnc
# --> at root (node 79): 1.214002, 95%CI = 0.537744 1.890261
# ace
ace.dimo.ace = ace(dimo, mamm.prun4.di.min, type="continuous")
ace.dimo.ace
# --> at root (node 79): 1.214002, 95%CI = 0.5313545 1.896650
# virtually the same
# --> center sexual dimorphism on ace.dimo.ace$ace[1] 


## number habitats
nohab<- as.numeric(d.spec[,"no.habitats"])
names(nohab)<- d.spec$Genus_species

nohab<- nohab[which(!is.na(nohab))]

# have to drop one Moschus_chrysogaster
nohab<- nohab[-92]

# reprun tree
mamm.prun5<- drop.tip(mamm, setdiff(mamm$tip.label,names(nohab)))
plot(mamm.prun5, type="fan", cex=0.5)
is.ultrametric(mamm.prun5)
is.rooted(mamm.prun5)
is.binary.tree(mamm.prun5) # FALSE
mamm.prun5.di <- multi2di(mamm.prun5)
is.binary.tree(mamm.prun5.di) # TRUE
mamm.prun5.di.min<- minBranchLength(mamm.prun5.di, 1e-10)
mamm.prun5.di.min$node.label<- NULL

# fastAnc
ace.nohab.fastAnc = fastAnc(mamm.prun5.di.min, nohab, CI = TRUE)
ace.nohab.fastAnc
# --> at root (node 98): 1.673256, 95%CI = -0.674268 4.020780
# ace
ace.nohab.ace = ace(nohab, mamm.prun5.di.min, type="continuous")
ace.nohab.ace
# virtually the same
# --> center no.habitats on ace.nohab.ace$ace[1], or just 1 


## breeding seasonality
seas<- d.spec[,"Seaonal/Non-seasonal Breeding"]
names(seas)<- d.spec$Genus_species

seas<- seas[which(!is.na(seas))]

# reprun tree
mamm.prun6<- drop.tip(mamm, setdiff(mamm$tip.label,names(seas)))
plot(mamm.prun6, type="fan", cex=0.5)
is.ultrametric(mamm.prun6)
is.rooted(mamm.prun6)
is.binary.tree(mamm.prun6) # FALSE
mamm.prun6.di <- multi2di(mamm.prun6)
is.binary.tree(mamm.prun6.di) # TRUE
mamm.prun6.di.min<- minBranchLength(mamm.prun6.di, 1e-10)
mamm.prun6.di.min$node.label<- NULL


# compare ER, ARD
ace.seas.ER = ace(seas, mamm.prun6.di.min, type="discrete", model="ER")
ace.seas.ARD = ace(seas, mamm.prun6.di.min, type="discrete", model="ARD")

anova(ace.seas.ER, ace.seas.ARD) # --> ARD preferred
#Likelihood Ratio Test Table
#  Log lik. Df Df change Resid. Dev Pr(>|Chi|)  
#1  -51.079  1                                  
#2  -48.872  2         1     4.4139    0.03565 *

ace.seas.ARD
#Scaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):
#Non-seasonal     Seasonal 
#   0.5017214    0.4982786
## --> stick to non-seasonal


		### phylogenetic mixed-effects models ###
		#########################################
	
## data prep

## fill in body size data
for(i in 1:nrow(d2)){
	if(is.na(d2[i,"Seaonal/Non-seasonal Breeding"])) d2[i,"Seaonal/Non-seasonal Breeding"]<- d2[i-1,"Seaonal/Non-seasonal Breeding"]
	if(is.na(d2[i,"M:F_BM"])) d2[i,"M:F_BM"]<- d2[i-1,"M:F_BM"]
	if(is.na(d2[i,"F_BM (kg)"])) d2[i,"F_BM (kg)"]<- d2[i-1,"F_BM (kg)"]
	if(is.na(d2[i,"Male_BM (kg)"])) d2[i,"Male_BM (kg)"]<- d2[i-1,"Male_BM (kg)"]
	if(is.na(d2[i,"Avg_BM (kg)"])) d2[i,"Avg_BM (kg)"]<- d2[i-1,"Avg_BM (kg)"]
	if(is.na(d2[i,"Dimorphic/Monomorphic"])) d2[i,"Dimorphic/Monomorphic"]<- d2[i-1,"Dimorphic/Monomorphic"]
	}

## exclude species with missing data and summaries of species
d3<- d2[d2$Continent!="overall",]

# check DV
summary(as.factor(d3[,"Social state - for R"]))
colnames(d3)[22]<- "Social.state.for.R"
d3$Social.state.for.R[which(d3$Social.state.for.R=="NA")]<- NA
	
# clean habitat --> if multiple, code as multiple
d3$Habitat[as.factor(d3$Habitat)=="Artificial (Terrestrial); Forest/Woodland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Artificial (Terrestrial); Native Grassland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Artificial (Terrestrial); Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Desert; Native Grassland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Desert; Shrubland; Wetland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Inland Rocky Areas"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Inland Rocky Areas; Native Grassland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Native Grassland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Native Grassland; Savanna"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Native Grassland; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Native Grassland; Wetland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Savanna"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woodland; Wetland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Inland Rocky Areas; Native Grassland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Inland Rocky Areas; Native Grassland; Savanna"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Inland Rocky Areas; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Native Grassland; Savanna"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Native Grassland; Savanna; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Native Grassland; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="savanna"]<- "Savanna"	
d3$Habitat[as.factor(d3$Habitat)=="Savanna; Shrubland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Shrubland; Wetland"]<- "Multiple"	
d3$Habitat[as.factor(d3$Habitat)=="Forest/Woosland"]<- "Forest/Woodland"	
d3$Habitat[as.factor(d3$Habitat)=="inland Rocky Areas"]<- "Inland Rocky Areas"	

# standardize predictors
d3$dimorph.s<- (d3[,"M:F_BM"]-as.numeric(ace.dimo.ace$ace[1]))/sd(d3[,"M:F_BM"], na.rm=TRUE)
d3$logBM.s<- (log(d3[,"F_BM (kg)"])-as.numeric(ace.bm.ace$ace[1]))/sd(log(d3[,"F_BM (kg)"]), na.rm=TRUE)
d3$n.studies.c<- as.numeric(d3$no.studies)-1
d3$n.habitats.c<- as.numeric(d3$no.habitats)-1 

d3$phylo<- d3$Genus_species
d3$Social.state.for.R<- relevel(as.factor(d3$Social.state.for.R), ref="Solitary")

## exclude pops with no social state
d3<- subset(d3, !is.na(d3$Social.state.for.R))

## check out missing cases
summary(d3$dimorph.c) # 2 NA
summary(d3$logBM.s) # 2 NA
summary(d3$n.studies.c)  # 2 NA
summary(as.factor(d3[,"Seaonal/Non-seasonal Breeding"]))  # 2 NA
d3<- rename(d3, c("Seaonal/Non-seasonal Breeding" = "SeasonBreed"))
summary(d3$n.habitats.c)  # 3 NA
summary(as.factor(d3$Continent))  # 2 NA
summary(as.factor(d3$Habitat))  # 3 NA

nrow(d3) # 246
# --> find missing habitat info?
d3[is.na(d3$no.habitats),] # still missing data on one deer population --> have to exclude
d3<- d3[!is.na(d3$no.habitats),]
	


	### DV = Social organization (categorical) ###
	##############################################


get_prior(Social.state.for.R~1+dimorph.s+logBM.s+n.habitats.c+SeasonBreed+n.studies.c+(1|phylo)+(1|Genus_species)+(1|Continent)+(1|Habitat), data=d3, family=categorical())

prior=c(prior(normal(0,50), class=Intercept), 
prior(cauchy(0,2), class=sd, dpar="muNonvariableGL"), 
prior(cauchy(0,2), class=sd, dpar="muPairliving"), 
prior(cauchy(0,2), class=sd, dpar="muSexspecificSolGroup"), 
prior(cauchy(0,2), class=sd, dpar="muVariableIVSO"), 
prior(normal(0,5), class=b))

m1a<- brm(Social.state.for.R~1+dimorph.s+logBM.s+n.habitats.c+SeasonBreed+n.studies.c+(1|phylo)+(1|Genus_species)+(1|Continent)+(1|Habitat), 
	data=d3, family=categorical(), cov_ranef = list(phylo = A), prior = prior, 
	chains = 2, cores = 2, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99, max_treedepth = 15))
summary(m1a, prob=0.95)		
plot(m1a)
# --> more likely variable if more dimorphic, more studies and small body size

# phylogenetic signal
VarPhy.var<- postm1a$sd_phylo__muVariableIVSO_Intercept
VarPhy.group<- postm1a$sd_phylo__muNonvariableGL_Intercept
VarPhy.pair<- postm1a$sd_phylo__muPairliving_Intercept
VarPhy.ss<- postm1a$sd_phylo__muSexspecificSolGroup_Intercept

VarCont.var<- postm1a$sd_Continent__muVariableIVSO_Intercept
VarCont.group<- postm1a$sd_Continent__muNonvariableGL_Intercept
VarCont.pair<- postm1a$sd_Continent__muPairliving_Intercept
VarCont.ss<- postm1a$sd_Continent__muSexspecificSolGroup_Intercept

VarHab.var<- postm1a$sd_Habitat__muVariableIVSO_Intercept
VarHab.group<- postm1a$sd_Habitat__muNonvariableGL_Intercept
VarHab.pair<- postm1a$sd_Habitat__muPairliving_Intercept
VarHab.ss<- postm1a$sd_Habitat__muSexspecificSolGroup_Intercept

VarSpec.var<- postm1a$sd_Genus_species__muVariableIVSO_Intercept
VarSpec.group<- postm1a$sd_Genus_species__muNonvariableGL_Intercept
VarSpec.pair<- postm1a$sd_Genus_species__muPairliving_Intercept
VarSpec.ss<- postm1a$sd_Genus_species__muSexspecificSolGroup_Intercept

VarDistro<- pi^2 / 3

lambda<- (VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss)/
	(VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss+
	VarCont.var+VarCont.group+VarCont.pair+VarCont.ss+
	VarHab.var+VarHab.group+VarHab.pair+VarHab.ss+
	VarSpec.var+VarSpec.group+VarSpec.pair+VarSpec.ss+VarDistro)
mean(lambda); median(lambda); HPDI(lambda, prob=0.95)
#[1] 0.2366697
#[1] 0.2131081
#     |0.95      0.95| 
#0.04800237 0.47639147


postm1a<- posterior_samples(m1a)
# probability of each state at the root (global intercept) -> code adapted from Koster & McElreath 2017 Behav Ecol Sociobiol Multinomial analysis of behavior: statistical methods
{    K <- 5
    ns <- nrow(postm1a)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm1a[,k] 
            } else {
                ptemp <- rep(0,ns)
            }
            return(ptemp)
        })
        ## The values are converted to probabilities using the softmax function
        ## which ensures that the predicted values across categories sum to
        ## 100% probabilities.
        for ( s in 1:ns ) p[[i]][s,] <- softmax2( p[[i]][s,] )
    }
}
p_mean.brms <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.brms <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m1a.brms<- cbind(p_mean.brms, p_HPDI.brms[c(1,3,5,7,9),], p_HPDI.brms[c(2,4,6,8,10),])
rownames(pred_probs_m1a.brms)<- c("Group-living", "Pair-living", "Sex-specific", "Variable", "Solitary")
colnames(pred_probs_m1a.brms)<- c("mean", "lwr95", "upr95")
round(pred_probs_m1a.brms,2)
#             mean lwr95 upr95
#Group-living 0.08  0.00  0.38
#Pair-living  0.00  0.00  0.00
#Sex-specific 0.03  0.00  0.12
#Variable     0.77  0.29  1.00
#Solitary     0.12  0.00  0.34

# pp's
sum(postm1a$b_muNonvariableGL_dimorph.s>0)/length(postm1a$b_muNonvariableGL_dimorph.s) # 0.845
sum(postm1a$b_muNonvariableGL_SeasonBreedSeasonal>0)/length(postm1a$b_muNonvariableGL_SeasonBreedSeasonal) # 0.8697

## habitat type
## forest/woodland vs savanna/grassland
open<- rowMeans(cbind(postm1a[,"r_Habitat__muVariableIVSO[Savanna,Intercept]"],postm1a[,"r_Habitat__muVariableIVSO[Native.Grassland,Intercept]"]))
mean(open); HPDI(open, prob=0.95)
openvsclosed<- postm1a[,"r_Habitat__muVariableIVSO[Forest/Woodland,Intercept]"] - open
sum(openvsclosed>0)/length(openvsclosed)

## group-living in open habitats
sum(rowMeans(cbind(postm1a[,"r_Habitat__muNonvariableGL[Savanna,Intercept]"],postm1a[,"r_Habitat__muNonvariableGL[Native.Grassland,Intercept]"]))>0)/length(openvsclosed)


### figures ###

## Figure 1: Phylogeny with ancestral social organization
# save phylo to file and read again to get correct order of tip labels
write.tree(mamm.prun2.di.min, "PhyForFigure.tre")
phy<- read.tree("PhyForFigure.tre")

# create vectors of SO per species
t<- table(d3$Genus_species, d3$Social.state.for.R)
t2<- t[ order(match(rownames(t), phy$tip.label)), ]


tiff("Figure 1.tif", compression="lzw", height=6.0, width=6.0, units="cm", res=600, pointsize=5)
par(mar = c(2, 1, 1, 0), xpd=NA)
plot(phy, cex=0.5, lwd=1.5, label.offset=5) # use d3$Social.state.for.R BUT this is at pop level -> how to represent IVSO?
axisPhylo() # units are 20my, see add.scale.bar(length=20)

points(x=71, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[1])
points(x=72, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[2])
points(x=73, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[3])
points(x=74, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[4])
points(x=75, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[5])

points(x=71, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[1])
points(x=72, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[2])
points(x=73, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[3])
points(x=74, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[4])
points(x=75, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[5])

for(i in 1:nrow(t2)){if(t[i,1]>0) points(x=71, y=100-i, pch=15, cex=0.7, col=wes_palettes$Zissou[1])}
for(i in 1:nrow(t2)){if(t[i,3]>0) points(x=72, y=100-i, pch=15, cex=0.7, col=wes_palettes$Zissou[2])}
for(i in 1:nrow(t2)){if(t[i,4]>0) points(x=73, y=100-i, pch=15, cex=0.7, col=wes_palettes$Zissou[3])}
for(i in 1:nrow(t2)){if(t[i,2]>0) points(x=74, y=100-i, pch=15, cex=0.7, col=wes_palettes$Zissou[4])}
for(i in 1:nrow(t2)){if(t[i,5]>0) points(x=75, y=100-i, pch=15, cex=0.7, col=wes_palettes$Zissou[5])}

## create inset plots
par(fig = c(0.05,0.4, 0.2, 0.6), mar = c(2, 2, 4, 2), new = T)

# SO probabilities at the root
plot(c(0.05,0.3,0.5,0.75,1)~c(1,2,3,4,5), col="white", xaxt="n", ylab="", xlab="", main="Probability at root")
axis(1, at=c(1,2,3,4,5), labels=FALSE) # c("Sol", "Pair", "SS", "Grp", "Var"), cex.axis=0.6
text(x=c(1,2,3,4,5), y=-0.375, xpd=TRUE, adj=0, srt=90, labels=c("Sol", "Pair", "SS", "Grp", "Var"))

points(1, pred_probs_m1a.brms[5,1], pch=16, col=wes_palettes$Zissou[1], cex=2)
arrows(1, pred_probs_m1a.brms[5,2], 1, pred_probs_m1a.brms[5,3], length=0, angle=90, col=wes_palettes$Zissou[1], lwd=2)
points(2, pred_probs_m1a.brms[2,1], pch=16, col=wes_palettes$Zissou[2], cex=2)
arrows(2, pred_probs_m1a.brms[2,2], 2, pred_probs_m1a.brms[2,3], length=0, angle=90, col=wes_palettes$Zissou[2], lwd=2)
points(3, pred_probs_m1a.brms[3,1], pch=16, col=wes_palettes$Zissou[3], cex=2)
arrows(3, pred_probs_m1a.brms[3,2], 3, pred_probs_m1a.brms[3,3], length=0, angle=90, col=wes_palettes$Zissou[3], lwd=2)
points(4, pred_probs_m1a.brms[1,1], pch=16, col=wes_palettes$Zissou[4], cex=2)
arrows(4, pred_probs_m1a.brms[1,2], 4, pred_probs_m1a.brms[1,3], length=0, angle=90, col=wes_palettes$Zissou[4], lwd=2)
points(5, pred_probs_m1a.brms[4,1], pch=16, col=wes_palettes$Zissou[5], cex=2)
arrows(5, pred_probs_m1a.brms[4,2], 5, pred_probs_m1a.brms[4,3], length=0, angle=90, col=wes_palettes$Zissou[5], lwd=2)

par(fig = c(0,1,0,1), new=FALSE)
dev.off()


## Figure 2: Transitions with significant predictors
tiff("Figure 2.tif", compression="lzw", height=6.0, width=6.0, units="cm", res=600, pointsize=5)
par(mfrow=c(1,2))
plot(n.habitats.c~logBM.s, col="white", data=d3, ylim=c(0,1), ylab="Pr (Variable SO)", xlab="Female body size [Z]")
logBM.s<- seq(min(d3$logBM.s),max(d3$logBM.s),length=100)
newdata<- as.data.frame(logBM.s)
newdata$dimorph.s<- 0
newdata$n.habitats.c<- 0
newdata$n.studies.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.SO.bm<- fitted(m1a, newdata, re_formula=NA, summary=FALSE)
pred.variable<- pred.SO.bm[,,5]
for (i in sample(1:nrow(pred.variable), 100)){
	lines(pred.variable[i,]~logBM.s, col=col.alpha(wes_palettes$Zissou[5], 0.2))
	}
lines(logBM.s, colMeans(pred.variable), lwd=2)
mtext("A)", adj=-0.15, line=2)

plot(n.habitats.c~logBM.s, col="white", data=d3, ylim=c(0,1), ylab="Pr (Pair-living)", xlab="Female body size [Z]")
logBM.s<- seq(min(d3$logBM.s),max(d3$logBM.s),length=100)
newdata<- as.data.frame(logBM.s)
newdata$dimorph.s<- 0
newdata$n.habitats.c<- 0
newdata$n.studies.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.pair<- pred.SO.bm[,,3]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~logBM.s, col=col.alpha(wes_palettes$Zissou[2], 0.2))
	}
lines(logBM.s, colMeans(pred.pair), lwd=2)
mtext("B)", adj=-0.15, line=2)
dev.off()



	### DV = # social organizations (numerical) ###
	###############################################


m1b<- brm(as.numeric(no.SO_allPop)~1+dimorph.s+logBM.s+n.habitats.c+SeasonBreed+n.studies.c+(1|phylo)+(1|Genus_species)+(1|Continent)+(1|Habitat), data=d3, 
	family=poisson(), cov_ranef = list(phylo = A),
	prior = c(prior(normal(0,50), class=Intercept), 
	prior(cauchy(0,2), class=sd), 
	prior(normal(0,5), class=b)),
	sample_prior = TRUE, chains = 2, cores = 2, 
	iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)
	)
summary(m1b, prob=0.95)		
plot(m1b)
# --> more if more studies

postm1b<- posterior_samples(m1b)
sum(postm1b[,"b_logBM.s"]<0)/length(postm1b[,"b_logBM.s"]) # 0.9175

## predicted # SO at ancestral state
exp(mean(postm1b[,"b_Intercept"]))
# 2.567788
exp(HPDI(postm1b[,"b_Intercept"],prob=0.95))
#   |0.95    0.95| 
# 1.914155 3.309647


## phylogenetic signal

# phylogenetic signal
VarDistro<- log(1+1/exp(postm1b$b_Intercept)) # see Nakagawa & Schielzeth 2013 appendix: b0 should be obtained from a model with centered or scaled variables, which is the case here

lambda<- postm1b$sd_phylo__Intercept/(postm1b$sd_phylo__Intercept+postm1b$sd_Habitat__Intercept+postm1b$sd_Continent__Intercept+postm1b$sd_Genus_species__Intercept+VarDistro)
mean(lambda); median(lambda); HPDI(lambda, prob=0.95)
#[1] 0.1286302
#[1] 0.1130974
#       |0.95        0.95| 
#1.243786e-05 2.944502e-01

## forest/woodland vs savanna/grassland
open<- rowMeans(cbind(postm1b[,"r_Habitat[Savanna,Intercept]"],postm1b[,"r_Habitat[Native.Grassland,Intercept]"]))
mean(open); HPDI(open, prob=0.95)
mean(postm1b[,"r_Habitat[Forest/Woodland,Intercept]"]); HPDI(postm1b[,"r_Habitat[Forest/Woodland,Intercept]"], prob=0.95)
openvsclosed<-  postm1b[,"r_Habitat[Forest/Woodland,Intercept]"] - rowMeans(cbind(postm1b[,"r_Habitat[Savanna,Intercept]"],postm1b[,"r_Habitat[Native.Grassland,Intercept]"]))
sum(openvsclosed>0)/length(openvsclosed)



