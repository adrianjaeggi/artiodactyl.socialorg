## load packages
library(ape)
library(wesanderson)
library(brms)
library(rethinking)
library(readxl)
library(paleotree)
library(phytools)
library(mice)
library(rstan)

#stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


	
	#### data prep ####
	###################

## read data file
d<- as.data.frame(read_excel("Artiodactyla_Dataset_corrected_11Mar.xlsx", sheet=1, skip=6)) # assuming working directory is set to the folder containing this file


# fix continent misspellings
d$Continent[as.factor(d$Continent)=="africa"]<- "Africa"	
d$Continent[as.factor(d$Continent)=="Europe (Introduced)"]<- "Europe"	
d$Continent[as.factor(d$Continent)=="North America (Introduced)"]<- "North America"	
d$Continent[as.factor(d$Continent)=="Oceania (Introduced)"]<- "Oceania"	
d$Continent[as.factor(d$Continent)=="Norther America"]<- "North America"	
	
# clean habitat --> if multiple, code as multiple
d$Habitat[as.factor(d$Habitat)=="Artificial (Terrestrial); Forest/Woodland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Artificial (Terrestrial); Native Grassland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Artificial (Terrestrial); Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Desert; Native Grassland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Desert; Shrubland; Wetland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Inland Rocky Areas"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Inland Rocky Areas; Native Grassland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Native Grassland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Native Grassland; Savanna"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Native Grassland; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Native Grassland; Wetland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Savanna"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woodland; Wetland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Inland Rocky Areas; Native Grassland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Inland Rocky Areas; Native Grassland; Savanna"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Inland Rocky Areas; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Native Grassland; Savanna"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Native Grassland; Savanna; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Native Grassland; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="savanna"]<- "Savanna"	
d$Habitat[as.factor(d$Habitat)=="Savanna; Shrubland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Shrubland; Wetland"]<- "Multiple"	
d$Habitat[as.factor(d$Habitat)=="Forest/Woosland"]<- "Forest/Woodland"	
d$Habitat[as.factor(d$Habitat)=="inland Rocky Areas"]<- "Inland Rocky Areas"	
## one NA

## prepare phylogeny
mamm<- read.tree("mammaltree.tre") # assuming working directory is set to the folder containing this file
plot(mamm)

## check species names match and trim trees to dataset
setdiff(d$Genus_species, mamm$tip.label) #    

## see iucnredlist.org for synonyms or subspecies names
# Beatragus_hunteri = Damaliscus_hunteri
# Capricornis_crispus = Naemorhedus_crispus
# Eudorcas_thomsonii = Gazella_thomsonii
# Nanger_granti = Gazella_granti
# Oryx_beisa = missing -> replace with dammah or leucoryx (both equally related to gazella)
# Lama_glama_guanicoe = Lama_glama
# Rucervus_duvaucelii = Cervus_duvaucelii
# Rusa_unicolor = Cervus_unicolor
# Moschus_leucogaster = Moschus_cupreus = Moschus_chrysogaster
d$Genus_species[d$Genus_species=="Beatragus_hunteri"]<- "Damaliscus_hunteri"
d$Genus_species[d$Genus_species=="Capricornis_crispus"]<- "Naemorhedus_crispus"
d$Genus_species[d$Genus_species=="Nanger_granti"]<- "Gazella_granti"
d$Genus_species[d$Genus_species=="Eudorcas_thomsonii"]<- "Gazella_thomsonii"
d$Genus_species[d$Genus_species=="Oryx_beisa"]<- "Oryx_dammah"
d$Genus_species[d$Genus_species=="Lama_glama_guanicoe"]<- "Lama_glama"
d$Genus_species[d$Genus_species=="Rucervus_duvaucelii"]<- "Cervus_duvaucelii"
d$Genus_species[d$Genus_species=="Rusa_unicolor"]<- "Cervus_unicolor"
d$Genus_species[d$Genus_species=="Moschus_leucogaster"]<- "Moschus_chrysogaster"
d$Genus_species[d$Genus_species=="Moschus_cupreus"]<- "Moschus_chrysogaster"
setdiff(d$Genus_species, mamm$tip.label) # none     
mamm.prun<- drop.tip(mamm, setdiff(mamm$tip.label,d$Genus_species))
plot(mamm.prun, type="fan", cex=0.5)
# check requirements for analysis
is.ultrametric(mamm.prun) # TRUE
is.rooted(mamm.prun) # TRUE
is.binary.tree(mamm.prun) # FALSE
mamm.prun.di <- multi2di(mamm.prun) # this forces the tree to be binary (as required for analysis)
is.binary.tree(mamm.prun.di) # TRUE
mamm.prun.di.min<- minBranchLength(mamm.prun.di, 1e-10) # makes the artificially introduced branches so small that it doesn't matter
mamm.prun.di.min$node.label<- NULL
plot(mamm.prun.di.min, type="fan", cex=0.5)

# convert to covariance matrix (see https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html)
A <- vcv.phylo(mamm.prun.di.min)





		### Ancestral state estimation for predictors ###
		#################################################
	
# load species-level data and change species names (see above)
d.spec<- as.data.frame(read_excel("Artiodacytla_Dataset_corrected_17 Dec.xlsx", sheet=2, skip=6)) # assuming working directory is set to the folder containing this file
d.spec$Genus_species[d.spec$Genus_species=="Beatragus_hunteri"]<- "Damaliscus_hunteri"
d.spec$Genus_species[d.spec$Genus_species=="Capricornis_crispus"]<- "Naemorhedus_crispus"
d.spec$Genus_species[d.spec$Genus_species=="Nanger_granti"]<- "Gazella_granti"
d.spec$Genus_species[d.spec$Genus_species=="Eudorcas_thomsonii"]<- "Gazella_thomsonii"
d.spec$Genus_species[d.spec$Genus_species=="Oryx_beisa"]<- "Oryx_dammah"
d.spec$Genus_species[d.spec$Genus_species=="Lama_glama_guanicoe"]<- "Lama_glama"
d.spec$Genus_species[d.spec$Genus_species=="Rucervus_duvaucelii"]<- "Cervus_duvaucelii"
d.spec$Genus_species[d.spec$Genus_species=="Rusa_unicolor"]<- "Cervus_unicolor"
d.spec$Genus_species[d.spec$Genus_species=="Moschus_leucogaster"]<- "Moschus_chrysogaster"
d.spec$Genus_species[d.spec$Genus_species=="Moschus_cupreus"]<- "Moschus_chrysogaster"

## body size
# create named list of trait values
bm<- log(d.spec[,"F_BM (kg)"])
names(bm)<- d.spec$Genus_species
bm<- bm[which(!is.na(bm))]

# reprune tree
mamm.prun2<- drop.tip(mamm, setdiff(mamm$tip.label,names(bm)))
plot(mamm.prun2, type="fan", cex=0.5)
is.ultrametric(mamm.prun2)
is.rooted(mamm.prun2)
is.binary.tree(mamm.prun2) # FALSE
mamm.prun2.di <- multi2di(mamm.prun2)
is.binary.tree(mamm.prun2.di) # TRUE
mamm.prun2.di.min<- minBranchLength(mamm.prun2.di, 1e-10)
mamm.prun2.di.min$node.label<- NULL

# fastAnc
ace.bm.fastAnc = fastAnc(mamm.prun2.di.min, bm, CI = TRUE)
ace.bm.fastAnc
# --> at root (node 79): 4.187368, 95%CI = 1.778210 6.596526
# ace
ace.bm.ace = ace(bm, mamm.prun2.di.min, type="continuous")
ace.bm.ace
# --> at root (node 79): 4.187368, 95%CI = 1.755541 6.619196
# virtually the same
# --> for phylogenetic mean, center body size on ace.bm.ace$ace[1] 


## sexual dimorphism M:F_BM
# create named list of trait values
dimo<- d.spec[,"M:F_BM"]
names(dimo)<- d.spec$Genus_species
dimo<- dimo[which(!is.na(dimo))]

# reprun tree
mamm.prun3<- drop.tip(mamm, setdiff(mamm$tip.label,names(dimo)))
plot(mamm.prun3, type="fan", cex=0.5)
is.ultrametric(mamm.prun3)
is.rooted(mamm.prun3)
is.binary.tree(mamm.prun3) # FALSE
mamm.prun3.di <- multi2di(mamm.prun3)
is.binary.tree(mamm.prun3.di) # TRUE
mamm.prun3.di.min<- minBranchLength(mamm.prun3.di, 1e-10)
mamm.prun3.di.min$node.label<- NULL

# fastAnc
ace.dimo.fastAnc = fastAnc(mamm.prun3.di.min, dimo, CI = TRUE)
ace.dimo.fastAnc
# --> at root (node 79): 1.214002, 95%CI = 0.537744 1.890261
# ace
ace.dimo.ace = ace(dimo, mamm.prun3.di.min, type="continuous")
ace.dimo.ace
# --> at root (node 79): 1.214002, 95%CI = 0.5313545 1.896650
# virtually the same
# --> for phylogenetic mean, center sexual dimorphism on ace.dimo.ace$ace[1] 


## breeding seasonality
# create named list of trait values
seas<- d.spec[,"Seaonal/Non-seasonal Breeding"]
names(seas)<- d.spec$Genus_species
seas<- seas[which(!is.na(seas))]

# reprun tree
mamm.prun4<- drop.tip(mamm, setdiff(mamm$tip.label,names(seas)))
plot(mamm.prun4, type="fan", cex=0.5)
is.ultrametric(mamm.prun4)
is.rooted(mamm.prun4)
is.binary.tree(mamm.prun4) # FALSE
mamm.prun4.di <- multi2di(mamm.prun4)
is.binary.tree(mamm.prun4.di) # TRUE
mamm.prun4.di.min<- minBranchLength(mamm.prun4.di, 1e-10)
mamm.prun4.di.min$node.label<- NULL

# compare ER, ARD
ace.seas.ER = ace(seas, mamm.prun4.di.min, type="discrete", model="ER")
ace.seas.ARD = ace(seas, mamm.prun4.di.min, type="discrete", model="ARD")

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



		##### Phylogenetic mixed-effects models #####
		#############################################
		

	## Model 1 ##
	#############
	
# check DV and relevel
summary(as.factor(d[,"Social state - for R"]))
colnames(d)[18]<- "Social.state.for.R"
d$Social.state.for.R<- relevel(as.factor(d$Social.state.for.R), ref="Solitary")

## use body size and dimorphism from fossil record
d$dimorph.s.foss<- (d[,"M:F_BM"]-1)/sd(d[,"M:F_BM"], na.rm=TRUE)
d$logBM.s.foss<- (log(d[,"F_BM (kg)"])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))

## center number of studies and habitats on 1
d$n.studies.c<- as.numeric(d$no.studies)-median(as.numeric(d$no.studies)) # median is 1
d$n.habitats.c<- as.numeric(d$no.habitats)-1 
d$phylo<- d$Genus_species

# use multiple imputation because brms cannot impute categorical predictors (seasonal-breeding); see https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html for details
# include all available body size data (e.g. sometimes avge mass and dimorphism are known, but not female mass
d.imp<- d[,c("Social.state.for.R", "Dimorphic/Monomorphic", "dimorph.s.foss", "logBM.s.foss", "Male_BM (kg)", "F_BM (kg)", "M:F_BM", "Avg_BM (kg)", "n.habitats.c", "SeasonBreed", "n.studies.c", "phylo", "Genus_species", "Continent", "Habitat")]
colnames(d.imp)[2]<- "dimorphcat"
colnames(d.imp)[5]<- "maleBM"
colnames(d.imp)[6]<- "femaleBM"
colnames(d.imp)[7]<- "dimorph"
colnames(d.imp)[8]<- "avgeBM"
# convert categorical variables to factors, otherwise they won't be imputed
d.imp$dimorphcat<- as.factor(d.imp$dimorphcat)
d.imp$SeasonBreed<- as.factor(d.imp$SeasonBreed)
d.imp$Habitat<- as.factor(d.imp$Habitat)

d.imp.mice<- mice(d.imp, m=10)

# extract dfs
list_nos <- seq(1:10)
fun_extract_dfs <- function(x) {
  complete(d.imp.mice, action = x) 
  }
list_d.imp <- lapply(list_nos, fun_extract_dfs)


prior=c(prior(normal(0,50), class=Intercept), 
prior(cauchy(0,2), class=sd, dpar="muGroupLiving"), 
prior(cauchy(0,2), class=sd, dpar="muPairLiving"), 
prior(cauchy(0,2), class=sd, dpar="muSexSpecific"), 
prior(cauchy(0,2), class=sd, dpar="muVariable"), 
prior(normal(0,5), class=b))

m1.imp<- brm_multiple(Social.state.for.R~1+dimorph.s.foss+logBM.s.foss+n.habitats.c+SeasonBreed+n.studies.c+(1|phylo)+(1|Genus_species)+(1|Continent)+(1|Habitat), 
	data=list_d.imp, family=categorical(), cov_ranef = list(phylo = A), prior = prior, 
	chains = 2, cores = 2, iter = 4000, warmup = 1000, thin=20, control = list(adapt_delta = 0.99, max_treedepth = 15))

saveRDS(m1.imp, "m1.imp.rds")
summary(m1.imp)
plot(m1.imp)

# rhat's high, this is expected because chains fit to slightly different data; plots look good though
# check rhat's within each dataset:
round(m1.imp$rhats, 2) # --> those look good!
max(apply(m1.imp$rhats, 1, max)) # 1.037 is the highest Rhat

postm1.imp<- posterior_samples(m1.imp)

# phylogenetic signal: sum of all SO variances / sum of all variance components + distribution-specific variance
VarPhy.var<- postm1.imp$sd_phylo__muVariable_Intercept^2
VarPhy.group<- postm1.imp$sd_phylo__muGroupLiving_Intercept^2
VarPhy.pair<- postm1.imp$sd_phylo__muPairLiving_Intercept^2
VarPhy.ss<- postm1.imp$sd_phylo__muSexSpecific_Intercept^2

VarCont.var<- postm1.imp$sd_Continent__muVariable_Intercept^2
VarCont.group<- postm1.imp$sd_Continent__muGroupLiving_Intercept^2
VarCont.pair<- postm1.imp$sd_Continent__muPairLiving_Intercept^2
VarCont.ss<- postm1.imp$sd_Continent__muSexSpecific_Intercept^2

VarHab.var<- postm1.imp$sd_Habitat__muVariable_Intercept^2
VarHab.group<- postm1.imp$sd_Habitat__muGroupLiving_Intercept^2
VarHab.pair<- postm1.imp$sd_Habitat__muPairLiving_Intercept^2
VarHab.ss<- postm1.imp$sd_Habitat__muSexSpecific_Intercept^2

VarSpec.var<- postm1.imp$sd_Genus_species__muVariable_Intercept^2
VarSpec.group<- postm1.imp$sd_Genus_species__muGroupLiving_Intercept^2
VarSpec.pair<- postm1.imp$sd_Genus_species__muPairLiving_Intercept^2
VarSpec.ss<- postm1.imp$sd_Genus_species__muSexSpecific_Intercept^2

VarDistro<- pi^2 / 3 # see Nakagawa & Schielzeth 2012 MEE

lambda<- (VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss)/
	(VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss+
	VarCont.var+VarCont.group+VarCont.pair+VarCont.ss+
	VarHab.var+VarHab.group+VarHab.pair+VarHab.ss+
	VarSpec.var+VarSpec.group+VarSpec.pair+VarSpec.ss+VarDistro)
mean(lambda); HPDI(lambda, prob=0.95); sum(lambda>0.01)/length(lambda)
#[1] 0.04628662
#       |0.95        0.95| 
#4.283857e-05 1.667173e-01
# [1] 0.709

# probability of each state at the root (global intercept) -> code adapted from Koster & McElreath 2017 Behav Ecol Sociobiol (could also use brms:fitted(), as further below)
{    K <- 5
    ns <- nrow(postm1.imp)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm1.imp[,k] 
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
p_mean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m1.imp<- cbind(p_mean, p_HPDI[c(1,3,5,7,9),], p_HPDI[c(2,4,6,8,10),])
rownames(pred_probs_m1.imp)<- c("Group-living", "Pair-living", "Sex-specific", "Variable", "Solitary")
colnames(pred_probs_m1.imp)<- c("mean", "lwr95", "upr95")
round(pred_probs_m1.imp,2) 
#             mean lwr95 upr95
#Group-living 0.00     0  0.00
#Pair-living  0.77     0  1.00
#Sex-specific 0.04     0  0.06
#Variable     0.17     0  0.98
#Solitary     0.03     0  0.15

# pp's
m1.imp.PP<- 0
for(i in 1:24){
	ifelse(mean(postm1.imp[,i])>0, m1.imp.PP[i]<- sum(postm1.imp[,i]>0)/nrow(postm1.imp), m1.imp.PP[i]<- sum(postm1.imp[,i]<0)/nrow(postm1.imp))
	}
names(m1.imp.PP)<- colnames(postm1.imp[,1:24])
m1.imp.PP
## strong (PP>0.9) associations for GL-body size, GL- seasonal, GL-n studies, PL-body size, SS-n.studies, variable-dimorph, variable-n.studies

## habitat type
## forest/woodland vs savanna/grassland
open<- rowMeans(cbind(postm1.imp[,"r_Habitat__muVariable[Savanna,Intercept]"],postm1.imp[,"r_Habitat__muVariable[Native.Grassland,Intercept]"]))
mean(open); HPDI(open, prob=0.95)
openvsclosed<- postm1.imp[,"r_Habitat__muVariable[Forest/Woodland,Intercept]"] - open
sum(openvsclosed>0)/length(openvsclosed) # 0.35

## group-living in open habitats
sum(rowMeans(cbind(postm1.imp[,"r_Habitat__muGroupLiving[Savanna,Intercept]"],postm1.imp[,"r_Habitat__muGroupLiving[Native.Grassland,Intercept]"]))>0)/length(openvsclosed) # 0.63

## pair-living in closed habitats
sum(postm1.imp[,"r_Habitat__muPairLiving[Forest/Woodland,Intercept]"]>0)/length(openvsclosed) # 0.30


# probability of each state at the phylogenetic mean of predictors (log body size = 4.187368, dimorphism = 1.214002)
bm.z<- (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))
dm.z<- (as.numeric(ace.dimo.ace$ace[1])-1)/sd(d[,"M:F_BM"], na.rm=TRUE)
{    K <- 5
    ns <- nrow(postm1.imp)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm1.imp[,k] + dm.z*postm1.imp[,k+(k*4)] + bm.z*postm1.imp[,k+(k*5)-(k-1)]
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
p_mean.phymean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.phymean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m1.imp.phymean<- cbind(p_mean.phymean, p_HPDI.phymean[c(1,3,5,7,9),], p_HPDI.phymean[c(2,4,6,8,10),])
rownames(pred_probs_m1.imp.phymean)<- c("Group-living", "Pair-living", "Sex-specific", "Variable", "Solitary")
colnames(pred_probs_m1.imp.phymean)<- c("mean", "lwr95", "upr95")
round(pred_probs_m1.imp.phymean,2)
#              mean lwr95 upr95
#Group-living 0.09  0.00  0.47
#Pair-living  0.00  0.00  0.00
#Sex-specific 0.02  0.00  0.01
#Variable     0.73  0.16  1.00
#Solitary     0.17  0.00  0.56


# same at lower 95% CI
bm.z<- (as.numeric(ace.bm.ace$CI95[1,1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))
dm.z<- (as.numeric(ace.dimo.ace$CI95[1,1])-1)/sd(d[,"M:F_BM"], na.rm=TRUE)
{    K <- 5
    ns <- nrow(postm1.imp)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm1.imp[,k] + dm.z*postm1.imp[,k+(k*4)] + bm.z*postm1.imp[,k+(k*5)-(k-1)]
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
p_mean.phymean.lwr95 <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.phymean.lwr95 <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m1.imp.phymean.lwr95<- cbind(p_mean.phymean.lwr95, p_HPDI.phymean.lwr95[c(1,3,5,7,9),], p_HPDI.phymean.lwr95[c(2,4,6,8,10),])
rownames(pred_probs_m1.imp.phymean.lwr95)<- c("Group-living", "Pair-living", "Sex-specific", "Variable", "Solitary")
colnames(pred_probs_m1.imp.phymean.lwr95)<- c("mean", "lwr95", "upr95")
round(pred_probs_m1.imp.phymean.lwr95,2)
#             mean lwr95 upr95
#Group-living 0.01     0  0.02
#Pair-living  0.53     0  1.00
#Sex-specific 0.04     0  0.11
#Variable     0.30     0  0.97
#Solitary     0.12     0  0.66

# and at upper 95% CI
bm.z<- (as.numeric(ace.bm.ace$CI95[1,2])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))
dm.z<- (as.numeric(ace.dimo.ace$CI95[1,2])-1)/sd(d[,"M:F_BM"], na.rm=TRUE)
{    K <- 5
    ns <- nrow(postm1.imp)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm1.imp[,k] + dm.z*postm1.imp[,k+(k*4)] + bm.z*postm1.imp[,k+(k*5)-(k-1)]
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
p_mean.phymean.upr95 <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI.phymean.upr95 <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m1.imp.phymean.upr95<- cbind(p_mean.phymean.upr95, p_HPDI.phymean.upr95[c(1,3,5,7,9),], p_HPDI.phymean.upr95[c(2,4,6,8,10),])
rownames(pred_probs_m1.imp.phymean.upr95)<- c("Group-living", "Pair-living", "Sex-specific", "Variable", "Solitary")
colnames(pred_probs_m1.imp.phymean.upr95)<- c("mean", "upr95", "upr95")
round(pred_probs_m1.imp.phymean.upr95,2)
#             mean upr95 upr95
#Group-living 0.43     0  0.99
#Pair-living  0.00     0  0.00
#Sex-specific 0.02     0  0.01
#Variable     0.48     0  0.98
#Solitary     0.07     0  0.33


## Figure 2: Phylogeny with social organization at root and phylo mean
# save phylo to file and read again to get correct order of tip labels
write.tree(mamm.prun.di.min, "PhyForFigure.tre")
phy<- read.tree("PhyForFigure.tre")
phy<- drop.tip(phy, setdiff(phy$tip.label,d$Genus_species))


# create vector of SO per species
t<- table(d$Genus_species, d$Social.state.for.R)
t2<- t[ order(match(rownames(t), phy$tip.label)), ]

tiff("Figure 2.tif", compression="lzw", height=6.0, width=6.0, units="cm", res=600, pointsize=5)
par(mar = c(2, 1, 1, 0), xpd=NA)
plot(phy, cex=0.5, lwd=1.5, label.offset=5)
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

for(i in 1:nrow(t2)){if(t2[i,"Solitary"]>0) points(x=71, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[1])}
for(i in 1:nrow(t2)){if(t2[i,"PairLiving"]>0) points(x=72, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[2])}
for(i in 1:nrow(t2)){if(t2[i,"SexSpecific"]>0) points(x=73, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[3])}
for(i in 1:nrow(t2)){if(t2[i,"GroupLiving"]>0) points(x=74, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[4])}
for(i in 1:nrow(t2)){if(t2[i,"Variable"]>0) points(x=75, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[5])}

## create inset plots
par(fig = c(0.05,0.4, 0.45, 0.75), mar = c(0, 2, 3, 2), new = T)

# SO probabilities at the root
plot(c(0.05,0.3,0.5,0.75,1)~c(1,2,3,4,5), col="white", xaxt="n", ylab="", xlab="", main="Ancestral state")
axis(1, at=c(1,2,3,4,5), labels=FALSE) # c("Sol", "Pair", "SS", "Grp", "Var"), cex.axis=0.6
# text(x=c(1,2,3,4,5), y=-0.375, xpd=TRUE, adj=0, srt=90, labels=c("Sol", "Pair", "SS", "Grp", "Var"))

points(1, pred_probs_m1.imp[5,1], pch=16, col=wes_palettes$Zissou[1], cex=2)
arrows(1, pred_probs_m1.imp[5,2], 1, pred_probs_m1.imp[5,3], length=0, angle=90, col=wes_palettes$Zissou[1], lwd=2)
points(2, pred_probs_m1.imp[2,1], pch=16, col=wes_palettes$Zissou[2], cex=2)
arrows(2, pred_probs_m1.imp[2,2], 2, pred_probs_m1.imp[2,3], length=0, angle=90, col=wes_palettes$Zissou[2], lwd=2)
points(3, pred_probs_m1.imp[3,1], pch=16, col=wes_palettes$Zissou[3], cex=2)
arrows(3, pred_probs_m1.imp[3,2], 3, pred_probs_m1.imp[3,3], length=0, angle=90, col=wes_palettes$Zissou[3], lwd=2)
points(4, pred_probs_m1.imp[1,1], pch=16, col=wes_palettes$Zissou[4], cex=2)
arrows(4, pred_probs_m1.imp[1,2], 4, pred_probs_m1.imp[1,3], length=0, angle=90, col=wes_palettes$Zissou[4], lwd=2)
points(5, pred_probs_m1.imp[4,1], pch=16, col=wes_palettes$Zissou[5], cex=2)
arrows(5, pred_probs_m1.imp[4,2], 5, pred_probs_m1.imp[4,3], length=0, angle=90, col=wes_palettes$Zissou[5], lwd=2)

par(fig = c(0.05,0.4, 0.1, 0.45), mar = c(2, 2, 3, 2), new = T)

plot(c(0.05,0.3,0.5,0.75,1)~c(1,2,3,4,5), col="white", xaxt="n", ylab="", xlab="", main="Phylogenetic mean")
axis(1, at=c(1,2,3,4,5), labels=FALSE) # c("Sol", "Pair", "SS", "Grp", "Var"), cex.axis=0.6
text(x=c(1,2,3,4,5), y=-0.375, xpd=TRUE, adj=0, srt=90, labels=c("Sol", "Pair", "SS", "Grp", "Var"))

points(1, pred_probs_m1.imp.phymean[5,1], pch=16, col=wes_palettes$Zissou[1], cex=2)
arrows(1, pred_probs_m1.imp.phymean[5,2], 1, pred_probs_m1.imp.phymean[5,3], length=0, angle=90, col=wes_palettes$Zissou[1], lwd=2)
points(2, pred_probs_m1.imp.phymean[2,1], pch=16, col=wes_palettes$Zissou[2], cex=2)
arrows(2, pred_probs_m1.imp.phymean[2,2], 2, pred_probs_m1.imp.phymean[2,3], length=0, angle=90, col=wes_palettes$Zissou[2], lwd=2)
points(3, pred_probs_m1.imp.phymean[3,1], pch=16, col=wes_palettes$Zissou[3], cex=2)
arrows(3, pred_probs_m1.imp.phymean[3,2], 3, pred_probs_m1.imp.phymean[3,3], length=0, angle=90, col=wes_palettes$Zissou[3], lwd=2)
points(4, pred_probs_m1.imp.phymean[1,1], pch=16, col=wes_palettes$Zissou[4], cex=2)
arrows(4, pred_probs_m1.imp.phymean[1,2], 4, pred_probs_m1.imp.phymean[1,3], length=0, angle=90, col=wes_palettes$Zissou[4], lwd=2)
points(5, pred_probs_m1.imp.phymean[4,1], pch=16, col=wes_palettes$Zissou[5], cex=2)
arrows(5, pred_probs_m1.imp.phymean[4,2], 5, pred_probs_m1.imp.phymean[4,3], length=0, angle=90, col=wes_palettes$Zissou[5], lwd=2)

par(fig = c(0,1,0,1), new=FALSE)
dev.off()



## Figure 3: Transitions with all predictors
## predictions for body size - everything else at baseline
logBM.s.foss<- seq(0,max(d$logBM.s.foss, na.rm=TRUE),length=100)
newdata<- as.data.frame(logBM.s.foss)
newdata$dimorph.s.foss<- 0
newdata$n.habitats.c<- 0
newdata$n.studies.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.SO.bm<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)


## predictions for dimorphism - for pair-living body size at baseline, for others at phylo mean
dimorph.s.foss<- seq(min(d$dimorph.s.foss, na.rm=TRUE),max(d$dimorph.s.foss, na.rm=TRUE),length=100)
newdata<- as.data.frame(dimorph.s.foss)
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
newdata$n.habitats.c<- 0
newdata$n.studies.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.SO.dm.anc<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.dm.phymean<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)


## predictions for n habitats - for pair-living body size at baseline, for others at phylo mean
n.habitats.c<- seq(0,max(d$n.habitats.c, na.rm=TRUE),length=100)
newdata<- as.data.frame(n.habitats.c)
newdata$dimorph.s.foss<- 0
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
newdata$n.studies.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.SO.hab.anc<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.hab.phymean<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)


## predictions for n studies - for pair-living body size at baseline, for others at phylo mean
n.studies.c<- seq(0,max(d$n.studies.c, na.rm=TRUE),length=100)
newdata<- as.data.frame(n.studies.c)
newdata$dimorph.s.foss<- 0
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
newdata$n.habitats.c<- 0
newdata$SeasonBreed<- "Non-seasonal"
pred.SO.stud.anc<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.stud.phymean<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)


## predictions for breeding seasonality - for pair-living body size at baseline, for others at phylo mean
SeasonBreed<- "Non-seasonal"
newdata<- as.data.frame(SeasonBreed)
newdata$dimorph.s.foss<- 0
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
newdata$n.habitats.c<- 0
newdata$n.studies.c<- 0
pred.SO.nonseas.anc<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.nonseas.phymean<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$SeasonBreed<- "Seasonal"
pred.SO.seas.phymean<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.seas.anc<- fitted(m1.imp, newdata, re_formula=NA, summary=FALSE)


## plot
tiff("Figure 3.tif", compression="lzw", height=10.0, width=10.0, units="cm", res=600, pointsize=5)
par(mfrow=c(5,5), mar = c(4, 4, 2, 2))

## body size
plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Solitary)", xlab="Female body size [Z]")
pred.sol<- pred.SO.bm[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.sol), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
mtext("A)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Pair-living)", xlab="Female body size [Z]")
pred.pair<- pred.SO.bm[,,3]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.pair), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muPairLiving_logBM.s.foss"]), 2))
mtext("B)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Sex-specific)", xlab="Female body size [Z]")
pred.ss<- pred.SO.bm[,,4]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.ss), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muSexSpecific_logBM.s.foss"]), 2))
mtext("C)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Group-living)", xlab="Female body size [Z]")
pred.group<- pred.SO.bm[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.group), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muGroupLiving_logBM.s.foss"]), 2))
mtext("D)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Variable)", xlab="Female body size [Z]")
pred.var<- pred.SO.bm[,,5]
for (i in sample(1:nrow(pred.var), 100)){
	lines(pred.var[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[5], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.var), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muVariable_logBM.s.foss"]), 2))
mtext("E)", adj=-0.5, line=0)


## dimorphism
plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Solitary)", xlab="Sexual dimorphism [Z]")
pred.sol<- pred.SO.dm.phymean[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.sol), lwd=1)
mtext("F)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Pair-living)", xlab="Sexual dimorphism [Z]")
pred.pair<- pred.SO.dm.anc[,,3]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.pair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muPairLiving_dimorph.s.foss"]), 2))
mtext("G)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Sex-specific)", xlab="Sexual dimorphism [Z]")
pred.ss<- pred.SO.dm.phymean[,,4]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.ss), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muSexSpecific_dimorph.s.foss"]), 2))
mtext("H)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Group-living)", xlab="Sexual dimorphism [Z]")
pred.group<- pred.SO.dm.phymean[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.group), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muGroupLiving_dimorph.s.foss"]), 2))
mtext("I)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Variable)", xlab="Sexual dimorphism [Z]")
pred.var<- pred.SO.dm.phymean[,,5]
for (i in sample(1:nrow(pred.var), 100)){
	lines(pred.var[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[5], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.var), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muVariable_dimorph.s.foss"]), 2))
mtext("J)", adj=-0.5, line=0)


## n habitats
plot(logBM.s.foss~n.habitats.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Solitary)", xlab="Number of habitats")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.sol<- pred.SO.hab.phymean[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~n.habitats.c, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(n.habitats.c, colMeans(pred.sol), lwd=1)
mtext("K)", adj=-0.5, line=0)

plot(logBM.s.foss~n.habitats.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Pair-living)", xlab="Number of habitats")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.pair<- pred.SO.hab.anc[,,3]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~n.habitats.c, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(n.habitats.c, colMeans(pred.pair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muPairLiving_n.habitats.c"]), 2))
mtext("L)", adj=-0.5, line=0)

plot(logBM.s.foss~n.habitats.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Sex-specific)", xlab="Number of habitats")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.ss<- pred.SO.hab.phymean[,,4]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~n.habitats.c, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(n.habitats.c, colMeans(pred.ss), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muSexSpecific_n.habitats.c"]), 2))
mtext("M)", adj=-0.5, line=0)

plot(logBM.s.foss~n.habitats.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Group-living)", xlab="Number of habitats")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.group<- pred.SO.hab.phymean[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~n.habitats.c, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(n.habitats.c, colMeans(pred.group), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muGroupLiving_n.habitats.c"]), 2))
mtext("N)", adj=-0.5, line=0)

plot(logBM.s.foss~n.habitats.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Variable)", xlab="Number of habitats")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.var<- pred.SO.hab.phymean[,,5]
for (i in sample(1:nrow(pred.var), 100)){
	lines(pred.var[i,]~n.habitats.c, col=col.alpha(wes_palettes$Zissou[5], 0.2), lwd=0.5)
	}
lines(n.habitats.c, colMeans(pred.var), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muVariable_n.habitats.c"]), 2))
mtext("O)", adj=-0.5, line=0)


## n studies
plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Solitary)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.sol<- pred.SO.stud.phymean[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.sol), lwd=1)
mtext("P)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Pair-living)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.pair<- pred.SO.stud.anc[,,3]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.pair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muPairLiving_n.studies.c"]), 2))
mtext("Q)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Sex-specific)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.ss<- pred.SO.stud.phymean[,,4]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.ss), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muSexSpecific_n.studies.c"]), 2))
mtext("R)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Group-living)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.group<- pred.SO.stud.phymean[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.group), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muGroupLiving_n.studies.c"]), 2))
mtext("S)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Variable)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.var<- pred.SO.stud.phymean[,,5]
for (i in sample(1:nrow(pred.var), 100)){
	lines(pred.var[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[5], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.var), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muVariable_n.studies.c"]), 2))
mtext("T)", adj=-0.5, line=0)


## breeding seasonality
plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xlim=c(-0.5,1.5), xaxt="n", ylab="Pr (Solitary)", xlab="Breeding seasonality")
axis(1, at=c(0,1), labels=c("Non-seas.", "Seas."))
points(0, mean(pred.SO.nonseas.phymean[,,1]), col=wes_palettes$Zissou[1], pch=16, cex=2)
arrows(0, HPDI(pred.SO.nonseas.phymean[,,1], prob=0.95)[1], 0, HPDI(pred.SO.nonseas.phymean[,,1], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[1])
points(1, mean(pred.SO.seas.phymean[,,1]), col=wes_palettes$Zissou[1], pch=16, cex=2)
arrows(1, HPDI(pred.SO.seas.phymean[,,1], prob=0.95)[1], 1, HPDI(pred.SO.seas.phymean[,,1], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[1])
mtext("U)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xlim=c(-0.5,1.5), xaxt="n", ylab="Pr (Pair-living)", xlab="Breeding seasonality")
axis(1, at=c(0,1), labels=c("Non-seas.", "Seas."))
points(0, mean(pred.SO.nonseas.anc[,,2]), col=wes_palettes$Zissou[2], pch=16, cex=2)
arrows(0, HPDI(pred.SO.nonseas.anc[,,2], prob=0.95)[1], 0, HPDI(pred.SO.nonseas.anc[,,2], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[2])
points(1, mean(pred.SO.seas.anc[,,2]), col=wes_palettes$Zissou[2], pch=16, cex=2)
arrows(1, HPDI(pred.SO.seas.anc[,,2], prob=0.95)[1], 1, HPDI(pred.SO.seas.anc[,,2], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[2])
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muPairLiving_SeasonBreedSeasonal"]), 2))
mtext("V)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xlim=c(-0.5,1.5), xaxt="n", ylab="Pr (Sex-specific)", xlab="Breeding seasonality")
axis(1, at=c(0,1), labels=c("Non-seas.", "Seas."))
points(0, mean(pred.SO.nonseas.phymean[,,3]), col=wes_palettes$Zissou[3], pch=16, cex=2)
arrows(0, HPDI(pred.SO.nonseas.phymean[,,3], prob=0.95)[1], 0, HPDI(pred.SO.nonseas.phymean[,,3], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[3])
points(1, mean(pred.SO.seas.phymean[,,3]), col=wes_palettes$Zissou[3], pch=16, cex=2)
arrows(1, HPDI(pred.SO.seas.phymean[,,3], prob=0.95)[1], 1, HPDI(pred.SO.seas.phymean[,,3], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[3])
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muSexSpecific_SeasonBreedSeasonal"]), 2))
mtext("W)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xlim=c(-0.5,1.5), xaxt="n", ylab="Pr (Group-living)", xlab="Breeding seasonality")
axis(1, at=c(0,1), labels=c("Non-seas.", "Seas."))
points(0, mean(pred.SO.nonseas.phymean[,,4]), col=wes_palettes$Zissou[4], pch=16, cex=2)
arrows(0, HPDI(pred.SO.nonseas.phymean[,,4], prob=0.95)[1], 0, HPDI(pred.SO.nonseas.phymean[,,4], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[4])
points(1, mean(pred.SO.seas.phymean[,,4]), col=wes_palettes$Zissou[4], pch=16, cex=2)
arrows(1, HPDI(pred.SO.seas.phymean[,,4], prob=0.95)[1], 1, HPDI(pred.SO.seas.phymean[,,4], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[4])
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muGroupLiving_SeasonBreedSeasonal"]), 2))
mtext("X)", adj=-0.5, line=0)

plot(logBM.s.foss~n.studies.c, col="white", data=d, ylim=c(0,1), xlim=c(-0.5,1.5), xaxt="n", ylab="Pr (Variable)", xlab="Breeding seasonality")
axis(1, at=c(0,1), labels=c("Non-seas.", "Seas."))
points(0, mean(pred.SO.nonseas.phymean[,,5]), col=wes_palettes$Zissou[5], pch=16, cex=2)
arrows(0, HPDI(pred.SO.nonseas.phymean[,,5], prob=0.95)[1], 0, HPDI(pred.SO.nonseas.phymean[,,5], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[5])
points(1, mean(pred.SO.seas.phymean[,,5]), col=wes_palettes$Zissou[5], pch=16, cex=2)
arrows(1, HPDI(pred.SO.seas.phymean[,,5], prob=0.95)[1], 1, HPDI(pred.SO.seas.phymean[,,5], prob=0.95)[2], length=0, angle=90, col=wes_palettes$Zissou[5])
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m1.imp.PP[names(m1.imp.PP)=="b_muVariable_SeasonBreedSeasonal"]), 2))
mtext("Y)", adj=-0.5, line=0)

dev.off()





	## Model 2 ##
	#############

# prepare DV
colnames(d)[19]<- "SO.new"
d$SO.new<- relevel(as.factor(d$SO.new), ref="Solitary")

# use multiple imputation because brms cannot impute categorical predictors (seasonal-breeding); see https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html for details
# include all available body size data (e.g. sometimes avge mass and dimorphism are known, but not female mass
d.imp2<- d[,c("SO.new", "Dimorphic/Monomorphic", "dimorph.s.foss", "logBM.s.foss", "Male_BM (kg)", "F_BM (kg)", "M:F_BM", "Avg_BM (kg)", "n.habitats.c", "SeasonBreed", "n.studies.c", "phylo", "Genus_species", "Continent", "Habitat")]
colnames(d.imp2)[2]<- "dimorphcat"
colnames(d.imp2)[5]<- "maleBM"
colnames(d.imp2)[6]<- "femaleBM"
colnames(d.imp2)[7]<- "dimorph"
colnames(d.imp2)[8]<- "avgeBM"
# convert categorical variables to factors, otherwise they won't be imputed
d.imp2$dimorphcat<- as.factor(d.imp2$dimorphcat)
d.imp2$SeasonBreed<- as.factor(d.imp2$SeasonBreed)
d.imp2$Habitat<- as.factor(d.imp2$Habitat)

d.imp2<- mice(d.imp2, m=10)

# extract dfs
list_nos <- seq(1:10)
fun_extract_dfs <- function(x) {
  complete(d.imp2, action = x) 
  }
list_d.imp2 <- lapply(list_nos, fun_extract_dfs)


prior=c(prior(normal(0,50), class=Intercept), 
prior(cauchy(0,2), class=sd, dpar="muNonVariableGroup"), 
prior(cauchy(0,2), class=sd, dpar="muPairGroup"), 
prior(cauchy(0,2), class=sd, dpar="muPairLiving"), 
prior(cauchy(0,2), class=sd, dpar="muSexSpecific"), 
prior(cauchy(0,2), class=sd, dpar="muSolGroup"), 
prior(cauchy(0,2), class=sd, dpar="muSolPair"), 
prior(cauchy(0,2), class=sd, dpar="muSolPairGroup"), 
prior(cauchy(0,2), class=sd, dpar="muVariableGroup"), 
prior(normal(0,5), class=b))

m2.imp<- brm_multiple(SO.new~dimorph.s.foss+logBM.s.foss+n.studies.c+(1|phylo)+(1|Genus_species), 
	data=list_d.imp2, family=categorical(), cov_ranef = list(phylo = A), prior = prior, 
	chains = 2, cores = 2, iter = 4000, warmup = 1000, thin=20, control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(m2.imp, "m2.imp.rds")
	
summary(m2.imp)
plot(m2.imp)

# rhat's high, this is expected because chains fit to slightly different data; plots look good though
# check rhat's within each dataset:
round(m2.imp$rhats, 2) # --> those look good!
max(apply(m2.imp$rhats, 1, max)) # 1.037 is the highest Rhat


postm2.imp<- posterior_samples(m2.imp)

# phylogenetic signal: sum of all SO variances / sum of all variance components + distribution-specific variance
VarPhy.group<- postm2.imp$sd_phylo__muNonVariableGroup_Intercept^2
VarPhy.pairgroup<- postm2.imp$sd_phylo__muPairGroup_Intercept^2
VarPhy.pair<- postm2.imp$sd_phylo__muPairLiving_Intercept^2
VarPhy.ss<- postm2.imp$sd_phylo__muSexSpecific_Intercept^2
VarPhy.solgroup<- postm2.imp$sd_phylo__muSolGroup_Intercept^2
VarPhy.solpair<- postm2.imp$sd_phylo__muSolPair_Intercept^2
VarPhy.solpairgroup<- postm2.imp$sd_phylo__muSolPairGroup_Intercept^2
VarPhy.var<- postm2.imp$sd_phylo__muVariableGroup_Intercept^2

VarSpec.group<- postm2.imp$sd_Genus_species__muNonVariableGroup_Intercept^2
VarSpec.pairgroup<- postm2.imp$sd_Genus_species__muPairGroup_Intercept^2
VarSpec.pair<- postm2.imp$sd_Genus_species__muPairLiving_Intercept^2
VarSpec.ss<- postm2.imp$sd_Genus_species__muSexSpecific_Intercept^2
VarSpec.solgroup<- postm2.imp$sd_Genus_species__muSolGroup_Intercept^2
VarSpec.solpair<- postm2.imp$sd_Genus_species__muSolPair_Intercept^2
VarSpec.solpairgroup<- postm2.imp$sd_Genus_species__muSolPairGroup_Intercept^2
VarSpec.var<- postm2.imp$sd_Genus_species__muVariableGroup_Intercept^2

VarDistro<- pi^2 / 3 # see Nakagawa & Schielzeth 2012 MEE

lambda<- (VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss+VarPhy.pairgroup+VarPhy.solgroup+VarPhy.solpair+VarPhy.solpairgroup)/
	(VarPhy.var+VarPhy.group+VarPhy.pair+VarPhy.ss+VarPhy.pairgroup+VarPhy.solgroup+VarPhy.solpair+VarPhy.solpairgroup+
	VarSpec.var+VarSpec.group+VarSpec.pair+VarSpec.ss+VarSpec.pairgroup+VarSpec.solgroup+VarSpec.solpair+VarSpec.solpairgroup)
mean(lambda); HPDI(lambda, prob=0.95); sum(lambda>0.01)/length(lambda)
#[1] 0.05964275
#      |0.95       0.95| 
#0.001059155 0.200408235
# [1] 0.8483333


# probability of each state at the root (global intercept) -> code adapted from Koster & McElreath 2017 Behav Ecol Sociobiol Multinomial analysis of behavior: statistical methods
{    K <- 9
    ns <- nrow(postm2.imp)
	n <- 1

    softmax2 <- function(x) {
        x <- max(x) - x
        exp(-x)/sum(exp(-x))
    }

    p <- list()

    for ( i in 1:n ) {
        p[[i]] <- sapply( 1:K , function(k) {
            if ( k < K ) {
                ptemp <- postm2.imp[,k] 
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
p_mean <- sapply( 1:length(p) , function(i) apply(p[[i]],2,mean) )
p_HPDI <- sapply( 1:length(p) , function(i) apply(p[[i]],2,HPDI, prob=0.95) )

pred_probs_m2.imp<- cbind(p_mean, p_HPDI[c(1,3,5,7,9,11,13,15,17),], p_HPDI[c(2,4,6,8,10,12,14,16,18),])
rownames(pred_probs_m2.imp)<- c("Group-living", "Pair-Group", "Pair-living", "Sex-specific", "Solitary-Group", "Solitary-Pair", "Solitary-Pair-Group", "Variable Group", "Solitary")
colnames(pred_probs_m2.imp)<- c("mean", "lwr95", "upr95")
round(pred_probs_m2.imp,2)
#                    mean lwr95 upr95
#Group-living        0.00     0  0.00
#Pair-Group          0.04     0  0.20
#Pair-living         0.54     0  1.00
#Sex-specific        0.01     0  0.00
#Solitary-Group      0.00     0  0.00
#Solitary-Pair       0.39     0  1.00
#Solitary-Pair-Group 0.02     0  0.06
#Variable Group      0.00     0  0.00
#Solitary            0.00     0  0.00


# pp's
m2.imp.PP<- 0
for(i in 1:32){
	ifelse(mean(postm2.imp[,i])>0, m2.imp.PP[i]<- sum(postm2.imp[,i]>0)/nrow(postm2.imp), m2.imp.PP[i]<- sum(postm2.imp[,i]<0)/nrow(postm2.imp))
	}
names(m2.imp.PP)<- colnames(postm2.imp[,1:32])
m2.imp.PP
## strong (PP>0.9) associations for GL-body size, GL-n studies (0.89), Pairgroup-body size (0.89), pairgroup-n studies, PL-body size, sexspecific-n studies, solgroup-dimorph,
# solgroup-body size, solgroup-n studies, solpair-body size, solpairgroup-body size, solpairgroup-n studies, variable-dimorphism, variable n studies


# probability of each state at the phylogenetic mean of predictors (log body size = 4.187368, dimorphism = 1.214002) -> using brms:fitted() this time
SeasonBreed<- "Non-seasonal"
newdata<- data.frame(SeasonBreed)
newdata$logBM.s.foss<- (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))
newdata$dimorph.s.foss<- (as.numeric(ace.dimo.ace$ace[1])-1)/sd(d[,"M:F_BM"], na.rm=TRUE)
newdata$n.studies.c<- 0
pred.phymean<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)
pred.phymean.sol.mean<- mean(pred.phymean[,,1])
pred.phymean.sol.HPDI<- HPDI(pred.phymean[,,1], prob=0.95)
pred.phymean.group.mean<- mean(pred.phymean[,,2])
pred.phymean.group.HPDI<- HPDI(pred.phymean[,,2], prob=0.95)
pred.phymean.pairgroup.mean<- mean(pred.phymean[,,3])
pred.phymean.pairgroup.HPDI<- HPDI(pred.phymean[,,3], prob=0.95)
pred.phymean.pair.mean<- mean(pred.phymean[,,4])
pred.phymean.pair.HPDI<- HPDI(pred.phymean[,,4], prob=0.95)
pred.phymean.ss.mean<- mean(pred.phymean[,,5])
pred.phymean.ss.HPDI<- HPDI(pred.phymean[,,5], prob=0.95)
pred.phymean.solgroup.mean<- mean(pred.phymean[,,6])
pred.phymean.solgroup.HPDI<- HPDI(pred.phymean[,,6], prob=0.95)
pred.phymean.solpair.mean<- mean(pred.phymean[,,7])
pred.phymean.solpair.HPDI<- HPDI(pred.phymean[,,7], prob=0.95)
pred.phymean.solpairgroup.mean<- mean(pred.phymean[,,8])
pred.phymean.solpairgroup.HPDI<- HPDI(pred.phymean[,,8], prob=0.95)
pred.phymean.vargroup.mean<- mean(pred.phymean[,,9])
pred.phymean.vargroup.HPDI<- HPDI(pred.phymean[,,9], prob=0.95)

pred.phymean.means<- rbind(pred.phymean.sol.mean, pred.phymean.group.mean, pred.phymean.pairgroup.mean, pred.phymean.pair.mean, pred.phymean.ss.mean, 
pred.phymean.solgroup.mean, pred.phymean.solpair.mean, pred.phymean.solpairgroup.mean, pred.phymean.vargroup.mean)

pred.phymean.HPDIs<- rbind(pred.phymean.sol.HPDI, pred.phymean.group.HPDI, pred.phymean.pairgroup.HPDI, pred.phymean.pair.HPDI, pred.phymean.ss.HPDI, 
pred.phymean.solgroup.HPDI, pred.phymean.solpair.HPDI, pred.phymean.solpairgroup.HPDI, pred.phymean.vargroup.HPDI)

round(cbind(pred.phymean.means, pred.phymean.HPDIs), 2)
#                                    |0.95 0.95|
#pred.phymean.sol.mean          0.25     0  0.50
#pred.phymean.group.mean        0.27     0  0.80
#pred.phymean.pairgroup.mean    0.01     0  0.04
#pred.phymean.pair.mean         0.00     0  0.00
#pred.phymean.ss.mean           0.01     0  0.02
#pred.phymean.solgroup.mean     0.16     0  0.38
#pred.phymean.solpair.mean      0.01     0  0.03
#pred.phymean.solpairgroup.mean 0.08     0  0.26
#pred.phymean.vargroup.mean     0.21     0  0.64


## Figure 2: Phylogeny with social organization at root and phylo mean
# save phylo to file and read again to get correct order of tip labels
write.tree(mamm.prun.di.min, "PhyForFigure.tre")
phy<- read.tree("PhyForFigure.tre")
phy<- drop.tip(phy, setdiff(phy$tip.label,d$Genus_species))


# create vector of SO per species
t<- table(d$Genus_species, d$SO.new)
t2<- t[ order(match(rownames(t), phy$tip.label)), ]

tiff("Figure S1.tif", compression="lzw", height=6.0, width=6.0, units="cm", res=600, pointsize=5)
par(mar = c(2, 1, 1, 0), xpd=NA)
plot(phy, cex=0.5, lwd=1.5, label.offset=9)

points(x=71, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[1])
points(x=72, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[2])
points(x=73, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[3])
points(x=74, y=100, pch=15, cex=0.7, col=wes_palettes$Zissou[4])
points(x=75, y=100, pch=15, cex=0.7, col=wes_palettes$Royal2[2])
points(x=76, y=100, pch=15, cex=0.7, col=wes_palettes$Royal2[3])
points(x=77, y=100, pch=15, cex=0.7, col=wes_palettes$Royal2[4])
points(x=78, y=100, pch=15, cex=0.7, col=wes_palettes$Royal2[5])
points(x=79, y=100, pch=15, cex=0.7, col=wes_palettes$Royal2[1])

points(x=71, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[1])
points(x=72, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[2])
points(x=73, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[3])
points(x=74, y=0, pch=15, cex=0.7, col=wes_palettes$Zissou[4])
points(x=75, y=0, pch=15, cex=0.7, col=wes_palettes$Royal2[2])
points(x=76, y=0, pch=15, cex=0.7, col=wes_palettes$Royal2[3])
points(x=77, y=0, pch=15, cex=0.7, col=wes_palettes$Royal2[4])
points(x=78, y=0, pch=15, cex=0.7, col=wes_palettes$Royal2[5])
points(x=79, y=0, pch=15, cex=0.7, col=wes_palettes$Royal2[1])

for(i in 1:nrow(t2)){if(t2[i,"Solitary"]>0) points(x=71, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[1])}
for(i in 1:nrow(t2)){if(t2[i,"PairLiving"]>0) points(x=72, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[2])}
for(i in 1:nrow(t2)){if(t2[i,"SexSpecific"]>0) points(x=73, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[3])}
for(i in 1:nrow(t2)){if(t2[i,"NonVariableGroup"]>0) points(x=74, y=0+i, pch=15, cex=0.7, col=wes_palettes$Zissou[4])}
for(i in 1:nrow(t2)){if(t2[i,"PairGroup"]>0) points(x=75, y=0+i, pch=15, cex=0.7, col=wes_palettes$Royal2[2])}
for(i in 1:nrow(t2)){if(t2[i,"SolGroup"]>0) points(x=76, y=0+i, pch=15, cex=0.7, col=wes_palettes$Royal2[3])}
for(i in 1:nrow(t2)){if(t2[i,"SolPair"]>0) points(x=77, y=0+i, pch=15, cex=0.7, col=wes_palettes$Royal2[4])}
for(i in 1:nrow(t2)){if(t2[i,"SolPairGroup"]>0) points(x=78, y=0+i, pch=15, cex=0.7, col=wes_palettes$Royal2[5])}
for(i in 1:nrow(t2)){if(t2[i,"VariableGroup"]>0) points(x=79, y=0+i, pch=15, cex=0.7, col=wes_palettes$Royal2[1])}

## create inset plots
par(fig = c(0.05,0.4, 0.45, 0.75), mar = c(0, 2, 3, 2), new = T)

# SO probabilities at the root
plot(c(0,0.05,0.1,0.3,0.5,0.75,0.8,0.9,1)~c(1,2,3,4,5,6,7,8,9), col="white", xaxt="n", ylab="", xlab="", main="Ancestral state")
axis(1, at=c(1,2,3,4,5,6,7,8,9), labels=FALSE) 

points(1, pred_probs_m2.imp[9,1], pch=16, col=wes_palettes$Zissou[1], cex=2)
arrows(1, pred_probs_m2.imp[9,2], 1, pred_probs_m2.imp[9,3], length=0, angle=90, col=wes_palettes$Zissou[1], lwd=2)
points(2, pred_probs_m2.imp[3,1], pch=16, col=wes_palettes$Zissou[2], cex=2)
arrows(2, pred_probs_m2.imp[3,2], 2, pred_probs_m2.imp[3,3], length=0, angle=90, col=wes_palettes$Zissou[2], lwd=2)
points(3, pred_probs_m2.imp[4,1], pch=16, col=wes_palettes$Zissou[3], cex=2)
arrows(3, pred_probs_m2.imp[4,2], 3, pred_probs_m2.imp[4,3], length=0, angle=90, col=wes_palettes$Zissou[3], lwd=2)
points(4, pred_probs_m2.imp[1,1], pch=16, col=wes_palettes$Zissou[4], cex=2)
arrows(4, pred_probs_m2.imp[1,2], 4, pred_probs_m2.imp[1,3], length=0, angle=90, col=wes_palettes$Zissou[4], lwd=2)
points(5, pred_probs_m2.imp[2,1], pch=16, col=wes_palettes$Royal2[2], cex=2)
arrows(5, pred_probs_m2.imp[2,2], 5, pred_probs_m2.imp[2,3], length=0, angle=90, col=wes_palettes$Royal2[2], lwd=2)
points(6, pred_probs_m2.imp[5,1], pch=16, col=wes_palettes$Royal2[3], cex=2)
arrows(6, pred_probs_m2.imp[5,2], 6, pred_probs_m2.imp[5,3], length=0, angle=90, col=wes_palettes$Royal2[3], lwd=2)
points(7, pred_probs_m2.imp[6,1], pch=16, col=wes_palettes$Royal2[4], cex=2)
arrows(7, pred_probs_m2.imp[6,2], 7, pred_probs_m2.imp[6,3], length=0, angle=90, col=wes_palettes$Royal2[4], lwd=2)
points(8, pred_probs_m2.imp[7,1], pch=16, col=wes_palettes$Royal2[5], cex=2)
arrows(8, pred_probs_m2.imp[7,2], 8, pred_probs_m2.imp[7,3], length=0, angle=90, col=wes_palettes$Royal2[5], lwd=2)
points(9, pred_probs_m2.imp[8,1], pch=16, col=wes_palettes$Royal2[1], cex=2)
arrows(9, pred_probs_m2.imp[8,2], 9, pred_probs_m2.imp[8,3], length=0, angle=90, col=wes_palettes$Royal2[1], lwd=2)

par(fig = c(0.05,0.4, 0, 0.45), mar = c(4, 2, 3, 2), new = T)

plot(c(0,0.05,0.1,0.3,0.5,0.75,0.8,0.9,1)~c(1,2,3,4,5,6,7,8,9), col="white", xaxt="n", ylab="", xlab="", main="Phylogenetic mean")
axis(1, at=c(1,2,3,4,5,6,7,8,9), labels=FALSE)
text(x=c(1,2,3,4,5,6,7,8,9), y=-0.5, xpd=TRUE, adj=0, srt=90, labels=c("Sol", "Pair", "SS", "Grp", "PrGp", "SlGp", "SlPr", "SPG", "VaGp"))

points(1, pred.phymean.sol.mean, pch=16, col=wes_palettes$Zissou[1], cex=2)
arrows(1, pred.phymean.sol.HPDI[1], 1, pred.phymean.sol.HPDI[2], length=0, angle=90, col=wes_palettes$Zissou[1], lwd=2)
points(2, pred.phymean.pair.mean, pch=16, col=wes_palettes$Zissou[2], cex=2)
arrows(2, pred.phymean.pair.HPDI[1], 2, pred.phymean.pair.HPDI[2], length=0, angle=90, col=wes_palettes$Zissou[2], lwd=2)
points(3, pred.phymean.ss.mean, pch=16, col=wes_palettes$Zissou[3], cex=2)
arrows(3, pred.phymean.ss.HPDI[1], 3, pred.phymean.ss.HPDI[2], length=0, angle=90, col=wes_palettes$Zissou[3], lwd=2)
points(4, pred.phymean.group.mean, pch=16, col=wes_palettes$Zissou[4], cex=2)
arrows(4, pred.phymean.group.HPDI[1], 4, pred.phymean.group.HPDI[2], length=0, angle=90, col=wes_palettes$Zissou[4], lwd=2)
points(5, pred.phymean.pairgroup.mean, pch=16, col=wes_palettes$Royal2[2], cex=2)
arrows(5, pred.phymean.pairgroup.HPDI[1], 5, pred.phymean.pairgroup.HPDI[2], length=0, angle=90, col=wes_palettes$Royal2[2], lwd=2)
points(6, pred.phymean.solgroup.mean, pch=16, col=wes_palettes$Royal2[3], cex=2)
arrows(6, pred.phymean.solgroup.HPDI[1], 6, pred.phymean.solgroup.HPDI[2], length=0, angle=90, col=wes_palettes$Royal2[3], lwd=2)
points(7, pred.phymean.solpair.mean, pch=16, col=wes_palettes$Royal2[4], cex=2)
arrows(7, pred.phymean.solpair.HPDI[1], 7, pred.phymean.solpair.HPDI[2], length=0, angle=90, col=wes_palettes$Royal2[4], lwd=2)
points(8, pred.phymean.solpairgroup.mean, pch=16, col=wes_palettes$Royal2[5], cex=2)
arrows(8, pred.phymean.solpairgroup.HPDI[1], 8, pred.phymean.solpairgroup.HPDI[2], length=0, angle=90, col=wes_palettes$Royal2[5], lwd=2)
points(9, pred.phymean.vargroup.mean, pch=16, col=wes_palettes$Royal2[1], cex=2)
arrows(9, pred.phymean.vargroup.HPDI[1], 9, pred.phymean.vargroup.HPDI[2], length=0, angle=90, col=wes_palettes$Royal2[1], lwd=2)

par(fig = c(0,1,0,1), new=FALSE)
dev.off()



## Figure S2: Transitions with all predictors
## predictions for body size - everything else at baseline
logBM.s.foss<- seq(0,max(d$logBM.s.foss, na.rm=TRUE),length=100)
newdata<- as.data.frame(logBM.s.foss)
newdata$dimorph.s.foss<- 0
newdata$n.studies.c<- 0
pred.SO.bm<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)

## predictions for dimorphism - for pair-living, pair-group, and sol-pair body size at baseline, for others at phylo mean
dimorph.s.foss<- seq(min(d$dimorph.s.foss, na.rm=TRUE),max(d$dimorph.s.foss, na.rm=TRUE),length=100)
newdata<- as.data.frame(dimorph.s.foss)
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
newdata$n.studies.c<- 0
pred.SO.dm.anc<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.dm.phymean<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)

## predictions for n studies - for pair-living, pair-group, and sol-pair body size at baseline, for others at phylo mean
n.studies.c<- seq(0,max(d$n.studies.c, na.rm=TRUE),length=100)
newdata<- as.data.frame(n.studies.c)
newdata$dimorph.s.foss<- 0
newdata$logBM.s.foss<- 0 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.stud.anc<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)

newdata$logBM.s.foss<- 4.478611 # at ancestral state - for phy mean use (as.numeric(ace.bm.ace$ace[1])-log(0.737)/sd(log(d[,"F_BM (kg)"]), na.rm=TRUE))=4.478611
pred.SO.stud.phymean<- fitted(m2.imp, newdata, re_formula=NA, summary=FALSE)


## plot
tiff("Figure S2.tif", compression="lzw", height=6.0, width=18.0, units="cm", res=600, pointsize=5)
par(mfrow=c(3,9), mar = c(4, 4, 2, 2))

## body size
plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Solitary)", xlab="Female body size [Z]")
pred.sol<- pred.SO.bm[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.sol), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
mtext("A)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Pair-living)", xlab="Female body size [Z]")
pred.pair<- pred.SO.bm[,,4]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.pair), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairLiving_logBM.s.foss"]), 2))
mtext("B)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Sex-specific)", xlab="Female body size [Z]")
pred.ss<- pred.SO.bm[,,5]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.ss), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSexSpecific_logBM.s.foss"]), 2))
mtext("C)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Group-living)", xlab="Female body size [Z]")
pred.group<- pred.SO.bm[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.group), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muNonVariableGroup_logBM.s.foss"]), 2))
mtext("D)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Pair-Group)", xlab="Female body size [Z]")
pred.pairgroup<- pred.SO.bm[,,3]
for (i in sample(1:nrow(pred.pairgroup), 100)){
	lines(pred.pairgroup[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Royal2[2], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.pairgroup), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairGroup_logBM.s.foss"]), 2))
mtext("E)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Solitary-Group)", xlab="Female body size [Z]")
pred.solgroup<- pred.SO.bm[,,6]
for (i in sample(1:nrow(pred.solgroup), 100)){
	lines(pred.solgroup[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Royal2[3], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.solgroup), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolGroup_logBM.s.foss"]), 2))
mtext("F)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Solitary-Pair)", xlab="Female body size [Z]")
pred.solpair<- pred.SO.bm[,,7]
for (i in sample(1:nrow(pred.solpair), 100)){
	lines(pred.solpair[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Royal2[4], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.solpair), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPair_logBM.s.foss"]), 2))
mtext("G)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Sol-Pair-Group)", xlab="Female body size [Z]")
pred.solpairgroup<- pred.SO.bm[,,8]
for (i in sample(1:nrow(pred.solpairgroup), 100)){
	lines(pred.solpairgroup[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Royal2[5], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.solpairgroup), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPairGroup_logBM.s.foss"]), 2))
mtext("H)", adj=-0.5, line=0)

plot(n.habitats.c~logBM.s.foss, col="white", data=d, ylim=c(0,1), xlim=c(0,7), ylab="Pr (Variable Group)", xlab="Female body size [Z]")
pred.vargroup<- pred.SO.bm[,,9]
for (i in sample(1:nrow(pred.vargroup), 100)){
	lines(pred.vargroup[i,]~logBM.s.foss, col=col.alpha(wes_palettes$Royal2[1], 0.2), lwd=0.5)
	}
lines(logBM.s.foss, colMeans(pred.vargroup), lwd=1)
abline(v=1.885, lty=2, lwd=0.5)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muVariableGroup_logBM.s.foss"]), 2))
mtext("I)", adj=-0.5, line=0)


## dimorphism
plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Solitary)", xlab="Sexual dimorphism [Z]")
pred.sol<- pred.SO.dm.phymean[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.sol), lwd=1)
mtext("J)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Pair-living)", xlab="Sexual dimorphism [Z]")
pred.pair<- pred.SO.dm.anc[,,4]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.pair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairLiving_dimorph.s.foss"]), 2))
mtext("K)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Sex-specific)", xlab="Sexual dimorphism [Z]")
pred.ss<- pred.SO.dm.phymean[,,5]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.ss), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSexSpecific_dimorph.s.foss"]), 2))
mtext("L)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Group-living)", xlab="Sexual dimorphism [Z]")
pred.group<- pred.SO.dm.phymean[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.group), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muNonVariableGroup_dimorph.s.foss"]), 2))
mtext("M)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Pair-Group)", xlab="Sexual dimorphism [Z]")
pred.pairgroup<- pred.SO.dm.anc[,,3]
for (i in sample(1:nrow(pred.pairgroup), 100)){
	lines(pred.pairgroup[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Royal2[2], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.pairgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairGroup_dimorph.s.foss"]), 2))
mtext("N)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Solitary-Group)", xlab="Sexual dimorphism [Z]")
pred.solgroup<- pred.SO.dm.phymean[,,6]
for (i in sample(1:nrow(pred.solgroup), 100)){
	lines(pred.solgroup[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Royal2[3], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.solgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolGroup_dimorph.s.foss"]), 2))
mtext("O)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Solitary-Pair)", xlab="Sexual dimorphism [Z]")
pred.solpair<- pred.SO.dm.anc[,,7]
for (i in sample(1:nrow(pred.solpair), 100)){
	lines(pred.solpair[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Royal2[4], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.solpair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPair_dimorph.s.foss"]), 2))
mtext("P)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Sol-Pair-Group)", xlab="Sexual dimorphism [Z]")
pred.solpairgroup<- pred.SO.dm.phymean[,,8]
for (i in sample(1:nrow(pred.solpairgroup), 100)){
	lines(pred.solpairgroup[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Royal2[5], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.solpairgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPairGroup_dimorph.s.foss"]), 2))
mtext("Q)", adj=-0.5, line=0)

plot(n.habitats.c~dimorph.s.foss, col="white", data=d, ylim=c(0,1), ylab="Pr (Variable Group)", xlab="Sexual dimorphism [Z]")
pred.vargroup<- pred.SO.dm.phymean[,,9]
for (i in sample(1:nrow(pred.vargroup), 100)){
	lines(pred.vargroup[i,]~dimorph.s.foss, col=col.alpha(wes_palettes$Royal2[1], 0.2), lwd=0.5)
	}
lines(dimorph.s.foss, colMeans(pred.vargroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muVariableGroup_dimorph.s.foss"]), 2))
mtext("R)", adj=-0.5, line=0)



## n studies
plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Solitary)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.sol<- pred.SO.stud.phymean[,,1]
for (i in sample(1:nrow(pred.sol), 100)){
	lines(pred.sol[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[1], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.sol), lwd=1)
mtext("S)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Pair-living)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.pair<- pred.SO.stud.anc[,,4]
for (i in sample(1:nrow(pred.pair), 100)){
	lines(pred.pair[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[2], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.pair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairLiving_n.studies.c"]), 2))
mtext("T)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Sex-specific)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.ss<- pred.SO.stud.phymean[,,5]
for (i in sample(1:nrow(pred.ss), 100)){
	lines(pred.ss[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[3], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.ss), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSexSpecific_n.studies.c"]), 2))
mtext("U)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Group-living)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.group<- pred.SO.stud.phymean[,,2]
for (i in sample(1:nrow(pred.group), 100)){
	lines(pred.group[i,]~n.studies.c, col=col.alpha(wes_palettes$Zissou[4], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.group), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muNonVariableGroup_n.studies.c"]), 2))
mtext("V)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Pair-Group)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.pairgroup<- pred.SO.stud.anc[,,3]
for (i in sample(1:nrow(pred.pairgroup), 100)){
	lines(pred.pairgroup[i,]~n.studies.c, col=col.alpha(wes_palettes$Royal2[2], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.pairgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muPairGroup_n.studies.c"]), 2))
mtext("W)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Solitary-Group)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.solgroup<- pred.SO.stud.phymean[,,6]
for (i in sample(1:nrow(pred.solgroup), 100)){
	lines(pred.solgroup[i,]~n.studies.c, col=col.alpha(wes_palettes$Royal2[3], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.solgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolGroup_n.studies.c"]), 2))
mtext("X)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Solitary-Pair)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.solpair<- pred.SO.stud.anc[,,7]
for (i in sample(1:nrow(pred.solpair), 100)){
	lines(pred.solpair[i,]~n.studies.c, col=col.alpha(wes_palettes$Royal2[4], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.solpair), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPair_n.studies.c"]), 2))
mtext("Y)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Sol-Pair-Group)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.solpairgroup<- pred.SO.stud.phymean[,,8]
for (i in sample(1:nrow(pred.solpairgroup), 100)){
	lines(pred.solpairgroup[i,]~n.studies.c, col=col.alpha(wes_palettes$Royal2[5], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.solpairgroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muSolPairGroup_n.studies.c"]), 2))
mtext("Z)", adj=-0.5, line=0)

plot(n.habitats.c~n.studies.c, col="white", data=d, ylim=c(0,1), xaxt="n", ylab="Pr (Variable Group)", xlab="Number of studies")
axis(1, at=c(0,1,2,3,4), labels=c("1", "2", "3", "4", "5"))
pred.vargroup<- pred.SO.stud.phymean[,,9]
for (i in sample(1:nrow(pred.vargroup), 100)){
	lines(pred.vargroup[i,]~n.studies.c, col=col.alpha(wes_palettes$Royal2[1], 0.2), lwd=0.5)
	}
lines(n.studies.c, colMeans(pred.vargroup), lwd=1)
legend("topright", x.intersp=0, box.lwd=0.5, legend=round(as.numeric(m2.imp.PP[names(m2.imp.PP)=="b_muVariableGroup_n.studies.c"]), 2))
mtext("Ä)", adj=-0.5, line=0)


dev.off()




