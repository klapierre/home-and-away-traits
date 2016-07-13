library(vegan)
library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)


setwd('C:\\Users\\Kim\\Dropbox\\NutNet\\traits - Firn\\home and away traits')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  



#####################
#####################
#read in data
traits <- read.csv('home_and_away.csv')%>%
  select(site_code, block, plot, year, Family, Taxon, Origin, trt.x, num_leaves, leaf_pct_N, leaf_pct_C, SLA, leaf_ppm_P, leaf_ppm_B, leaf_ppm_Na, leaf_ppm_Ca, leaf_ppm_Fe, leaf_ppm_K, leaf_ppm_Mg, leaf_ppm_Mn, leaf_ppm_Cu, leaf_ppm_Zn, leaf_ppm_Sr, country, region, RAIN_PET, MAT, MAP, N, P, K, Exclose)

# traitsComplete <- traits[complete.cases(traits),]


#####################
#subset out NPK experiment only
traitsNPK <- traits%>%
  filter(Exclose==0)

#check data
boxplot(traitsNPK$leaf_pct_N)
boxplot(traitsNPK$leaf_pct_C)
# boxplot(traitsNPK$leaf_area_mm2)
# boxplot(traitsNPK$leaf_dry_mass_g)
boxplot(traitsNPK$SLA)
boxplot(traitsNPK$leaf_ppm_P)
boxplot(traitsNPK$leaf_ppm_B)
boxplot(traitsNPK$leaf_ppm_Na)
boxplot(traitsNPK$leaf_ppm_Ca)
boxplot(traitsNPK$leaf_ppm_Fe)
boxplot(traitsNPK$leaf_ppm_K)
boxplot(traitsNPK$leaf_ppm_Mg)
boxplot(traitsNPK$leaf_ppm_Mn)
boxplot(traitsNPK$leaf_ppm_Cu)
boxplot(traitsNPK$leaf_ppm_Zn)
boxplot(traitsNPK$leaf_ppm_Sr)

#generate table of predictor variables
designNPK <- traitsNPK%>%
  select(site_code, N, P, K, Origin)

#generate traits table
traitsNPKdata <- traitsNPK%>%
  select(-site_code:-num_leaves, -country:-Exclose)

rdaModel <- rda(traitsNPKdata ~ Origin, designNPK, na=na.omit, subset=complete.cases(traitsNPKdata))

anova(rdaModel, by="term", permu=999)
anova(rdaModel, by="axis")
plot(rdaModel, scaling=3, type="text")
screeplot(rdaModel, type="lines")
summary(rdaModel)


###only look at spp that cross treatments (Anthox. and Alopec.)
traitsAllTrt <- subset(traitsNPK, subset=(Taxon=='ALOPECURUS PRATENSIS'|Taxon=='ANTHOXANTHUM ODORATUM'))

#generate table of predictor variables
designAllTrt <- traitsAllTrt%>%
  select(site_code, Origin, Taxon, N, P, K)

#generate traits table
traitsAllTrtData <- traitsAllTrt%>%
  select(-site_code:-num_leaves, -country:-Exclose, -leaf_ppm_B:-leaf_ppm_Fe, -leaf_ppm_Mg:-leaf_ppm_Sr)

rdaModel <- rda(traitsAllTrtData ~ Origin*Taxon*N*P*K, designAllTrt, na=na.omit, subset=complete.cases(traitsAllTrtData))

anova(rdaModel, by="term", permu=999)
anova(rdaModel, by="axis")
plot(rdaModel, scaling=3, type="text")
screeplot(rdaModel, type="lines")
summary(rdaModel)


#SLA figures
aloSLA <- ggplot(data=barGraphStats(data=subset(traitsAllTrt, Taxon=='ALOPECURUS PRATENSIS'), variable='SLA', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Specific Leaf Area', breaks=seq(0,20000,5000)) +
  scale_x_discrete(name='', limits=c('Control', 'N', 'P', 'K', 'NP', 'NK', 'PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=20000, label='(a) A. pratensis', size=10, hjust='left')

antSLA <- ggplot(data=barGraphStats(data=subset(traitsAllTrt, Taxon=='ANTHOXANTHUM ODORATUM'), variable='SLA', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Specific Leaf Area', breaks=seq(0,2500,500)) +
  scale_x_discrete(name='', limits=c('Control', 'N', 'P', 'K', 'NP', 'NK', 'PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=2500, label='(b) A. odoratum', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,1)))
print(aloSLA, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(antSLA, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1500 x 2000



#tissue N figures
aloN <- ggplot(data=barGraphStats(data=subset(traitsAllTrt, Taxon=='ALOPECURUS PRATENSIS'), variable='leaf_pct_N', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Leaf Percent N', breaks=seq(0,6,2)) +
  scale_x_discrete(name='', limits=c('Control', 'N', 'P', 'K', 'NP', 'NK', 'PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=6, label='(a) A. pratensis', size=10, hjust='left')

antN <- ggplot(data=barGraphStats(data=subset(traitsAllTrt, Taxon=='ANTHOXANTHUM ODORATUM'), variable='leaf_pct_N', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Leaf Percent N', breaks=seq(0,6,2)) +
  scale_x_discrete(name='', limits=c('Control', 'N', 'P', 'K', 'NP', 'NK', 'PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=6, label='(b) A. odoratum', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,1)))
print(aloN, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(antN, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1500 x 2000





###check out poa -- doesn't change in abundance home and away
traitsPoa <- subset(traitsNPK, subset=(Taxon=='POA PRATENSIS'))

poaSLA <- ggplot(data=barGraphStats(data=traitsPoa, variable='SLA', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Specific Leaf Area', breaks=seq(0,30000,5000)) +
  scale_x_discrete(name='', limits=c('PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=30000, label='(a)', size=10, hjust='left')

poaN <- ggplot(data=barGraphStats(data=traitsPoa, variable='leaf_pct_N', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Leaf Percent N', breaks=seq(0,5,1)) +
  scale_x_discrete(name='', limits=c('PK', 'NPK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=5, label='(b)', size=10, hjust='left') +
  theme(axis.title.y=element_text(margin=margin(r=50)))
  
pushViewport(viewport(layout=grid.layout(2,1)))
print(poaSLA, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(poaN, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1000 x 1200
  
  

  
###check out holcus -- higher in abundance home than away
traitsHol <- subset(traitsNPK, subset=(Taxon=='HOLCUS LANATUS'))

holSLA <- ggplot(data=barGraphStats(data=traitsHol, variable='SLA', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Specific Leaf Area', breaks=seq(0,30000,5000)) +
  scale_x_discrete(name='', limits=c('N', 'P', 'PK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=30000, label='(a)', size=10, hjust='left')

holN <- ggplot(data=barGraphStats(data=traitsHol, variable='leaf_pct_N', byFactorNames=c('Origin', 'trt.x')), aes(x=trt.x, y=mean, shape=Origin)) +
  geom_point(position=position_dodge(0.3), size=6) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.3)) +
  scale_y_continuous(name='Leaf Percent N', breaks=seq(0,5,1)) +
  scale_x_discrete(name='', limits=c('N', 'P', 'PK')) +
  scale_shape_discrete(labels=c('Introduced', 'Native')) +
  annotate('text', x=0.5, y=5, label='(b)', size=10, hjust='left') +
  theme(axis.title.y=element_text(margin=margin(r=50)))

pushViewport(viewport(layout=grid.layout(2,1)))
print(holSLA, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(holN, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1000 x 1200
  