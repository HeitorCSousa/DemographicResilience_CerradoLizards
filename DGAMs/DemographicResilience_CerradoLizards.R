###########################################
#Demographic resilience of Cerrado lizards#
###########################################
rm(list = ls())

#Set working directory
setwd("~/Documents/GitHub/DemographicResilience_CerradoLizards/DGAMs")

# Load packages -----------------------------------------------------------

library(brms)
library(sjPlot)
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(reshape2)
library(lubridate)
library(tidybayes)
library(viridis)
library(modelr)
library(ggplot2)
library(forecast)
library(mvgam)
library(rgl)
library(plot3D)
library(rstan)
library(psych)
library(GGally)
library(ggfortify)
library(emmeans)
library(R2ucare)
library(tidyverse)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# options(brms.backend = "cmdstanr")


# Read and organize data --------------------------------------------------


res.lh.param.Cn <- readRDS("res_lh_param_Cn.rds")
res.lh.param.Ma <- readRDS("res_lh_param_Ma.rds")
res.lh.param.Ti <- readRDS("res_lh_param_Ti.rds")

Cnigropunctatum.data <- readRDS("Cnigropunctatum_data.rds")
Matticolus.data <- readRDS("Matticolus_data.rds")
Titambere.data <- readRDS("Titambere_data.rds")
dim(Titambere.data$amb)

env.array.df <- function(env.array){
  env.df <- data.frame(plot.int = rep(1:5,170),
                       t = rep(1:170, each=5),
                       date = rep(seq(as.Date('2005-11-01'),
                                   as.Date('2019-12-31'),'month'),each=5),
                       tmed2m = as.numeric(env.array[,,1]),
                       RHmax= as.numeric(env.array[,,2]),
                       sol= as.numeric(env.array[,,3]),
                       tmed0cm= as.numeric(env.array[,,4]),
                       tmin0cm=  as.numeric(env.array[,,5]),
                       precip= as.numeric(env.array[,,6]),
                       perf=  as.numeric(env.array[,,7]),
                       ha_90=  as.numeric(env.array[,,8]),
                       fire=  as.numeric(env.array[,,9]),
                       TSLF=  as.numeric(env.array[,,10]))

  env.df <- env.df[order(env.df$plot, env.df$t),]
  return(env.df)
}

env.df.Cn <- env.array.df(Cnigropunctatum.data$amb)
env.df.Ma <- env.array.df(Matticolus.data$amb)
env.df.Ti <- env.array.df(Titambere.data$amb)
env.df.cerrado.liz <- rbind(env.df.Cn,env.df.Ma,env.df.Ti)

res.cerrado.liz <- rbind(res.lh.param.Cn, 
                         res.lh.param.Ma, 
                         res.lh.param.Ti)

res.env.cerrado.liz <- cbind(res.cerrado.liz, env.df.cerrado.liz)

res.env.cerrado.liz$year <- year(res.env.cerrado.liz$date)
res.env.cerrado.liz$month <- month(res.env.cerrado.liz$date)
res.env.cerrado.liz$plot <- factor(res.env.cerrado.liz$plot, levels = c("C", "Q", "EB", "MB", "LB"), ordered = T)

res.env.cerrado.liz <- res.env.cerrado.liz[res.env.cerrado.liz$t < 170,]

#Check distribution
hist(res.env.cerrado.liz$fst.amp)
hist(res.env.cerrado.liz$fst.att)
hist(res.env.cerrado.liz$recovery.time)

hist(log10(res.env.cerrado.liz$fst.amp))
hist(log10(res.env.cerrado.liz$fst.att))
hist(log10(res.env.cerrado.liz$recovery.time))

descdist(res.env.cerrado.liz$fst.amp)
descdist(log10(res.env.cerrado.liz$fst.amp))

descdist(res.env.cerrado.liz$fst.att)
descdist(log10(res.env.cerrado.liz$fst.att))

descdist(res.env.cerrado.liz$recovery.time)
descdist(log10(res.env.cerrado.liz$recovery.time))

#Log10 compensation and 1 - resistance+1
res.env.cerrado.liz$fst.amp <- log10(res.env.cerrado.liz$fst.amp)

fit.fst.amp <- fitDist(res.env.cerrado.liz$fst.amp
                       , k = 2, type = "realAll", trace = FALSE, try.gamlss = TRUE)
fit.fst.amp$fits
summary(fit.fst.amp)

fit.fst.att <- fitDist(res.env.cerrado.liz$fst.att
                       , k = 2, type = "realAll", trace = FALSE, try.gamlss = TRUE)
fit.fst.att$fits
summary(fit.fst.att)

fit.recov.time <- fitDist(res.env.cerrado.liz$recovery.time
                          , k = 2, type = "realAll", trace = FALSE, try.gamlss = TRUE)
fit.recov.time$fits
summary(fit.recov.time)

summary(res.env.cerrado.liz)


# Check correlations ------------------------------------------------------

#Life-history
cor(res.env.cerrado.liz[,3:12], use="na.or.complete", method="spearman")
pairs.panels(res.env.cerrado.liz[,3:12])
quartz(height = 8, width = 12)
ggpairs(res.env.cerrado.liz[,3:12], 
        lower = list(continuous = wrap("smooth_loess", 
                                       colour = rgb(1,0,0,0.1))),
        diag = list(continuous = wrap("densityDiag", 
                                      fill = rgb(1,0,0,0.1), 
                                      colour = rgb(1,0,0,1))))

#Resilience
cor(res.env.cerrado.liz[,13:16],use="na.or.complete")
pairs.panels(res.env.cerrado.liz[,13:16])
quartz(height = 8, width = 12)
ggpairs(res.env.cerrado.liz[,13:16], 
        lower = list(continuous = wrap("smooth_loess", 
                                       colour = rgb(1,0,0,0.1))),
        diag = list(continuous = wrap("densityDiag", 
                                      fill = rgb(1,0,0,0.1), 
                                      colour = rgb(1,0,0,1))))

#Resilience + Life-history
pairs.panels(res.env.cerrado.liz[,3:16])
quartz(height = 8, width = 12)
ggpairs(res.env.cerrado.liz[,3:16], 
        lower = list(continuous = wrap("smooth_loess", 
                                       colour = rgb(1,0,0,0.1))),
        diag = list(continuous = wrap("densityDiag", 
                                      fill = rgb(1,0,0,0.1), 
                                      colour = rgb(1,0,0,1))))

#Resilience + Environment
quartz(height = 8, width = 12)
ggpairs(res.env.cerrado.liz[,c(1,2,13:17,20:31)], 
        lower = list(continuous = wrap("smooth_loess", 
                                       colour = rgb(1,0,0,0.1))),
        diag = list(continuous = wrap("densityDiag", 
                                      fill = rgb(1,0,0,0.1), 
                                      colour = rgb(1,0,0,1))))


# PCA ---------------------------------------------------------------------

#Resilience components
pca.res<-prcomp(res.env.cerrado.liz[,c(13,14,16)],scale=T, center=T)
quartz(height = 8, width = 10)
autoplot(pca.res, data=res.env.cerrado.liz, colour = "plot.int",shape ="species", label=F,alpha=.4,size=2,
         loadings=T, loadings.label = TRUE,scale=0)+
  scale_colour_gradientn(colours=turbo(5))

#Life history components
pca.lh<-prcomp(na.omit(res.env.cerrado.liz[,c(3:12)]),scale=T,center=T)

quartz(height = 8, width = 10)
autoplot(pca.lh, data=na.omit(res.env.cerrado.liz), colour = "plot.int",shape ="species", label=F,size=1.5,
         loadings=T, loadings.label = TRUE,scale=0)+
  scale_colour_gradientn(colours=turbo(5))

quartz(height = 8, width = 10)
autoplot(pca.lh, data=na.omit(res.env.cerrado.liz), colour = "plot",shape ="species", 
         frame = T, frame.color = "plot",
         label=F,size=1.5,
         loadings=T, loadings.label = TRUE,scale=0)+
  scale_colour_manual(values=turbo(5, alpha = 0.5), name = "Fire severity")+
  scale_fill_manual(values=turbo(5, alpha = 0.5), name = "Fire severity")

quartz(height = 8, width = 10)
autoplot(pca.lh, data=na.omit(res.env.cerrado.liz), colour = "species",
         frame = T, frame.color = "species",
         label=F, size=2,
         loadings=T, loadings.label = TRUE,scale=0)+
  scale_colour_manual(values=alpha.col(c("black", "blue", "brown"), 0.4), name = "Species")+
  scale_fill_manual(values=alpha.col(c("black", "blue", "brown"), 0.4), name = "Species")

#Life history and resilience components
pca.lh.res<-prcomp(na.omit(res.lh.env.cerrado.liz[,c(3:16)]),scale=T, center=T)

quartz(height = 8, width = 10)
autoplot(pca.lh.res, data=res.lh.env.cerrado.liz, colour = "species",
         label=F, frame = T, frame.colour = "species",
         loadings=T, loadings.label = TRUE,scale=0)

quartz(height = 8, width = 10)
  autoplot(pca.lh.res, data=res.lh.env.cerrado.liz, colour = "plot.int",
           label=F, shape = "species", size = 2,
           loadings=T, loadings.label = TRUE,scale=0)+
  scale_colour_gradientn(colours=turbo(5))


# Scale life-history predictors of interest-------------------------------------

hist(res.env.cerrado.liz$life.expect)
hist(scale(log10(res.env.cerrado.liz$life.expect)))
shapiro.test(scale(log10(res.env.cerrado.liz$life.expect)))

hist(res.env.cerrado.liz$gen.time)
hist(scale(log10(res.env.cerrado.liz$gen.time)))
shapiro.test(scale(log10(res.env.cerrado.liz$gen.time)))

hist(res.env.cerrado.liz$repro.value)
hist(scale(log10(res.env.cerrado.liz$repro.value)))
shapiro.test(scale(log10(res.env.cerrado.liz$repro.value)))

hist(res.env.cerrado.liz$semel)
hist(scale(log10(res.env.cerrado.liz$semel)))
shapiro.test(scale(log10(res.env.cerrado.liz$semel)))

hist(res.env.cerrado.liz$shape.surv)
hist(scale(res.env.cerrado.liz$shape.surv))

#Colour palletes
pal.plot <- turbo(5,alpha = .6, begin = 1, end=0,direction=-1)
pal.TSLF <- cividis(length(levels(as.factor(res.env.cerrado.liz$TSLF))),alpha=.5)

quartz(height = 9, width = 6)
ggplot(res.env.cerrado.liz, aes(x = as.factor(plot.int), y = gen.time, colour = plot.int, fill = plot.int, alpha = 0.6)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Generation time") +
  coord_cartesian( clip = "off")

quartz(height = 9, width = 6)
ggplot(res.env.cerrado.liz, aes(x = as.factor(plot.int), y = repro.value, 
                                colour = plot.int, fill = plot.int, alpha = 0.6)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Reproductive output") +
  coord_cartesian( clip = "off")


quartz(height = 9, width = 6)
ggplot(res.env.cerrado.liz, aes(x = as.factor(plot.int), y = semel, 
                                colour = plot.int, fill = plot.int, alpha = 0.6)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Degree of iteroparity") +
  coord_cartesian( clip = "off")

plot(repro.value ~ semel, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Reproductive value", xlab="Degree of iteroparity",bty="n")

plot(gen.time ~ semel, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Generation time", xlab="Degree of iteroparity",bty="n")

plot(gen.time ~ scale(log10(longev)), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Generation time", xlab="Longevity",bty="n")

plot(gen.time ~ life.expect, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Generation time", xlab="Life expectancy",bty="n")

plot(gen.time ~ repro.value, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Generation time", xlab="Reproductive value",bty="n")

quartz(height = 9, width = 6)
ggplot(res.env.cerrado.liz, aes(x = gen.time, y = repro.value, 
                                colour = plot.int)) + 
  geom_point() +
  scale_colour_gradientn(colours=turbo(5, alpha = .5), name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Generation time", y="Reproductive output") +
  coord_cartesian( clip = "off")


#Amplification
#Life-history traits
#pal.TSLF[as.numeric(as.factor(res.env.cerrado.liz$TSLF))]
  
plot(fst.amp ~ gen.time, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation", xlab="Generation time",bty="n")

plot(fst.amp ~ life.expect, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation", xlab="Mean life expectancy",bty="n")

plot(fst.amp ~ semel, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation", xlab="Degree of iteroparity",bty="n")

plot(fst.amp ~ gen.time, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="LB",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="LB"]))+14,
      xlab="Generation time",bty="n")#Change the plot to check if the pattern is the same

plot(fst.amp ~ gen.time, data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],pch=16,
     col=pal.plot[res.env.cerrado.liz$plot.int],xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.amp ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],pch=16,
     col=pal.plot[res.env.cerrado.liz$plot.int],xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.amp ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],pch=16,
     col=pal.plot[res.env.cerrado.liz$plot.int],xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.amp ~ repro.value, data=res.env.cerrado.liz, 
     pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],
     ylab="Compensation", xlab="Reproductive value",bty="n")

plot(fst.amp ~ repro.value, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="C",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="LB"]))+14,
     xlab="Reproductive value",bty="n")#Change the plot to check if the pattern is the same

plot(fst.amp ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same

plot(fst.amp ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same

plot(fst.amp ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same


boxplot(fst.amp ~ species, data=res.env.cerrado.liz,ylab="Compensation",xlab="Species")
boxplot(fst.amp ~ plot.int, data=res.env.cerrado.liz,ylab="Compensation",xlab="Fire regime")
boxplot(fst.amp ~ plot.int*species, data=res.env.cerrado.liz,ylab="Compensation",xlab="Fire regime")
boxplot(fst.amp ~ year, data=res.env.cerrado.liz,ylab="Compensation",xlab="Year")
boxplot(fst.amp ~ month, data=res.env.cerrado.liz,ylab="Compensation",xlab="Month")

ggplot(res.env.cerrado.liz, 
       aes(x = as.factor(plot.int), y = fst.amp)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Compensation") 

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",], 
       aes(x = as.factor(plot.int), y = fst.amp)) + 
  ggdist::stat_halfeye(
    adjust = 5, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Compensation")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",], 
       aes(x = as.factor(plot.int), y = fst.amp)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Compensation")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",], 
       aes(x = as.factor(plot.int), y = fst.amp)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Compensation")

#DA476AFF
#Weather and microclimate
plot(fst.amp ~ tmed2m, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.amp ~ RHmax, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.amp ~ sol, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.amp ~ tmed0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.amp ~ tmin0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.amp ~ precip, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Fire
plot(fst.amp ~ jitter(plot.int), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation",xlab="Fire regime")
plot(fst.amp ~ TSLF, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation",xlab="Time since last fire")
plot(fst.amp ~ jitter(fire), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation",xlab="Fire occurrence")

#Ecophysiology
plot(fst.amp ~ ha_90, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation",xlab="Hours of activity")
plot(fst.amp ~ perf, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Compensation",xlab="Locomotor performance")

#Attenuation
#Life-history traits
plot(fst.att ~ gen.time, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Resistance",xlab="Generation time",bty="n")

plot(fst.att ~ gen.time, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="LB",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="C"]))+14,
     xlab="Generation time",bty="n")#Change the plot to check if the pattern is the same

plot(fst.att ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.att ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.att ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n")#Change the species to check if the pattern is the same

plot(fst.att ~ repro.value, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],ylab="Resistance",xlab="Reproductive value",bty="n")

plot(fst.att ~ repro.value, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="LB",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="C"]))+14,
     xlab="Reproductive value",bty="n")#Change the plot to check if the pattern is the same

plot(fst.att ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same

plot(fst.att ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same

plot(fst.att ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n")#Change the species to check if the pattern is the same

# plot(fst.att ~ shape.surv, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
#      col=pal.plot[res.env.cerrado.liz$plot.int])
# plot(fst.att ~ shape.rep, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
#      col=pal.plot[res.env.cerrado.liz$plot.int])


boxplot(fst.att ~ species, data=res.env.cerrado.liz,ylab="Resistance",xlab="Species")
boxplot(fst.att ~ plot.int, data=res.env.cerrado.liz,ylab="Resistance",xlab="Fire regime")
boxplot(fst.att ~ plot.int*species, data=res.env.cerrado.liz,ylab="Resistance",xlab="Fire regime")

boxplot(fst.att ~ year, data=res.env.cerrado.liz,ylab="Compensation",xlab="Year")
boxplot(fst.att ~ month, data=res.env.cerrado.liz,ylab="Compensation",xlab="Month")

quartz(8,8)
ggplot(res.env.cerrado.liz, 
       aes(x = as.factor(plot.int), y = fst.att)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Resistance")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",], 
       aes(x = as.factor(plot.int), y = fst.att)) + 
  ggdist::stat_halfeye(
    adjust = 5, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Resistance")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",], 
       aes(x = as.factor(plot.int), y = fst.att)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Resistance")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",], 
       aes(x = as.factor(plot.int), y = fst.att)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Resistance")

#Weather and microclimate
plot(fst.att ~ tmed2m, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ RHmax, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ sol, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ tmed0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ tmin0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ precip, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Fire
plot(fst.att ~ jitter(plot.int), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ TSLF, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ jitter(fire), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Ecophysiology
plot(fst.att ~ ha_90, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(fst.att ~ perf, data=res.env.cerrado.liz,pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Recovery time
#Life-history traits
plot(recovery.time ~ gen.time, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],bty="n",xlab="Generation time",ylab="Recovery time (months)")

plot(recovery.time~ gen.time, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="C",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="C"]))+14,
     xlab="Generation time",bty="n",ylab="Recovery time (months)")#Change the plot to check if the pattern is the same

plot(recovery.time ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],
     pch=16,col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

plot(recovery.time ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],
     pch=16,col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

plot(recovery.time ~ gen.time, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],
     pch=16,col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Generation time",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

plot(recovery.time ~ repro.value, 
     data=res.env.cerrado.liz, 
     pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int],xlab="Reproductive value",ylab="Recovery time (months)")

plot(recovery.time~ repro.value, data=res.env.cerrado.liz[res.env.cerrado.liz$plot=="C",], pch=as.numeric(as.factor(res.env.cerrado.liz$species[res.env.cerrado.liz$plot=="C"]))+14,
     xlab="Reproductive output",bty="n",ylab="Recovery time (months)")#Change the plot to check if the pattern is the same

plot(recovery.time ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

plot(recovery.time ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

plot(recovery.time ~ repro.value, 
     data=res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",],
     pch=16, col=pal.plot[res.env.cerrado.liz$plot.int],
     xlab="Reproductive value",bty="n",ylab="Recovery time (months)")#Change the species to check if the pattern is the same

# plot(recovery.time ~ shape.surv, data=res.env.cerrado.liz,pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
#      col=pal.plot[res.env.cerrado.liz$plot.int])
# plot(recovery.time ~ shape.rep, data=res.env.cerrado.liz,pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
#      col=pal.plot[res.env.cerrado.liz$plot.int])
boxplot(recovery.time ~ species, data=res.env.cerrado.liz)
boxplot(recovery.time ~ plot.int*species, data=res.env.cerrado.liz)
boxplot(recovery.time ~ year, data=res.env.cerrado.liz)
boxplot(recovery.time ~ month, data=res.env.cerrado.liz)

ggplot(res.env.cerrado.liz, 
       aes(x = as.factor(plot.int), y = recovery.time)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Recovery time")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="C_nigropunctatum",], 
       aes(x = as.factor(plot.int), y = recovery.time)) + 
  ggdist::stat_halfeye(
    adjust = 5, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Recovery time (months)")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="M_atticolus",], 
       aes(x = as.factor(plot.int), y = recovery.time)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Recovery time (months)")

quartz(8,8)
ggplot(res.env.cerrado.liz[res.env.cerrado.liz$species=="T_itambere",], 
       aes(x = as.factor(plot.int), y = recovery.time)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = .9, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian( clip = "off")+
  labs(x="Fire severity", y="Recovery time (months)")


#Weather and microclimate
plot(recovery.time ~ tmed2m, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ RHmax, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ sol, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ tmed0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ tmin0cm, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ precip, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Fire
plot(recovery.time ~ jitter(plot.int), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ TSLF, data=res.env.cerrado.liz,pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ jitter(fire), data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])

#Ecophysiology
plot(recovery.time ~ ha_90, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])
plot(recovery.time ~ perf, data=res.env.cerrado.liz, pch=as.numeric(as.factor(res.env.cerrado.liz$species))+14,
     col=pal.plot[res.env.cerrado.liz$plot.int])


###########
#Modelling#
###########

#Scale predictors
res.env.cerrado.liz$gen.time <- scale(log10(res.env.cerrado.liz$gen.time))
res.env.cerrado.liz$repro.value <- scale(log10(res.env.cerrado.liz$repro.value))
res.env.cerrado.liz$sp.plot <- interaction(res.env.cerrado.liz$species, res.env.cerrado.liz$plot.int)

res.env.cerrado.liz$fst.amp.lag1 <- lag(res.env.cerrado.liz$fst.amp, n = 1)
res.env.cerrado.liz$fst.att.lag1 <- lag(res.env.cerrado.liz$fst.att, n = 1)
res.env.cerrado.liz$recov.lag1 <- lag(res.env.cerrado.liz$recovery.time, n = 1)

plot(fst.amp ~ fst.amp.lag1, data = res.env.cerrado.liz,
     bty = "n", pch = 19, col = rgb(0,0,0,.5))

plot(fst.att ~ fst.att.lag1, data = res.env.cerrado.liz,
     bty = "n", pch = 19, col = rgb(0,0,0,.5))

plot(recovery.time ~ recov.lag1, data = res.env.cerrado.liz,
     bty = "n", pch = 19, col = rgb(0,0,0,.5))

######################################################
#First-step amplification (reactivity) = compensation#
######################################################
load("brm_res_CerradoLizards.RData")

#Species and plot



ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.5"]))

ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="M_atticolus.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="M_atticolus.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="M_atticolus.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="M_atticolus.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="M_atticolus.5"]))

ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="T_itambere.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="T_itambere.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="T_itambere.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="T_itambere.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$sp.plot=="T_itambere.5"]))

forecast::auto.arima(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum"]),
                     xreg = res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species=="C_nigropunctatum"],
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::auto.arima(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::auto.arima(ts(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

res.env.cerrado.liz$plot.int <- as.factor(res.env.cerrado.liz$plot.int)

(priors<-get_prior(fst.amp ~ -1 + species:plot.int + (1 | year/month) + arma(time = t, gr = sp.plot),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:18]<-"normal(0,1)"

m1.amp.plot.species <- brm(fst.amp ~ -1 + species:plot.int + (1 | year/month),
                           data = res.env.cerrado.liz , family = "gamma",
                           iter = 5000,
                           prior = priors,
                           #control = list(adapt_delta = 0.99),
                           cores = getOption("mc.cores", 4),
                           algorithm = "laplace",
                           backend = "cmdstanr")

summary(m1.amp.plot.species)
plot(m1.amp.plot.species)

bayes_R2(m1.amp.plot.species)

pp_check(m1.amp.plot.species)

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5])


forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1])


forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5])


forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5])


## MVGAM

data.mvgam <- res.env.cerrado.liz[,c("fst.amp", "fst.att", "recovery.time",
                                     "fst.amp.lag1", "fst.att.lag1", "recov.lag1",
                                        "gen.time", "repro.value", 
                                        "species", "plot.int", "t", "sp.plot", "year", "month")]

names(data.mvgam) <- c("fst.amp", "fst.att", "recovery.time",
                       "fst.amp.lag1", "fst.att.lag1", "recov.lag1",
                       "gen.time", "repro.value", 
                       "species", "plot.int", 
                       "time", "series", "year", "month")

data.mvgam <- na.omit(data.mvgam)

data.mvgam <- data.mvgam[data.mvgam$time!=1,]
summary(data.mvgam)

plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 1)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 2)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 3)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 4)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 5)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 6)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 7)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 8)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 9)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 10)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 11)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 12)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 13)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 14)
plot_mvgam_series(data = data.mvgam, y = "fst.amp", series = 15)

data.mvgam$species <- as.factor(data.mvgam$species)
data.mvgam$gen.time <- as.numeric(data.mvgam$gen.time)
data.mvgam$repro.value <- as.numeric(data.mvgam$repro.value)


(priors <- get_mvgam_priors(fst.amp ~ species:plot.int, 
                          data = na.omit(data.mvgam),
                          trend_model = "AR1MA",
                          family = Gamma()))

priors$prior[c(1:16)] <- priors$example_change[c(1:16)]
priors$prior[17] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[18] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"


m1.amp.plot.species.mvgam <- mvgam(fst.amp ~ -1 + species:plot.int, 
                                   data = data.mvgam,
                                   trend_model = "AR1MA",
                                   family = Gamma(),
                                   drift = T,
                                   priors = priors,
                                   noncentred = T,
                                   # adapt_delta = 0.9,
                                   burnin = 3000,
                                   samples = 2000)
summary(m1.amp.plot.species.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
    ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.amp.plot.species.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.amp.plot.species.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.amp.plot.species.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.amp.plot.species.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.amp.plot.species.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.amp.plot.species.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.amp.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.amp.plot.species.mvgam, type = "dens_overlay", variable = "betas")

mcmc_plot(m1.amp.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.amp.plot.species.mvgam, type = "dens_overlay", variable = "betas")

avg_comparisons(m1.amp.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species")


plot(conditional_effects(m1.amp.plot.species.mvgam))


#Forecasts
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 1)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 2)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 3)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 4)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 5)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 6)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 7)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 8)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 9)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 10)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 11)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 12)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 13)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 14)
plot(m1.amp.plot.species.mvgam, type = "forecast", series = 15)


#Residuals
#Residuals
plot(m1.amp.plot.species.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.amp.plot.species.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.plot.species.mvgam$resids$T_itambere.5) ~ 1)


# Posterior checks
ppc(m1.amp.plot.species.mvgam, type = 'hist')
ppc(m1.amp.plot.species.mvgam, type = 'density')
ppc(m1.amp.plot.species.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.amp.plot.species.mvgam)
pp_check(m1.amp.plot.species.mvgam, type = "ecdf_overlay")

pp_check(m1.amp.plot.species.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.amp.plot.species.mvgam, type = 'freqpoly')
pp_check(m1.amp.plot.species.mvgam, type = 'dens_overlay')
pp_check(m1.amp.plot.species.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.amp.plot.species.mvgam, type = 'dens_overlay')
pp_check(m1.amp.plot.species.mvgam, type = 'dens_overlay_grouped', group = "species")


#Smooths
p1 <- plot_predictions(m1.amp.plot.species.mvgam, condition = c("species", "plot.int"), draw = F)

res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
p1$plot.int <- as.integer(p1$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, fst.amp))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
                       ) + 
  ## add justified jitter from the {gghalves} package

  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                                    ymin = conf.low, 
                                                     ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
                  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3, scales = "free_y")+
  labs(x="Fire severity", y="Compensation")

options(marginaleffects_posterior_interval = "hdi")
options(marginaleffects_posterior_center = mean)
avg_comparisons(m1.amp.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species", type = "response")

plot_comparisons(m1.amp.plot.species.mvgam,
                 variables = list(plot.int = "pairwise"),
                 by = "species", type = "response")



#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, fst.amp))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
  ) + 
  ## add justified jitter from the {gghalves} package
  
  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                                    ymin = conf.low, 
                                                    ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3)+
  labs(x="Fire severity", y="Compensation")

#Generation time
(priors<-get_prior(fst.amp ~ gen.time*species + (1|plot/species) + (1|year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:16]<-"normal(0,1)"
m1.amp.gen.time <- brm(fst.amp ~ gen.time*species + (1|plot/species) + (1|year/month),
                       data = res.env.cerrado.liz , family = "gamma",
                       iter = 5000,
                       prior = priors,
                       control = list(adapt_delta = 0.999),
                       cores = getOption("mc.cores", 4))

summary(m1.amp.gen.time)
plot(m1.amp.gen.time)
plot(conditional_effects(m1.amp.gen.time), points=T)

#MVGAM
(priors <- get_mvgam_priors(fst.amp ~ gen.time * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m1.amp.gen.time.mvgam <- mvgam(fst.amp ~ -1 + gen.time*species + t2(plot.int, bs = "re", by = species), 
                                   data = data.mvgam,
                                   trend_model = "AR1MA",
                                   family = Gamma(),
                                   priors = priors,
                                   drift = T,
                                   noncentred = T,
                                   # adapt_delta = 0.9,
                                   burnin = 3000,
                                   samples = 2000)
summary(m1.amp.gen.time.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.amp.gen.time.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.amp.gen.time.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.amp.gen.time.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.amp.gen.time.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.amp.gen.time.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.amp.gen.time.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.amp.gen.time.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.amp.gen.time.mvgam, type = "dens_overlay", variable = "betas")

coef(m1.amp.gen.time.mvgam)
plot(conditional_effects(m1.amp.gen.time.mvgam$mgcv_model))
ggeffects::ggpredict(m1.amp.gen.time.mvgam$mgcv_model, terms = "gen.time[all]")

plot_model(m1.amp.gen.time.mvgam$mgcv_model, type = "pred", terms = "gen.time")

#Forecasts
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 1)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 2)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 3)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 4)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 5)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 6)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 7)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 8)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 9)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 10)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 11)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 12)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 13)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 14)
plot(m1.amp.gen.time.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.amp.gen.time.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.amp.gen.time.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.amp.gen.time.mvgam$resids$T_itambere.5) ~ 1)


# Posterior checks
ppc(m1.amp.gen.time.mvgam, type = 'hist')
ppc(m1.amp.gen.time.mvgam, type = 'density')
ppc(m1.amp.gen.time.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.amp.gen.time.mvgam)
pp_check(m1.amp.gen.time.mvgam, type = "ecdf_overlay")

pp_check(m1.amp.gen.time.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.amp.gen.time.mvgam, type = 'freqpoly')
pp_check(m1.amp.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.amp.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.amp.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.amp.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.amp.gen.time.mvgam, type = "pterms")


#Smooths
p1 <- plot_predictions(m1.amp.gen.time.mvgam, condition = c("gen.time", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,fst.amp))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=fst.amp,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Compensation") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,fst.amp))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=fst.amp,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Compensation") 


#Reproductive value
(priors<-get_prior(fst.amp ~ repro.value*species + (1|plot/species) + (1|year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:16]<-"normal(0,1)"

m1.amp.repro.value <- brm(fst.amp ~ repro.value*species + (1|plot/species) + (1|year/month),
                        data = res.env.cerrado.liz , family = "gamma",
                        iter = 5000,
                        prior = priors,
                        control = list(adapt_delta = 0.99),
                        cores = getOption("mc.cores", 4))

summary(m1.amp.repro.value)
plot(m1.amp.repro.value)
plot(conditional_effects(m1.amp.repro.value),  points = T)

#MVGAM
(priors <- get_mvgam_priors(fst.amp ~ repro.value * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m1.amp.repro.value.mvgam <- mvgam(fst.amp ~ -1 + repro.value*species + t2(plot.int, bs = "re", by = species), 
                               data = data.mvgam,
                               trend_model = "AR1MA",
                               family = Gamma(),
                               priors = priors,
                               drift = T,
                               noncentred = T,
                               adapt_delta = 0.9,
                               burnin = 4000,
                               samples = 4000)
summary(m1.amp.repro.value.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.amp.repro.value.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.amp.repro.value.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.amp.repro.value.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.amp.repro.value.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.amp.repro.value.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.amp.repro.value.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.amp.repro.value.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.amp.repro.value.mvgam, type = "dens_overlay", variable = "betas")

#Forecasts
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 1)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 2)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 3)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 4)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 5)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 6)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 7)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 8)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 9)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 10)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 11)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 12)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 13)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 14)
plot(m1.amp.repro.value.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.amp.repro.value.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.1)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.1)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.1)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.1)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.1)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.1)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.1)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.2)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.2)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.2)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.2)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.2)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.2)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.2)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.3)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.3)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.3)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.3)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.3)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.3)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.3)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.4)[-168] )
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.4)[-168] )

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.4)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.4)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.4)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.4)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.4)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.5)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$C_nigropunctatum.5)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.5)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$M_atticolus.5)[-168])

plot(m1.amp.repro.value.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.5)[-168])
lmtest::bgtest(colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.5)[-1] ~ colMeans(m1.amp.repro.value.mvgam$resids$T_itambere.5)[-168])

# Posterior checks
ppc(m1.amp.repro.value.mvgam, type = 'hist')
ppc(m1.amp.repro.value.mvgam, type = 'density')
ppc(m1.amp.repro.value.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.amp.repro.value.mvgam)
pp_check(m1.amp.repro.value.mvgam, type = "ecdf_overlay")

pp_check(m1.amp.repro.value.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.amp.repro.value.mvgam, type = 'freqpoly')
pp_check(m1.amp.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.amp.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.amp.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.amp.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.amp.repro.value.mvgam, type = "pterms")


#Smooths
p1 <- plot_predictions(m1.amp.repro.value.mvgam, condition = c("repro.value", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)

quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,fst.amp))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=fst.amp,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Compensation") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,fst.amp))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=fst.amp,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Compensation") 

#Compare with null model
(priors <- get_prior(fst.amp ~ species + (1|plot/species) + (1|year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:13]<-"normal(0,1)"

m0.amp <- brm(fst.amp ~ species + (1|plot/species) + (1|year/month),
              data = res.env.cerrado.liz , family = "gamma",
              iter = 5000, prior = priors,
              control = list(adapt_delta = 0.99),
              cores = getOption("mc.cores", 4))

summary(m0.amp)

loo(m0.amp, m1.amp.gen.time, m1.amp.repro.value)

#MVGAM
(priors <- get_mvgam_priors(fst.amp ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:3)] <- priors$example_change[c(1:3)]

priors$prior[5] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[6] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m0.amp.mvgam <- mvgam(fst.amp ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                                  data = data.mvgam,
                                  trend_model = "AR1MA",
                                  family = Gamma(),
                                  priors = priors,
                                  drift = T,
                                  noncentred = T,
                                  # adapt_delta = 0.9,
                                  burnin = 3000,
                                  samples = 2000)
summary(m0.amp.mvgam)
loo(m0.amp, m1.amp.plot.species, m1.amp.gen.time, m1.amp.repro.value)

#MVGAM
ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m0.amp.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m0.amp.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m0.amp.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m0.amp.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m0.amp.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m0.amp.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m0.amp.mvgam, type = "trace", variable = "betas")
mcmc_plot(m0.amp.mvgam, type = "dens_overlay", variable = "betas")

loo_compare(m0.amp.mvgam, m1.amp.plot.species.mvgam, m1.amp.gen.time.mvgam, m1.amp.repro.value.mvgam)

#####################################
#First-step attenuation = resistance#
#####################################
res.env.cerrado.liz$plot.int <- as.factor(res.env.cerrado.liz$plot.int)
(priors<-get_prior(bf(fst.att ~ species:plot.int + (1|year/month),
                      phi ~ species),
                   data = res.env.cerrado.liz , family = "beta"))

priors$prior[1:26]<-"normal(0,1)"
priors$prior[23:26]<-"normal(0,5)"


m1.att.plot.species <- brm(bf(fst.att ~ species:plot.int + (1|year/month),
                              phi ~ species),
                           data = res.env.cerrado.liz , family = "beta",
                           iter = 5000, prior = priors,
                           #control = list(adapt_delta = 0.95),
                           cores = getOption("mc.cores", 4))

summary(m1.att.plot.species)
plot(m1.att.plot.species)
plot(conditional_effects(m1.att.plot.species))

#MVGAM
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.5"]))

ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="M_atticolus.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="M_atticolus.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="M_atticolus.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="M_atticolus.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="M_atticolus.5"]))

ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="T_itambere.1"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="T_itambere.2"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="T_itambere.3"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="T_itambere.4"]))
ar.mle(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$sp.plot=="T_itambere.5"]))

forecast::auto.arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::Arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum"]), order = c(1,0,1))

forecast::auto.arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::Arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus"]), order = c(1,0,1))


forecast::auto.arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::Arima(ts(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere"]), order = c(1,0,1))

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5])


forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1])


forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5])


forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4])

forecast::ggAcf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5])
forecast::ggPacf(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5])



plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 1)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 2)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 3)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 4)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 5)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 6)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 7)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 8)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 9)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 10)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 11)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 12)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 13)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 14)
plot_mvgam_series(data = data.mvgam, y = "fst.att", series = 15)


(priors <- get_mvgam_priors(fst.att ~ species:plot.int, 
                            data = na.omit(data.mvgam),
                            trend_model = "AR1MA",
                            family = betar()))

priors$prior[c(1:16)] <- priors$example_change[c(1:16)]
priors$prior[17] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[18] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m1.att.plot.species.mvgam <- mvgam(fst.att ~ -1 + species:plot.int, 
                                   data = data.mvgam,
                                   trend_model = "AR1MA",
                                   family = betar(),
                                   priors = priors,
                                   drift = T,
                                   noncentred = T,
                                   adapt_delta = 0.99,
                                   burnin = 6000,
                                   samples = 2000)
summary(m1.att.plot.species.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.att.plot.species.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.att.plot.species.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.att.plot.species.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.att.plot.species.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.att.plot.species.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.att.plot.species.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.att.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.att.plot.species.mvgam, type = "dens_overlay", variable = "betas")

mcmc_plot(m1.att.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.att.plot.species.mvgam, type = "dens_overlay", variable = "betas")

avg_comparisons(m1.att.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species")

plot(conditional_effects(m1.att.plot.species.mvgam))


#Forecasts
plot(m1.att.plot.species.mvgam, type = "forecast", series = 1)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 2)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 3)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 4)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 5)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 6)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 7)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 8)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 9)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 10)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 11)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 12)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 13)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 14)
plot(m1.att.plot.species.mvgam, type = "forecast", series = 15)


#Residuals
#Residuals
plot(m1.att.plot.species.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.att.plot.species.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.plot.species.mvgam$resids$T_itambere.5) ~ 1)

# Posterior checks
ppc(m1.att.plot.species.mvgam, type = 'hist')
ppc(m1.att.plot.species.mvgam, type = 'density')
ppc(m1.att.plot.species.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.att.plot.species.mvgam)
pp_check(m1.att.plot.species.mvgam, type = "ecdf_overlay")

pp_check(m1.att.plot.species.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.att.plot.species.mvgam, type = 'freqpoly')
pp_check(m1.att.plot.species.mvgam, type = 'dens_overlay')
pp_check(m1.att.plot.species.mvgam, type = 'dens_overlay_grouped', group = "species")


#Smooths
p1 <- plot_predictions(m1.att.plot.species.mvgam, condition = c("species", "plot.int"), draw = F)

res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
p1$plot.int <- as.integer(p1$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, fst.att))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
  ) + 
  ## add justified jitter from the {gghalves} package
  
  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                 ymin = conf.low, 
                                 ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3, scales = "free_y")+
  labs(x="Fire severity", y="Resistance")

options(marginaleffects_posterior_interval = "hdi")
options(marginaleffects_posterior_center = mean)
avg_comparisons(m1.att.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species", type = "response")

plot_comparisons(m1.att.plot.species.mvgam,
                 variables = list(plot.int = "pairwise"),
                 by = "species", type = "response")



#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, fst.att))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
  ) + 
  ## add justified jitter from the {gghalves} package
  
  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                 ymin = conf.low, 
                                 ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3)+
  labs(x="Fire severity", y="Resistance")



#Generation time
(priors<-get_prior(bf(fst.att ~ gen.time*species + (1|plot/species) + (1|year/month),
                      phi ~ species),
                   data = res.env.cerrado.liz , family = "beta"))

priors$prior[1:20]<-"normal(0,1)"
priors$prior[17:20]<-"normal(0,5)"

m1.att.gen.time <- brm(bf(fst.att ~ gen.time*species + (1|plot/species) + (1|year/month),
                          phi ~ species),
                       data = res.env.cerrado.liz , family = "beta",
                       prior = priors,
                       iter = 5000,
                       control = list(adapt_delta = 0.99),
                       cores = getOption("mc.cores", 4))


summary(m1.att.gen.time)
plot(m1.att.gen.time)

plot(conditional_effects(m1.att.gen.time),  points = T)

#MVGAM
(priors <- get_mvgam_priors(fst.att ~ gen.time * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = betar()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m1.att.gen.time.mvgam <- mvgam(fst.att ~ -1 + gen.time*species + t2(plot.int, bs = "re", by = species), 
                               data = data.mvgam,
                               trend_model = "AR1MA",
                               family = betar(),
                               priors = priors,
                               drift = T,
                               noncentred = T,
                               adapt_delta = 0.99,
                               burnin = 6000,
                               samples = 2000)
summary(m1.att.gen.time.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.att.gen.time.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.att.gen.time.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.att.gen.time.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.att.gen.time.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.att.gen.time.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.att.gen.time.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.att.gen.time.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.att.gen.time.mvgam, type = "dens_overlay", variable = "betas")

plot_model(m1.att.gen.time.mvgam$mgcv_model, type = "pred", terms = "gen.time")

#Forecasts
plot(m1.att.gen.time.mvgam, type = "forecast", series = 1)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 2)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 3)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 4)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 5)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 6)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 7)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 8)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 9)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 10)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 11)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 12)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 13)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 14)
plot(m1.att.gen.time.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.att.gen.time.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.att.gen.time.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.gen.time.mvgam$resids$T_itambere.5) ~ 1)


# Posterior checks
ppc(m1.att.gen.time.mvgam, type = 'hist')
ppc(m1.att.gen.time.mvgam, type = 'density')
ppc(m1.att.gen.time.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.att.gen.time.mvgam)
pp_check(m1.att.gen.time.mvgam, type = "ecdf_overlay")

pp_check(m1.att.gen.time.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.att.gen.time.mvgam, type = 'freqpoly')
pp_check(m1.att.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.att.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.att.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.att.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.att.gen.time.mvgam, type = "pterms")


#Smooths
p1 <- plot_predictions(m1.att.gen.time.mvgam, condition = c("gen.time", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,fst.att))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=fst.att,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Resistance") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,fst.att))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=fst.att,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Resistance") 


#Reproductive value
(priors<-get_prior(bf(fst.att ~ repro.value*species + (1|plot/species) + (1|year/month),
                      phi ~ species),
                   data = res.env.cerrado.liz , family = "beta"))

priors$prior[1:20]<-"normal(0,1)"
priors$prior[17:20]<-"normal(0,5)"

m1.att.repro.value <- brm(bf(fst.att ~ repro.value*species + (1|plot/species) + (1|year/month),
                             phi ~ species),
                        data = res.env.cerrado.liz , family = "beta",
                        prior = priors,
                        iter = 10000,
                        control = list(adapt_delta = 0.99),
                        cores = getOption("mc.cores", 4),
                        silent = 0)

summary(m1.att.repro.value)
plot(m1.att.repro.value)
plot(conditional_effects(m1.att.repro.value),  points = T)

#MVGAM
(priors <- get_mvgam_priors(fst.att ~ repro.value * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = betar()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
priors$prior[10] <- "phi ~ gamma(10, 0.005);"

# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"

c(prior(normal(1000,200)),  # brms default
  prior(gamma(10, 0.005))) %>%  # Solomon's alternative
  parse_dist() %>% 
  
  ggplot(aes(xdist = .dist_obj, y = prior)) + 
  stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
  scale_x_continuous(expression(italic(p)(phi))) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  labs(title = "Beta ANOVA") +
  coord_cartesian(xlim = c(0, 5000))

m1.att.repro.value.mvgam <- mvgam(fst.att ~ -1 + repro.value*species + t2(plot.int, bs = "re", by = species), 
                                  data = data.mvgam,
                                  trend_model = "AR1MA",
                                  family = betar(),
                                  priors = priors,
                                  drift = T,
                                  noncentred = T,
                                  #adapt_delta = 0.999,
                                  burnin = 6000,
                                  samples = 4000)
summary(m1.att.repro.value.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.att.repro.value.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.att.repro.value.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.att.repro.value.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.att.repro.value.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.att.repro.value.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.att.repro.value.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.att.repro.value.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.att.repro.value.mvgam, type = "dens_overlay", variable = "betas")

#Forecasts
plot(m1.att.repro.value.mvgam, type = "forecast", series = 1)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 2)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 3)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 4)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 5)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 6)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 7)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 8)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 9)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 10)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 11)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 12)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 13)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 14)
plot(m1.att.repro.value.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.att.repro.value.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.att.repro.value.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.att.repro.value.mvgam$resids$T_itambere.5) ~ 1)

# Posterior checks
ppc(m1.att.repro.value.mvgam, type = 'hist')
ppc(m1.att.repro.value.mvgam, type = 'density')
ppc(m1.att.repro.value.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.att.repro.value.mvgam)
pp_check(m1.att.repro.value.mvgam, type = "ecdf_overlay")

pp_check(m1.att.repro.value.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.att.repro.value.mvgam, type = 'freqpoly')
pp_check(m1.att.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.att.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.att.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.att.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.att.repro.value.mvgam, type = "pterms")


#Smooths
p1 <- plot_predictions(m1.att.repro.value.mvgam, condition = c("repro.value", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)

quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,fst.att))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=fst.att,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Resistance") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,fst.att))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=fst.att,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Resistance") 

#Compare with null model
(priors <- get_prior(fst.att ~ species + (1|plot/species) + (1|year/month),
                     data = res.env.cerrado.liz , family = "beta"))

priors$prior[1:13]<-"normal(0,1)"

m0.att <- brm(fst.att ~ species + (1|plot/species) + (1|year/month),
              data = res.env.cerrado.liz , family = "beta",
              iter = 5000, prior = priors,
              control = list(adapt_delta = 0.99),
              cores = getOption("mc.cores", 4))

summary(m0.att)

loo(m0.att, m1.att.gen.time, m1.att.repro.value)

#MVGAM
(priors <- get_mvgam_priors(fst.att ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = betar()))

priors$prior[c(1:3)] <- priors$example_change[c(1:3)]

priors$prior[5] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[6] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m0.att.mvgam <- mvgam(fst.att ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                      data = data.mvgam,
                      trend_model = "AR1MA",
                      family = betar(),
                      priors = priors,
                      drift = T,
                      noncentred = T,
                      adapt_delta = 0.99,
                      burnin = 6000,
                      samples = 2000)
summary(m0.att.mvgam)

#MVGAM
ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m0.att.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m0.att.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m0.att.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m0.att.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m0.att.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m0.att.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m0.att.mvgam, type = "trace", variable = "betas")
mcmc_plot(m0.att.mvgam, type = "dens_overlay", variable = "betas")

loo_compare(m0.att.mvgam, m1.att.plot.species.mvgam, m1.att.gen.time.mvgam, m1.att.repro.value.mvgam)

###############
#Recovery time#
###############
#Species and plot

ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.1"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.2"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.3"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.4"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="C_nigropunctatum.5"]))

ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="M_atticolus.1"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="M_atticolus.2"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="M_atticolus.3"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="M_atticolus.4"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="M_atticolus.5"]))

ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="T_itambere.1"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="T_itambere.2"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="T_itambere.3"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="T_itambere.4"]))
ar.mle(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$sp.plot=="T_itambere.5"]))

forecast::auto.arima(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species=="C_nigropunctatum"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::auto.arima(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species=="M_atticolus"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

forecast::auto.arima(ts(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species=="T_itambere"]),
                     stepwise = F, approximation = F, parallel = T, num.cores = 4, max.order = 12)

res.env.cerrado.liz$plot.int <- as.factor(res.env.cerrado.liz$plot.int)

(priors<-get_prior(recovery.time ~ -1 + species:plot.int + (1 | year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:18]<-"normal(0,1)"

m1.recov.t.plot.species <- brm(recovery.time ~ -1 + species:plot.int + (1 | year/month),
                           data = res.env.cerrado.liz , family = "gamma",
                           iter = 5000,
                           prior = priors,
                           #control = list(adapt_delta = 0.99),
                           cores = getOption("mc.cores", 4),
                           backend = "cmdstanr")

summary(m1.recov.t.plot.species)
plot(m1.recov.t.plot.species)

bayes_R2(m1.recov.t.plot.species)

pp_check(m1.recov.t.plot.species)

res.env.cerrado.liz %>%
  add_residual_draws(m1.recov.t.plot.species) %>%
  ggplot(aes(x = t, y = .residual)) +
  stat_lineribbon()+
  facet_wrap(~species+plot.int, scales = "free", ncol = 5)

resids <- residuals(m1.recov.t.plot.species)

forecast::ggAcf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==1,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==2,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==3,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==4,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="C_nigropunctatum" & res.env.cerrado.liz$plot.int==5,1])


forecast::ggAcf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==1,1])


forecast::ggAcf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==2,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==3,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==4,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="M_atticolus" & res.env.cerrado.liz$plot.int==5,1])


forecast::ggAcf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==1,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==2,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==3,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==4,1])

forecast::ggAcf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5,1])
forecast::ggPacf(resids[res.env.cerrado.liz$species=="T_itambere" & res.env.cerrado.liz$plot.int==5,1])


## MVGAM

plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 1)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 2)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 3)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 4)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 5)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 6)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 7)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 8)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 9)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 10)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 11)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 12)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 13)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 14)
plot_mvgam_series(data = data.mvgam, y = "recovery.time", series = 15)


(priors <- get_mvgam_priors(recovery.time ~ species:plot.int, 
                            data = na.omit(data.mvgam),
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:16)] <- priors$example_change[c(1:16)]
priors$prior[17] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[18] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"


m1.recov.t.plot.species.mvgam <- mvgam(recovery.time ~ -1 + species:plot.int, 
                                   data = data.mvgam,
                                   trend_model = "AR1MA",
                                   family = Gamma(),
                                   drift = T,
                                   priors = priors,
                                   noncentred = T,
                                   # adapt_delta = 0.9,
                                   burnin = 3000,
                                   samples = 2000)
summary(m1.recov.t.plot.species.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.recov.t.plot.species.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.recov.t.plot.species.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.recov.t.plot.species.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.recov.t.plot.species.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.recov.t.plot.species.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.recov.t.plot.species.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.recov.t.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.recov.t.plot.species.mvgam, type = "dens_overlay", variable = "betas")

mcmc_plot(m1.recov.t.plot.species.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.recov.t.plot.species.mvgam, type = "dens_overlay", variable = "betas")

avg_comparisons(m1.recov.t.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species")


plot(conditional_effects(m1.recov.t.plot.species.mvgam))


#Forecasts
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 1)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 2)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 3)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 4)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 5)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 6)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 7)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 8)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 9)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 10)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 11)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 12)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 13)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 14)
plot(m1.recov.t.plot.species.mvgam, type = "forecast", series = 15)


#Residuals
#Residuals
plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.recov.t.plot.species.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.plot.species.mvgam$resids$T_itambere.5) ~ 1)


# Posterior checks
ppc(m1.recov.t.plot.species.mvgam, type = 'hist')
ppc(m1.recov.t.plot.species.mvgam, type = 'density')
ppc(m1.recov.t.plot.species.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.recov.t.plot.species.mvgam)
pp_check(m1.recov.t.plot.species.mvgam, type = "ecdf_overlay")

pp_check(m1.recov.t.plot.species.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.recov.t.plot.species.mvgam, type = 'freqpoly')
pp_check(m1.recov.t.plot.species.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.plot.species.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.recov.t.plot.species.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.plot.species.mvgam, type = 'dens_overlay_grouped', group = "species")


#Smooths
p1 <- plot_predictions(m1.recov.t.plot.species.mvgam, condition = c("species", "plot.int"), draw = F)

res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
p1$plot.int <- as.integer(p1$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, recovery.time))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
  ) + 
  ## add justified jitter from the {gghalves} package
  
  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                 ymin = conf.low, 
                                 ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3, scales = "free_y")+
  labs(x="Fire severity", y="Recovery time (months)")

options(marginaleffects_posterior_interval = "hdi")
options(marginaleffects_posterior_center = mean)
avg_comparisons(m1.recov.t.plot.species.mvgam,
                variables = list(plot.int = "pairwise"),
                by = "species", type = "response")

plot_comparisons(m1.recov.t.plot.species.mvgam,
                 variables = list(plot.int = "pairwise"),
                 by = "species", type = "response")

pairs(emmeans(m1.recov.t.plot.species, ~ plot.int|species))


#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(plot.int, recovery.time))+
  geom_point(aes(x = jitter(plot.int), col = plot.int),alpha = .1, size = 3) +
  ggdist::stat_halfeye(aes(fill = plot.int),
                       alpha = 0.5,
                       adjust = 1, 
                       width = .5, 
                       .width = 0, 
                       justification = -.5, 
                       point_colour = NA
  ) + 
  ## add justified jitter from the {gghalves} package
  
  # Use geom_pointrange
  geom_pointrange(data = p1, aes(x = plot.int, y = estimate,
                                 ymin = conf.low, 
                                 ymax= conf.high,   fill = plot.int),
                  shape = 21, colour = "black"
  ) +
  scale_colour_gradientn(colours=turbo(5), name="Fire severity")+
  scale_fill_gradientn(colours=turbo(5), name="Fire severity")+
  facet_wrap(~ species, nrow = 3)+
  labs(x="Fire severity", y="Recovery time (months)")


#Generation time
(priors<-get_prior(recovery.time ~ gen.time*species + (1|plot/species) + (1|year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:16]<-"normal(0,1)"
m1.recov.t.gen.time <- brm(recovery.time ~ gen.time*species + (1|plot/species) + (1|year/month),
                       data = res.env.cerrado.liz , family = "gamma",
                       iter = 5000,
                       prior = priors,
                       control = list(adapt_delta = 0.999),
                       cores = getOption("mc.cores", 4))

summary(m1.recov.t.gen.time)
plot(m1.recov.t.gen.time)
plot(conditional_effects(m1.recov.t.gen.time), points=T)

#MVGAM
(priors <- get_mvgam_priors(recovery.time ~ gen.time * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m1.recov.t.gen.time.mvgam <- mvgam(recovery.time ~ -1 + gen.time*species + t2(plot.int, bs = "re", by = species), 
                               data = data.mvgam,
                               trend_model = "AR1MA",
                               family = Gamma(),
                               priors = priors,
                               drift = T,
                               noncentred = T,
                               # adapt_delta = 0.9,
                               burnin = 3000,
                               samples = 2000)
summary(m1.recov.t.gen.time.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.recov.t.gen.time.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.recov.t.gen.time.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.recov.t.gen.time.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.recov.t.gen.time.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.recov.t.gen.time.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.recov.t.gen.time.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.recov.t.gen.time.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.recov.t.gen.time.mvgam, type = "dens_overlay", variable = "betas")

coef(m1.recov.t.gen.time.mvgam)
plot(conditional_effects(m1.recov.t.gen.time.mvgam))

#Forecasts
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 1)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 2)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 3)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 4)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 5)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 6)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 7)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 8)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 9)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 10)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 11)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 12)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 13)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 14)
plot(m1.recov.t.gen.time.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.1) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.1) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.1) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.1) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.2) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.2) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.2) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.2) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.3) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.3) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.3) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.3) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.4) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.4) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.4) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.4) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$C_nigropunctatum.5) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$M_atticolus.5) ~ 1)

plot(m1.recov.t.gen.time.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.5) ~ 1)
lmtest::bgtest(colMeans(m1.recov.t.gen.time.mvgam$resids$T_itambere.5) ~ 1)


# Posterior checks
ppc(m1.recov.t.gen.time.mvgam, type = 'hist')
ppc(m1.recov.t.gen.time.mvgam, type = 'density')
ppc(m1.recov.t.gen.time.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.recov.t.gen.time.mvgam)
pp_check(m1.recov.t.gen.time.mvgam, type = "ecdf_overlay")

pp_check(m1.recov.t.gen.time.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.recov.t.gen.time.mvgam, type = 'freqpoly')
pp_check(m1.recov.t.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.recov.t.gen.time.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.gen.time.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.recov.t.gen.time.mvgam, type = "pterms")


#Smooths
p1 <- plot_predictions(m1.recov.t.gen.time.mvgam, condition = c("gen.time", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,recovery.time))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=recovery.time,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Compensation") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(gen.time,recovery.time))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=gen.time,y=recovery.time,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=gen.time))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low,ymax=conf.high,x=gen.time,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  # coord_cartesian(ylim=c(0,1))+
  labs(x="Generation time", y="Recovery time (months)") 


#Reproductive value
(priors<-get_prior(recovery.time ~ repro.value*species + (1|plot/species) + (1|year/month),
                   data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:16]<-"normal(0,1)"

m1.recov.t.repro.value <- brm(recovery.time ~ repro.value*species + (1|plot/species) + (1|year/month),
                          data = res.env.cerrado.liz , family = "gamma",
                          iter = 5000,
                          prior = priors,
                          control = list(adapt_delta = 0.99),
                          cores = getOption("mc.cores", 4))

summary(m1.recov.t.repro.value)
plot(m1.recov.t.repro.value)
plot(conditional_effects(m1.recov.t.repro.value),  points = T)

#MVGAM
(priors <- get_mvgam_priors(recovery.time ~ repro.value * species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:6)] <- priors$example_change[c(1:6)]

priors$prior[8] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[9] <- "sigma ~ normal(0, 1);"
priors$prior[10] <- "shape ~ gamma(10, 0.005);"


m1.recov.t.repro.value.mvgam <- mvgam(recovery.time ~ -1 + repro.value*species + t2(plot.int, bs = "re", by = species), 
                                  data = data.mvgam,
                                  trend_model = "AR1MA",
                                  family = Gamma(),
                                  priors = priors,
                                  drift = T,
                                  noncentred = T,
                                  adapt_delta = 0.9,
                                  burnin = 10000,
                                  samples = 4000)

summary(m1.recov.t.repro.value.mvgam)

ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m1.recov.t.repro.value.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m1.recov.t.repro.value.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m1.recov.t.repro.value.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m1.recov.t.repro.value.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m1.recov.t.repro.value.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m1.recov.t.repro.value.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m1.recov.t.repro.value.mvgam, type = "trace", variable = "betas")
mcmc_plot(m1.recov.t.repro.value.mvgam, type = "dens_overlay", variable = "betas")

#Forecasts
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 1)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 2)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 3)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 4)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 5)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 6)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 7)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 8)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 9)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 10)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 11)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 12)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 13)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 14)
plot(m1.recov.t.repro.value.mvgam, type = "forecast", series = 15)


#Residuals
plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 1)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.1)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.1)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 2)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.1)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.1)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 3)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.1)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.1)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.1)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 4)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.2)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.2)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 5)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.2)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.2)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 6)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.2)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.2)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.2)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 7)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.3)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.3)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 8)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.3)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.3)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 9)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.3)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.3)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.3)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 10)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.4)[-168] )
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.4)[-168] )

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 11)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.4)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.4)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 12)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.4)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.4)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.4)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 13)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.5)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$C_nigropunctatum.5)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 14)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.5)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$M_atticolus.5)[-168])

plot(m1.recov.t.repro.value.mvgam, type = "residuals", series = 15)
lmtest::dwtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.5)[-168])
lmtest::bgtest(colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.5)[-1] ~ colMeans(m1.recov.t.repro.value.mvgam$resids$T_itambere.5)[-168])

# Posterior checks
ppc(m1.recov.t.repro.value.mvgam, type = 'hist')
ppc(m1.recov.t.repro.value.mvgam, type = 'density')
ppc(m1.recov.t.repro.value.mvgam, type = 'cdf')

# Many more options are available with pp_check()
pp_check(m1.recov.t.repro.value.mvgam)
pp_check(m1.recov.t.repro.value.mvgam, type = "ecdf_overlay")

pp_check(m1.recov.t.repro.value.mvgam, type = "ecdf_overlay_grouped", group = "species")
pp_check(m1.recov.t.repro.value.mvgam, type = 'freqpoly')
pp_check(m1.recov.t.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")
pp_check(m1.recov.t.repro.value.mvgam, type = 'dens_overlay')
pp_check(m1.recov.t.repro.value.mvgam, type = 'dens_overlay_grouped', group = "species")

plot(m1.recov.t.repro.value.mvgam, type = "pterms")

plot(conditional_effects(m1.recov.t.repro.value.mvgam), points = 0.5)

#Smooths
p1 <- plot_predictions(m1.recov.t.repro.value.mvgam, condition = c("repro.value", "species"), type = "response", draw = F)


res.env.cerrado.liz$plot.int <- as.integer(res.env.cerrado.liz$plot.int)

quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,recovery.time))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=recovery.time,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Recovery time (months)") 

#Fixed labs
quartz(width=6,height=9)
ggplot(res.env.cerrado.liz,aes(repro.value,recovery.time))+
  geom_point(data=res.env.cerrado.liz,
             mapping=aes(x=repro.value,y=recovery.time,color=plot.int),
             size=5, alpha=.4)+
  scale_colour_gradientn(colours=turbo(5, alpha = .4), name="Fire severity")+
  geom_line(data=p1,
            aes(y=estimate,x=repro.value))+
  scale_linetype_manual(values = c("dotted", "dashed", "twodash"))+
  geom_ribbon(data=p1,
              aes(ymin=conf.low, ymax=conf.high,x=repro.value,fill=species),
              alpha=.2)+
  scale_fill_manual(values = c("black", "blue", "brown"))+
  facet_wrap(~ species, nrow = 3) +
  #coord_cartesian(ylim=c(0,.5))+
  labs(x="Reproductive output", y="Recovery time (months)") 

#Compare with null model
(priors <- get_prior(recovery.time ~ species + (1|plot/species) + (1|year/month),
                     data = res.env.cerrado.liz , family = "gamma"))

priors$prior[1:13]<-"normal(0,1)"
priors$prior[14] <- "gamma(10, 0.005)"

m0.recov.t <- brm(recovery.time ~ species + (1|plot/species) + (1|year/month),
              data = res.env.cerrado.liz , family = "gamma",
              iter = 5000, prior = priors,
              control = list(adapt_delta = 0.99),
              cores = getOption("mc.cores", 4))

summary(m0.recov.t)

loo(m0.recov.t, m1.recov.t.gen.time, m1.recov.t.repro.value)

#MVGAM
(priors <- get_mvgam_priors(recovery.time ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                            data = data.mvgam,
                            trend_model = "AR1MA",
                            family = Gamma()))

priors$prior[c(1:3)] <- priors$example_change[c(1:3)]

priors$prior[5] <- "ar1 ~ normal(0.8,0.1);"
priors$prior[6] <- "sigma ~ normal(0, 1);"
# priors$prior[21] <- "sigma_obs ~ normal(0, 0.1);"



m0.recov.t.mvgam <- mvgam(recovery.time ~ -1 + species + t2(plot.int, bs = "re", by = species), 
                      data = data.mvgam,
                      trend_model = "AR1MA",
                      family = Gamma(),
                      priors = priors,
                      drift = T,
                      noncentred = T,
                      # adapt_delta = 0.9,
                      burnin = 3000,
                      samples = 2000)
summary(m0.recov.t.mvgam)
loo(m0.recov.t, m1.recov.t.plot.species, m1.recov.t.gen.time, m1.recov.t.repro.value)

#MVGAM
ar1_pars <- rep(NA, 15)
for(i in 1:15){
  ar1_pars[i] <- paste0('ar1[', i, ']')
}


mcmc_plot(m0.recov.t.mvgam, type = "trace", variable = ar1_pars)
mcmc_plot(m0.recov.t.mvgam, type = "dens_overlay", variable = ar1_pars)

theta_pars <- rep(NA, 15)
for(i in 1:15){
  theta_pars[i] <- paste0('theta[', i, ']')
}

mcmc_plot(m0.recov.t.mvgam, type = "trace", variable = theta_pars)
mcmc_plot(m0.recov.t.mvgam, type = "dens_overlay", variable = theta_pars)

sigma_pars <- rep(NA, 15)
for(i in 1:15){
  sigma_pars[i] <- paste0('sigma[', i, ']')
}

mcmc_plot(m0.recov.t.mvgam, type = "trace", variable = sigma_pars)
mcmc_plot(m0.recov.t.mvgam, type = "dens_overlay", variable = sigma_pars)


mcmc_plot(m0.recov.t.mvgam, type = "trace", variable = "betas")
mcmc_plot(m0.recov.t.mvgam, type = "dens_overlay", variable = "betas")

loo_compare(m0.recov.t.mvgam, m1.recov.t.plot.species.mvgam, m1.recov.t.gen.time.mvgam, m1.recov.t.repro.value.mvgam)



mean(res.env.cerrado.liz$recovery.time)
round(tapply(res.env.cerrado.liz$recovery.time, list(res.env.cerrado.liz$species), mean),3)
round(tapply(res.env.cerrado.liz$recovery.time, list(res.env.cerrado.liz$species), sd), 3)



save(list = c("m0.amp", "m1.amp.plot.species", "m1.amp.gen.time", "m1.amp.repro.value",
       "m0.att", "m1.att.plot.species", "m1.att.gen.time", "m1.att.repro.value",
       "m0.recov.t", "m1.recov.t.plot.species", "m1.recov.t.gen.time", "m1.recov.t.repro.value",
       "m0.amp.mvgam", "m1.amp.plot.species.mvgam", "m1.amp.gen.time.mvgam", "m1.amp.repro.value.mvgam",
       "m0.att.mvgam", "m1.att.plot.species.mvgam", "m1.att.gen.time.mvgam", "m1.att.repro.value.mvgam",
       "m0.recov.t.mvgam", "m1.recov.t.plot.species.mvgam", "m1.recov.t.gen.time.mvgam", "m1.recov.t.repro.value.mvgam"),
     file = "brm_res_CerradoLizards.RData")

##########
#3D plots#
##########
r3dDefaults$useFreeType <- FALSE
open3d()
plot3d(res.env.cerrado.liz$fst.amp,
       res.env.cerrado.liz$fst.att,
       res.env.cerrado.liz$recovery.time,
       xlab = "Compensation",
       ylab = "Resistance",
       zlab = "Recovery time",
       col = c("black", "blue", "brown")[as.numeric(as.factor(res.env.cerrado.liz$species))],
       alpha = 0.1,
       radius = 0.05,
       type = "s",
       aspect = c(1, 1, 1))

ellips.Cn <- ellipse3d(cov(cbind(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "C_nigropunctatum"],
                                 res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "C_nigropunctatum"],
                                 res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "C_nigropunctatum"])),
                       centre=c(mean(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "C_nigropunctatum"]),
                                mean(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "C_nigropunctatum"]),
                                mean(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "C_nigropunctatum"])),
                       level = 0.95)

ellips.Ma <- ellipse3d(cov(cbind(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "M_atticolus"],
                                 res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "M_atticolus"],
                                 res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "M_atticolus"])),
                       centre=c(mean(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "M_atticolus"]),
                                mean(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "M_atticolus"]),
                                mean(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "M_atticolus"])),
                       level = 0.95)

ellips.Ti <- ellipse3d(cov(cbind(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "T_itambere"],
                                 res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "T_itambere"],
                                 res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "T_itambere"])),
                       centre=c(mean(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "T_itambere"]),
                                mean(res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "T_itambere"]),
                                mean(res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "T_itambere"])),
                       level = 0.95)

plot3d(ellips.Cn, col = "black", alpha = 0.2, add = TRUE, box = FALSE, type = "wire")
plot3d(ellips.Ma, col = "blue", alpha = 0.2, add = TRUE, box = FALSE, type = "wire")
plot3d(ellips.Ti, col = "brown", alpha = 0.2, add = TRUE, box = FALSE, type = "wire")
grid3d(c("x+", "y+", "z"))

rgl.snapshot("resilience_sp_3d.png",fmt="png")

close3d()

open3d()
plot3d(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "C_nigropunctatum"],
       res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "C_nigropunctatum"],
       res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "C_nigropunctatum"],
       xlab = "Compensation",
       ylab = "Resistance",
       zlab = "Recovery time",
       col = pal.plot[res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "C_nigropunctatum"]],
       alpha = 0.1,
       radius = 0.005,
       type = "s",
       aspect = c(1, 1, 1))

grid3d(c("x", "y", "z"))

rgl.snapshot("resilience_Cnigro_3d.png",fmt="png")

close3d()

open3d()
plot3d(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "M_atticolus"],
       res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "M_atticolus"],
       res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "M_atticolus"],
       xlab = "Compensation",
       ylab = "Resistance",
       zlab = "Recovery time",
       col = pal.plot[res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "M_atticolus"]],
       alpha = 0.1,
       radius = 0.025,
       type = "s",
       aspect = c(1, 1, 1))

grid3d(c("x", "y", "z"))

rgl.snapshot("resilience_Matticolus_3d.png",fmt="png")

close3d()

open3d()
plot3d(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "T_itambere"],
       res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "T_itambere"],
       res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "T_itambere"],
       xlab = "Compensation",
       ylab = "Resistance",
       zlab = "Recovery time",
       col = pal.plot[res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "T_itambere"]],
       alpha = 0.1,
       radius = 0.025,
       type = "s",
       aspect = c(1, 1, 1))

grid3d(c("x", "y", "z"))

rgl.snapshot("resilience_Titambere_3d.png",fmt="png")

close3d()

#######################
#Using plot 3d package#
#######################

quartz(width = 10, height = 8)
scatter3D(res.env.cerrado.liz$fst.amp,
          res.env.cerrado.liz$fst.att,
          res.env.cerrado.liz$recovery.time, 
          bty = "g", pch = 20, cex = 2,
          xlab = "Compensation",
          ylab = "Resistance",
          zlab = "Recovery time",
          colvar = as.integer(as.factor(res.env.cerrado.liz$species)), 
          col = c(alpha.col("black", 0.2), 
                  alpha.col("blue", 0.2), 
                  alpha.col("brown", 0.2)),
          ticktype = "detailed",
          phi = 20, theta = 30,
          colkey = list(at = c(1, 2, 3), side = 4,
                        addlines = T, length = 0.5, width = 0.5,
                        labels = c("C_nigropunctatum", "M_atticolus", "T_itambere"))
         )

# Add small dots on basal plane and on the depth plane
scatter3D_fancy <- function(x, y, z, colvar, col, colkey, ...)
{
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = 20, 
              cex = 1, add = TRUE, col = col, colkey = colkey)
    
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = 20, 
              cex = 1, add = TRUE, col = col, colkey = colkey)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst, colkey = colkey, col = col) 
}

quartz(width = 10, height = 8)
scatter3D_fancy(res.env.cerrado.liz$fst.amp,
                res.env.cerrado.liz$fst.att,
                res.env.cerrado.liz$recovery.time, 
                bty = "g", pch = 20, cex = 2,
                xlab = "Compensation",
                ylab = "Resistance",
                zlab = "Recovery time",
                colvar = as.integer(as.factor(res.env.cerrado.liz$species)), 
                col = c(alpha.col("black", 0.2), 
                        alpha.col("blue", 0.2), 
                        alpha.col("brown", 0.2)),
                ticktype = "detailed",
                phi = 20, theta = 30,
                colkey = list(at = c(1, 2, 3), side = 4,
                              addlines = T, length = 0.5, width = 0.5,
                              labels = c("C_nigropunctatum", "M_atticolus", "T_itambere")))

quartz(width = 10, height = 8)
scatter3D(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "C_nigropunctatum"],
          res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "C_nigropunctatum"],
          res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "C_nigropunctatum"], 
          bty = "g", pch = 20,
          xlab = "Compensation",
          ylab = "Resistance",
          zlab = "Recovery time",
          colvar = res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "C_nigropunctatum"], 
          col = pal.plot,
          ticktype = "detailed",
          NAcol = "",
          alpha = 0.2,
          theta = 300, phi = 20
)

quartz(width = 10, height = 8)
scatter3D(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "M_atticolus"],
          res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "M_atticolus"],
          res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "M_atticolus"], 
          bty = "g", pch = 20,
          xlab = "Compensation",
          ylab = "Resistance",
          zlab = "Recovery time",
          colvar = res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "M_atticolus"], 
          col = pal.plot,
          ticktype = "detailed",
          NAcol = "",
          alpha = 0.2,
          theta = 300, phi = 20
)

quartz(width = 10, height = 8)
scatter3D(res.env.cerrado.liz$fst.amp[res.env.cerrado.liz$species == "T_itambere"],
          res.env.cerrado.liz$fst.att[res.env.cerrado.liz$species == "T_itambere"],
          res.env.cerrado.liz$recovery.time[res.env.cerrado.liz$species == "T_itambere"], 
          bty = "g", pch = 20,
          xlab = "Compensation",
          ylab = "Resistance",
          zlab = "Recovery time",
          colvar = res.env.cerrado.liz$plot.int[res.env.cerrado.liz$species == "T_itambere"], 
          col = pal.plot,
          ticktype = "detailed",
          NAcol = "",
          alpha = 0.2,
          theta = 300, phi = 20
)


# Heat map ----------------------------------------------------------------

ggplot(res.env.cerrado.liz, 
       aes(x=fst.att, y=fst.amp, color=recovery.time))+
  geom_point(aes(x=fst.att, y=fst.amp, color = recovery.time, shape = species), size = 3, alpha =0.5)+
  scale_color_viridis_c(name = "Recovery time")+
  # scale_fill_gradientn(colours=viridis(10), name = "Recovery time")+
  # geom_point(data = res.env.cerrado.liz, 
  #            aes(x = fst.att, y = fst.amp, group = species, z = 0),
  #            alpha = 0.5, size = 3) +
  # scale_color_viridis_c(name = "species")+
  # scale_alpha(range = c(1, 0.1)) +
  labs(x = "Resistance", y = "Compensation")

ggplot(res.env.cerrado.liz)+
  geom_point(aes(x=fst.amp, y=recovery.time, color = fst.att, shape = species), size = 3, alpha =0.5)+
  scale_color_viridis_c(name = "Resistance")+
  # scale_fill_gradientn(colours=viridis(10), name = "Recovery time")+
  # geom_point(data = res.env.cerrado.liz, 
  #            aes(x = fst.att, y = fst.amp, group = species, z = 0),
  #            alpha = 0.5, size = 3) +
  # scale_color_viridis_c(name = "species")+
  # scale_alpha(range = c(1, 0.1)) +
  # facet_wrap(~species, nrow = 3, scales = "free")+
  labs(x = "Compensation", y = "Recovery time")

quartz(width = 10, height = 8)

ggplot(res.env.cerrado.liz)+
  geom_point(aes(x=fst.att, y=recovery.time, color = fst.amp, shape = species), size = 3, alpha =0.5)+
  scale_color_viridis_c(name = "Compensation")+
  # scale_fill_gradientn(colours=viridis(10), name = "Recovery time")+
  # geom_point(data = res.env.cerrado.liz, 
  #            aes(x = fst.att, y = fst.amp, group = species, z = 0),
  #            alpha = 0.5, size = 3) +
  # scale_color_viridis_c(name = "species")+
  # scale_alpha(range = c(1, 0.1)) +
  # facet_wrap(~species, nrow = 3, scales = "free")+
  labs(x = "Resistance", y = "Recovery time")

# Plot of captures through time -------------------------------------------

plot(colSums(Cnigropunctatum.data$y), type = "l")  


#Copeoglossum nigropunctatum
cap.plot <- matrix(nrow = 5, ncol = 170)

mplot <- data.frame(C = as.numeric(Cnigropunctatum.data$plot==1),
                    Q = as.numeric(Cnigropunctatum.data$plot==2),
                    BP = as.numeric(Cnigropunctatum.data$plot==3),
                    BM = as.numeric(Cnigropunctatum.data$plot==4),
                    BT = as.numeric(Cnigropunctatum.data$plot==5))

for(i in 1:5){
  cap.plot[i,] <- colSums(Cnigropunctatum.data$y[mplot[,i]==1,])
  
}

cap.plot[2,] <- colSums(Cnigropunctatum.data$y[mplot[,2]==1,])



cap.plot.df <- as.data.frame(t(cap.plot))

names(cap.plot.df) <- c("C", "Q", "EB", "MB", "LB")

cap.plot.df$total <- rowSums(cap.plot.df)
cap.plot.df$t <- 1:170
cap.plot.df$time <- seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/12/01"), by = "month")
cap.plot.df$year <- lubridate::year(cap.plot.df$time)
cap.plot.df$sp <- rep("C_nigropunctatum", 170)

#Micrablepharus atticolus
cap.Ma.plot <- matrix(nrow = 5, ncol = 170)

mplot.Ma <- data.frame(C = as.numeric(Matticolus.data$plot==1),
                    Q = as.numeric(Matticolus.data$plot==2),
                    BP = as.numeric(Matticolus.data$plot==3),
                    BM = as.numeric(Matticolus.data$plot==4),
                    BT = as.numeric(Matticolus.data$plot==5))

for(i in 1:5){
  cap.Ma.plot[i,] <- colSums(Matticolus.data$y[mplot.Ma[,i]==1,])
  
}

cap.Ma.plot.df <- as.data.frame(t(cap.Ma.plot))

names(cap.Ma.plot.df) <- c("C", "Q", "EB", "MB", "LB")

cap.Ma.plot.df$total <- rowSums(cap.Ma.plot.df)
cap.Ma.plot.df$t <- 1:170
cap.Ma.plot.df$time <- seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/12/01"), by = "month")
cap.Ma.plot.df$year <- lubridate::year(cap.Ma.plot.df$time)
cap.Ma.plot.df$sp <- rep("M_atticolus", 170)


#Tropidurus itambere
cap.Ti.plot <- matrix(nrow = 5, ncol = 170)

mplot.Ti <- data.frame(C = as.numeric(Titambere.data$plot==1),
                       Q = as.numeric(Titambere.data$plot==2),
                       BP = as.numeric(Titambere.data$plot==3),
                       BM = as.numeric(Titambere.data$plot==4),
                       BT = as.numeric(Titambere.data$plot==5))

for(i in 1:5){
  cap.Ti.plot[i,] <- colSums(Titambere.data$y[mplot.Ti[,i]==1,])
  
}

(cap.Ti.plot.df <- as.data.frame(t(cap.Ti.plot)))

names(cap.Ti.plot.df) <- c("C", "Q", "EB", "MB", "LB")

cap.Ti.plot.df$total <- rowSums(cap.Ti.plot.df)
cap.Ti.plot.df$t <- 1:170
cap.Ti.plot.df$time <- seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/12/01"), by = "month")
cap.Ti.plot.df$year <- lubridate::year(cap.Ti.plot.df$time)
cap.Ti.plot.df$sp <- rep("T_itambere", 170)

#All spp
cap.spp.plot.df <- rbind(cap.plot.df, cap.Ma.plot.df, cap.Ti.plot.df)
cap.spp.plot.year.df <- cap.spp.plot.df %>% group_by(sp, year) %>% 
  summarise(C = sum(C),
            Q = sum(Q),
            EB = sum(EB),
            MB = sum(MB),
            LB = sum(LB),
            total = sum(total))


ggplot(cap.spp.plot.year.df, aes(x = year, y = total)) +
  geom_col(aes(x = year, y = total), position = "dodge") +
  facet_wrap(~ sp, nrow = 3)

cap.spp.plot.year.long.df <- pivot_longer(data = cap.spp.plot.year.df,
                                          cols = `C`:`LB`,
                                          names_to = "plot",
                                          values_to = "value")

cap.spp.plot.year.long.df$plot <- factor(cap.spp.plot.year.long.df$plot, 
                                         levels = c("C", "Q", "EB", "MB", "LB"))
quartz(height = 8, width = 10)
ggplot(cap.spp.plot.year.long.df, aes(x = year, y = value, fill = plot)) +
  geom_col(aes(x = year, y = value)) +
  facet_wrap(~ sp, nrow = 3) +
  scale_fill_manual(values = turbo(5)) +
  labs(x="Year", y="Number of captures")

  
capts.spp <- data.frame(Cnigro = colSums(Cnigropunctatum.data$y),
                        Matticolus = colSums(Matticolus.data$y),
                        Titambere = colSums(Titambere.data$y),
                        t = 1:170)

capts.spp$time <- seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/12/01"), by = "month")

capts.spp$year <- lubridate::year(capts.spp$time)

ggplot(capts.spp, aes(x = time, y = Cnigro)) +
  geom_line(aes(x = time, y = Cnigro), colour = "black", alpha = 0.5, linewidth = 1.2) +
  geom_line(aes(x = time, y = Matticolus), colour = "blue", alpha = 0.5, linewidth = 1.2) +
  geom_line(aes(x = time, y = Titambere), colour = "brown", alpha = 0.5, linewidth = 1.2)
  
capts.spp.year <- capts.spp %>% group_by(year) %>% 
  summarise(Cnigro = sum(Cnigro),
            Matticolus = sum(Matticolus),
            Titambere = sum(Titambere))

quartz(height = 8, width = 10)
ggplot(capts.spp.year, aes(x = year, y = Cnigro)) +
  geom_line(aes(x = year, y = Cnigro), colour = "black", alpha = 0.5, linewidth = 1.2) +
  geom_line(aes(x = year, y = Matticolus), colour = "blue", alpha = 0.5, linewidth = 1.2) +
  geom_line(aes(x = year, y = Titambere), colour = "brown", alpha = 0.5, linewidth = 1.2) +
  labs(x="Year", y="Number of captures")

ggplot(capts.spp.year, aes(x = year, y = Cnigro)) +
  geom_col(aes(x = year, y = Cnigro), colour = "black", fill = "black", alpha = 0.5, linewidth = 1.2) +
  geom_col(aes(x = year, y = Matticolus), colour = "blue",fill = "blue", alpha = 0.5, linewidth = 1.2) +
  geom_col(aes(x = year, y = Titambere), colour = "brown", fill = "brown", alpha = 0.5, linewidth = 1.2)

quartz(height = 8, width = 10)
barplot(t(as.matrix(capts.spp.year[,-1])), col = scales::alpha(c("black", "blue", "brown"),0.5), 
        legend.text = T, args.legend = list(x="topleft"), names.arg = capts.spp.year$year,
        xlab = "Year", ylab = "Number of captures")

quartz(height = 8, width = 10)
barplot(t(as.matrix(capts.spp.year[,-1])), col = scales::alpha(c("black", "blue", "brown"),0.5), 
        legend.text = T, args.legend = list(x="topleft"), names.arg = capts.spp.year$year,
        xlab = "Year", ylab = "Number of captures", beside = T)

quartz(height = 8, width = 10)
barplot(prop.table(t(as.matrix(capts.spp.year[,-1])),margin = 2), col = scales::alpha(c("black", "blue", "brown"),0.5), 
        names.arg = capts.spp.year$year,
        xlab = "Year", ylab = "Proportion of captures", beside = F)

plot(1:15,capts.spp.year$Cnigro, col = scales::alpha(c("black", 0.5)))

  
capts.plot.spp <- data.frame(Freq = c(table(Cnigropunctatum.data$plot),
                                        table(Matticolus.data$plot),
                                        table(Titambere.data$plot)),
                             species = as.factor(c(rep("C_nigropunctatum", 5),
                                         rep("M_atticolus", 5),
                                         rep("Titambere", 5))),
                             plot = factor(rep(c("C", "Q", "EB", "MB", "LB"), 3),
                                           levels = c("C", "Q", "EB", "MB", "LB"))
                             )
capts.plot.spp$plot.int <- as.integer(capts.plot.spp$plot)

quartz(height = 9, width = 6)
ggplot(capts.plot.spp, aes(x = plot, y = Freq, colour = plot.int, fill = plot.int)) +
  geom_col() +
  scale_colour_gradientn(colours = turbo(5), name = "Fire severity")+
  scale_fill_gradientn(colours = turbo(5), name = "Fire severity")+
  facet_wrap(~ species, nrow = 3) +
  # coord_cartesian(ylim=c(0,5))+
  labs(x="Plot", y="Individuals")


# Survival and recruitment estimates x time -------------------------------
#C. nigropunctatum
pradel.Cnigro.df <- readRDS("results_pradel_Cnigro_df.rds")

f.pradel.Cnigro <- as.data.frame(pradel.Cnigro.df[grep(pattern = "f", 
                                                x = row.names(pradel.Cnigro.df))[1:845],])
tail(f.pradel.Cnigro)

phi.pradel.Cnigro <- as.data.frame(pradel.Cnigro.df[grep(pattern = "phi", 
                                                  x = row.names(pradel.Cnigro.df))[1:845],])
tail(phi.pradel.Cnigro)

f.pradel.Cnigro$plot <- as.factor(rep(1:5,169))
f.pradel.Cnigro$time <- rep(1:169, each = 5)

phi.pradel.Cnigro$plot <- rep(1:5,169)
phi.pradel.Cnigro$time <- rep(1:169, each = 5)

f.pradel.Cnigro$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
f.pradel.Cnigro$year <- lubridate::year(f.pradel.Cnigro$date)
f.pradel.Cnigro$month <- lubridate::month(f.pradel.Cnigro$date)
f.pradel.Cnigro$sp <- rep("C_nigropunctatum", 845)

phi.pradel.Cnigro$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
phi.pradel.Cnigro$year <- lubridate::year(phi.pradel.Cnigro$date)

phi.pradel.Cnigro$sp <- rep("C_nigropunctatum", 845)

quartz(width=8,height=9)
ggplot(f.pradel.Cnigro, aes(x = as.factor(month), y = mean, colour = plot))+
  geom_boxplot(aes(x = as.factor(month), y = mean, colour = plot), alpha=0.5) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25)
  #geom_path(data=f.pradel[f.ipm2$plot==3,], aes(x = time, y= mean, colour = plot), linetype = "dashed" )+
  # ylim(c(0,0.3))+
  scale_color_manual(values=turbo(5))+
  labs(y = "Recruitment per capita", x = "Month")


svl.Cnigro <- readRDS("Cnigropunctatum.RECOR.imp.rds")

ggplot(svl.Cnigro, aes(x = as.Date(data), y = crc))+
  geom_point(alpha=0.25)+
  geom_line(aes(y = 50), linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Time")

quartz(width=8,height=9)
ggplot(svl.Cnigro, aes(x = as.factor(month), y = crc))+
  geom_boxplot(aes(x = as.factor(month), y = crc), alpha=0.5)+
  geom_abline(intercept = 50, slope = 0, linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Month")

ggplot(f.pradel.Cnigro, aes(x = date, y = mean, colour = plot))+
  geom_line(aes(x = date, y = mean, colour = plot), alpha=0.95, linewidth = 2) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25, color = NA)+
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = plot), alpha = 0.25, color = NA)+
  facet_wrap(~ plot, nrow = 5)+
  scale_color_manual(values=turbo(5))+
  scale_fill_manual(values = turbo(5))

#M. atticolus
pradel.Matticolus.df <- read.csv("results.ipm2.Matticolus.df_100000iters.csv")

f.pradel.Matticolus <- as.data.frame(pradel.Matticolus.df[grep(pattern = "f", 
                                                               x = pradel.Matticolus.df$X)[1:845],])

tail(f.pradel.Matticolus)

phi.pradel.Matticolus <- as.data.frame(pradel.Matticolus.df[grep(pattern = "phi", 
                                                         x = pradel.Matticolus.df$X)[1:845],])
tail(phi.pradel.Matticolus)

f.pradel.Matticolus$plot <- as.factor(rep(1:5,169))
f.pradel.Matticolus$time <- rep(1:169, each = 5)

phi.pradel.Matticolus$plot <- rep(1:5,169)
phi.pradel.Matticolus$time <- rep(1:169, each = 5)

f.pradel.Matticolus$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
f.pradel.Matticolus$year <- lubridate::year(f.pradel.Matticolus$date)
f.pradel.Matticolus$month <- lubridate::month(f.pradel.Matticolus$date)
f.pradel.Matticolus$sp <- rep("M_atticolus", 845)

phi.pradel.Matticolus$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
phi.pradel.Matticolus$year <- lubridate::year(phi.pradel.Matticolus$date)

phi.pradel.Matticolus$sp <- rep("M_atticolus", 845)

quartz(width=8,height=9)
ggplot(f.pradel.Matticolus, aes(x = as.factor(month), y = mean, colour = plot))+
  geom_boxplot(aes(x = as.factor(month), y = mean, colour = plot), alpha=0.5) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25)
  #geom_path(data=f.pradel[f.ipm2$plot==3,], aes(x = time, y= mean, colour = plot), linetype = "dashed" )+
  # ylim(c(0,0.3))+
  scale_color_manual(values=turbo(5))+
  labs(y = "Recruitment per capita", x = "Month")


svl.Matticolus <- readRDS("Matticolus.RECOR.imp.rds")

ggplot(svl.Matticolus, aes(x = as.Date(data), y = crc))+
  geom_point(alpha=0.25)+
  geom_line(aes(y = 30), linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Time")

quartz(width=8,height=9)
ggplot(svl.Matticolus, aes(x = as.factor(month), y = crc))+
  geom_boxplot(aes(x = as.factor(month), y = crc), alpha=0.5)+
  geom_abline(intercept = 30, slope = 0, linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Month")

ggplot(f.pradel.Matticolus, aes(x = date, y = mean, colour = plot))+
  geom_line(aes(x = date, y = mean, colour = plot), alpha=0.95, linewidth = 2) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25, color = NA)+
  geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, fill = plot), alpha = 0.25, color = NA)+
  facet_wrap(~ plot, nrow = 5)+
  scale_color_manual(values=turbo(5))+
  scale_fill_manual(values = turbo(5))

#T. itambere
pradel.Titambere.df <- read.csv("results.pradel.itambere.df_400000_noecophys.csv")

f.pradel.Titambere <- as.data.frame(pradel.Titambere.df[grep(pattern = "f", 
                                                               x = pradel.Titambere.df$X)[1:845],])

tail(f.pradel.Titambere)

phi.pradel.Titambere <- as.data.frame(pradel.Titambere.df[grep(pattern = "phi", 
                                                                 x = pradel.Titambere.df$X)[1:845],])
tail(phi.pradel.Titambere)

f.pradel.Titambere$plot <- as.factor(rep(1:5,169))
f.pradel.Titambere$time <- rep(1:169, each = 5)

phi.pradel.Titambere$plot <- rep(1:5,169)
phi.pradel.Titambere$time <- rep(1:169, each = 5)

f.pradel.Titambere$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
f.pradel.Titambere$year <- lubridate::year(f.pradel.Titambere$date)
f.pradel.Titambere$month <- lubridate::month(f.pradel.Titambere$date)
f.pradel.Titambere$sp <- rep("T_itambere", 845)

phi.pradel.Titambere$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
phi.pradel.Titambere$year <- lubridate::year(phi.pradel.Titambere$date)

phi.pradel.Titambere$sp <- rep("T_itambere", 845)

quartz(width=8,height=9)
ggplot(f.pradel.Titambere, aes(x = as.factor(month), y = mean, colour = plot))+
  geom_boxplot(aes(x = as.factor(month), y = mean, colour = plot), alpha=0.5) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25)
  #geom_path(data=f.pradel[f.ipm2$plot==3,], aes(x = time, y= mean, colour = plot), linetype = "dashed" )+
  # ylim(c(0,0.3))+
  scale_color_manual(values=turbo(5))+
  labs(y = "Recruitment per capita", x = "Month")

svl.Titambere <- readRDS("Titambere.RECOR.imp.rds")

ggplot(svl.Titambere, aes(x = as.Date(data), y = crc))+
  geom_point(alpha=0.25)+
  geom_line(aes(y = 50), linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Time")

quartz(width=8,height=9)
ggplot(svl.Titambere, aes(x = as.factor(month), y = crc))+
  geom_boxplot(aes(x = as.factor(month), y = crc), alpha=0.5)+
  geom_abline(intercept = 50, slope = 0, linetype = 2)+
  labs(y = "Snout-vent length (mm)", x = "Month")

ggplot(f.pradel.Titambere, aes(x = date, y = mean, colour = plot))+
  geom_line(aes(x = date, y = mean, colour = plot), alpha=0.95, linewidth = 2) +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean +sd, fill = plot), alpha = 0.25, color = NA)+
  geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, fill = plot), alpha = 0.25, color = NA)+
  facet_wrap(~ plot, nrow = 5)+
  scale_color_manual(values=turbo(5))+
  scale_fill_manual(values = turbo(5))


# Population growth x Resilience ------------------------------------------
rho.pradel.Cnigro <- as.data.frame(pradel.Cnigro.df[grep(pattern = "rho", 
                                                       x = row.names(pradel.Cnigro.df))[1:845],])
tail(rho.pradel.Cnigro)

rho.pradel.Cnigro$plot <- factor(rep(c("C", "Q", "EB", "MB", "LB"),169), levels = c("C", "Q", "EB", "MB", "LB"), ordered = T)
rho.pradel.Cnigro$t <- rep(1:169, each = 5)
rho.pradel.Cnigro$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
rho.pradel.Cnigro$year <- lubridate::year(rho.pradel.Cnigro$date)
rho.pradel.Cnigro$month <- lubridate::month(rho.pradel.Cnigro$date)
rho.pradel.Cnigro$species <- rep("C_nigropunctatum", 845)

rho.pradel.Matticolus <- as.data.frame(pradel.Matticolus.df[grep(pattern = "rho", 
                                                               x = pradel.Matticolus.df$X)[1:845],])

tail(rho.pradel.Matticolus)

rho.pradel.Matticolus$plot <- factor(rep(c("C", "Q", "EB", "MB", "LB"),169), levels = c("C", "Q", "EB", "MB", "LB"), ordered = T)
rho.pradel.Matticolus$t <- rep(1:169, each = 5)
rho.pradel.Matticolus$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
rho.pradel.Matticolus$year <- lubridate::year(rho.pradel.Matticolus$date)
rho.pradel.Matticolus$month <- lubridate::month(rho.pradel.Matticolus$date)
rho.pradel.Matticolus$species <- rep("M_atticolus", 845)

rho.pradel.Titambere <- as.data.frame(pradel.Titambere.df[grep(pattern = "rho", 
                                                             x = pradel.Titambere.df$X)[1:845],])

tail(rho.pradel.Titambere)

rho.pradel.Titambere$plot <- factor(rep(c("C", "Q", "EB", "MB", "LB"),169), levels = c("C", "Q", "EB", "MB", "LB"), ordered = T)
rho.pradel.Titambere$t <- rep(1:169, each = 5)
rho.pradel.Titambere$date <- rep(seq.Date(from = as.Date("2005/11/01"), to = as.Date("2019/11/01"), by = "month"), each = 5)
rho.pradel.Titambere$year <- lubridate::year(rho.pradel.Titambere$date)
rho.pradel.Titambere$month <- lubridate::month(rho.pradel.Titambere$date)
rho.pradel.Titambere$species <- rep("T_itambere", 845)

names(rho.pradel.Cnigro) <- names(rho.pradel.Matticolus[,-1])
rho.df <- rbind(rho.pradel.Cnigro, rho.pradel.Matticolus[,-1], rho.pradel.Titambere[,-1], deparse.level = 2)

rho.res.df <- full_join(res.env.cerrado.liz, rho.df, by = c("plot", "t", "date", "year", "month", "species"))

quartz(width=6,height=9)
ggplot(rho.res.df, aes(y = fst.amp, x = log(mean), color = plot))+
  geom_point(alpha = 0.5, size = 3) +
  # geom_smooth(aes(fill = plot), method = "loess")+
  facet_wrap("species", nrow = 3, scale = "free_y") +
  scale_colour_manual(values = turbo(5))+
  # scale_fill_manual(values = turbo(5))+
  labs(x = "Realized population growth (PJS)", y = "Compensation (IPM)")

quartz(width=6,height=9)
ggplot(rho.res.df, aes(y = fst.att, x = log(mean), color = plot))+
  # geom_smooth(aes(fill = plot), method = "loess")+
  geom_point(alpha = 0.5, size = 3) +
  facet_wrap("species", nrow = 3, scale = "free_y") +
  scale_colour_manual(values = turbo(5))+
  # scale_fill_manual(values = turbo(5))+
  labs(x = "Population growth (PJS)", y = "Resistance (IPM)")

quartz(width=6,height=9)
ggplot(rho.res.df, aes(y = recovery.time, x = log(mean), color = plot))+
  geom_point(alpha = 0.5, size = 3) +
  facet_wrap("species",  nrow = 3, scale = "free_y") +
  # geom_smooth(aes(fill = plot), method = "loess")+
  scale_colour_manual(values = turbo(5))+
  # scale_fill_manual(values = turbo(5))+
  labs(x = "Population growth (PJS)", y = "Recovery time (IPM)")

#Pop. growth x Time
quartz(width=6,height=9)
ggplot(rho.res.df, aes(x = date, y = log(mean), colour = plot, fill = plot)) + 
  geom_line(alpha = 0.5, linewidth = 1.5)+
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Population growth (PJS)") +
  coord_cartesian( clip = "off")


#Pop. growth x plot x species

cv <- function(x){sd(x)/mean(x)}

round(tapply(log(rho.res.df$mean), INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = mean),3)
round(tapply(log(rho.res.df$mean), INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = median),3)

round(tapply(log(rho.res.df$mean), INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = sd),3)
round(tapply(log(rho.res.df$mean), INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = cv),3)

round(tapply(log(rho.res.df$mean[rho.res.df$mean<1]), 
       INDEX = list(rho.res.df$plot[rho.res.df$mean<1], 
                    rho.res.df$species[rho.res.df$mean<1]), FUN = mean),3)

round(tapply(log(rho.res.df$mean[rho.res.df$mean<1]), 
       INDEX = list(rho.res.df$plot[rho.res.df$mean<1], 
                    rho.res.df$species[rho.res.df$mean<1]), FUN = sd),3)

round(tapply(log(rho.res.df$mean[rho.res.df$mean>1]), 
       INDEX = list(rho.res.df$plot[rho.res.df$mean>1], 
                    rho.res.df$species[rho.res.df$mean>1]), FUN = mean),3)

round(tapply(log(rho.res.df$mean[rho.res.df$mean>1]), 
       INDEX = list(rho.res.df$plot[rho.res.df$mean>1], 
                    rho.res.df$species[rho.res.df$mean>1]), FUN = sd),3)

exp(tapply(log(rho.res.df$X97.5.), INDEX = list(rho.res.df$species, rho.res.df$plot), FUN = mean))
exp(tapply(log(rho.res.df$X2.5.), INDEX = list(rho.res.df$species, rho.res.df$plot), FUN = mean))
round(tapply(rho.res.df$mean>1, INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = sum)/169,2)
round(tapply(rho.res.df$mean<1, INDEX = list(rho.res.df$plot, rho.res.df$species), FUN = sum)/169,2)


#Boxplots
quartz(height = 9, width = 6)
ggplot(rho.res.df, aes(x = plot, y = log(mean), colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1,
    width = 1,
    .width = 0,
    justification = -.1,
    point_colour = NA,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = .1, 
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  # scale_y_continuous(limits = c(0.5,2))+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Realized population growth (PJS)") +
  coord_cartesian( clip = "off")

quartz(height = 9, width = 6)
ggplot(rho.res.df[rho.res.df$mean<1,], aes(x = plot, y = mean, colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1,
    width = 1,
    .width = 0,
    justification = -.1,
    point_colour = NA,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = .1, 
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  # scale_y_continuous(limits = c(0.5,2))+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Resistance (PJS)") +
  coord_cartesian( clip = "off")

quartz(height = 9, width = 6)
ggplot(rho.res.df, aes(x = plot, y = fst.att, colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA, 
    alpha = 0.6
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA, 
    alpha = 0.6
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Resistance (IPM)") +
  coord_cartesian( clip = "off")

quartz(height = 9, width = 6)
ggplot(rho.res.df[rho.res.df$mean>1,], aes(x = plot, y = log10(mean), colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1,
    width = 1,
    .width = 0,
    justification = -.1,
    point_colour = NA,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = .1, 
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  # scale_y_continuous(limits = c(0.5,2))+
  facet_wrap(~species, nrow = 3, scales = "free")+
  labs(x="Fire severity", y="Compensation (PJS)") +
  coord_cartesian( clip = "off")

ggplot(rho.res.df, aes(x = plot, y = fst.amp, colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA, 
    alpha = 0.6
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA, 
    alpha = 0.6
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3)+
  labs(x="Fire severity", y="Compensation (IPM)") +
  coord_cartesian( clip = "off")

ggplot(rho.res.df, aes(x = plot, y = recovery.time, colour = plot, fill = plot)) + 
  ggdist::stat_halfeye(
    adjust = 1, 
    width = .5, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA, 
    alpha = 0.6
  ) + 
  geom_boxplot(
    width = .05, 
    outlier.shape = NA, 
    alpha = 0.6
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .5, 
    ## add some transparency
    alpha = .3
  ) +
  scale_colour_manual(values=turbo(5), name="Fire severity")+
  scale_fill_manual(values=turbo(5),  name="Fire severity")+
  facet_wrap(~species, nrow = 3, scale = "free")+
  labs(x="Fire severity", y="Recovery time (IPM)") +
  coord_cartesian( clip = "off")

#Plots of environmental effects on vital rates (survival and recruitment) --------
#Organize data

#Copeoglossum nigropunctatum
pradel.Cnigro.df <- readRDS("results_pradel_Cnigro_df.rds")


pradel.Cnigro.env.df <- rbind(
  #Environmental slopes for survival
  pradel.Cnigro.df['betaphiJS[1]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[2]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[3]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[4]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[5]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[6]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[7]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[8]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[9]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaphiJS[10]',c('mean', '2.5%', '97.5%')],
  
  pradel.Cnigro.df['sigma.phiJS',c('mean', '2.5%', '97.5%')],
  
  #Environmental slopes for reproduction
  pradel.Cnigro.df['betaf[1]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[2]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[3]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[4]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[5]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[6]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[7]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[8]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[9]',c('mean', '2.5%', '97.5%')],
  pradel.Cnigro.df['betaf[10]',c('mean', '2.5%', '97.5%')],
  
  pradel.Cnigro.df['sigma.f',c('mean', '2.5%', '97.5%')]
  
)
pradel.Cnigro.env.df$vitalrate <- factor(c(rep("Survival", 11), 
                                    rep("Recruitment", 11)),
                                  levels = c("Survival", "Recruitment"))
pradel.Cnigro.env.df$cov <- factor(c(rep(c("tmed2m",
                                    "RHmax",
                                    "sol",
                                    "tmed0cm",
                                    "tmin0cm",
                                    "precip",
                                    "perf",
                                    "ha_90",
                                    "fire",
                                    "TSLF",
                                    "random"),
                                    2)
                                  ),
                                  levels = c("tmed2m",
                                               "sol",
                                               "precip",
                                               "tmed0cm",
                                               "tmin0cm",
                                               "RHmax",
                                               "perf",
                                               "ha_90",
                                               "fire",
                                               "TSLF",
                                               "random"))

names(pradel.Cnigro.env.df) <- c("mean", "lower", "upper", "vitalrate", "cov")

#Micrablepharus atticolus
pradel.Matticolus.df <- read.csv("results.ipm2.Matticolus.df_100000iters.csv")


pradel.Matticolus.env.df <- rbind(
  #Environmental slopes for survival
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[1]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[2]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[3]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[4]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[5]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[6]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[7]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[8]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[9]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaphiJS[10]',c('mean', 'X2.5.', 'X97.5.')],
  
  pradel.Matticolus.df[pradel.Matticolus.df=='sigma.phiJS',c('mean', 'X2.5.', 'X97.5.')],
  
  #Environmental slopes for reproduction
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[1]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[2]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[3]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[4]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[5]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[6]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[7]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[8]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[9]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Matticolus.df[pradel.Matticolus.df=='betaf[10]',c('mean', 'X2.5.', 'X97.5.')],
  
  pradel.Matticolus.df[pradel.Matticolus.df=='sigma.f',c('mean', 'X2.5.', 'X97.5.')]
  
)

pradel.Matticolus.env.df <- na.omit(pradel.Matticolus.env.df)

pradel.Matticolus.env.df$vitalrate <- factor(c(rep("Survival", 11), 
                                           rep("Recruitment", 11)),
                                         levels = c("Survival", "Recruitment"))

pradel.Matticolus.env.df$cov <- factor(c(rep(c("tmed2m",
                                           "RHmax",
                                           "sol",
                                           "tmed0cm",
                                           "tmin0cm",
                                           "precip",
                                           "perf",
                                           "ha_90",
                                           "fire",
                                           "TSLF",
                                           "random"),
                                           2)
                                         ),
                                       levels = c("tmed2m",
                                                  "sol",
                                                  "precip",
                                                  "tmed0cm",
                                                  "tmin0cm",
                                                  "RHmax",
                                                  "perf",
                                                  "ha_90",
                                                  "fire",
                                                  "TSLF",
                                                  "random")
                                       )

names(pradel.Matticolus.env.df) <- c("mean", "lower", "upper", "vitalrate", "cov")

#Tropidurus itambere
pradel.Titambere.df <- read.csv("results.pradel.itambere.df_400000_noecophys.csv")


pradel.Titambere.env.df <- rbind(
  #Environmental slopes for survival
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[1]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[2]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[3]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[4]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[5]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[6]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[7]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[8]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[9]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaphiJS[10]',c('mean', 'X2.5.', 'X97.5.')],
  
  pradel.Titambere.df[pradel.Titambere.df=='sigma.phiJS',c('mean', 'X2.5.', 'X97.5.')],
  
  #Environmental slopes for reproduction
  pradel.Titambere.df[pradel.Titambere.df=='betaf[1]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[2]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[3]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[4]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[5]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[6]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[7]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[8]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[9]',c('mean', 'X2.5.', 'X97.5.')],
  pradel.Titambere.df[pradel.Titambere.df=='betaf[10]',c('mean', 'X2.5.', 'X97.5.')],
  
  pradel.Titambere.df[pradel.Titambere.df=='sigma.f',c('mean', 'X2.5.', 'X97.5.')]
  
)

pradel.Titambere.env.df <- na.omit(pradel.Titambere.env.df)

pradel.Titambere.env.df$vitalrate <- factor(c(rep("Survival", 11), 
                                               rep("Recruitment", 11)),
                                             levels = c("Survival", "Recruitment"))

pradel.Titambere.env.df$cov <- factor(c(rep(c("tmed2m",
                                               "RHmax",
                                               "sol",
                                               "tmed0cm",
                                               "tmin0cm",
                                               "precip",
                                               "perf",
                                               "ha_90",
                                               "fire",
                                               "TSLF",
                                               "random"),
                                             2)
                                        ),
                                      levels = c("tmed2m",
                                                 "sol",
                                                 "precip",
                                                 "tmed0cm",
                                                 "tmin0cm",
                                                 "RHmax",
                                                 "perf",
                                                 "ha_90",
                                                 "fire",
                                                 "TSLF",
                                                 "random")
                                      )

names(pradel.Titambere.env.df) <- c("mean", "lower", "upper", "vitalrate", "cov")


#Join the three species
pradel.Cnigro.env.df$species <- c("C_nigropunctatum")
pradel.Matticolus.env.df$species <- c("M_atticolus")
pradel.Titambere.env.df$species <- c("T_itambere")

pradel.env.df <- rbind(pradel.Cnigro.env.df,
                       pradel.Matticolus.env.df,
                       pradel.Titambere.env.df)

quartz(width=6,height=9)
ggplot(pradel.env.df, aes(x = cov, y = mean, colour = vitalrate)) +
  geom_pointrange(aes(ymin = lower, ymax = upper, colour = vitalrate), size = 1.25, linewidth = 2,
                  position = position_dodge(width = 1)) +
  facet_wrap(~ species, nrow = 3) + 
  scale_color_manual(values = viridis::cividis(2), name = "Vital rate")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  labs(x = "Environmental covariate", y = "Estimate")

# GOF tests ---------------------------------------------------------------

#C. nigropunctatum
# Cnigropunctatum.data$y

overall_CJS(as.matrix(Cnigropunctatum.data$y), rep(1, nrow(Cnigropunctatum.data$y)))
overall_CJS(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==1,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==1])
overall_CJS(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==2,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==2])
overall_CJS(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==3,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==3])
overall_CJS(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==4,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==4])
overall_CJS(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==5,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==5])

#Test 3.SR - transients
#H0 = Newly encountered individuals have the same chance to be later reobserved  as  recaptured  (previously  encountered)  individuals
test3sr(as.matrix(Cnigropunctatum.data$y), rep(1, nrow(Cnigropunctatum.data$y)))$test3sr
test3sr(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==1,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==1])$test3sr
test3sr(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==2,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==2])$test3sr
test3sr(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==3,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==3])$test3sr
test3sr(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==4,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==4])$test3sr
test3sr(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==5,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==5])$test3sr

#Test 2.CT - trap dependence
#H0 = Missed individuals have the same chance to be recaptured at the next occasion as currently captured individuals
test2ct(as.matrix(Cnigropunctatum.data$y), rep(1, nrow(Cnigropunctatum.data$y)))$test2ct
test2ct(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==1,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==1])$test2ct
test2ct(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==2,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==2])$test2ct
test2ct(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==3,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==3])$test2ct
test2ct(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==4,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==4])$test2ct
test2ct(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==5,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==5])$test2ct

#Test 3.SM
#H0 = Among those individuals seen again, when they were seen does not differ among previously and newly marked individuals
test3sm(as.matrix(Cnigropunctatum.data$y), rep(1, nrow(Cnigropunctatum.data$y)))$test3sm
test3sm(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==1,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==1])$test3sm
test3sm(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==2,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==2])$test3sm
test3sm(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==3,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==3])$test3sm
test3sm(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==4,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==4])$test3sm
test3sm(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==5,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==5])$test3sm

#Test 2.CL
#H0 = There is no difference in the timing of reencounters between the individuals encountered and not encountered at occasion i, con-ditional on presence at both occasions i and i +2
test2cl(as.matrix(Cnigropunctatum.data$y), rep(1, nrow(Cnigropunctatum.data$y)))$test2cl
test2cl(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==1,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==1])$test2cl
test2cl(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==2,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==2])$test2cl
test2cl(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==3,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==3])$test2cl
test2cl(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==4,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==4])$test2cl
test2cl(as.matrix(Cnigropunctatum.data$y)[Cnigropunctatum.data$plot==5,], rep(1, nrow(Cnigropunctatum.data$y))[Cnigropunctatum.data$plot==5])$test2cl

#M. atticolus
# Matticolus.data$y

overall_CJS(as.matrix(Matticolus.data$y), rep(1, nrow(Matticolus.data$y)))
overall_CJS(as.matrix(Matticolus.data$y)[Matticolus.data$plot==1,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==1])
overall_CJS(as.matrix(Matticolus.data$y)[Matticolus.data$plot==2,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==2])
overall_CJS(as.matrix(Matticolus.data$y)[Matticolus.data$plot==3,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==3])
overall_CJS(as.matrix(Matticolus.data$y)[Matticolus.data$plot==4,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==4])
overall_CJS(as.matrix(Matticolus.data$y)[Matticolus.data$plot==5,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==5])

#Test 3.SR - transients
#H0 = Newly encountered individuals have the same chance to be later reobserved  as  recaptured  (previously  encountered)  individuals
test3sr(as.matrix(Matticolus.data$y), rep(1, nrow(Matticolus.data$y)))$test3sr
test3sr(as.matrix(Matticolus.data$y)[Matticolus.data$plot==1,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==1])$test3sr
test3sr(as.matrix(Matticolus.data$y)[Matticolus.data$plot==2,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==2])$test3sr
test3sr(as.matrix(Matticolus.data$y)[Matticolus.data$plot==3,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==3])$test3sr
test3sr(as.matrix(Matticolus.data$y)[Matticolus.data$plot==4,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==4])$test3sr
test3sr(as.matrix(Matticolus.data$y)[Matticolus.data$plot==5,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==5])$test3sr

#Test 2.CT - trap dependence
#H0 = Missed individuals have the same chance to be recaptured at the next occasion as currently captured individuals
test2ct(as.matrix(Matticolus.data$y), rep(1, nrow(Matticolus.data$y)))$test2ct
test2ct(as.matrix(Matticolus.data$y)[Matticolus.data$plot==1,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==1])$test2ct
test2ct(as.matrix(Matticolus.data$y)[Matticolus.data$plot==2,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==2])$test2ct
test2ct(as.matrix(Matticolus.data$y)[Matticolus.data$plot==3,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==3])$test2ct
test2ct(as.matrix(Matticolus.data$y)[Matticolus.data$plot==4,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==4])$test2ct
test2ct(as.matrix(Matticolus.data$y)[Matticolus.data$plot==5,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==5])$test2ct

#Test 3.SM
#H0 = Among those individuals seen again, when they were seen does not differ among previously and newly marked individuals
test3sm(as.matrix(Matticolus.data$y), rep(1, nrow(Matticolus.data$y)))$test3sm
test3sm(as.matrix(Matticolus.data$y)[Matticolus.data$plot==1,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==1])$test3sm
test3sm(as.matrix(Matticolus.data$y)[Matticolus.data$plot==2,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==2])$test3sm
test3sm(as.matrix(Matticolus.data$y)[Matticolus.data$plot==3,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==3])$test3sm
test3sm(as.matrix(Matticolus.data$y)[Matticolus.data$plot==4,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==4])$test3sm
test3sm(as.matrix(Matticolus.data$y)[Matticolus.data$plot==5,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==5])$test3sm

#Test 2.CL
#H0 = There is no difference in the timing of reencounters between the individuals encountered and not encountered at occasion i, conditional on presence at both occasions i and i +2
test2cl(as.matrix(Matticolus.data$y), rep(1, nrow(Matticolus.data$y)))$test2cl
test2cl(as.matrix(Matticolus.data$y)[Matticolus.data$plot==1,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==1])$test2cl
test2cl(as.matrix(Matticolus.data$y)[Matticolus.data$plot==2,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==2])$test2cl
test2cl(as.matrix(Matticolus.data$y)[Matticolus.data$plot==3,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==3])$test2cl
test2cl(as.matrix(Matticolus.data$y)[Matticolus.data$plot==4,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==4])$test2cl
test2cl(as.matrix(Matticolus.data$y)[Matticolus.data$plot==5,], rep(1, nrow(Matticolus.data$y))[Matticolus.data$plot==5])$test2cl

#T. itambere
# Titambere.data$y

overall_CJS(as.matrix(Titambere.data$y), rep(1, nrow(Titambere.data$y)))
overall_CJS(as.matrix(Titambere.data$y)[Titambere.data$plot==1,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==1])
overall_CJS(as.matrix(Titambere.data$y)[Titambere.data$plot==2,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==2])
overall_CJS(as.matrix(Titambere.data$y)[Titambere.data$plot==3,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==3])
overall_CJS(as.matrix(Titambere.data$y)[Titambere.data$plot==4,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==4])
overall_CJS(as.matrix(Titambere.data$y)[Titambere.data$plot==5,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==5])

#Test 3.SR - transients
#H0 = Newly encountered individuals have the same chance to be later reobserved  as  recaptured  (previously  encountered)  individuals
test3sr(as.matrix(Titambere.data$y), rep(1, nrow(Titambere.data$y)))$test3sr
test3sr(as.matrix(Titambere.data$y)[Titambere.data$plot==1,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==1])$test3sr
test3sr(as.matrix(Titambere.data$y)[Titambere.data$plot==2,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==2])$test3sr
test3sr(as.matrix(Titambere.data$y)[Titambere.data$plot==3,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==3])$test3sr
test3sr(as.matrix(Titambere.data$y)[Titambere.data$plot==4,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==4])$test3sr
test3sr(as.matrix(Titambere.data$y)[Titambere.data$plot==5,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==5])$test3sr

#Test 2.CT - trap dependence
#H0 = Missed individuals have the same chance to be recaptured at the next occasion as currently captured individuals
test2ct(as.matrix(Titambere.data$y), rep(1, nrow(Titambere.data$y)))$test2ct
test2ct(as.matrix(Titambere.data$y)[Titambere.data$plot==1,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==1])$test2ct
test2ct(as.matrix(Titambere.data$y)[Titambere.data$plot==2,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==2])$test2ct
test2ct(as.matrix(Titambere.data$y)[Titambere.data$plot==3,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==3])$test2ct
test2ct(as.matrix(Titambere.data$y)[Titambere.data$plot==4,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==4])$test2ct
test2ct(as.matrix(Titambere.data$y)[Titambere.data$plot==5,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==5])$test2ct

#Test 3.SM
#H0 = Among those individuals seen again, when they were seen does not differ among previously and newly marked individuals
test3sm(as.matrix(Titambere.data$y), rep(1, nrow(Titambere.data$y)))$test3sm
test3sm(as.matrix(Titambere.data$y)[Titambere.data$plot==1,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==1])$test3sm
test3sm(as.matrix(Titambere.data$y)[Titambere.data$plot==2,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==2])$test3sm
test3sm(as.matrix(Titambere.data$y)[Titambere.data$plot==3,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==3])$test3sm
test3sm(as.matrix(Titambere.data$y)[Titambere.data$plot==4,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==4])$test3sm
test3sm(as.matrix(Titambere.data$y)[Titambere.data$plot==5,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==5])$test3sm

#Test 2.CL
#H0 = There is no difference in the timing of reencounters between the individuals encountered and not encountered at occasion i, conditional on presence at both occasions i and i +2
test2cl(as.matrix(Titambere.data$y), rep(1, nrow(Titambere.data$y)))$test2cl
test2cl(as.matrix(Titambere.data$y)[Titambere.data$plot==1,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==1])$test2cl
test2cl(as.matrix(Titambere.data$y)[Titambere.data$plot==2,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==2])$test2cl
test2cl(as.matrix(Titambere.data$y)[Titambere.data$plot==3,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==3])$test2cl
test2cl(as.matrix(Titambere.data$y)[Titambere.data$plot==4,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==4])$test2cl
test2cl(as.matrix(Titambere.data$y)[Titambere.data$plot==5,], rep(1, nrow(Titambere.data$y))[Titambere.data$plot==5])$test2cl
