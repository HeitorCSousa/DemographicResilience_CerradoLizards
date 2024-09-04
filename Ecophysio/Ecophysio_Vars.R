Sys.setenv(TZ = "America/Sao_Paulo")#set time zone
#set working directory
setwd("~/Documents/GitHub/DemographicResilience_CerradoLizards/Ecophysio/")


# Load packages -----------------------------------------------------------

# devtools::install_github("ilyamaclean/microclima")
# remotes::install_github("everydayduffy/mcera5")
# devtools::install_github("mrke/NicheMapR")
# install.packages("/Users/heito/Downloads/NicheMapR_3.1.0.tgz", repos = NULL, 
#                  type = .Platform$pkgType)

library(raster)
library(NicheMapR)
library(microclima)
library(ecmwfr)
library(mcera5)
library(lubridate)
library(dplyr)
library(tidync)
library(ncdf4)
library(lifecycle)
library(ggplot2)
library(gamm4)
library(MODISTools)
library(gdalUtils)
library(lubridate)
library(tidyverse)
library(lme4)
library(mgcv)
library(sjPlot)
library(viridis)
library(Mapinguari)
library(usdm)
library(corrplot)

# Generate ecophysiological-related variables -----------------------------




# Performance
physio_BSB<- read.table("Ecophysio_BSB.txt", h=T)
head(physio_BSB)

perf_df <-
  data.frame(
    temp = c(physio_BSB$sprint_1_tb, 
             physio_BSB$sprint_2_tb, 
             physio_BSB$sprint_3_tb, 
             physio_BSB$CTmin, 
             physio_BSB$CTmax),
    perf = c(physio_BSB$sprint_1_speed, 
             physio_BSB$sprint_2_speed, 
             physio_BSB$sprint_3_speed, 
             rep(0, nrow(physio_BSB)*2)),
    sp = rep(physio_BSB$Species,5),
    id = rep(physio_BSB$ID,5),
    crc = rep(physio_BSB$SVL,5))

completos <- complete.cases(perf_df[, c("temp","perf")]) #remove NAs
perf_df <- droplevels(perf_df[completos, ])

perf_Cnigropunctatum_df <- perf_df[perf_df$sp=="C_nigropunctatum",]
perf_Matticolus_df <- perf_df[perf_df$sp=="M_atticolus",]
perf_Titambere_df <- perf_df[perf_df$sp=="T_itambere",]

#How many individuals?
length(unique(perf_Cnigropunctatum_df$id))
length(unique(perf_Matticolus_df$id))
length(unique(perf_Titambere_df$id))

summary(perf_Cnigropunctatum_df)
summary(perf_Matticolus_df)
summary(perf_Titambere_df)

quartz(8,8)
plot(perf_Cnigropunctatum_df$temp, perf_Cnigropunctatum_df$perf, bty="n",
     xlab="Temperature (°C)", ylab = "Speed (cm/s)",xlim = c(0,50),
     main = expression(italic("Copeoglossum nigropunctatum")))

quartz(8,8)
plot(perf_Cnigropunctatum_df$crc, perf_Cnigropunctatum_df$perf, bty="n",
     xlab="Size (mm)", ylab = "Speed (cm/s)",xlim = c(65,105),
     main = expression(italic("Copeoglossum nigropunctatum")))


quartz(8,8)
plot(perf_Matticolus_df$temp, perf_Matticolus_df$perf, bty="n",
     xlab="Temperature (°C)", ylab = "Speed (cm/s)",xlim = c(0,50),
     main = expression(italic("Micrablepharus atticolus")))

quartz(8,8)
plot(perf_Matticolus_df$crc, perf_Matticolus_df$perf, bty="n",
     xlab="Size (mm)", ylab = "Speed (cm/s)",xlim = c(30,50),
     main = expression(italic("Micrablepharus atticolus")))


quartz(8,8)
plot(perf_Titambere_df$temp, perf_Titambere_df$perf, bty="n",
     xlab="Temperature (°C)", ylab = "Speed (cm/s)",xlim = c(0,50),
     main = expression(italic("Tropidurus itambere")))

quartz(8,8)
plot(perf_Titambere_df$crc, perf_Titambere_df$perf, bty="n",
     xlab="Size (mm)", ylab = "Speed (cm/s)",xlim = c(40,100),
     main = expression(italic("Tropidurus itambere")))


tpc_Cnigropunctatum <- 
  gamm4::gamm4(perf ~ s(temp, k = 5) + crc, random = ~(1|id), 
               data = perf_Cnigropunctatum_df)

summary(tpc_Cnigropunctatum$gam)

tpc_Cnigropunctatum <- 
  gamm4::gamm4(perf ~ s(temp, k = 5), random = ~(1|id), 
               data = perf_Cnigropunctatum_df)

summary(tpc_Cnigropunctatum$gam)

tpc_Matticolus <- 
  gamm4::gamm4(perf ~ s(temp, k = 5) + crc, random = ~(1|id), 
               data = perf_Matticolus_df)

summary(tpc_Matticolus$gam)

tpc_Matticolus <- 
  gamm4::gamm4(perf ~ s(temp, k = 5), random = ~(1|id), 
               data = perf_Matticolus_df)

summary(tpc_Matticolus$gam)

tpc_Titambere <- 
  gamm4::gamm4(perf ~ s(temp, k = 5) + crc, random = ~(1|id), 
               data = perf_Titambere_df)

summary(tpc_Titambere$gam)

tpc_Titambere <- 
  gamm4::gamm4(perf ~ s(temp, k = 5), random = ~(1|id), 
               data = perf_Titambere_df)

summary(tpc_Titambere$gam)

preddata_tpc_Cnigropunctatum<-data.frame(temp=seq(10,50,0.1))

preddata_tpc_Matticolus<-data.frame(temp=seq(10,50,0.1))

preddata_tpc_Titambere<-data.frame(temp=seq(10,50,0.1))

# Faz a predicao 
pred_tpc_Cnigropunctatum <- predict(tpc_Cnigropunctatum$gam,
                                    preddata_tpc_Cnigropunctatum,
                                    se.fit=T)

# Anexa as predicoes e erros na tabela de dados
preddata_tpc_Cnigropunctatum$predicted <- pred_tpc_Cnigropunctatum$fit
preddata_tpc_Cnigropunctatum$se <- pred_tpc_Cnigropunctatum$se.fit
preddata_tpc_Cnigropunctatum$lower <- pred_tpc_Cnigropunctatum$fit-pred_tpc_Cnigropunctatum$se.fit
preddata_tpc_Cnigropunctatum$upper <- pred_tpc_Cnigropunctatum$fit + pred_tpc_Cnigropunctatum$se.fit

quartz(8,8)

ggplot(preddata_tpc_Cnigropunctatum, aes(x=temp, y=predicted)) + geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)+
  labs(title='Copeoglossum nigropunctatum',y='Predicted Speed (cm/s)',x='Temperature (°C)') +
  lims(x=c(10,40),y=c(0,15)) + theme(plot.title = element_text(hjust = 0.5))




#Matticolus
## Faz a predicao 
pred_tpc_Matticolus <- predict(tpc_Matticolus$gam,
                               preddata_tpc_Matticolus,
                               se.fit=T)

## Anexa as predicoes e erros na tabela de dados
preddata_tpc_Matticolus$predicted <- pred_tpc_Matticolus$fit
preddata_tpc_Matticolus$se <- pred_tpc_Matticolus$se.fit
preddata_tpc_Matticolus$lower <- pred_tpc_Matticolus$fit-pred_tpc_Matticolus$se.fit
preddata_tpc_Matticolus$upper <- pred_tpc_Matticolus$fit + pred_tpc_Matticolus$se.fit

quartz(8,8)

ggplot(preddata_tpc_Matticolus, aes(x=temp, y=predicted)) + geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)+
  labs(title='Micrablepharus atticolus',y='Predicted Speed (cm/s)',x='Temperature (°C)') +
  lims(x=c(10,50),y=c(0,10)) + theme(plot.title = element_text(hjust = 0.5))

#Titambere
## Faz a predicao 
pred_tpc_Titambere <- predict(tpc_Titambere$gam,
                              preddata_tpc_Titambere,
                              se.fit=T)

## Anexa as predicoes e erros na tabela de dados
preddata_tpc_Titambere$predicted <- pred_tpc_Titambere$fit
preddata_tpc_Titambere$se <- pred_tpc_Titambere$se.fit
preddata_tpc_Titambere$lower <- pred_tpc_Titambere$fit-pred_tpc_Titambere$se.fit
preddata_tpc_Titambere$upper <- pred_tpc_Titambere$fit + pred_tpc_Titambere$se.fit

quartz(8,8)

ggplot(preddata_tpc_Titambere, aes(x=temp, y=predicted)) + geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)+
  labs(title='Tropidurus itambere',y='Predicted Speed (cm/s)',x='Temperature (°C)') +
  lims(x=c(10,45),y=c(0,15)) + theme(plot.title = element_text(hjust = 0.5))



#################
#Leaf Area Index#
#################
#install.packages("MODISTools", dep=T)

mt_products()
mt_bands(product = "MCD15A2H")


pts.arms <- read.table("Pontos_Armadilhas.txt", h=T)
pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="5"]


C_lai <- mt_subset(product = "MCD15A2H",
                   lat = pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="47"],
                   lon = pts.arms$long[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="47"],
                   band = "Lai_500m",
                   start = "2005-11-01",
                   end = "2019-12-31",
                   site_name = "C",
                   internal = T,
                   progress = T)

tail(C_lai)

Q_lai <- mt_subset(product = "MCD15A2H",
                   lat = pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="15"],
                   lon = pts.arms$long[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="15"],
                   band = "Lai_500m",
                   start = "2005-11-01",
                   end = "2019-12-31",
                   site_name = "Q",
                   internal = T,
                   progress = T)

tail(Q_lai)

BP_lai <- mt_subset(product = "MCD15A2H",
                    lat = pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="25"],
                    lon = pts.arms$long[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="25"],
                    band = "Lai_500m",
                    start = "2005-11-01",
                    end = "2019-12-31",
                    site_name = "BP",
                    internal = T,
                    progress = T)

tail(BP_lai)

BM_lai <- mt_subset(product = "MCD15A2H",
                    lat = pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="35"],
                    lon = pts.arms$long[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="35"],
                    band = "Lai_500m",
                    start = "2005-11-01",
                    end = "2019-12-31",
                    site_name = "BM",
                    internal = T,
                    progress = T)

tail(BM_lai)

BT_lai <- mt_subset(product = "MCD15A2H",
                    lat = pts.arms$lat[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="5"],
                    lon = pts.arms$long[pts.arms$local=="RECOR_FOGO" & pts.arms$armadilha=="5"],
                    band = "Lai_500m",
                    start = "2005-11-01",
                    end = "2019-12-31",
                    site_name = "BT",
                    internal = T,
                    progress = T)

tail(BT_lai)

all.equal(Q_lai$value, BP_lai$value)
all.equal(Q_lai$value, BT_lai$value)
all.equal(BM_lai$value, BP_lai$value)

sapply(list(C_lai$value,
            Q_lai$value,
            BP_lai$value,
            BM_lai$value,
            BT_lai$value), summary)*0.1

sapply(list(C_lai$value,
            Q_lai$value,
            BP_lai$value,
            BM_lai$value,
            BT_lai$value), sd)*0.1

lai.fire <- rbind.data.frame(C_lai,
                             Q_lai,
                             BP_lai,
                             BM_lai,
                             BT_lai)

lai.fire$plot <- c(rep("C", nrow(C_lai)),
                   rep("Q", nrow(Q_lai)),
                   rep("BP", nrow(BP_lai)),
                   rep("BM", nrow(BM_lai)),
                   rep("BT", nrow(BT_lai)))

lai.fire$year<- as.integer(substr(lai.fire$calendar_date,1,4))
lai.fire$month<- as.integer(substr(lai.fire$calendar_date,6,7))

summary(lai.fire)

lai.fire.month <-
  lai.fire %>%
  dplyr::group_by(plot, year, month) %>% 
  dplyr::summarise(lai = mean(value/10))
summary(lai.fire.month)
tapply(lai.fire.month$lai,lai.fire.month$plot,FUN=summary)

quartz(8,8)
interaction.plot(lai.fire.month$month, lai.fire.month$plot, lai.fire.month$lai)

quartz(8,8)
boxplot(lai ~ plot, data = lai.fire.month)

quartz(8,8)
interaction.plot(lai.fire.month$year, lai.fire.month$plot, lai.fire.month$lai)

saveRDS(lai.fire.month,"lai_fire_month.rds")

##############
#Microclimate#
##############


######################
#Mechanistic approach#
######################


# get ERA5 climate data with package mcera5 (just do once for region and time of interest)

# assign your credentials (register here: https://cds.climate.copernicus.eu/user/register)

uid <- "86619"
cds_api_key <- "151ac044-f8cf-4ee0-aa00-64da02cf0295"

ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")

pts.arms <- read.table("Pontos_Armadilhas.txt", h=T)

pts.arms$lat[pts.arms$local=="RECOR_FOGO"]
# bounding coordinates (in WGS84 / EPSG:4326)

(xmn <- min(pts.arms$long[pts.arms$local=="RECOR_FOGO"]))
(xmx <- max(pts.arms$long[pts.arms$local=="RECOR_FOGO"]))
(ymn <- min(pts.arms$lat[pts.arms$local=="RECOR_FOGO"]))
(ymx <- max(pts.arms$lat[pts.arms$local=="RECOR_FOGO"]))

# temporal extent
st_time <- lubridate::ymd("2005:11:01")
en_time <- lubridate::ymd("2019:12:31")

# filename and location for downloaded .nc files
file_prefix <- "era5"
op <- "/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis/Spatial_Data"

# build a request (covering multiple years)
req <- build_era5_request(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx,
                          start_time = st_time, #if gives error, update with new start time to download the rest of the years, lubridate::ymd("2019:01:01"),
                          end_time = en_time,
                          outfile_name = file_prefix)
str(req)
request_era5(request = req, uid = uid, out_path = op)


# list the path of the .nc file for a given year
my_nc<-list.files(pattern = "era5", recursive = TRUE)
my_nc
# for a single point (make sure it's within the bounds of your .nc file)
(x <- mean(c(xmn,xmx)))
(y <- mean(c(ymn,ymx)))

st_time_list <- c(lubridate::ymd("2005:11:01"),
                  lubridate::ymd("2006:01:01"),
                  lubridate::ymd("2007:01:01"),
                  lubridate::ymd("2008:01:01"),
                  lubridate::ymd("2009:01:01"),
                  lubridate::ymd("2010:01:01"),
                  lubridate::ymd("2011:01:01"),
                  lubridate::ymd("2012:01:01"),
                  lubridate::ymd("2013:01:01"),
                  lubridate::ymd("2014:01:01"),
                  lubridate::ymd("2015:01:01"),
                  lubridate::ymd("2016:01:01"),
                  lubridate::ymd("2017:01:01"),
                  lubridate::ymd("2018:01:01"),
                  lubridate::ymd("2019:01:01"))

en_time_list <- c(lubridate::ymd("2005:12:31"),
                  lubridate::ymd("2006:12:31"),
                  lubridate::ymd("2007:12:31"),
                  lubridate::ymd("2008:12:31"),
                  lubridate::ymd("2009:12:31"),
                  lubridate::ymd("2010:12:31"),
                  lubridate::ymd("2011:12:31"),
                  lubridate::ymd("2012:12:31"),
                  lubridate::ymd("2013:12:31"),
                  lubridate::ymd("2014:12:31"),
                  lubridate::ymd("2015:12:31"),
                  lubridate::ymd("2016:12:31"),
                  lubridate::ymd("2017:12:31"),
                  lubridate::ymd("2018:12:31"),
                  lubridate::ymd("2019:12:31"))

length(my_nc)
length(st_time_list)
length(en_time_list)
point_out <- list()
point_out_precip <- list()

for(i in 1:length(my_nc)){
  # gather all hourly variables
  point_out[[i]] <- extract_clim(nc = my_nc[i], long = x, lat = y,
                            start_time = st_time_list[i],  end_time = en_time_list[i])
  
  # gather daily precipitation (we specify to convert precipitation from hourly
  # to daily, which is already the default behavior)
  point_out_precip[[i]] <- extract_precip(nc = my_nc[i], long = x, lat = y,
                                     start_time = st_time_list[i],  
                                     end_time = en_time_list[i],
                                     convert_daily = TRUE)
  
}

point_out.df <- point_out %>% reduce(rbind.data.frame) 
summary(point_out.df)

point_out_precip.df <- point_out_precip %>% reduce(c) 
summary(point_out_precip.df)
length(point_out_precip.df)==nrow(point_out.df)

saveRDS(point_out.df,"microclima_vars_df.rds")
saveRDS(point_out_precip.df, "microclima_precip.rds")

point_out <-readRDS("microclima_vars_df.rds")

point_out_precip <- readRDS("microclima_precip.rds")
point_out_precip.df <- data.frame(precip = point_out_precip,
                                  obs_time = seq(as.Date('2005-11-01'),
                                                 as.Date('2019-12-31'),'days'))

nrow(filter(point_out,obs_time >= "2005-11-01" & obs_time < "2006-01-01"))
summary(filter(point_out_precip.df,obs_time >= "2005-11-01" & obs_time <= "2005-12-31")$precip)

#We need some canopy cover estimate to predict microclimate data in different fire regimes
#We use the leaf area index we downloaded previously and relate to canopy cover measures from 2009 and 2010

lai.fire.month <- readRDS("lai_fire_month.rds") #MacOS
head(lai.fire.month)

lai.month <- lai.fire.month %>% 
  group_by(year, month) %>%  
  summarise(lai = mean(lai))


#Canopy cover measures
MH.df <- read.table("MHbase.txt",h=T)
head(MH.df)

#windows(10,10)#Windows
quartz(8,8)#MacOS
with(MH.df,{
  interaction.plot(Campanha,Parcela,Ceu/25)
})

tapply(MH.df$Ceu/25, list(MH.df$Parcela,MH.df$Campanha), mean)*100
summary(lm(Ceu/25 ~ Parcela , data=MH.df))
tapply(lai.fire.month$lai,list(lai.fire.month$plot,lai.fire.month$month),mean)



(prop.lai.fire <- scale(tapply(MH.df$Ceu/25, list(MH.df$Parcela), mean),scale=F))

lai.month$lai.C <- lai.month$lai + (lai.month$lai *  prop.lai.fire[4])

lai.month$lai.Q <- lai.month$lai + (lai.month$lai * prop.lai.fire[5])

lai.month$lai.BP <-  lai.month$lai + (lai.month$lai *  prop.lai.fire[2])

lai.month$lai.BM <- lai.month$lai + (lai.month$lai *  prop.lai.fire[1])

lai.month$lai.BT <- lai.month$lai + (lai.month$lai *  prop.lai.fire[3])


# windows(10,10)
quartz(8,8)
plot(lai.month$lai,ylim = c(0.8,3),type="l")
lines(lai.month$lai.C, col="darkgreen")
lines(lai.month$lai.Q, col="blue")
lines(lai.month$lai.BP, col="yellow")
lines(lai.month$lai.BM, col="orange")
lines(lai.month$lai.BT, col="red")

MH.df$year <- c(rep(2009,50*3),rep(2010,50*9))
MH.df$month <- c(rep(10,50),rep(11,50),rep(12,50),rep(1,50),rep(2,50),rep(3,50),rep(4,50),
                 rep(5,50),rep(6,50),rep(7,50),rep(8,50),rep(9,50))
MH.df$CeuAb <- 25 - MH.df$Ceu 
y <- cbind(MH.df$Ceu, MH.df$CeuAb)
y <- round(y,0)
MH.lai.df <- left_join(MH.df, lai.month,by=c("month","year"))

MH.lai.df$Parcela <- factor(MH.lai.df$Parcela, levels = c("C","Q", "BP", "BM", "BT"))

glm1 <- glm(y ~ lai + Parcela, data=MH.lai.df, family="binomial")
summary(glm1)

quartz(10,10)
plot(Ceu/25~lai,data=MH.lai.df,col = as.numeric(as.factor(MH.lai.df$Parcela)))
xx <- seq(0,3,.1)
lines(xx,predict(glm1,newdata = data.frame(Parcela = rep("BM",31),lai =xx),type="response"),col=1)
lines(xx,predict(glm1,newdata = data.frame(Parcela = rep("BP",31),lai =xx),type="response"),col=2)
lines(xx,predict(glm1,newdata = data.frame(Parcela = rep("BT",31),lai =xx),type="response"),col=3)
lines(xx,predict(glm1,newdata = data.frame(Parcela = rep("C",31),lai =xx),type="response"),col=4)
lines(xx,predict(glm1,newdata = data.frame(Parcela = rep("Q",31),lai =xx),type="response"),col=5)

gam1 <- gam(y ~ s(lai) + Parcela, data=MH.lai.df, family="binomial")
summary(gam1)

quartz(8,8)
plot(gam1)
vis.gam(gam1,theta=60)


quartz(8,8)
plot_model(glm1,type = "pred",terms = c("lai","Parcela"),colors = "warm")



lai.fire.hour<-data.frame(obs_time = seq(as.POSIXlt('2005-11-01'),
               as.POSIXlt('2020-01-01'),'hour')[-124177])
lai.fire.hour$year <- year(lai.fire.hour$obs_time)
lai.fire.hour$month <- month(lai.fire.hour$obs_time)

lai.fire.hour <-left_join(lai.fire.hour,lai.month)

lai.fire.hour$lai.C <- lai.fire.hour$lai + (lai.fire.hour$lai *  prop.lai.fire[4])

lai.fire.hour$lai.Q <- lai.fire.hour$lai + (lai.fire.hour$lai * prop.lai.fire[5])

lai.fire.hour$lai.BP <-  lai.fire.hour$lai + (lai.fire.hour$lai *  prop.lai.fire[2])

lai.fire.hour$lai.BM <- lai.fire.hour$lai + (lai.fire.hour$lai *  prop.lai.fire[1])

lai.fire.hour$lai.BT <- lai.fire.hour$lai + (lai.fire.hour$lai *  prop.lai.fire[3])


#Canopy

canopy.C <- predict(glm1,data.frame(lai = lai.fire.hour$lai,
                                    Parcela = rep("C",nrow(lai.fire.hour))),
                    type="response")
canopy.Q <- predict(glm1,data.frame(lai = lai.fire.hour$lai,
                                    Parcela = rep("Q",nrow(lai.fire.hour))),
                    type="response")
canopy.BP <- predict(glm1,data.frame(lai = lai.fire.hour$lai,
                                     Parcela = rep("BP",nrow(lai.fire.hour))),
                     type="response")
canopy.BM <- predict(glm1,data.frame(lai = lai.fire.hour$lai,
                                     Parcela = rep("BM",nrow(lai.fire.hour))),
                     type="response")
canopy.BT <- predict(glm1,data.frame(lai = lai.fire.hour$lai,
                                     Parcela = rep("BT",nrow(lai.fire.hour))),
                     type="response")


canopy.fire <- cbind.data.frame(canopy.C, canopy.Q, canopy.BP, canopy.BM, canopy.BT)

canopy.fire$obs_time<- seq(as.POSIXlt('2005-11-01'),
                                         as.POSIXlt('2020-01-01'),'hour')[-124177]
canopy.fire$year <- year(canopy.fire$obs_time)
canopy.fire$month <- month(canopy.fire$obs_time)
canopy.fire$day <- day(canopy.fire$obs_time)

canopy.fire.day <- canopy.fire %>%
  dplyr::group_by(year, month, day) %>%
  dplyr::summarise(canopy.C = mean(canopy.C),
                   canopy.Q = mean(canopy.Q),
                   canopy.BP = mean(canopy.BP),
                   canopy.BM = mean(canopy.BM),
                   canopy.BT = mean(canopy.BT))

saveRDS(canopy.fire.day,"canopy.fire.day.rds")

#Using NichemapR micro_era5 function#
#####################################
# run micro_era5 for a location (make sure it's within the bounds of your .nc files)
#Ran in HPC Oxford cluster

dstart <- "01/11/2005"
dfinish <- "31/12/2019"
(x <- mean(c(xmn,xmx)))
(y <- mean(c(ymn,ymx)))
loc <- c(x, y)

micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish,
                  spatial = '~/Documents/GitHub/DemographicResilience_CerradoLizards/Ecophysio/Spatial_Data/era5')
 
saveRDS(micro,'micro_era5.rds')

micro.C<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish,
                    minshade = canopy.fire.day$canopy.C*100,
                    spatial = 'Spatial_Data/era5')

saveRDS(micro.C,'micro_C_era5.rds')

micro.Q<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish,
                    minshade = canopy.fire.day$canopy.Q*100,runshade = 0,
                    spatial = 'Spatial_Data/era5')

saveRDS(micro.Q,'micro_Q_era5.rds')

micro.BP<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, minshade = canopy.fire.day$canopy.BP*100,runshade = 0,
                     spatial = 'Spatial_Data/era5')

saveRDS(micro.BP,'micro_BP_era5.rds')

micro.BM<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, minshade = canopy.fire.day$canopy.BM*100,runshade = 0,
                     spatial = 'Spatial_Data/era5')

saveRDS(micro.BM,'micro_BM_era5.rds')

micro.BT<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, minshade = canopy.fire.day$canopy.BT*100,runshade = 0,
                     spatial = 'Spatial_Data/era5')

saveRDS(micro.BT,'micro_BT_era5.rds')

micro<-readRDS('micro_era5.rds')
micro.C<-readRDS('micro_C_era5.rds')
micro.Q<-readRDS('micro_Q_era5.rds')
micro.BP<-readRDS('micro_BP_era5.rds')
micro.BM<-readRDS('micro_BM_era5.rds')
micro.BT<-readRDS('micro_BT_era5.rds')


metmin<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
metmax <-as.data.frame(micro$shadmet)#above ground microclimatic conditions, max shade
soilmin<-as.data.frame(micro$soil) # soil temperatures, minimum shade
soilmax <- as.data.frame(micro$shadsoil)# soil temperatures, max shade
met.C<-as.data.frame(micro.C$met) # above ground microclimatic conditions, C fire regime
soil.C <- as.data.frame(micro.C$soil)# soil temperatures, C fire regime
met.Q<-as.data.frame(micro.Q$met) # above ground microclimatic conditions, Q fire regime
soil.Q <- as.data.frame(micro.Q$soil)# soil temperatures, Q fire regime
met.BP<-as.data.frame(micro.BP$met) # above ground microclimatic conditions, BP fire regime
soil.BP <- as.data.frame(micro.BP$soil)# soil temperatures, BP fire regime
met.BM<-as.data.frame(micro.BM$met) # above ground microclimatic conditions, BM fire regime
soil.BM <- as.data.frame(micro.BM$soil)# soil temperatures, BM fire regime
met.BT<-as.data.frame(micro.BT$met) # above ground microclimatic conditions, BT fire regime
soil.BT <- as.data.frame(micro.BT$soil)# soil temperatures, BT fire regime


metmin <- metmin[,-c(7:12,15:19)]
metmax <- metmax[,-c(7:12,15:19)]
met.C <- met.C[,-c(7:12,15:19)]
met.Q <- met.Q[,-c(7:12,15:19)]
met.BP <- met.BP[,-c(7:12,15:19)]
met.BM <- met.BM[,-c(7:12,15:19)]
met.BT <- met.BT[,-c(7:12,15:19)]

soilmin <- soilmin[-c(7:12)]
soilmax <- soilmax[-c(7:12)]
soil.C <- soil.C[-c(7:12)]
soil.Q <- soil.Q[-c(7:12)]
soil.BP <- soil.BP[-c(7:12)]
soil.BM <- soil.BM[-c(7:12)]
soil.BT <- soil.BT[-c(7:12)]


sapply(list(metmin,metmax,met.C,met.Q,met.BP,met.BM,met.BT),summary)
sapply(list(soilmin,soilmax,soil.C,soil.Q,soil.BP,soil.BM,soil.BT),summary)

# append dates
tzone<-paste("Etc/GMT",0,sep="")
dates<-seq(as.POSIXlt('2005-11-01'),
           as.POSIXlt('2020-01-01'),'hour')[-124177]


micro.0 <- data.frame(dates,metmin,soilmin,canopy = rep("min",nrow(metmin)))
micro.90 <- data.frame(dates,metmax,soilmax,canopy = rep("max",nrow(metmax)))
micro.fire.C <- data.frame(dates,met.C,soil.C,canopy = rep("C",nrow(met.C)))
micro.fire.Q <- data.frame(dates,met.Q,soil.Q,canopy = rep("Q",nrow(met.Q)))
micro.fire.BP <- data.frame(dates,met.BP,soil.BP,canopy = rep("BP",nrow(met.BP)))
micro.fire.BM <- data.frame(dates,met.BM,soil.BM,canopy = rep("BM",nrow(met.BM)))
micro.fire.BT <- data.frame(dates,met.BT,soil.BT,canopy = rep("BT",nrow(met.BT)))

micro.fire <- rbind(micro.90, micro.fire.C, micro.fire.Q,
                    micro.fire.BP, micro.fire.BM, micro.fire.BT, micro.0)

micro.fire <- micro.fire[,-c(10,11)]

quartz(8,8)

interaction.plot(micro.fire$TIME, micro.fire$canopy, micro.fire$TALOC)

micro.fire$year <- year(micro.fire$dates)
micro.fire$month <- month(micro.fire$dates)
micro.fire$hour <- hour(micro.fire$dates)
summary(micro.fire)


quartz(8,8)
boxplot(TALOC ~ canopy, data=micro.fire)


#setwd("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/Microclima/")
climate.df <- readRDS("climate.df.rds")

microclim.BSB <- readRDS("microclimate_BSB.rds")
summary(microclim.BSB,2)


quartz(8,8)
par(mfrow=c(3,1))
interaction.plot(microclim.BSB$hour, microclim.BSB$model_type,microclim.BSB$temp)
interaction.plot(microclim.BSB$hour, microclim.BSB$microhabitat,microclim.BSB$temp)
interaction.plot(microclim.BSB$hour, microclim.BSB$plot,microclim.BSB$temp)
summary(microclim.BSB$temp)


## Predicoes para o periodo de estudo

pred_temp <-
  data.frame(temp = micro.fire$TALOC)

micro.fire$Cnigro_perf <-
  predict.gam(tpc_Cnigropunctatum$gam, newdata = pred_temp)

micro.fire$Matticolus_perf <-
  predict.gam(tpc_Matticolus$gam, newdata = pred_temp)

micro.fire$Titambere_perf <-
  predict.gam(tpc_Titambere$gam, newdata = pred_temp)


summary(micro.fire)

# Horas de atividade




tpref_BSB <- read.table("Tpref_BSB.txt", h=T)
str(tpref_BSB)

tpref_Cnigropunctatum <- dplyr::filter(tpref_BSB, Species == "C_nigropunctatum")
tpref_Matticolus <- dplyr::filter(tpref_BSB, Species == "M_atticolus")
tpref_Titambere <- dplyr::filter(tpref_BSB, Species == "T_itambere")


hvtFUN <- function(temp.envr,temp.lab, quantiles, radiation) {
  vtmin <- quantile(temp.lab, quantiles[1], na.rm = TRUE)
  vtmax <- quantile(temp.lab, quantiles[2], na.rm = TRUE)
  hv <- ifelse(temp.envr > vtmin & temp.envr < vtmax, 1, 0)
  hv[radiation == 0] <- 0
  hv
}

micro.fire$Cnigro_ha90 <- hvtFUN(micro.fire$TALOC,
                                 tpref_Cnigropunctatum$Temp,
                                 c(0.05,0.95),
                                 micro.fire$SOLR)

micro.fire$Matticolus_ha90 <- hvtFUN(micro.fire$TALOC,
                                     tpref_Matticolus$Temp,
                                     c(0.05,0.95),
                                     micro.fire$SOLR)

micro.fire$Titambere_ha90 <- hvtFUN(micro.fire$TALOC,
                                    tpref_Titambere$Temp,
                                    c(0.05,0.95),
                                    micro.fire$SOLR)



quartz(8,8)
interaction.plot(micro.fire$month,
                 micro.fire$canopy,
                 micro.fire$Cnigro_ha90, mean)

quartz(8,8)
interaction.plot(micro.fire$month,
                 micro.fire$canopy,
                 micro.fire$Matticolus_ha90, mean)

quartz(8,8)
interaction.plot(micro.fire$month,
                 micro.fire$canopy,
                 micro.fire$Titambere_ha90, mean)


saveRDS(micro.fire, "climate_ecophysio.rds")

climate.ecophysio <- readRDS("climate_ecophysio.rds")
summary(climate.ecophysio)

climate.ecophysio <- dplyr::filter(climate.ecophysio, canopy != "min")
climate.ecophysio <- dplyr::filter(climate.ecophysio, canopy != "max")

climate.ecophysio$day <- day(climate.ecophysio$dates)

climate.ecophysio.day <-
  climate.ecophysio %>%
  dplyr::group_by(canopy, year, month, day) %>% 
  dplyr::summarise(tmed = mean(TALOC),
                   tmin = min(TALOC),
                   tmax = max(TALOC),
                   tmed2m = mean(TAREF),
                   tmin2m = min(TAREF),
                   tmax2m = max(TAREF),
                   RH = mean(RHLOC),
                   RHmin = min(RHLOC),
                   RHmax = max(RHLOC),
                   sol = mean(SOLR),
                   tmed0cm = mean(D0cm),
                   tmin0cm = min(D0cm),
                   tmax0cm = max(D0cm),
                   Cnigro_perf = mean(Cnigro_perf),
                   Matticolus_perf = mean(Matticolus_perf),
                   Titambere_perf = mean(Titambere_perf),
                   Cnigro_ha90 = sum(Cnigro_ha90),
                   Matticolus_ha90 = sum(Matticolus_ha90),
                   Titambere_ha90 = sum(Titambere_ha90)
                   )

climate.ecophysio.month <-
  climate.ecophysio.day %>%
  dplyr::group_by(canopy, year, month) %>% 
  dplyr::summarise(tmed = mean(tmed),
                   tmin = mean(tmin),
                   tmax = mean(tmax),
                   tmed2m = mean(tmed2m),
                   tmin2m = mean(tmin2m),
                   tmax2m = mean(tmax2m),
                   RH = mean(RH),
                   RHmin = mean(RHmin),
                   RHmax = mean(RHmax),
                   sol = mean(sol),
                   tmed0cm = mean(tmed0cm),
                   tmin0cm = mean(tmin0cm),
                   tmax0cm = mean(tmax0cm),
                   Cnigro_perf = mean(Cnigro_perf),
                   Matticolus_perf = mean(Matticolus_perf),
                   Titambere_perf = mean(Titambere_perf),
                   Cnigro_ha90 = sum(Cnigro_ha90),
                   Matticolus_ha90 = sum(Matticolus_ha90),
                   Titambere_ha90 = sum(Titambere_ha90)
  )

point_out_precip.df$month <- month(point_out_precip.df$obs_time)
point_out_precip.df$year <- year(point_out_precip.df$obs_time)

precip.month.df <-
  point_out_precip.df %>%
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(precip = sum(precip))

nrow(precip.month.df);nrow(climate.ecophysio.month)
climate.ecophysio.month <- left_join(climate.ecophysio.month, precip.month.df,
                                     by = c("month","year"))
summary(climate.ecophysio.month)

vifcor(as.data.frame(climate.ecophysio.month[,c(4:16,23)]),th=.8) #take off tmin, tmax0cm, RH, tmin2m, tmed, tmax, and RHmin

m <- cor(climate.ecophysio.month[,c(4:16,23)])
quartz(8, 8)
corrplot.mixed(m, order = 'AOE')

climate.ecophysio.month <- subset(climate.ecophysio.month,
                                  select = -c(tmin,tmax0cm, RH, tmin2m,
                                              tmed,tmax,tmax2m,RHmin))

saveRDS(climate.ecophysio.month,"climate_ecophysio_month.rds")

