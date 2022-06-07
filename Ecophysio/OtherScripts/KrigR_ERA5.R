Sys.setenv(TZ = "America/Sao_Paulo")#set time zone
setwd("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/")
# 



##################################
#Download climate data from ERA 5#
##################################

#Load package KrigR 
#devtools::install_github("https://github.com/ErikKusch/KrigR")
library(KrigR)

#Check if other packages are instalee
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

package_vec <- c(
  "tidyr", # for turning rasters into ggplot-dataframes
  "ggplot2", # for plotting
  "viridis", # colour palettes
  "cowplot", # gridding multiple plots
  "rosm", # obtaining satellite maps
  "gimms" # to get some pre-existing data to match in our downscaling
)
sapply(package_vec, install.load.package)

#rasterOptions(tmpdir=file.path("D:/Documentos/Raster")) 


#Before we can proceed, you need to open up an account at the CDS and create 
#an API key by following this link: https://cds.climate.copernicus.eu/api-how-to. 
#This is required so that you may issue download requests through the KrigR package. 

#Register the user ID and API Key as characters 

API_User <- "86619"
API_Key <- "151ac044-f8cf-4ee0-aa00-64da02cf0295"

#A data directory for all of our individual Kriging processes
#A shapefile directory (located within our data directory) for 
#all of the shapefiles that we will use

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data
Dir.Shapes <- file.path(Dir.Data, "Shapes") # folder path for shapefiles
Dirs <- sapply(c(Dir.Data, Dir.Shapes), function(x) if (!dir.exists(x)) dir.create(x))

#Read and plot points coordinates
arms.pts <- read.table("Pontos_Armadilhas.txt",h=T)
arms.pts

(arms.pts.RECOR <- subset(arms.pts,local=="RECOR_FOGO"))
coordinates(arms.pts.RECOR) <- c("long","lat")
proj4string(arms.pts.RECOR) <- CRS("+proj=longlat +datum=WGS84")
arms.pts.RECOR
arms.pts.RECOR<-st_as_sf(arms.pts.RECOR)

#Read and plot UCs area
UCsMask <- readOGR(Dir.Shapes,"UCs_Armadilhas_Repteis", verbose = FALSE) # read
UCsMask
RECORMask <- subset(UCsMask, WDPA_PID=="478492")

quartz(8,8)
ggplot() +
  geom_sf(data = st_as_sf(RECORMask), colour = "darkred") +
  geom_sf(data=arms.pts.RECOR)+
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")


##############
#Climate Data#
##############
Extent <- extent(arms.pts.RECOR)
quartz(8,8)
bmaps.plot(bbox = RECORMask, type = "AerialWithLabels", quiet = TRUE, progress = "none", colour = "darkred")

Dir.RECORExt <- file.path(Dir.Data, "RECOR_Extent")
#dir.create(Dir.RECORExt)

# 
# RECOR_Raw_temp <- download_ERA(
#   Variable = "2m_temperature",
#   DataSet = "era5",
#   Type = "monthly_averaged_reanalysis_by_hour_of_day",
#   DateStart = "2005-11-01",
#   DateStop = "2019-12-31",
#   TResolution = "hour",
#   Extent = RECORMask,
#   Dir = Dir.RECORExt,
#   API_User = API_User,
#   API_Key = API_Key
# )
# 
# RECOR_Raw_temp
# 
# Cerrado_Raw_skin_temp <- download_ERA(
#   Variable = "skin_temperature",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_skin_temp
# 
# Cerrado_Raw_temp_lv1 <- download_ERA(
#   Variable = "soil_temperature_level_1",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_temp_lv1
# 
# Cerrado_Raw_temp_lv2 <- download_ERA(
#   Variable = "soil_temperature_level_1",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_temp_lv2
# 
# Cerrado_Raw_sol <- download_ERA(
#   Variable = "surface_solar_radiation_downwards",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_sol
# 
# Cerrado_Raw_water_lv1 <- download_ERA(
#   Variable = "volumetric_soil_water_layer_1",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_water_lv1
# 
# Cerrado_Raw_water_lv2 <- download_ERA(
#   Variable = "volumetric_soil_water_layer_1",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_water_lv2
# 
# Cerrado_Raw_water_veg <- download_ERA(
#   Variable = "skin_reservoir_content",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_water_veg
# 
# Cerrado_Raw_precip <- download_ERA(
#   Variable = "total_precipitation",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_precip
# 
# Cerrado_Raw_dewpoint_temp <- download_ERA(
#   Variable = "2m_dewpoint_temperature",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_dewpoint_temp
# 
# Cerrado_Raw_surf_press <- download_ERA(
#   Variable = "surface_pressure",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_surf_press
# 
# 
# Cerrado_Raw_wind_u <- download_ERA(
#   Variable = "10m_u_component_of_wind",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_wind_u
# 
# Cerrado_Raw_wind_v <- download_ERA(
#   Variable = "10m_v_component_of_wind",
#   DataSet = "era5-land",
#   DateStart = "1982-01-01",
#   DateStop = "2020-12-31",
#   TResolution = "month",
#   Extent = CerradoMask,
#   Dir = Dir.CerradoExt,
#   API_User = API_User,
#   API_Key = API_Key,
#   SingularDL = TRUE
# )
# 
# Cerrado_Raw_wind_v

#setwd("/Volumes/GoogleDrive/My Drive/Doutorado/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/Data/RECOR_Extent")

# files <- list.files(pattern = ".nc", recursive = TRUE)
# files
# nc.climate <- list()
# library(ncdf4)
# 
# #Read files
# for(i in 1:length(files)){
#   nc.climate[[i]] <- nc_open(files[i])
# }
# 
# 
# # #Save metadata files
#   for(i in 1:length(files)){
#     {
#   sink(paste0(files[i],"_metadata.txt"))
#   print(nc.climate[[i]])
#   sink()
# }
# }
# 
# lon <- ncvar_get(nc.climate[[1]], "longitude")
# lat <- ncvar_get(nc.climate[[1]], "latitude", verbose = F)
# t <- ncvar_get(nc.climate[[1]], "time")
# dim(t)
# 
# arrays.climate<-list()
# for(i in 1:length(files)){
#   arrays.climate[[i]] <- ncvar_get(nc.climate[[i]])
# }
# 
# dim(arrays.climate[[1]])
# 
# 
# 
# library(tidyverse)
# climate.df <- arrays.climate %>% reduce(cbind.data.frame)
# # replace netCDF fill values with NA's
# climate.df[climate.df == -32767] <- NA
# names(climate.df) <- files
# climate.df$hour <- rep(0:23,12*15)
# climate.df$month <- rep(c(rep(1,24),rep(2,24),rep(3,24),rep(4,24),rep(5,24),rep(6,24),
#                       rep(7,24),rep(8,24),rep(9,24),rep(10,24),rep(11,24),rep(12,24)),15)
# climate.df$year <- c(rep(2005,24*12), rep(2006,24*12), rep(2007,24*12),
#                      rep(2008,24*12), rep(2009,24*12), rep(2010,24*12),
#                      rep(2011,24*12), rep(2012,24*12), rep(2013,24*12),
#                      rep(2014,24*12), rep(2015,24*12), rep(2016,24*12),
#                      rep(2017,24*12), rep(2018,24*12), rep(2019,24*12))
# 
# climate.df$hour.GMT_3 <- rep(c(21:23,0:20),12*15)
# climate.df$month.GMT_3 <- rep(c(rep(12,3),rep(1,24),rep(2,24),rep(3,24),rep(4,24),rep(5,24),rep(6,24),
#                           rep(7,24),rep(8,24),rep(9,24),rep(10,24),rep(11,24),rep(12,21)),15)
# climate.df$year.GMT_3 <- c(rep(2004,3),rep(2005,24*12), rep(2006,24*12), rep(2007,24*12),
#                      rep(2008,24*12), rep(2009,24*12), rep(2010,24*12),
#                      rep(2011,24*12), rep(2012,24*12), rep(2013,24*12),
#                      rep(2014,24*12), rep(2015,24*12), rep(2016,24*12),
#                      rep(2017,24*12), rep(2018,24*12), rep(2019,24*11),rep(2019,21))
# 
# 
# 
# summary(climate.df)
# saveRDS(climate.df,"/Volumes/GoogleDrive/My Drive/Doutorado/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/Microclima/climate.df.rds")
