
##############
#Microclimate#
##############


library(lubridate)

#Le dados microclimaticos
# setwd("G:/Meu Drive/Doutorado/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/Microclima/HOBOs Cerrado/Toper")
# files<-list.files(pattern = ".txt", recursive = TRUE)
# files
# 
# dados.list <- list()
# 
# for(i in 1:length(files)){
#   dados.list[[i]] <- read.table(files[i],h=F,skip=2,sep=";")[,1:5]
#   dados.list[[i]]$day <- substr(dados.list[[i]]$V2,1,2)
#   dados.list[[i]]$month <- substr(dados.list[[i]]$V2,4,5)
#   dados.list[[i]]$year <- substr(dados.list[[i]]$V2,7,10)
#   dados.list[[i]]$V2 <- as.POSIXlt(dados.list[[i]]$V2,
#                                    format="%d/%m/%Y")
#   dados.list[[i]]$hour <- substr(dados.list[[i]]$V3,1,2)
#   dados.list[[i]]$file <- rep(files[i],nrow(dados.list[[i]]))
#   file.name <- str_split(dados.list[[i]]$file, '/', simplify = TRUE)
#   dados.list[[i]]$ID <- file.name[,ncol(file.name)]
#   if (max(dados.list[[i]]$month) > 12){
#     dados.list[[i]]$day1 <- dados.list[[i]]$month
#     dados.list[[i]]$month <- dados.list[[i]]$day
#     dados.list[[i]]$day <- dados.list[[i]]$day1
#     dados.list[[i]] <- dados.list[[i]][,-9]
#   }
# }
# 
# 
# library(reshape)
# library(tidyverse)
# dados.merged <- dados.list %>% reduce(rbind.data.frame)
# dados.merged$day <- as.numeric (dados.merged$day)
# dados.merged$month <- as.numeric (dados.merged$month)
# dados.merged$year <- as.numeric (dados.merged$year)
# dados.merged$hour <- as.numeric (dados.merged$hour)



setwd("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/Microclima")
# microclima.2013 <- read_delim("Microclima_IBGE_fogo_2013.csv", delim = ";")
# microclima.2014 <- read_delim("Microclima_IBGE_fogo_2014.csv", delim = ";")
# 
# microclima.2013_14 <- rbind.data.frame(microclima.2013, microclima.2014)
# levels(as.factor(microclima.2013_14$filename))
# 
# microclima.2013_14$plot <- rep(NA, nrow(microclima.2013_14))
# microclima.2013_14$plot[grep(pattern = "C",microclima.2013_14$filename)] <- "C"
# microclima.2013_14$plot[grep(pattern = "Q",microclima.2013_14$filename)] <- "Q"
# microclima.2013_14$plot[grep(pattern = "BP",microclima.2013_14$filename)] <- "BP"
# microclima.2013_14$plot[grep(pattern = "BM",microclima.2013_14$filename)] <- "BM"
# microclima.2013_14$plot[grep(pattern = "BT",microclima.2013_14$filename)] <- "BT"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo1.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo2.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo3.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo4.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo5.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo6.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo7.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo8.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo9.csv"] <-  "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo10.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo11.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo12.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo13.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo14.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo15.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo16.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo17.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo18.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo19.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo20.csv"] <- "C"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo21.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo22.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo23.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo23.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo24.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo25.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo26.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo27.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo28.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo29.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo30.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo31.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo32.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo33.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo34.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo35.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo36.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo37.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo38.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo39.csv"] <- "BM"
# microclima.2013_14$plot[microclima.2013_14$filename == "Fogo40.csv"] <- "BM"
# 
# 
# levels(as.factor(microclima.2013_14$plot))
# 
# summary(microclima.2013_14)
# 
# quartz(8,8)
# interaction.plot(microclima.2013_14$date, microclima.2013_14$filename, microclima.2013_14$temp)
# 
# quartz(8,8)
# boxplot(temp ~ plot, data = microclima.2013_14)
# 
# quartz(8,8)
# boxplot(t_air_mean ~ plot, data = microclima.2013_14[microclima.2013_14$month==6,])
# 
# quartz(8,8)
# boxplot(t_air_max ~ plot, data = microclima.2013_14)
# 
# quartz(8,8)
# boxplot(t_air_min ~ plot, data = microclima.2013_14)
# 
# quartz(8,8)
# plot(temp ~ t_air_mean, data = microclima.2013_14)
# 
# saveRDS(microclima.2013_14,
#         "microclimate_BSB.rds")

