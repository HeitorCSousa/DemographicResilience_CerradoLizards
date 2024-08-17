
# Hierarchical models to estimate vital rates and their relationship with the environmental variables-------------------------------------------------------------------------

## Copeoglossum nigropunctatum ----

#Load packages
library(runjags)
library(parallel)
library(rjags)
library(reshape)
library(ggplot2)
library(missForest)
library(lubridate)
library(ggmcmc)
library(MCMCvis)
library(viridis)


setwd("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis")

### Read data -----------------------------

#Read fecundity data
fecundidade.Cn <- read.csv("fecundidade_Cn.csv", h=T)
summary(fecundidade.Cn)

#Probability of reproduction

Cnigropunctatum.RECOR.imp<-readRDS("Cnigropunctatum_RECOR_imp.rds") #Read data
Cnigropunctatum.RECOR.females<-subset(Cnigropunctatum.RECOR.imp,sexo=="F")# subset only females
summary(as.factor(Cnigropunctatum.RECOR.females$ovos))

#We only recorded females bearing eggs/embryos. To become binary, we need to determine the females without eggs/embryos
Cnigropunctatum.RECOR.females$ovos[is.na(Cnigropunctatum.RECOR.females$ovos)]<-"n" 
table(Cnigropunctatum.RECOR.females$ovos) #number of observations (females with eggs/embryos)
head(Cnigropunctatum.RECOR.females)
Cnigropunctatum.RECOR.females$ovos <- as.numeric(as.factor(Cnigropunctatum.RECOR.females$ovos))-1

Cnigropunctatum.RECOR.females <- Cnigropunctatum.RECOR.females[,c("campanha","mes","ano","massa","crc","ovos")]#subset only columns of interest
Cnigropunctatum.RECOR.females <- na.omit(Cnigropunctatum.RECOR.females) #remove other lines than females
Cnigropunctatum.RECOR.females <- Cnigropunctatum.RECOR.females[Cnigropunctatum.RECOR.females$ano < 2020,] #Subset time frame
summary(Cnigropunctatum.RECOR.females)

## Clean demographic data ----------
Cnigropunctatum.RECOR.imp<-readRDS("Cnigropunctatum_RECOR_imp.rds") #read data

head(Cnigropunctatum.RECOR.imp)
tail(Cnigropunctatum.RECOR.imp)
dados.demografia <- Cnigropunctatum.RECOR.imp[Cnigropunctatum.RECOR.imp$ano < 2020, ] #remove data later than 2019
tail(dados.demografia)
str(dados.demografia)

#Remove dead animals
dados.demografia$morto[is.na(dados.demografia$morto)] <- "n"
dados.demografia <- dados.demografia[dados.demografia$morto=="n", ]
head(dados.demografia)
table(dados.demografia$campanha)
table(dados.demografia$identidade)

completos <- complete.cases(dados.demografia[, c("identidade", "parcela")]) #remove NAs
nigro.planilha <- droplevels(dados.demografia[completos, ])
head(nigro.planilha)
str(nigro.planilha)
table(nigro.planilha$campanha)
table(nigro.planilha$identidade)


## Prepare input file and run monthly data ("campanha")


# Create ID
IDENT <- paste(nigro.planilha$parcela, nigro.planilha$identidade, nigro.planilha$ciclo, sep="")
head(IDENT)
nigro.planilha <- data.frame(nigro.planilha, IDENT)
rm(IDENT)
str(nigro.planilha)

# Subset variables of interest
table1 <- nigro.planilha[nigro.planilha$recaptura!="(s)", c("IDENT", "campanha", "sexo", "crc","massa", "recaptura", "parcela")]
str(table1)

# Identifies captures without IDs
table2 <- complete.cases(table1[, c("IDENT")])

# Remove captures without ID and SVL
table3 <- table1[table2, ]

# Order the data by ID and time (campanha)
table4 <- table3[order(table3$IDENT, table3$campanha), ]
str(table4)
summary(table4)

# Convert IDENT from factor to character
table5 <- droplevels(table4)
table5$IDENT <- as.character(table5$IDENT)
table5$parcela <- as.character(table5$parcela)
str(table5)


# Calculate recapture frequencies
recap.table <- data.frame(table(table5$IDENT))
names(recap.table) <- c("identidade", "capturas")
recap.table
table(recap.table$capturas)
943+(2*182)+(3*64)+(4*15)+(5*2)+(6*1) #two captures=1, 3 captures=2, 4 captures=3 ...

# Filter data.frame records to use in the analysis

head(table5)

Age<-c(rep(NA,nrow(table5)))
Age

datA<-data.frame(table5$campanha,table5$sexo,table5$IDENT,table5$crc,table5$massa,Age,table5$parcela)
names(datA)<-c("Year","Sex","TrueID","CRC","Massa","Age","Parcela")
datA$Year<-datA$Year+2000
datA$Age[datA$CRC<=40]<-0

head(datA)
tail(datA)
str(datA)

##Creates variables for the growth model -------

###  del is the time period since the first capture of an individual (0 in the first capture)

del<-c()   ### months since the first capture

for(i in 1:nrow(datA)){
  del[i]<-datA$Year[i]-min(datA$Year[datA$TrueID==datA$TrueID[i]])
}

plot<-cast(datA, TrueID~., value="Parcela", fun.aggregate=function(x) tail(x,1))  ###determine the plor for each individual
plot<-as.character(plot[,2])
plot
plot[plot=="BM"]<-4 # Ordering by fire severity
plot[plot=="BP"]<-3 # Ordering by fire severity
plot[plot=="BT"]<-5 # Ordering by fire severity
plot[plot=="C"]<-1 # Ordering by fire severity
plot[plot=="Q"]<-2 # Ordering by fire severity
plot<-as.numeric(plot)
plot

sex <- cast(datA, TrueID~., value="Sex", fun.aggregate=function(x) tail(x,1))  ###determine the sex for each individual - filter the value of last capture
sex <- factor(sex[,2], levels = c("F", "M"))
sex
sex <- as.integer(sex)
sex
sex <- sex - 1

ind = as.numeric(factor(datA$TrueID)) #ID
y = datA$CRC #SVL
n = max(ind)  ### number of individuals
m = nrow(datA)### number of observations

age<- c()  ## age at first capture
for (a in 1:n){ age[a] <- datA$Age[ind==a][1]}

year <- c()
for (a in 1:n){ year[a] <- datA$Year[ind==a][1]}

head(datA)
tail(datA)
str(datA)

## Create capture histories for the survival model-------
known.states.cjs<-function(ch){
  state<-ch
  for (i in 1:dim(ch)[1]){
    n1<-min(which(ch[i,]==1))
    n2<-max(which(ch[i,]==1))
    state[i,n1:n2]<-1
    state[i,n1]<-NA
  }
  state[state==0]<-NA
  return(state)
}

cjs.init.z<-function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2<-max(which(ch[i,]==1))
    ch[i,f[i]:n2]<-NA
  }
  for (i in 1:dim(ch)[1])
  { ch[i,1:f[i]]<-NA
  }
  return(ch)
}

#Capture histories
eh <- cast(datA,TrueID ~ Year, fun.aggregate = function(x) as.numeric(length(x) >0),value="CRC");eh <- eh[,2:ncol(eh)]
eh.all <- seq(min(datA$Year), max(datA$Year)) #preencher todos os anos ignorados
missing <- eh.all[!(eh.all %in% names(eh))]
col=matrix(0,nrow=nrow(eh),ncol=length(missing))
colnames(col) <- missing
eh <- cbind(eh, col)
eh <- eh[,sort(colnames(eh))]
head(eh)
(f <- apply(eh,1,function(x) which(x==1)[1]))#time since first capture
(nind <- nrow(eh)) #Number of individuals
(n.occasions <- ncol(eh)) #Number of occasions
m #Number of observations
n #Number of individuals

#Create data.frame of plots
mplot <- data.frame(C = as.numeric(plot==1),
                    Q = as.numeric(plot==2),
                    BP = as.numeric(plot==3),
                    BM = as.numeric(plot==4),
                    BT = as.numeric(plot==5))

# Create matrix X indicating SVL (crc)
x <- cast(datA,
          TrueID ~ Year,
          fun.aggregate = function(x) mean(x),
          value = "CRC",
          fill = NA)

x <- x[, 2:ncol(x)]
x.all <- seq(min(datA$Year), max(datA$Year)) #fill all the ignored months
missing <- x.all[!(x.all %in% names(x))]
col=matrix(NA,nrow=nrow(x),ncol=length(missing))
colnames(col) <- missing
x <- cbind(x, col)
x <- x[,sort(colnames(x))]
#x[is.na(x)] <- 0
head(x)

# Create matrix X indicating mean svl
# xpop <- x
# xpop[is.na(xpop)] <- 0
# head(xpop)

xpopC  <- matrix(0,nrow = nind, ncol = n.occasions)
xpopQ  <- matrix(0,nrow = nind, ncol = n.occasions)
xpopBP <- matrix(0,nrow = nind, ncol = n.occasions)
xpopBM <- matrix(0,nrow = nind, ncol = n.occasions)
xpopBT <- matrix(0,nrow = nind, ncol = n.occasions)

for(i in 1:nind){
  for (t in 1:(n.occasions)){
    xpopC[i,t]  <- x[i,t] * mplot[i,1]
    xpopQ[i,t]  <- x[i,t] * mplot[i,2]
    xpopBP[i,t] <- x[i,t] * mplot[i,3]
    xpopBM[i,t] <- x[i,t] * mplot[i,4]
    xpopBT[i,t] <- x[i,t] * mplot[i,5]
  }
}

xpopC[xpopC==0] <- NA
xpopQ[xpopQ==0] <- NA
xpopBP[xpopBP==0] <- NA
xpopBM[xpopBM==0] <- NA
xpopBT[xpopBT==0] <- NA


(xpopmeanC <- colMeans(xpopC,na.rm=T))
(xpopmeanQ <- colMeans(xpopQ,na.rm=T))
(xpopmeanBP <- colMeans(xpopBP,na.rm=T))
(xpopmeanBM <- colMeans(xpopBM,na.rm=T))
(xpopmeanBT <- colMeans(xpopBT,na.rm=T))

xpop <- c(xpopmeanC,xpopmeanQ,xpopmeanBP,xpopmeanBM,xpopmeanBT)


# Environmental variables -------------------------------------------------

#Read data
# var_amb <- readRDS("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis/Ecophysio/climate.ecophysio.month.rds")
var_amb <- readRDS("climate_ecophysio_month.rds")
var_amb$canopy <- factor(var_amb$canopy,levels = c("C","Q","BP","BM","BT"))
var_amb <- var_amb[order(var_amb$canopy),]
var_amb$canopy
var_amb$xpop <- xpop

#Impute missing data
xpop.imp <- missForest::missForest(as.data.frame(var_amb,
                                                 verbose = T,
                                                 variablewise = T,
                                                 maxiter = 100,
                                                 ntree = 1000,
                                                 parallelize = "variables"))

names(xpop.imp)
xpop.imp$OOBerror

xpop.fire <- rbind(xpop.imp$ximp$xpop[xpop.imp$ximp$canopy=="C"] ,
                   xpop.imp$ximp$xpop[xpop.imp$ximp$canopy=="Q"] ,
                   xpop.imp$ximp$xpop[xpop.imp$ximp$canopy=="BP"],
                   xpop.imp$ximp$xpop[xpop.imp$ximp$canopy=="BM"],
                   xpop.imp$ximp$xpop[xpop.imp$ximp$canopy=="BT"])
xpop.fire

# Derive data for the Pradel-Jolly-Seber (PJS) Model
# e = index of oldest observation
get.first <- function (x) min(which (x!=0))
e <- apply (eh, 1, get.first)
e <- data.frame(plot=plot,e=e)


# l = index of last observation
get.last <- function(x) max(which(x!=0))
l <- apply (eh,1, get.last )
l <- data.frame(plot=plot,l=l)


# u = number of observed animals for the first time in i
u1 <- data.frame(table(e$e[e$plot==1]))
u2 <- data.frame(table(e$e[e$plot==2]))
u3 <- data.frame(table(e$e[e$plot==3]))
u4 <- data.frame(table(e$e[e$plot==4]))
u5 <- data.frame(table(e$e[e$plot==5]))

u.all <- seq(1, n.occasions) # Fill the ignored months
missing <- u.all[!(u.all %in% u1$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u1) <- names(df)
u1 <- rbind(u1, df)
u1$e<-as.numeric(levels(u1$e))[u1$e]
u1<-u1[order(u1$e), ]
u1<-u1$Freq
u1

missing <- u.all[!(u.all %in% u2$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u2) <- names(df)
u2 <- rbind(u2, df)
u2$e<-as.numeric(levels(u2$e))[u2$e]
u2<-u2[order(u2$e), ]
u2<-u2$Freq
u2

missing <- u.all[!(u.all %in% u3$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u3) <- names(df)
u3 <- rbind(u3, df)
u3$e<-as.numeric(levels(u3$e))[u3$e]
u3<-u3[order(u3$e), ]
u3<-u3$Freq
u3

missing <- u.all[!(u.all %in% u4$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u4) <- names(df)
u4 <- rbind(u4, df)
u4$e<-as.numeric(levels(u4$e))[u4$e]
u4<-u4[order(u4$e), ]
u4<-u4$Freq
u4

missing <- u.all[!(u.all %in% u5$Var1)]
df<-data.frame(e=as.factor(missing),Freq=rep(0,length(missing)))
names(u5) <- names(df)
u5 <- rbind(u5, df)
u5$e<-as.numeric(levels(u5$e))[u5$e]
u5<-u5[order(u5$e), ]
u5<-u5$Freq
u5

u <- rbind(u1,u2,u3,u4,u5)

# n = number of observed animals in i
n1 <- colSums(eh[plot==1,])
n2 <- colSums(eh[plot==2,])
n3 <- colSums(eh[plot==3,])
n4 <- colSums(eh[plot==4,])
n5 <- colSums(eh[plot==5,])

n<- rbind(n1,n2,n3,n4,n5)

colSums(n) == colSums(eh)

# v = number of observed animals for the last time in i
v1 <- data.frame(table(l$l[l$plot==1]))
v2 <- data.frame(table(l$l[l$plot==2]))
v3 <- data.frame(table(l$l[l$plot==3]))
v4 <- data.frame(table(l$l[l$plot==4]))
v5 <- data.frame(table(l$l[l$plot==5]))

v.all <- seq(1, n.occasions) # Fill the ignored months
missing <- v.all[!(v.all %in% v1$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v1) <- names(df)
v1 <- rbind(v1, df)
v1$l<-as.numeric(levels(v1$l))[v1$l]
v1<-v1[order(v1$l), ]
v1<-v1$Freq
v1

missing <- v.all[!(v.all %in% v2$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v2) <- names(df)
v2 <- rbind(v2, df)
v2$l<-as.numeric(levels(v2$l))[v2$l]
v2<-v2[order(v2$l), ]
v2<-v2$Freq
v2

missing <- v.all[!(v.all %in% v3$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v3) <- names(df)
v3 <- rbind(v3, df)
v3$l<-as.numeric(levels(v3$l))[v3$l]
v3<-v3[order(v3$l), ]
v3<-v3$Freq
v3

missing <- v.all[!(v.all %in% v4$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v4) <- names(df)
v4 <- rbind(v4, df)
v4$l<-as.numeric(levels(v4$l))[v4$l]
v4<-v4[order(v4$l), ]
v4<-v4$Freq
v4

missing <- v.all[!(v.all %in% v5$Var1)]
df<-data.frame(l=as.factor(missing),Freq=rep(0,length(missing)))
names(v5) <- names(df)
v5 <- rbind(v5, df)
v5$l<-as.numeric(levels(v5$l))[v5$l]
v5<-v5[order(v5$l), ]
v5<-v5$Freq
v5

v <- rbind(v1,v2,v3,v4,v5)
v

# d = number of removed animals from the population at moment i (in our case, only zeroes)
d <- rep(0,dim(eh)[2])
d <- rbind(d,d,d,d,d)
d

# Time covariable
time <- c(1:(dim(eh)[2]-1))

# standardize time
stand_time <- (time - mean(time))/sd(time)

#Calculate Time since last fire (TSLF)

BP<-c(rep(0,5),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,138))
BP
length(BP)
BM<-c(rep(0,7),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,136))
length(BM)

BT<-c(rep(0,8),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,23),1,
      rep(0,35),1,
      rep(0,99))
length(BT)

C<-c(rep(0,236),1,rep(0,99))
length(C)

Q<-c(rep(0,19),1,
     rep(0,23),1,
     rep(0,47),1,
     rep(0,47),1,
     rep(0,47),1,
     rep(0,148))
length(Q)

fire.time <- seq.Date(as.Date("1992/01/01"),as.Date("2019/12/31"), by="month")
fire.time.df <- data.frame(time = fire.time,
                           C = C,
                           Q = Q,
                           BP = BP,
                           BM = BM,
                           BT = BT)
last.fire <- data.frame(time = as.Date("1971/12/31"),
                        C = 1,
                        Q = 1,
                        BP = 1,
                        BM = 1,
                        BT = 1)

fire.time.df <- rbind(last.fire, fire.time.df)

# make an index of the latest events
last_event_index_C <- cumsum(fire.time.df$C) + 1

# shift it by one to the right
last_event_index_C <- c(1, last_event_index_C[1:length(last_event_index_C) - 1])

# get the dates of the events and index the vector with the last_event_index,
# added an NA as the first date because there was no event
TSLF_C <- c(as.Date(NA), fire.time.df[which(fire.time.df$C==1), "time"])[last_event_index_C]

# substract the event's date with the date of the last event
fire.time.df$TSLF_C <- (fire.time.df$time - TSLF_C)/30

# make an index of the latest events
last_event_index_Q <- cumsum(fire.time.df$Q) + 1

# shift it by one to the right
last_event_index_Q <- c(1, last_event_index_Q[1:length(last_event_index_Q) - 1])

# get the dates of the events and index the vector with the last_event_index,
# added an NA as the first date because there was no event
TSLF_Q <- c(as.Date(NA), fire.time.df[which(fire.time.df$Q==1), "time"])[last_event_index_Q]

# substract the event's date with the date of the last event
fire.time.df$TSLF_Q <- (fire.time.df$time - TSLF_Q)/30

# make an index of the latest events
last_event_index_BP <- cumsum(fire.time.df$BP) + 1

# shift it by one to the right
last_event_index_BP <- c(1, last_event_index_BP[1:length(last_event_index_BP) - 1])

# get the dates of the events and index the vector with the last_event_index,
# added an NA as the first date because there was no event
TSLF_BP <- c(as.Date(NA), fire.time.df[which(fire.time.df$BP==1), "time"])[last_event_index_BP]

# substract the event's date with the date of the last event
fire.time.df$TSLF_BP <- (fire.time.df$time - TSLF_BP)/30

# make an index of the latest events
last_event_index_BM <- cumsum(fire.time.df$BM) + 1

# shift it by one to the right
last_event_index_BM <- c(1, last_event_index_BM[1:length(last_event_index_BM) - 1])

# get the dates of the events and index the vector with the last_event_index,
# added an NA as the first date because there was no event
TSLF_BM <- c(as.Date(NA), fire.time.df[which(fire.time.df$BM==1), "time"])[last_event_index_BM]

# substract the event's date with the date of the last event
fire.time.df$TSLF_BM <- (fire.time.df$time - TSLF_BM)/30

# make an index of the latest events
last_event_index_BT <- cumsum(fire.time.df$BT) + 1

# shift it by one to the right
last_event_index_BT <- c(1, last_event_index_BT[1:length(last_event_index_BT) - 1])

# get the dates of the events and index the vector with the last_event_index,
# added an NA as the first date because there was no event
TSLF_BT <- c(as.Date(NA), fire.time.df[which(fire.time.df$BT==1), "time"])[last_event_index_BT]

# substract the event's date with the date of the last event
fire.time.df$TSLF_BT <- (fire.time.df$time - TSLF_BT)/30

fire.time.df <- fire.time.df[fire.time.df$time >= "2005-11-01",]
nrow(fire.time.df)

#Standardize variables and create array
(mean.fire <- mean(c(fire.time.df$C,fire.time.df$Q,fire.time.df$BP,
                     fire.time.df$BM,fire.time.df$BT)))

(sd.fire <- sd(c(fire.time.df$C,fire.time.df$Q,fire.time.df$BP,
                 fire.time.df$BM,fire.time.df$BT)))

(mean.TSLF <- as.numeric(mean(c(fire.time.df$TSLF_C,fire.time.df$TSLF_Q,fire.time.df$TSLF_BP,
                                fire.time.df$TSLF_BM,fire.time.df$TSLF_BT))))
(sd.TSLF <- as.numeric(sd(c(fire.time.df$TSLF_C,fire.time.df$TSLF_Q,fire.time.df$TSLF_BP,
                            fire.time.df$TSLF_BM,fire.time.df$TSLF_BT))))

amb <- array(c(rbind((var_amb$tmed2m[var_amb$canopy=="C"]-mean(var_amb$tmed2m))/sd(var_amb$tmed2m),
                     (var_amb$tmed2m[var_amb$canopy=="Q"]-mean(var_amb$tmed2m))/sd(var_amb$tmed2m),
                     (var_amb$tmed2m[var_amb$canopy=="BP"]-mean(var_amb$tmed2m))/sd(var_amb$tmed2m),
                     (var_amb$tmed2m[var_amb$canopy=="BM"]-mean(var_amb$tmed2m))/sd(var_amb$tmed2m),
                     (var_amb$tmed2m[var_amb$canopy=="BT"]-mean(var_amb$tmed2m))/sd(var_amb$tmed2m)),
               rbind((var_amb$RHmax[var_amb$canopy=="C"]-mean(var_amb$RHmax))/sd(var_amb$RHmax),
                     (var_amb$RHmax[var_amb$canopy=="Q"]-mean(var_amb$RHmax))/sd(var_amb$RHmax),
                     (var_amb$RHmax[var_amb$canopy=="BP"]-mean(var_amb$RHmax))/sd(var_amb$RHmax),
                     (var_amb$RHmax[var_amb$canopy=="BM"]-mean(var_amb$RHmax))/sd(var_amb$RHmax),
                     (var_amb$RHmax[var_amb$canopy=="BT"]-mean(var_amb$RHmax))/sd(var_amb$RHmax)),
               rbind((var_amb$sol[var_amb$canopy=="C"]-mean(var_amb$sol))/sd(var_amb$sol),
                     (var_amb$sol[var_amb$canopy=="Q"]-mean(var_amb$sol))/sd(var_amb$sol),
                     (var_amb$sol[var_amb$canopy=="BP"]-mean(var_amb$sol))/sd(var_amb$sol),
                     (var_amb$sol[var_amb$canopy=="BM"]-mean(var_amb$sol))/sd(var_amb$sol),
                     (var_amb$sol[var_amb$canopy=="BT"]-mean(var_amb$sol))/sd(var_amb$sol)),
               rbind((var_amb$tmed0cm[var_amb$canopy=="C"]-mean(var_amb$tmed0cm))/sd(var_amb$tmed0cm),
                     (var_amb$tmed0cm[var_amb$canopy=="Q"]-mean(var_amb$tmed0cm))/sd(var_amb$tmed0cm),
                     (var_amb$tmed0cm[var_amb$canopy=="BP"]-mean(var_amb$tmed0cm))/sd(var_amb$tmed0cm),
                     (var_amb$tmed0cm[var_amb$canopy=="BM"]-mean(var_amb$tmed0cm))/sd(var_amb$tmed0cm),
                     (var_amb$tmed0cm[var_amb$canopy=="BT"]-mean(var_amb$tmed0cm))/sd(var_amb$tmed0cm)),
               rbind((var_amb$tmin0cm[var_amb$canopy=="C"]-mean(var_amb$tmin0cm))/sd(var_amb$tmin0cm),
                     (var_amb$tmin0cm[var_amb$canopy=="Q"]-mean(var_amb$tmin0cm))/sd(var_amb$tmin0cm),
                     (var_amb$tmin0cm[var_amb$canopy=="BP"]-mean(var_amb$tmin0cm))/sd(var_amb$tmin0cm),
                     (var_amb$tmin0cm[var_amb$canopy=="BM"]-mean(var_amb$tmin0cm))/sd(var_amb$tmin0cm),
                     (var_amb$tmin0cm[var_amb$canopy=="BT"]-mean(var_amb$tmin0cm))/sd(var_amb$tmin0cm)),
               rbind((var_amb$precip[var_amb$canopy=="C"]-mean(var_amb$precip))/sd(var_amb$precip),
                     (var_amb$precip[var_amb$canopy=="Q"]-mean(var_amb$precip))/sd(var_amb$precip),
                     (var_amb$precip[var_amb$canopy=="BP"]-mean(var_amb$precip))/sd(var_amb$precip),
                     (var_amb$precip[var_amb$canopy=="BM"]-mean(var_amb$precip))/sd(var_amb$precip),
                     (var_amb$precip[var_amb$canopy=="BT"]-mean(var_amb$precip))/sd(var_amb$precip)),
               rbind((var_amb$Cnigro_perf[var_amb$canopy=="C"]-mean(var_amb$Cnigro_perf))/sd(var_amb$Cnigro_perf),
                     (var_amb$Cnigro_perf[var_amb$canopy=="Q"]-mean(var_amb$Cnigro_perf))/sd(var_amb$Cnigro_perf),
                     (var_amb$Cnigro_perf[var_amb$canopy=="BP"]-mean(var_amb$Cnigro_perf))/sd(var_amb$Cnigro_perf),
                     (var_amb$Cnigro_perf[var_amb$canopy=="BM"]-mean(var_amb$Cnigro_perf))/sd(var_amb$Cnigro_perf),
                     (var_amb$Cnigro_perf[var_amb$canopy=="BT"]-mean(var_amb$Cnigro_perf))/sd(var_amb$Cnigro_perf)),
               rbind((var_amb$Cnigro_ha90[var_amb$canopy=="C"]-mean(var_amb$Cnigro_ha90))/sd(var_amb$Cnigro_ha90),
                     (var_amb$Cnigro_ha90[var_amb$canopy=="Q"]-mean(var_amb$Cnigro_ha90))/sd(var_amb$Cnigro_ha90),
                     (var_amb$Cnigro_ha90[var_amb$canopy=="BP"]-mean(var_amb$Cnigro_ha90))/sd(var_amb$Cnigro_ha90),
                     (var_amb$Cnigro_ha90[var_amb$canopy=="BM"]-mean(var_amb$Cnigro_ha90))/sd(var_amb$Cnigro_ha90),
                     (var_amb$Cnigro_ha90[var_amb$canopy=="BT"]-mean(var_amb$Cnigro_ha90))/sd(var_amb$Cnigro_ha90)),
               rbind(fire.time.df$C,
                     fire.time.df$Q ,
                     fire.time.df$BP,
                     fire.time.df$BM,
                     fire.time.df$BT),
               rbind((fire.time.df$TSLF_C - mean.TSLF)/sd.TSLF,
                     (fire.time.df$TSLF_Q - mean.TSLF)/sd.TSLF,
                     (fire.time.df$TSLF_BP - mean.TSLF)/sd.TSLF,
                     (fire.time.df$TSLF_BM - mean.TSLF)/sd.TSLF,
                     (fire.time.df$TSLF_BT - mean.TSLF)/sd.TSLF)),
             dim = c(5,170,10))
dim(amb)
str(amb)


# Data list for JAGS
bugs.data <- list(u = u, n = n, v = v, d = d, first = f, nind = dim(eh)[1], n.occasions = dim (eh)[2],
                  y = eh, amb = amb, x = as.matrix(x), xpop = xpop.fire, z = known.states.cjs(eh),
                  mu.L0 = mean(datA$CRC[datA$CRC<=40],na.rm=T),
                  tau.L0 = var(datA$CRC[datA$CRC<=40],na.rm=T),
                  # mu.LI = max(datA$CRC,na.rm=T),
                  AFC = as.numeric(age),
                  #mplot = mplot,
                  sex = sex,
                  plot = plot,
                  neggs = fecundidade.Cn$Nembr,
                  xfec = fecundidade.Cn$CRC,
                  n.fec = length(fecundidade.Cn$CRC),
                  eggs = Cnigropunctatum.RECOR.females$ovos,
                  xprep = as.numeric(Cnigropunctatum.RECOR.females$crc),
                  n.probrep = length(Cnigropunctatum.RECOR.females$ovos))

saveRDS(bugs.data, "Cnigropunctatum.data.rds")

bugs.data <- readRDS("Cnigropunctatum.data.rds")

#Mark-recapture statistics
table(rowSums(bugs.data$y))#Number of captures for each individual
(sum(bugs.data$y)-bugs.data$nind)/bugs.data$nind#Mean recapture rate
table(rowSums(bugs.data$z, na.rm = T))#time period between first and last captures
mean(table(rowSums(bugs.data$z, na.rm = T))[-1])#Mean time period between first and last captures
sd(table(rowSums(bugs.data$z, na.rm = T))[-1])#SD time period between first and last captures


# Initial values
inits <- function (){ list(tauphib = 1, betaTphi = rep(0,10), varphi = rep(0,10),
                           sigma.phiJS = runif(1, 1, 2),
                           sigma.f = runif(1, 0.5, 1),
                           taufb = 1, betaTf = rep(0,10), varf = rep(0,10),
                           #alpha.pJS = runif(1, -0.5, 0.5),
                           taupb = 1,  betaTp = rep(0,10), varp = rep(0,10),
                           sigma.pJS = runif(1, 0.1, 0.5)
)}

# Define parameters to monitour
parameters <- c("phiJS", "alpha.phiJS", "sigma.phiJS",
                "betaphiJS","varphi",
                "f", "alpha.f", "sigma.f",
                "betaf","varf",
                "pJS", "alpha.pJS", "sigma.pJS",
                "betap","varp",
                "rho",
                "beta.phi", "beta2.phi","beta.p", "beta2.p",
                "p.AFC","r.AFC","var.AFC","mn.AFC","mu.K",
                "mu.LI",
                "alpha.fec","beta1.fec",
                "alpha.prep", "beta1.prep","beta2.prep"
)


# Specify model in BUGS language
sink("vitalrates-nigro-crc.jags")
cat("

data{
for(j in 1:5){
C[j]<-10000
zeros[j]<-0
}}

model {

  #########################
  ## SURVIVAL/GROWTH MODEL#
  #########################

  
   for(i in 1:nind){
    AFC[i] ~ dnegbin(p.AFC , r.AFC)T(0,50)  ###Change trucationfor different species. AFC is age of first capture - known in many cases (for turtles (newborns) and estimated when not known
    L0[i] ~ dnorm(mu.L0, tau.L0)  ### draw values for intial size
    LI[i] ~ dnorm(mu.LI,tau.LI)T(0,) ### asymptotic size - taubeta allows for individual variation, while mean size is plot dependent
    newLI[i] ~ dnorm(mu.LI,tau.LI)T(0,)
    LLoldLI[i] <-logdensity.norm(LI[i],mu.LI,tau.LI)
    LLnewLI[i] <-logdensity.norm(newLI[i],mu.LI,tau.LI)
    logit(K[i]) <- K.L[i]  ## mean growth rate is plot dependent with variation defined by xi*theta
    K.L[i] ~ dnorm(mu.K[plot[i]],tau.K)
  }
  
  ### Priors for the growth model
  for(i in 1:5){
    mu.K1[i] ~ dunif(0.5,1)
    mu.K[i] <- log(mu.K1[i]) - log(1-mu.K1[i])
  }
  r.AFC ~ dgamma(0.01,0.01)
  p.AFC <- r.AFC/(r.AFC+mn.AFC)
  mn.AFC ~ dgamma(0.01,0.01)
  var.AFC <- r.AFC*(1-p.AFC)/(p.AFC*p.AFC)
  sd.sample ~ dt(0,0.0004,3)T(0,) ### t priors as in Schofield et al. 2013
  mu.LI ~ dnorm(110,20) #Change for different species
  sd.LI ~ dt(0,0.0004,3)T(0,)
  sd.K ~ dt(0,0.0004,3)T(0,)
  tau.sample <- 1/(sd.sample^2)
  tau.LI <- 1/(sd.LI^2)
  tau.K <- 1/(sd.K^2)    
  
# Priors and constraints
for (i in 1:nind){
  for (t in first[i]:(n.occasions-1)){
    logit(phi[i,t]) <- logit.phiJS[plot[i],t] + beta.phi*x[i,t] + beta2.phi*pow(x[i,t],2)
    logit(p[i,t]) <-  logit.pJS[plot[i],t] + beta.p*x[i,t] + beta2.p*pow(x[i,t],2)
                      # Growth
  } #t
} #i

#### PRIORS
beta.phi ~ dnorm(0, 0.01)           # Prior for slope parameter
beta2.phi ~ dnorm(0, 0.01) 
beta.p ~ dnorm(0, 0.01)           # Prior for slope parameter
beta2.p ~ dnorm(0, 0.01)

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,first[i]] <- 1
   for (t in (first[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      x[i,t-1] ~ dnorm(L[i,t-1],tau.sample)T(0,)
      L[i,t-1] <- (L0[i] + (LI[i]-L0[i])*(1-K[i]^(AFC[i]+(t-1)-first[i]))) 
      } #t
   } #i
   
   ################
   #Number of eggs#
   ################
   # Priors
alpha.fec ~ dunif(-10, 10)
beta1.fec ~ dunif(-1, 1)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n.fec){
   neggs[i] ~ dpois(fecundity[i])          # 1. Distribution for random part
   log(fecundity[i]) <- log.fecundity[i]  # 2. Link function
   log.fecundity[i] <- alpha.fec + beta1.fec * xfec[i]   # 3. Linear predictor
   } #i
   
   ############################
   #Probability of reproducing#
   ############################
   # Priors
alpha.prep ~ dunif(-10, 10)
beta1.prep ~ dunif(-1, 1)
beta2.prep ~ dunif(-1, 1)


# Likelihood: Note key components of a GLM on one line each
for (i in 1:n.probrep){
   eggs[i] ~ dbern(probrep[i])          # 1. Distribution for random part
   logit(probrep[i]) <- alpha.prep + beta1.prep * xprep[i] + beta2.prep*pow(xprep[i],2) # 3. Linear predictor
   } #i

   #################
   #Pradel JS model#
   #################
   
###########PRIORS#######################
for(j in 1:5){
gamma[j, 1]<-0
phiJS[j, n.occasions]<-0
}

for(t in 1:n.occasions){
for(j in 1:5){
muJS[j,t]~dunif(0,1)
}}

#logit constraint for survival probability(phiJS)
for(j in 1:5){
alpha.phiJS[j] ~ dnorm(0.5,0.01)
mean.phiJS[j] <- 1/(1+exp(-alpha.phiJS[j]))#alpha.phiJS on prob scale
for(t in 1:(n.occasions-1)){
phiJS[j, t] <- 1/(1+exp(-logit.phiJS[j, t]))
logit.phiJS[j, t] <- alpha.phiJS[j] + eps.phiJS[j,t] + inprod(amb[j,t,],betaphiJS)
eps.phiJS[j,t] ~ dnorm(0,tau.phiJS)
}
}

for(j in 1:10){
varphi[j]~dbern(0.5)
betaTphi[j]~dnorm(0,tauphib)
betaphiJS[j]<-varphi[j]*betaTphi[j]
}

#environmental parameters
tauphib~dgamma(1,0.001)

#log constraint for recruitment rate(f)

for(j in 1:5){
alpha.f[j] ~ dnorm(-0.5,0.01)
mean.f[j] <- exp(alpha.f[j])#alpha.f on prob scale
for(t in 1:(n.occasions-1)){
f[j,t] <- exp(log.f[j,t])
log.f[j,t]<- alpha.f[j] + eps.f[j,t]+ inprod(amb[j,t,],betaf)
eps.f[j,t] ~ dnorm(0,tau.f)
}
}

for(j in 1:10){
varf[j]~dbern(0.5)
betaTf[j]~dnorm(0,taufb)
betaf[j]<-varf[j]*betaTf[j]
}

#environmental parameters
taufb~dgamma(1,0.001)

for(j in 1:5){
alpha.pJS[j] ~ dnorm(-0.5,0.01)
mean.pJS[j] <- 1/(1+exp(-alpha.pJS[j])) #alpha.pJS on prob scale
for(t in 1:n.occasions){
#logit constraint for detectability (p)
pJS[j,t] <- 1/(1+exp(-logit.pJS[j,t]))
logit.pJS[j,t] <- alpha.pJS[j] + eps.pJS[j,t] + inprod(amb[j,t,],betap)
eps.pJS[j,t] ~ dnorm(0,tau.pJS)
}}

for(j in 1:10){
varp[j]~dbern(0.5)
betaTp[j]~dnorm(0,taupb)
betap[j]<-varp[j]*betaTp[j]
}

#environmental parameters
taupb~dgamma(1,0.001)

#temporal random variation
tau.f<-1/(sigma.f*sigma.f)
sigma.f~dunif(0,2)

tau.phiJS<-1/(sigma.phiJS*sigma.phiJS)
sigma.phiJS~dunif(0,2)

tau.pJS<-1/(sigma.pJS*sigma.pJS)
sigma.pJS~dunif(0,2)

###########LIKELIHOOD(ZERO-TRICK)######
for(j in 1:5){
zeros[j]~dpois(zero.mean[j])
zero.mean[j]<--LJS[j]+C[j]
LJS[j]<-sum(l.num[j, 1:n.occasions])-l.denom[j]

#####log-likelihood for the first occasion
l.num[j,1]<-(u[j,1]*log(xi[j,1]))+(n[j,1]*log(pJS[j,1]))+(secondexpo[j,1]*log(1-pJS[j,1]))+
(thirdexpo[j,1]*log(phiJS[j,1]))+(fourthexpo[j,1]*log(muJS[j,1]))+
(d[j,1]*log(1-muJS[j,1]))+(fifthexpo[j,1]*log(1-(pJS[j,1]*(1-muJS[j,1]))))+
(sixthexpo[j,1]*log(chi[j,1]))
xi[j,1]<-1
secondexpo_a[j,1]<-sum(u[j, 1:1])
secondexpo_b[j,1]<-0
secondexpo[j,1]<-secondexpo_a[j,1]-secondexpo_b[j,1]-n[j,1]
thirdexpo[j,1]<-sum(v[j,2:n.occasions])
fourthexpo[j,1]<-n[j,1]-d[j,1]
fifthexpo[j,1]<-sum(u[j,2:n.occasions])
sixthexpo[j,1]<-v[j,1]-d[j,1]

#####log-likelihood for the last occasion
l.num[j,n.occasions]<-(u[j,n.occasions]*log(xi[j,n.occasions]))+(firstexpo[j,n.occasions]*(log(phiJS[j,n.occasions-1])-log(phiJS[j,n.occasions-1]+f[j,n.occasions-1])))+
(n[j,n.occasions]*log(pJS[j,n.occasions]))+(secondexpo[j,n.occasions]*log(1-pJS[j,n.occasions]))+
(fourthexpo[j,n.occasions]*log(muJS[j,n.occasions]))+(d[j,n.occasions]*log(1-muJS[j,n.occasions]))+
(fifthexpo[j,n.occasions]*log(1-(pJS[j,n.occasions]*(1-muJS[j,n.occasions]))))+
(sixthexpo[j,n.occasions]*log(chi[j,n.occasions]))
chi[j,n.occasions]<-1

firstexpo[j,n.occasions]<-sum(u[j,1:(n.occasions-1)])
secondexpo_a[j,n.occasions]<-sum(u[j,1:n.occasions])
secondexpo_b[j,n.occasions]<-sum(v[j,1:(n.occasions-1)])
secondexpo[j,n.occasions]<-secondexpo_a[j,n.occasions]-secondexpo_b[j,n.occasions]-n[j,n.occasions]
fourthexpo[j,n.occasions]<-n[j,n.occasions]-d[j,n.occasions]
fifthexpo[j,n.occasions]<-0
sixthexpo[j,n.occasions]<-v[j,n.occasions]-d[j,n.occasions]
}

#####likelihood from occasion 2 to n.occasions-1
for(j in 1:5){
for(i in 2:(n.occasions-1)){
l.num[j,i]<-(u[j,i]*log(xi[j,i]))+(firstexpo[j,i]*(log(phiJS[j,i-1])-log(phiJS[j,i-1]+f[j,i-1])))+
(n[j,i]*log(pJS[j,i]))+(secondexpo[j,i]*log(1-pJS[j,i]))+
(thirdexpo[j,i]*log(phiJS[j,i]))+(fourthexpo[j,i]*log(muJS[j,i]))+
(d[j,i]*log(1-muJS[j,i]))+(fifthexpo[j,i]*log(1-(pJS[j,i]*(1-muJS[j,i]))))+
(sixthexpo[j,i]*log(chi[j,i]))

#first exponent
firstexpo[j,i]<-sum(u[j,1:(i-1)])

#second exponent
secondexpo_a[j,i]<-sum(u[j,1:i])
secondexpo_b[j,i]<-sum(v[j,1:(i-1)])
secondexpo[j,i]<-secondexpo_a[j,i]-secondexpo_b[j,i]-n[j,i]

#third exponent
thirdexpo[j,i]<-sum(v[j,(i+1):n.occasions])

#fourth exponent
fourthexpo[j,i]<-n[j,i]-d[j,i]

#fifth exponent
fifthexpo[j,i]<-sum(u[j,(i+1):n.occasions])

#sixth exponent
sixthexpo[j,i]<-v[j,i]-d[j,i]
}
}

#####likelihood denominator
#1st product
PROD1.1[1]<-1
PROD1.2[1]<-1
PROD1.3[1]<-1
PROD1.4[1]<-1
PROD1.5[1]<-1

for(j in 1:(n.occasions-1)){
PROD1_tmp1[1,j]<-0
PROD1_tmp2[1,j]<-0
PROD1_tmp3[1,j]<-0
PROD1_tmp4[1,j]<-0
PROD1_tmp5[1,j]<-0
}

#fill part of PROD1_tmp
for(i in 2:(n.occasions-1)){
for(j in i:(n.occasions-1)){
PROD1_tmp1[i,j]<-0
PROD1_tmp2[i,j]<-0
PROD1_tmp3[i,j]<-0
PROD1_tmp4[i,j]<-0
PROD1_tmp5[i,j]<-0
}
}

for(i in 2:n.occasions){
for(j in 1:(i-1)){
PROD1_tmp1[i,j]<-phiJS[1,j]*(1-(pJS[1,j]*(1-muJS[1,j])))
PROD1_tmp2[i,j]<-phiJS[2,j]*(1-(pJS[2,j]*(1-muJS[2,j])))
PROD1_tmp3[i,j]<-phiJS[3,j]*(1-(pJS[3,j]*(1-muJS[3,j])))
PROD1_tmp4[i,j]<-phiJS[4,j]*(1-(pJS[4,j]*(1-muJS[4,j])))
PROD1_tmp5[i,j]<-phiJS[5,j]*(1-(pJS[5,j]*(1-muJS[5,j])))
}
}


PROD1.1[2]<-PROD1_tmp1[2,1]
PROD1.2[2]<-PROD1_tmp2[2,1]
PROD1.3[2]<-PROD1_tmp3[2,1]
PROD1.4[2]<-PROD1_tmp4[2,1]
PROD1.5[2]<-PROD1_tmp5[2,1]

for(i in 3:n.occasions){
PROD1.1[i]<-prod(PROD1_tmp1[i,1:(i-1)])
PROD1.2[i]<-prod(PROD1_tmp2[i,1:(i-1)])
PROD1.3[i]<-prod(PROD1_tmp3[i,1:(i-1)])
PROD1.4[i]<-prod(PROD1_tmp4[i,1:(i-1)])
PROD1.5[i]<-prod(PROD1_tmp5[i,1:(i-1)])
}

#2nd product
PROD2.1[n.occasions]<-1
PROD2.2[n.occasions]<-1
PROD2.3[n.occasions]<-1
PROD2.4[n.occasions]<-1
PROD2.5[n.occasions]<-1

for(i in 1:(n.occasions-1)){
for(j in (i+1):n.occasions){
PROD2_tmp1[i,j]<-gamma[1,j]
PROD2_tmp2[i,j]<-gamma[2,j]
PROD2_tmp3[i,j]<-gamma[3,j]
PROD2_tmp4[i,j]<-gamma[4,j]
PROD2_tmp5[i,j]<-gamma[5,j]
}
}

#fill part of PROD2_tmp
for(i in 1:(n.occasions-1)){
for(j in 1:i){
PROD2_tmp1[i,j]<-0
PROD2_tmp2[i,j]<-0
PROD2_tmp3[i,j]<-0
PROD2_tmp4[i,j]<-0
PROD2_tmp5[i,j]<-0
}
}

PROD2.1[n.occasions-1]<-PROD2_tmp1[(n.occasions-1),n.occasions]
PROD2.2[n.occasions-1]<-PROD2_tmp2[(n.occasions-1),n.occasions]
PROD2.3[n.occasions-1]<-PROD2_tmp3[(n.occasions-1),n.occasions]
PROD2.4[n.occasions-1]<-PROD2_tmp4[(n.occasions-1),n.occasions]
PROD2.5[n.occasions-1]<-PROD2_tmp5[(n.occasions-1),n.occasions]

for(i in 1:(n.occasions-2)){
PROD2.1[i]<-prod(PROD2_tmp1[i,(i+1):n.occasions])
PROD2.2[i]<-prod(PROD2_tmp2[i,(i+1):n.occasions])
PROD2.3[i]<-prod(PROD2_tmp3[i,(i+1):n.occasions])
PROD2.4[i]<-prod(PROD2_tmp4[i,(i+1):n.occasions])
PROD2.5[i]<-prod(PROD2_tmp5[i,(i+1):n.occasions])

}
for(i in 1:n.occasions){
denom_base_tmp1[i]<-xi[1,i]*PROD1.1[i]*PROD2.1[i]*pJS[1,i]
denom_base_tmp2[i]<-xi[2,i]*PROD1.2[i]*PROD2.2[i]*pJS[2,i]
denom_base_tmp3[i]<-xi[3,i]*PROD1.3[i]*PROD2.3[i]*pJS[3,i]
denom_base_tmp4[i]<-xi[4,i]*PROD1.4[i]*PROD2.4[i]*pJS[4,i]
denom_base_tmp5[i]<-xi[5,i]*PROD1.5[i]*PROD2.5[i]*pJS[5,i]

}

denom_base1 <- sum(denom_base_tmp1[])
denom_base2 <- sum(denom_base_tmp2[])
denom_base3 <- sum(denom_base_tmp3[])
denom_base4 <- sum(denom_base_tmp4[])
denom_base5 <- sum(denom_base_tmp5[])

denom_expo1 <- sum(u[1,1:n.occasions])
denom_expo2 <- sum(u[2,1:n.occasions])
denom_expo3 <- sum(u[3,1:n.occasions])
denom_expo4 <- sum(u[4,1:n.occasions])
denom_expo5 <- sum(u[5,1:n.occasions])

l.denom[1] <- denom_expo1 * log(denom_base1)
l.denom[2] <- denom_expo2 * log(denom_base2)
l.denom[3] <- denom_expo3 * log(denom_base3)
l.denom[4] <- denom_expo4 * log(denom_base4)
l.denom[5] <- denom_expo5 * log(denom_base5)


#################Define xi and chi
for(i in 2:n.occasions){
for(j in 1:5){
xi.tmp[j,i]<-(1-gamma[j,i])+
(gamma[j,i]*((1-pJS[j,i-1])/(1-(pJS[j,i-1]*(1-muJS[j,i-1]))))*xi[j,i-1])
xi[j,i]<-max(xi.tmp[j,i],0.00001)
}
}

for(i in 1:(n.occasions-1)){
for(j in 1:5){
chi[j,i]<-(1-phiJS[j,i])+(phiJS[j,i]*(1-pJS[j,i+1])*chi[j,i+1])
}
}

#################Gamma and rho as derived parameter
for(i in 2:n.occasions){
for(j in 1:5){
rho[j,i]<-phiJS[j,i-1]+f[j,i-1]
}
}

for(i in 2:n.occasions){
for(j in 1:5){
gamma[j,i]<-phiJS[j,i-1]/(phiJS[j,i-1]+f[j,i-1])
}
}

}
",fill = TRUE)
sink()

# MCMC settings
ni <- 10000
nt <- 10
nb <- 50000
nc <- 4
na <- 50000

#Call JAGS from R
bugs.data$y <- as.matrix(bugs.data$y)
bugs.data$z <- as.matrix(bugs.data$z)
runjags.options(jagspath = "/usr/local/bin/jags")
vitalrates.Cnigro <- run.jags(data=bugs.data, inits=inits, monitor=parameters, model="vitalrates-nigro-crc.jags",
                                n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                                method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                       summarise = FALSE,
                         modules = c("glm"))

results.vitalrates.Cnigro <- results.jags(vitalrates.Cnigro)
results.vitalrates.Cnigro <- add.summary(results.vitalrates.Cnigro)

results.vitalrates.Cnigro.df <- summary(results.vitalrates.Cnigro)
View(results.vitalrates.Cnigro.df)

saveRDS(results.vitalrates.Cnigro.df,"results_vitalrates_Cnigro_df.rds")
saveRDS(results.vitalrates.Cnigro, "results_vitalrates_Cnigro.rds")
summary(vitalrates.Cnigro)

results.vitalrates.Cnigro.df <- readRDS("results_vitalrates_Cnigro_df.rds")
write.csv(results.vitalrates.Cnigro.df, "results_vitalrates_Cnigro_df.csv")

S <- ggs(results.vitalrates.Cnigro$mcmc[,3459:3460])
ggs_density(S,family="beta.phi")
ggs_density(S,family="beta2.phi")
ggs_traceplot(S,family="beta.phi")
ggs_traceplot(S,family="beta2.phi")

#Some parameters have not converged. Next, we run each model separately

#PJS model-----------------------------------------------


# Initial values
inits <- function (){ list(tauphib = 1, betaTphi = rep(0,10), varphi = rep(0,10),
                           sigma.phiJS = runif(1, 1, 2),
                           sigma.f = runif(1, 0.5, 1),
                           taufb = 1, betaTf = rep(0,10), varf = rep(0,10),
                           #alpha.pJS = runif(1, -0.5, 0.5),
                           taupb = 1,  betaTp = rep(0,10), varp = rep(0,10),
                           sigma.pJS = runif(1, 0.1, 0.5)
)}

# Define the parameters to monitour
parameters <- c("phiJS", "alpha.phiJS", "sigma.phiJS",
                "betaphiJS","varphi",
                "f", "alpha.f", "sigma.f",
                "betaf","varf",
                "pJS", "alpha.pJS", "sigma.pJS",
                "betap","varp",
                "rho")


# Specify model in BUGS language
sink("pradel-nigro-crc.jags")
cat("

data{
for(j in 1:5){
C[j]<-10000
zeros[j]<-0
}}

model {

 
   #################
   #Pradel JS model#
   #################
   
###########PRIORS#######################
for(j in 1:5){
gamma[j, 1]<-0
phiJS[j, n.occasions]<-0
}

for(t in 1:n.occasions){
for(j in 1:5){
muJS[j,t]~dunif(0,1)
}}

#logit constraint for survival probability(phiJS)
for(j in 1:5){
alpha.phiJS[j] ~ dnorm(0.5,0.01)
mean.phiJS[j] <- 1/(1+exp(-alpha.phiJS[j]))#alpha.phiJS on prob scale
for(t in 1:(n.occasions-1)){
phiJS[j, t] <- 1/(1+exp(-logit.phiJS[j, t]))
logit.phiJS[j, t] <- alpha.phiJS[j] + eps.phiJS[j,t] + inprod(amb[j,t,],betaphiJS)
eps.phiJS[j,t] ~ dnorm(0,tau.phiJS)
}
}

for(j in 1:10){
varphi[j]~dbern(0.5)
betaTphi[j]~dnorm(0,tauphib)
betaphiJS[j]<-varphi[j]*betaTphi[j]
}

#environmental parameters
tauphib~dgamma(1,0.001)

#log constraint for recruitment rate(f)

for(j in 1:5){
alpha.f[j] ~ dnorm(-0.5,0.01)
mean.f[j] <- exp(alpha.f[j])#alpha.f on prob scale
for(t in 1:(n.occasions-1)){
f[j,t] <- exp(log.f[j,t])
log.f[j,t]<- alpha.f[j] + eps.f[j,t]+ inprod(amb[j,t,],betaf)
eps.f[j,t] ~ dnorm(0,tau.f)
}
}

for(j in 1:10){
varf[j]~dbern(0.5)
betaTf[j]~dnorm(0,taufb)
betaf[j]<-varf[j]*betaTf[j]
}

#environmental parameters
taufb~dgamma(1,0.001)

for(j in 1:5){
alpha.pJS[j] ~ dnorm(-0.5,0.01)
mean.pJS[j] <- 1/(1+exp(-alpha.pJS[j])) #alpha.pJS on prob scale
for(t in 1:n.occasions){
#logit constraint for detectability (p)
pJS[j,t] <- 1/(1+exp(-logit.pJS[j,t]))
logit.pJS[j,t] <- alpha.pJS[j] + eps.pJS[j,t] + inprod(amb[j,t,],betap)
eps.pJS[j,t] ~ dnorm(0,tau.pJS)
}}

for(j in 1:10){
varp[j]~dbern(0.5)
betaTp[j]~dnorm(0,taupb)
betap[j]<-varp[j]*betaTp[j]
}

#environmental parameters
taupb~dgamma(1,0.001)

#temporal random variation
tau.f<-1/(sigma.f*sigma.f)
sigma.f~dunif(0,2)

tau.phiJS<-1/(sigma.phiJS*sigma.phiJS)
sigma.phiJS~dunif(0,2)

tau.pJS<-1/(sigma.pJS*sigma.pJS)
sigma.pJS~dunif(0,2)

###########LIKELIHOOD(ZERO-TRICK)######
for(j in 1:5){
zeros[j]~dpois(zero.mean[j])
zero.mean[j]<--LJS[j]+C[j]
LJS[j]<-sum(l.num[j, 1:n.occasions])-l.denom[j]

#####log-likelihood for the first occasion
l.num[j,1]<-(u[j,1]*log(xi[j,1]))+(n[j,1]*log(pJS[j,1]))+(secondexpo[j,1]*log(1-pJS[j,1]))+
(thirdexpo[j,1]*log(phiJS[j,1]))+(fourthexpo[j,1]*log(muJS[j,1]))+
(d[j,1]*log(1-muJS[j,1]))+(fifthexpo[j,1]*log(1-(pJS[j,1]*(1-muJS[j,1]))))+
(sixthexpo[j,1]*log(chi[j,1]))
xi[j,1]<-1
secondexpo_a[j,1]<-sum(u[j, 1:1])
secondexpo_b[j,1]<-0
secondexpo[j,1]<-secondexpo_a[j,1]-secondexpo_b[j,1]-n[j,1]
thirdexpo[j,1]<-sum(v[j,2:n.occasions])
fourthexpo[j,1]<-n[j,1]-d[j,1]
fifthexpo[j,1]<-sum(u[j,2:n.occasions])
sixthexpo[j,1]<-v[j,1]-d[j,1]

#####log-likelihood for the last occasion
l.num[j,n.occasions]<-(u[j,n.occasions]*log(xi[j,n.occasions]))+(firstexpo[j,n.occasions]*(log(phiJS[j,n.occasions-1])-log(phiJS[j,n.occasions-1]+f[j,n.occasions-1])))+
(n[j,n.occasions]*log(pJS[j,n.occasions]))+(secondexpo[j,n.occasions]*log(1-pJS[j,n.occasions]))+
(fourthexpo[j,n.occasions]*log(muJS[j,n.occasions]))+(d[j,n.occasions]*log(1-muJS[j,n.occasions]))+
(fifthexpo[j,n.occasions]*log(1-(pJS[j,n.occasions]*(1-muJS[j,n.occasions]))))+
(sixthexpo[j,n.occasions]*log(chi[j,n.occasions]))
chi[j,n.occasions]<-1

firstexpo[j,n.occasions]<-sum(u[j,1:(n.occasions-1)])
secondexpo_a[j,n.occasions]<-sum(u[j,1:n.occasions])
secondexpo_b[j,n.occasions]<-sum(v[j,1:(n.occasions-1)])
secondexpo[j,n.occasions]<-secondexpo_a[j,n.occasions]-secondexpo_b[j,n.occasions]-n[j,n.occasions]
fourthexpo[j,n.occasions]<-n[j,n.occasions]-d[j,n.occasions]
fifthexpo[j,n.occasions]<-0
sixthexpo[j,n.occasions]<-v[j,n.occasions]-d[j,n.occasions]
}

#####likelihood from occasion 2 to n.occasions-1
for(j in 1:5){
for(i in 2:(n.occasions-1)){
l.num[j,i]<-(u[j,i]*log(xi[j,i]))+(firstexpo[j,i]*(log(phiJS[j,i-1])-log(phiJS[j,i-1]+f[j,i-1])))+
(n[j,i]*log(pJS[j,i]))+(secondexpo[j,i]*log(1-pJS[j,i]))+
(thirdexpo[j,i]*log(phiJS[j,i]))+(fourthexpo[j,i]*log(muJS[j,i]))+
(d[j,i]*log(1-muJS[j,i]))+(fifthexpo[j,i]*log(1-(pJS[j,i]*(1-muJS[j,i]))))+
(sixthexpo[j,i]*log(chi[j,i]))

#first exponent
firstexpo[j,i]<-sum(u[j,1:(i-1)])

#second exponent
secondexpo_a[j,i]<-sum(u[j,1:i])
secondexpo_b[j,i]<-sum(v[j,1:(i-1)])
secondexpo[j,i]<-secondexpo_a[j,i]-secondexpo_b[j,i]-n[j,i]

#third exponent
thirdexpo[j,i]<-sum(v[j,(i+1):n.occasions])

#fourth exponent
fourthexpo[j,i]<-n[j,i]-d[j,i]

#fifth exponent
fifthexpo[j,i]<-sum(u[j,(i+1):n.occasions])

#sixth exponent
sixthexpo[j,i]<-v[j,i]-d[j,i]
}
}

#####likelihood denominator
#1st product
PROD1.1[1]<-1
PROD1.2[1]<-1
PROD1.3[1]<-1
PROD1.4[1]<-1
PROD1.5[1]<-1

for(j in 1:(n.occasions-1)){
PROD1_tmp1[1,j]<-0
PROD1_tmp2[1,j]<-0
PROD1_tmp3[1,j]<-0
PROD1_tmp4[1,j]<-0
PROD1_tmp5[1,j]<-0
}

#fill part of PROD1_tmp
for(i in 2:(n.occasions-1)){
for(j in i:(n.occasions-1)){
PROD1_tmp1[i,j]<-0
PROD1_tmp2[i,j]<-0
PROD1_tmp3[i,j]<-0
PROD1_tmp4[i,j]<-0
PROD1_tmp5[i,j]<-0
}
}

for(i in 2:n.occasions){
for(j in 1:(i-1)){
PROD1_tmp1[i,j]<-phiJS[1,j]*(1-(pJS[1,j]*(1-muJS[1,j])))
PROD1_tmp2[i,j]<-phiJS[2,j]*(1-(pJS[2,j]*(1-muJS[2,j])))
PROD1_tmp3[i,j]<-phiJS[3,j]*(1-(pJS[3,j]*(1-muJS[3,j])))
PROD1_tmp4[i,j]<-phiJS[4,j]*(1-(pJS[4,j]*(1-muJS[4,j])))
PROD1_tmp5[i,j]<-phiJS[5,j]*(1-(pJS[5,j]*(1-muJS[5,j])))
}
}


PROD1.1[2]<-PROD1_tmp1[2,1]
PROD1.2[2]<-PROD1_tmp2[2,1]
PROD1.3[2]<-PROD1_tmp3[2,1]
PROD1.4[2]<-PROD1_tmp4[2,1]
PROD1.5[2]<-PROD1_tmp5[2,1]

for(i in 3:n.occasions){
PROD1.1[i]<-prod(PROD1_tmp1[i,1:(i-1)])
PROD1.2[i]<-prod(PROD1_tmp2[i,1:(i-1)])
PROD1.3[i]<-prod(PROD1_tmp3[i,1:(i-1)])
PROD1.4[i]<-prod(PROD1_tmp4[i,1:(i-1)])
PROD1.5[i]<-prod(PROD1_tmp5[i,1:(i-1)])
}

#2nd product
PROD2.1[n.occasions]<-1
PROD2.2[n.occasions]<-1
PROD2.3[n.occasions]<-1
PROD2.4[n.occasions]<-1
PROD2.5[n.occasions]<-1

for(i in 1:(n.occasions-1)){
for(j in (i+1):n.occasions){
PROD2_tmp1[i,j]<-gamma[1,j]
PROD2_tmp2[i,j]<-gamma[2,j]
PROD2_tmp3[i,j]<-gamma[3,j]
PROD2_tmp4[i,j]<-gamma[4,j]
PROD2_tmp5[i,j]<-gamma[5,j]
}
}

#fill part of PROD2_tmp
for(i in 1:(n.occasions-1)){
for(j in 1:i){
PROD2_tmp1[i,j]<-0
PROD2_tmp2[i,j]<-0
PROD2_tmp3[i,j]<-0
PROD2_tmp4[i,j]<-0
PROD2_tmp5[i,j]<-0
}
}

PROD2.1[n.occasions-1]<-PROD2_tmp1[(n.occasions-1),n.occasions]
PROD2.2[n.occasions-1]<-PROD2_tmp2[(n.occasions-1),n.occasions]
PROD2.3[n.occasions-1]<-PROD2_tmp3[(n.occasions-1),n.occasions]
PROD2.4[n.occasions-1]<-PROD2_tmp4[(n.occasions-1),n.occasions]
PROD2.5[n.occasions-1]<-PROD2_tmp5[(n.occasions-1),n.occasions]

for(i in 1:(n.occasions-2)){
PROD2.1[i]<-prod(PROD2_tmp1[i,(i+1):n.occasions])
PROD2.2[i]<-prod(PROD2_tmp2[i,(i+1):n.occasions])
PROD2.3[i]<-prod(PROD2_tmp3[i,(i+1):n.occasions])
PROD2.4[i]<-prod(PROD2_tmp4[i,(i+1):n.occasions])
PROD2.5[i]<-prod(PROD2_tmp5[i,(i+1):n.occasions])

}
for(i in 1:n.occasions){
denom_base_tmp1[i]<-xi[1,i]*PROD1.1[i]*PROD2.1[i]*pJS[1,i]
denom_base_tmp2[i]<-xi[2,i]*PROD1.2[i]*PROD2.2[i]*pJS[2,i]
denom_base_tmp3[i]<-xi[3,i]*PROD1.3[i]*PROD2.3[i]*pJS[3,i]
denom_base_tmp4[i]<-xi[4,i]*PROD1.4[i]*PROD2.4[i]*pJS[4,i]
denom_base_tmp5[i]<-xi[5,i]*PROD1.5[i]*PROD2.5[i]*pJS[5,i]

}

denom_base1 <- sum(denom_base_tmp1[])
denom_base2 <- sum(denom_base_tmp2[])
denom_base3 <- sum(denom_base_tmp3[])
denom_base4 <- sum(denom_base_tmp4[])
denom_base5 <- sum(denom_base_tmp5[])

denom_expo1 <- sum(u[1,1:n.occasions])
denom_expo2 <- sum(u[2,1:n.occasions])
denom_expo3 <- sum(u[3,1:n.occasions])
denom_expo4 <- sum(u[4,1:n.occasions])
denom_expo5 <- sum(u[5,1:n.occasions])

l.denom[1] <- denom_expo1 * log(denom_base1)
l.denom[2] <- denom_expo2 * log(denom_base2)
l.denom[3] <- denom_expo3 * log(denom_base3)
l.denom[4] <- denom_expo4 * log(denom_base4)
l.denom[5] <- denom_expo5 * log(denom_base5)


#################Define xi and chi
for(i in 2:n.occasions){
for(j in 1:5){
xi.tmp[j,i]<-(1-gamma[j,i])+
(gamma[j,i]*((1-pJS[j,i-1])/(1-(pJS[j,i-1]*(1-muJS[j,i-1]))))*xi[j,i-1])
xi[j,i]<-max(xi.tmp[j,i],0.00001)
}
}

for(i in 1:(n.occasions-1)){
for(j in 1:5){
chi[j,i]<-(1-phiJS[j,i])+(phiJS[j,i]*(1-pJS[j,i+1])*chi[j,i+1])
}
}

#################Gamma and rho as derived parameter
for(i in 2:n.occasions){
for(j in 1:5){
rho[j,i]<-phiJS[j,i-1]+f[j,i-1]
}
}

for(i in 2:n.occasions){
for(j in 1:5){
gamma[j,i]<-phiJS[j,i-1]/(phiJS[j,i-1]+f[j,i-1])
}
}

}
",fill = TRUE)
sink()

# MCMC settings
ni <- 100000
nt <- 1
nb <- 200000
nc <- 4
na <- 50000

#Call JAGS from R
bugs.data$y <- as.matrix(bugs.data$y)
bugs.data$z <- as.matrix(bugs.data$z)
runjags.options(jagspath = "/usr/local/bin/jags")
pradel.Cnigro<- run.jags(data=bugs.data, inits=inits, monitor=parameters, model="pradel-nigro-crc.jags",
                       n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                       method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                       summarise = FALSE,
                       modules = c("glm"))

results.pradel.Cnigro <- results.jags(pradel.Cnigro)

results.pradel.Cnigro.df <- MCMCsummary(results.pradel.Cnigro$mcmc)
View(results.pradel.Cnigro.df)

saveRDS(results.pradel.Cnigro.df,"results_pradel_Cnigro_df.rds")
write.csv(results.pradel.Cnigro.df, "results_pradel_Cnigro_df.csv")
saveRDS(results.pradel.Cnigro, "results_pradel_Cnigro.rds")


#CJS model--------------------------------------

# Initial values
inits <- function (){}

# Define the parameters to monitour
parameters <- c("alpha.phi","beta.phi", "beta2.phi","alpha.p","beta.p", "beta2.p",
                "p.AFC","r.AFC","var.AFC","mn.AFC",
                "mu.K","mu.LI"
)


# Specify model in BUGS language
sink("cjs-nigro-crc.jags")
cat("

model {

  #########################
  ## SURVIVAL/GROWTH MODEL#
  #########################

  
   for(i in 1:nind){
    AFC[i] ~ dnegbin(p.AFC , r.AFC)T(0,50)  ###Change trucationfor different species. AFC is age of first capture - known in many cases (for newborns) and estimated when not known
    L0[i] ~ dnorm(mu.L0, tau.L0)  ### draw values for intial size
    LI[i] ~ dnorm(mu.LI,tau.LI)T(0,) ### asymptotic size - taubeta allows for individual variation, while mean size is plot dependent
    newLI[i] ~ dnorm(mu.LI,tau.LI)T(0,)
    LLoldLI[i] <-logdensity.norm(LI[i],mu.LI,tau.LI)
    LLnewLI[i] <-logdensity.norm(newLI[i],mu.LI,tau.LI)
    logit(K[i]) <- K.L[i]  ## mean growth rate is plot dependent with variation defined by xi*theta
    K.L[i] ~ dnorm(mu.K[plot[i]],tau.K)
  }
  
  ### Priors for the growth model
  for(i in 1:5){
    mu.K1[i] ~ dunif(0.5,1)
    mu.K[i] <- log(mu.K1[i]) - log(1-mu.K1[i])
  }
  r.AFC ~ dgamma(0.01,0.01)
  p.AFC <- r.AFC/(r.AFC+mn.AFC)
  mn.AFC ~ dgamma(0.01,0.01)
  var.AFC <- r.AFC*(1-p.AFC)/(p.AFC*p.AFC)
  sd.sample ~ dt(0,0.0004,3)T(0,) ### t priors as in Schofield et al. 2013
  mu.LI ~ dnorm(110,20) #Change for different species
  sd.LI ~ dt(0,0.0004,3)T(0,)
  sd.K ~ dt(0,0.0004,3)T(0,)
  tau.sample <- 1/(sd.sample^2)
  tau.LI <- 1/(sd.LI^2)
  tau.K <- 1/(sd.K^2)    

  
# Priors and constraints
for (i in 1:nind){
  for (t in first[i]:(n.occasions-1)){
    logit(phi[i,t]) <- alpha.phi + beta.phi*x[i,t] + beta2.phi*pow(x[i,t],2)
    logit(p[i,t]) <-  alpha.p + beta.p*x[i,t] + beta2.p*pow(x[i,t],2)
                      # Growth
  } #t
} #i

#### PRIORS
alpha.phi ~ dnorm(0, 0.01)
beta.phi ~ dnorm(0, 0.01)           # Prior for slope parameter
alpha.p ~ dnorm(0, 0.01)
beta.p ~ dnorm(0, 0.01)           # Prior for slope parameter



# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,first[i]] <- 1
   for (t in (first[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      x[i,t-1] ~ dnorm(L[i,t-1],tau.sample)T(0,)
      L[i,t-1] <- (L0[i] + (LI[i]-L0[i])*(1-K[i]^(AFC[i]+(t-1)-first[i]))) 
      } #t
   } #i
   


}
",fill = TRUE)
sink()

# MCMC settings
ni <- 20000
nt <- 1
nb <- 100000
nc <- 4
na <- 50000

#Call JAGS from R
bugs.data$y <- as.matrix(bugs.data$y)
bugs.data$z <- as.matrix(bugs.data$z)

pradel.Cnigro.df <- read.csv("results_pradel_Cnigro_df.csv")

p.pradel <- pradel.Cnigro.df[grep(pattern = "f", x = pradel.Cnigro.df$X)[1:850],]

phi.pradel <- pradel.Cnigro.df[grep(pattern = "phi", x = pradel.Cnigro.df$X)[1:850],]


p.pradel$plot <- rep(1:5,170)
p.pradel$time <- rep(1:170, each = 5)

phi.pradel$plot <- rep(1:5,170)
phi.pradel$time <- rep(1:170, each = 5)


phi.pradel$plot <- as.factor(phi.pradel$plot)

ggplot(phi.pradel, aes(x = time, y = mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  ylim(c(0.75,1))+
  scale_color_manual(values=turbo(5))

bugs.data$phiJS <- matrix(phi.pradel$mean,nrow = 5, ncol = 170)
bugs.data$pJS <- matrix(p.pradel$mean,nrow = 5, ncol = 170)
runjags.options(jagspath = "/usr/local/bin/jags")

cjs.Cnigro<- run.jags(data=bugs.data, inits=inits, monitor=parameters, model="cjs-nigro-crc.jags",
                       n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                       method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                       summarise = TRUE,
                       modules = c("glm"))

results.cjs.Cnigro <- results.jags("runjagsfiles_4")#runjagsfiles_4
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

cjs.Cnigro.extend <- extend.jags(results.cjs.Cnigro,
                                       adapt = na,thin = nt, sample = ni, 
                                       method = "bgparallel", jags.refresh = 30,
                                       keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro <- results.jags(cjs.Cnigro.extend,  combine = T)#runjagsfiles
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

ni = 100000
cjs.Cnigro.extend <- extend.jags(results.cjs.Cnigro,
                                 adapt = na,thin = nt, sample = ni, 
                                 method = "bgparallel", jags.refresh = 30,
                                 keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro <- results.jags(cjs.Cnigro.extend,  combine = T)#runjagsfiles_5
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

na = 60000
cjs.Cnigro.extend2 <- extend.jags(results.cjs.Cnigro,
                                 adapt = na,thin = nt, sample = ni, 
                                 method = "bgparallel", jags.refresh = 30,
                                 keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro <- results.jags(cjs.Cnigro.extend2,  combine = T)#runjagsfiles_6
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

na = 50000
cjs.Cnigro.extend3 <- extend.jags(results.cjs.Cnigro,
                                  adapt = na,thin = nt, sample = ni, 
                                  method = "bgparallel", jags.refresh = 30,
                                  keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro <- results.jags(cjs.Cnigro.extend3,  combine = T)#runjagsfiles_7
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

na = 50000
ni = 160000
cjs.Cnigro.extend4 <- extend.jags(results.cjs.Cnigro,
                                  adapt = na,thin = nt, sample = ni, 
                                  method = "bgparallel", jags.refresh = 30,
                                  keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro <- results.jags(cjs.Cnigro.extend4,  combine = T)#runjagsfiles_8
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro)
View(results.cjs.Cnigro.df)

saveRDS(results.cjs.Cnigro.df,"results_cjs_Cnigro_df.rds")
saveRDS(results.cjs.Cnigro, "results_cjs_Cnigro.rds")

write.csv(results.cjs.Cnigro.df,"results_cjs_Cnigro_df.csv")

S <- ggs(results.cjs.Cnigro2$mcmc)

quartz(8,12)
ggs_density(S,family="alpha.p")
ggs_density(S,family="beta.p")
ggs_density(S,family="beta2.p")

ggs_traceplot(S,family="beta.p")
ggs_traceplot(S,family="beta2.p")
ggs_traceplot(S,family="alpha.p")

#CJS model per sex -------------------------

setwd("/Volumes/Macintosh HD/Users/heito/Documents/IPMs/Cap2_LizardsDemography_Cerrado")

Cnigropunctatum.RECOR.imp<-readRDS("Cnigropunctatum_RECOR_imp.rds") #read data

head(Cnigropunctatum.RECOR.imp)
tail(Cnigropunctatum.RECOR.imp)
dados.demografia <- Cnigropunctatum.RECOR.imp[Cnigropunctatum.RECOR.imp$ano < 2020, ] #remove data later than 2019
dados.demografia <- Cnigropunctatum.RECOR.imp[Cnigropunctatum.RECOR.imp$campanha > 52, ] #remove data when sex was not registered
tail(dados.demografia)
str(dados.demografia)

#Remove dead animals
dados.demografia$morto[is.na(dados.demografia$morto)] <- "n"
dados.demografia <- dados.demografia[dados.demografia$morto=="n", ]
head(dados.demografia)
table(dados.demografia$campanha)
table(dados.demografia$identidade)

completos <- complete.cases(dados.demografia[, c("identidade", "parcela")]) #remove NAs
nigro.planilha <- droplevels(dados.demografia[completos, ])
head(nigro.planilha)
str(nigro.planilha)
table(nigro.planilha$campanha)
table(nigro.planilha$identidade)


## Prepare input file and run monthly data ("campanha")


# Create ID
IDENT <- paste(nigro.planilha$parcela, nigro.planilha$identidade, nigro.planilha$ciclo, sep="")
head(IDENT)
nigro.planilha <- data.frame(nigro.planilha, IDENT)
rm(IDENT)
str(nigro.planilha)

# Subset variables of interest
table1 <- nigro.planilha[nigro.planilha$recaptura!="(s)", c("IDENT", "campanha", "sexo", "crc","massa", "recaptura", "parcela")]
str(table1)

# Identifies captures without IDs
table2 <- complete.cases(table1[, c("IDENT")])

# Remove captures without ID and SVL
table3 <- table1[table2, ]

# Order the data by ID and time (campanha)
table4 <- table3[order(table3$IDENT, table3$campanha), ]
str(table4)
summary(table4)

# Convert IDENT from factor to character
table5 <- droplevels(table4)
table5$IDENT <- as.character(table5$IDENT)
table5$parcela <- as.character(table5$parcela)
str(table5)


# Calculate recapture frequencies
recap.table <- data.frame(table(table5$IDENT))
names(recap.table) <- c("identidade", "capturas")
recap.table
table(recap.table$capturas)
872+(2*169)+(3*52)+(4*10)+(5*3)+(6*1)+(7*1) #two captures=1, 3 captures=2, 4 captures=3 ...

# Filter data.frame records to use in the analysis

head(table5)

Age<-c(rep(NA,nrow(table5)))
Age

datA<-data.frame(table5$campanha,table5$sexo,table5$IDENT,table5$crc,table5$massa,Age,table5$parcela)
names(datA)<-c("Year","Sex","TrueID","CRC","Massa","Age","Parcela")
datA$Year<-datA$Year+2000
datA$Age[datA$CRC<=40]<-0

head(datA)
tail(datA)
str(datA)

##Creates variables for the growth model -------

###  del is the time period since the first capture of an individual (0 in the first capture)

del<-c()   ### months since the first capture

for(i in 1:nrow(datA)){
  del[i]<-datA$Year[i]-min(datA$Year[datA$TrueID==datA$TrueID[i]])
}

plot<-cast(datA, TrueID~., value="Parcela", fun.aggregate=function(x) tail(x,1))  ###determine the plor for each individual
plot<-as.character(plot[,2])
plot
plot[plot=="BM"]<-4 # Ordering by fire severity
plot[plot=="BP"]<-3 # Ordering by fire severity
plot[plot=="BT"]<-5 # Ordering by fire severity
plot[plot=="C"]<-1 # Ordering by fire severity
plot[plot=="Q"]<-2 # Ordering by fire severity
plot<-as.numeric(plot)
plot

sex <- cast(datA, TrueID~., value="Sex", fun.aggregate=function(x) tail(x,1))  ###determine the sex for each individual - filter the value of last capture
sex <- factor(sex[,2], levels = c("F", "M"))
sex
sex <- as.integer(sex)
sex
sex <- sex - 1

ind = as.numeric(factor(datA$TrueID)) #ID
y = datA$CRC #SVL
n = max(ind)  ### number of individuals
m = nrow(datA)### number of observations

age<- c()  ## age at first capture
for (a in 1:n){ age[a] <- datA$Age[ind==a][1]}

year <- c()
for (a in 1:n){ year[a] <- datA$Year[ind==a][1]}

head(datA)
tail(datA)
str(datA)

## Create capture histories for the survival model-------
known.states.cjs<-function(ch){
  state<-ch
  for (i in 1:dim(ch)[1]){
    n1<-min(which(ch[i,]==1))
    n2<-max(which(ch[i,]==1))
    state[i,n1:n2]<-1
    state[i,n1]<-NA
  }
  state[state==0]<-NA
  return(state)
}

cjs.init.z<-function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2<-max(which(ch[i,]==1))
    ch[i,f[i]:n2]<-NA
  }
  for (i in 1:dim(ch)[1])
  { ch[i,1:f[i]]<-NA
  }
  return(ch)
}

#Capture histories
eh <- cast(datA,TrueID ~ Year, fun.aggregate = function(x) as.numeric(length(x) >0),value="CRC");eh <- eh[,2:ncol(eh)]
eh.all <- seq(min(datA$Year), max(datA$Year)) #preencher todos os anos ignorados
missing <- eh.all[!(eh.all %in% names(eh))]
col=matrix(0,nrow=nrow(eh),ncol=length(missing))
colnames(col) <- missing
eh <- cbind(eh, col)
eh <- eh[,sort(colnames(eh))]
head(eh)
(f <- apply(eh,1,function(x) which(x==1)[1]))#time since first capture
(nind <- nrow(eh)) #Number of individuals
(n.occasions <- ncol(eh)) #Number of occasions
m #Number of observations
n #Number of individuals

#Create data.frame of plots
mplot <- data.frame(C = as.numeric(plot==1),
                    Q = as.numeric(plot==2),
                    BP = as.numeric(plot==3),
                    BM = as.numeric(plot==4),
                    BT = as.numeric(plot==5))

# Create matrix X indicating SVL (crc)
x <- cast(datA,
          TrueID ~ Year,
          fun.aggregate = function(x) mean(x),
          value = "CRC",
          fill = NA)

x <- x[, 2:ncol(x)]
x.all <- seq(min(datA$Year), max(datA$Year)) #fill all the ignored months
missing <- x.all[!(x.all %in% names(x))]
col=matrix(NA,nrow=nrow(x),ncol=length(missing))
colnames(col) <- missing
x <- cbind(x, col)
x <- x[,sort(colnames(x))]
#x[is.na(x)] <- 0
head(x)

bugs.data.sex <- list(first = f, nind = dim(eh)[1], n.occasions = dim (eh)[2],
                  y = eh, x = as.matrix(x), z = known.states.cjs(eh),
                  mu.L0 = mean(datA$CRC[datA$CRC<=40],na.rm=T),
                  tau.L0 = var(datA$CRC[datA$CRC<=40],na.rm=T),
                  # mu.LI = max(datA$CRC,na.rm=T),
                  AFC = as.numeric(age),
                  #mplot = mplot,
                  sex = sex,
                  plot = plot)

# Initial values
inits <- function (){}

# Define the parameters to monitour
parameters <- c("alpha.phi","beta.phi", "beta2.phi","alpha.p","beta.p", "beta2.p",
                "p.AFC","r.AFC","var.AFC","mn.AFC",
                "mu.K","mu.LI"
)


# Specify model in BUGS language
sink("cjs-sex-nigro-crc.jags")
cat("

model {

  #########################
  ## SURVIVAL/GROWTH MODEL#
  #########################

  
   for(i in 1:nind){
    AFC[i] ~ dnegbin(p.AFC[sex[i]+1], r.AFC[sex[i]+1])T(0,50)  ###Change trucationfor different species. AFC is age of first capture - known in many cases (for newborns) and estimated when not known
    sex[i] ~ dbern(psi)   ### allows for sex to be predicted if unknown for an individual
    L0[i] ~ dnorm(mu.L0, tau.L0)  ### draw values for intial size
    LI[i] ~ dnorm(mu.LI[sex[i]+1],tau.LI)T(0,) ### asymptotic size - taubeta allows for individual variation, while mean size is sex dependent
    newLI[i] ~ dnorm(mu.LI[sex[i]+1],tau.LI)T(0,)
    LLoldLI[i] <-logdensity.norm(LI[i],mu.LI[sex[i]+1],tau.LI)
    LLnewLI[i] <-logdensity.norm(newLI[i],mu.LI[sex[i]+1],tau.LI)
    logit(K[i]) <- K.L[i]  ## mean growth rate is sex dependent with variation defined by xi*theta
    K.L[i] ~ dnorm(mu.K[sex[i]+1],tau.K)
  }
  
  ### Priors for the growth model
  for(i in 1:2){
    mu.K1[i] ~ dunif(0.5,1)
    mu.K[i] <- log(mu.K1[i]) - log(1-mu.K1[i])
  r.AFC[i] ~ dgamma(0.01,0.01)
  p.AFC[i] <- r.AFC[i]/(r.AFC[i]+mn.AFC[i])
  mn.AFC[i] ~ dgamma(0.01,0.01)
  var.AFC[i] <- r.AFC[i]*(1-p.AFC[i])/(p.AFC[i]*p.AFC[i])
  mu.LI[i] ~ dnorm(110,20) #Change for different species
  }
  
  sd.sample ~ dt(0,0.0004,3)T(0,) ### t priors as in Schofield et al. 2013
  sd.LI ~ dt(0,0.0004,3)T(0,)
  sd.K ~ dt(0,0.0004,3)T(0,)
  tau.sample <- 1/(sd.sample^2)
  tau.LI <- 1/(sd.LI^2)
  tau.K <- 1/(sd.K^2)    
  psi ~ dbeta(1,1)
  
# Priors and constraints
for (i in 1:nind){
  for (t in first[i]:(n.occasions-1)){
    logit(phi[i,t]) <- alpha.phi[sex[i]+1] + beta.phi[sex[i]+1]*x[i,t] + beta2.phi[sex[i]+1]*pow(x[i,t],2)
    logit(p[i,t]) <-  alpha.p[sex[i]+1] + beta.p[sex[i]+1]*x[i,t] + beta2.p[sex[i]+1]*pow(x[i,t],2)
                      # Growth
  } #t
} #i

#### PRIORS
for(i in 1:2){
alpha.phi[i] ~ dnorm(0, 0.01)
beta.phi[i] ~ dnorm(0, 0.01)           # Prior for slope parameter
beta2.phi[i] ~ dnorm(0, 0.01)           # Prior for slope parameter

alpha.p[i] ~ dnorm(0, 0.01)
beta.p[i] ~ dnorm(0, 0.01)           # Prior for slope parameter
beta2.p[i] ~ dnorm(0, 0.01)           # Prior for slope parameter

}


# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,first[i]] <- 1
   for (t in (first[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      x[i,t-1] ~ dnorm(L[i,t-1],tau.sample)T(0,)
      L[i,t-1] <- (L0[i] + (LI[i]-L0[i])*(1-K[i]^(AFC[i]+(t-1)-first[i]))) 
      } #t
   } #i
   


}
",fill = TRUE)
sink()

# MCMC settings
ni <- 100000
nt <- 1
nb <- 200000
nc <- 4
na <- 50000

#Call JAGS from R
bugs.data.sex$y <- as.matrix(bugs.data.sex$y)
bugs.data.sex$z <- as.matrix(bugs.data.sex$z)

runjags.options(jagspath = "/usr/local/bin/jags")

cjs.Cnigro.sex <- run.jags(data=bugs.data.sex, inits=inits, monitor=parameters, model="cjs-sex-nigro-crc.jags",
                      n.chains = nc, adapt = na,thin = nt, sample = ni, burnin = nb,
                      method = "bgparallel", jags.refresh = 30,keep.jags.files = TRUE,
                      summarise = TRUE,
                      modules = c("glm"))

results.cjs.Cnigro.sex <- results.jags(cjs.Cnigro.sex)#runjagsfiles
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro.sex)
View(results.cjs.Cnigro.df)

ni = 200000

cjs.Cnigro.extend.sex <- extend.jags(results.cjs.Cnigro.sex,
                                 adapt = na,thin = nt, sample = ni, 
                                 method = "bgparallel", jags.refresh = 30,
                                 keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro.sex <- results.jags(cjs.Cnigro.extend.sex,  combine = T)#runjagsfiles_1
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro.sex)
View(results.cjs.Cnigro.df)

cjs.Cnigro.extend.sex2 <- extend.jags(results.cjs.Cnigro.sex,
                                 adapt = na,thin = nt, sample = ni, 
                                 method = "bgparallel", jags.refresh = 30,
                                 keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro.sex <- results.jags(cjs.Cnigro.extend.sex2,  combine = T)#runjagsfiles_2
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro.sex)
View(results.cjs.Cnigro.df)

cjs.Cnigro.extend.sex3 <- extend.jags(results.cjs.Cnigro.sex,
                                  adapt = na,thin = nt, sample = ni, 
                                  method = "bgparallel", jags.refresh = 30,
                                  keep.jags.files = TRUE,jags = "/usr/local/bin/jags")

results.cjs.Cnigro.sex <- results.jags(cjs.Cnigro.extend.sex3,  combine = T)#runjagsfiles_6
# results.cjs.Cnigro <- add.summary(results.cjs.Cnigro)

results.cjs.Cnigro.df <- summary(results.cjs.Cnigro.sex)
View(results.cjs.Cnigro.df)

#Diagnostic plots

S <- ggs(results.cjs.Cnigro.sex$mcmc)

quartz(8,12)
ggs_density(S,family="alpha.phi")
ggs_density(S,family="alpha.p")

ggs_density(S,family="beta.p")
ggs_density(S,family="beta2.p")

ggs_density(S,family="mu.K")

ggs_traceplot(S,family="beta2.phi")
ggs_traceplot(S,family="beta.phi")
ggs_traceplot(S,family="beta.p")
ggs_traceplot(S,family="beta2.p")
ggs_traceplot(S,family="mu.K")

saveRDS(results.cjs.Cnigro.df,"results_cjs_Cnigro_sex_df.rds")
saveRDS(results.cjs.Cnigro.sex, "results_cjs_Cnigro_sex.rds")

write.csv(results.cjs.Cnigro.df,"results_cjs_Cnigro_sex_df.csv")

#Function plots 
results.cjs.Cnigro.sex.df <- read.csv("results_cjs_Cnigro_sex_df.csv")
results.cjs.Cnigro.df <- read.csv("results_cjs_Cnigro_df.csv")

#Function to estimate size from age
age_to_size <- function(x,mu.L0,mu.LI,K) mu.L0 + (mu.LI-mu.L0)*(1-plogis(K)^x)

#Function to estimate age from size
size_to_age <- function(x,mu.L0,mu.LI,K) log(1-((x - mu.L0)/(mu.LI - mu.L0)))/log(inv_logit(K))

#Function to estimate size in t1 from size in t0
sizet0_t1 <- function(x,mu.L0,mu.LI,K) age_to_size(size_to_age(x,mu.L0,mu.LI,K)+1,mu.L0,mu.LI,K)

#Growth parameters
(mu.L0 <- bugs.data.sex$mu.L0)
(mu.LI <- results.cjs.Cnigro.sex.df$Mean[grep(pattern = "mu.LI", 
                                          x = results.cjs.Cnigro.sex.df$X)])

(K <- results.cjs.Cnigro.sex.df$Mean[grep(pattern = "mu.K", 
                                           x = results.cjs.Cnigro.sex.df$X)])

(K.lw <- results.cjs.Cnigro.sex.df$Lower95[grep(pattern = "mu.K", 
                                                 x = results.cjs.Cnigro.sex.df$X)])

(K.up <- results.cjs.Cnigro.sex.df$Upper95[grep(pattern = "mu.K", 
                                                 x = results.cjs.Cnigro.sex.df$X)])

xx <- 0:120 #time in months

#Create data frame
growth.df <- data.frame(sex = as.factor(rep(c("F", "M"), each = length(xx))),
                        mean = c(age_to_size(xx, mu.L0, mu.LI, K[1]),
                                 age_to_size(xx, mu.L0, mu.LI, K[2])),
                        lower = c(age_to_size(xx, mu.L0, mu.LI, K.lw[1]),
                                  age_to_size(xx, mu.L0, mu.LI, K.lw[2])),
                        upper = c(age_to_size(xx, mu.L0, mu.LI, K.up[1]),
                                  age_to_size(xx, mu.L0, mu.LI, K.up[2])),
                        time = rep(xx, 2),
                        class = "sex"
)

#Plot
quartz(width = 10, height = 8)

ggplot(data = growth.df, aes(x = time, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Time (months)", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(2), name = "Sex") +
  scale_fill_manual(values = turbo(2), name = "Sex")

#Results ignoring sex
(mu.LI.mn <- results.cjs.Cnigro.df$Mean[grep(pattern = "mu.LI", 
                                              x = results.cjs.Cnigro.df$X)])

(K.mn <- mean(results.cjs.Cnigro.df$Mean[grep(pattern = "mu.K", 
                                          x = results.cjs.Cnigro.df$X)]))

(K.lw.mn <- mean(results.cjs.Cnigro.df$Lower95[grep(pattern = "mu.K", 
                                                x = results.cjs.Cnigro.df$X)]))

(K.up.mn <- mean(results.cjs.Cnigro.df$Upper95[grep(pattern = "mu.K", 
                                                x = results.cjs.Cnigro.df$X)]))

#Create data frame
growth.nosex.df <- data.frame(sex = as.factor("Both"),
                              mean = c(age_to_size(xx, mu.L0, mu.LI.mn, K.mn[1])),
                              lower = c(age_to_size(xx, mu.L0, mu.LI.mn, K.lw.mn[1])),
                              upper = c(age_to_size(xx, mu.L0, mu.LI.mn, K.up.mn[1])),
                              time = rep(xx),
                              class = "nosex"
)

#Plot
quartz(width = 10, height = 8)

ggplot(data = growth.nosex.df, aes(x = time, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Time (months)", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(2), name = "Sex") +
  scale_fill_manual(values = turbo(2), name = "Sex")

quartz(width = 12, height = 6)

ggplot(data = rbind(growth.df,growth.nosex.df), aes(x = time, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Time (months)", y = "Snout-vent length (mm)") +
  scale_colour_manual(values = turbo(3), name = "Sex") +
  scale_fill_manual(values = turbo(3), name = "Sex") +
  facet_wrap("class")

#Survival parameters
(alpha.phi <- results.cjs.Cnigro.sex.df$Mean[grep(pattern = "alpha.phi", 
                                              x = results.cjs.Cnigro.sex.df$X)])

(alpha.phi.lw <- alpha.phi - results.cjs.Cnigro.sex.df$SD[grep(pattern = "alpha.phi", 
                                                  x = results.cjs.Cnigro.sex.df$X)])

(alpha.phi.up <- alpha.phi + results.cjs.Cnigro.sex.df$SD[grep(pattern = "alpha.phi", 
                                                        x = results.cjs.Cnigro.sex.df$X)])

(beta.phi <- results.cjs.Cnigro.sex.df$Mean[grep(pattern = "beta.phi", 
                                                  x = results.cjs.Cnigro.sex.df$X)])

(beta.phi.lw <- beta.phi - results.cjs.Cnigro.sex.df$SD[grep(pattern = "beta.phi", 
                                                 x = results.cjs.Cnigro.sex.df$X)])

(beta.phi.up <- beta.phi + results.cjs.Cnigro.sex.df$SD[grep(pattern = "beta.phi", 
                                                       x = results.cjs.Cnigro.sex.df$X)])

(beta2.phi <- results.cjs.Cnigro.sex.df$Mean[grep(pattern = "beta2.phi", 
                                                  x = results.cjs.Cnigro.sex.df$X)])

(beta2.phi.lw <- beta2.phi - results.cjs.Cnigro.sex.df$SD[grep(pattern = "beta2.phi", 
                                                  x = results.cjs.Cnigro.sex.df$X)])

(beta2.phi.up <- beta2.phi + results.cjs.Cnigro.sex.df$SD[grep(pattern = "beta2.phi", 
                                                               x = results.cjs.Cnigro.sex.df$X)])
svl.xx <- seq(mu.L0, max(mu.LI), by = 0.1)

#Create data frame
surv.df <- data.frame(sex = as.factor(rep(c("F", "M"), each = length(svl.xx))),
                        mean = c(plogis(alpha.phi[1] + beta.phi[1]*svl.xx + beta2.phi[1]*svl.xx^2),
                                 plogis(alpha.phi[2] + beta.phi[2]*svl.xx + beta2.phi[2]*svl.xx^2)),
                        lower = c(plogis(alpha.phi.up[1] + beta.phi.lw[1]*svl.xx + beta2.phi.lw[1]*svl.xx^2),
                                  plogis(alpha.phi.up[2] + beta.phi.lw[2]*svl.xx + beta2.phi.lw[2]*svl.xx^2)),
                        upper = c(plogis(alpha.phi.lw[1] + beta.phi.up[1]*svl.xx + beta2.phi.up[1]*svl.xx^2),
                                  plogis(alpha.phi.lw[2] + beta.phi.up[2]*svl.xx + beta2.phi.up[2]*svl.xx^2)),
                        svl = rep(svl.xx, 2),
                      class = "sex"
)

#Plot
quartz(width = 10, height = 8)

ggplot(data = surv.df, aes(x = svl, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Snout-vent length (mm)", y = "Survival") +
  scale_colour_manual(values = turbo(2), name = "Sex") +
  scale_fill_manual(values = turbo(2), name = "Sex")

#Results ignoring sex (more data available)
#Survival parameters
(alpha.phi.mn <- results.cjs.Cnigro.df$Mean[grep(pattern = "alpha.phi", 
                                                  x = results.cjs.Cnigro.df$X)])

(alpha.phi.lw.mn <- alpha.phi.mn - results.cjs.Cnigro.df$SD[grep(pattern = "alpha.phi", 
                                                        x = results.cjs.Cnigro.df$X)])

(alpha.phi.up.mn <- alpha.phi.mn + results.cjs.Cnigro.df$SD[grep(pattern = "alpha.phi", 
                                                        x = results.cjs.Cnigro.df$X)])

(beta.phi.mn <- results.cjs.Cnigro.df$Mean[grep(pattern = "beta.phi", 
                                                 x = results.cjs.Cnigro.df$X)])

(beta.phi.lw.mn <- beta.phi.mn - results.cjs.Cnigro.df$SD[grep(pattern = "beta.phi", 
                                                       x = results.cjs.Cnigro.df$X)])

(beta.phi.up.mn <- beta.phi.mn + results.cjs.Cnigro.df$SD[grep(pattern = "beta.phi", 
                                                       x = results.cjs.Cnigro.df$X)])

(beta2.phi.mn <- results.cjs.Cnigro.df$Mean[grep(pattern = "beta2.phi", 
                                                  x = results.cjs.Cnigro.df$X)])

(beta2.phi.lw.mn <- beta2.phi.mn - results.cjs.Cnigro.df$SD[grep(pattern = "beta2.phi", 
                                                        x = results.cjs.Cnigro.df$X)])

(beta2.phi.up.mn <- beta2.phi.mn + results.cjs.Cnigro.df$SD[grep(pattern = "beta2.phi", 
                                                        x = results.cjs.Cnigro.df$X)])
svl.xx <- seq(mu.L0, max(mu.LI), by = 0.1)

#Create data frame
surv.nosex.df <- data.frame(sex = as.factor("Both"),
                      mean = c(plogis(alpha.phi.mn[1] + beta.phi.mn[1]*svl.xx + beta2.phi.mn[1]*svl.xx^2)),
                      lower = c(plogis(alpha.phi.up.mn[1] + beta.phi.lw.mn[1]*svl.xx + beta2.phi.lw.mn[1]*svl.xx^2)),
                      upper = c(plogis(alpha.phi.lw.mn[1] + beta.phi.up.mn[1]*svl.xx + beta2.phi.up.mn[1]*svl.xx^2)),
                      svl = rep(svl.xx, 1),
                      class = "nosex"
)

#Plot
quartz(width = 10, height = 8)

ggplot(data = surv.nosex.df, aes(x = svl, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Snout-vent length (mm)", y = "Survival") +
  scale_colour_manual(values = turbo(1), name = "Sex") +
  scale_fill_manual(values = turbo(1), name = "Sex")

quartz(width = 12, height = 6)

ggplot(data = rbind(surv.df,surv.nosex.df), aes(x = svl, y = mean, colour = sex))+
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
  labs(x = "Snout-vent length (mm)", y = "Survival") +
  scale_colour_manual(values = turbo(3), name = "Sex") +
  scale_fill_manual(values = turbo(3), name = "Sex")+
  facet_wrap("class")




