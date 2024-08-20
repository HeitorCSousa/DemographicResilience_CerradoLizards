rm(list = ls())

#Set working directory
setwd("~/Documents/GitHub/DemographicResilience_CerradoLizards/IPMs")

#Load packages
library(jagsUI)
library(rjags)
library(ipmr)
library(Rage)
library(popdemo)
library(tidyverse)
library(BayesPostEst)
library(MCMCvis)
library(ggplot2)
library(viridis)
library(psych)

#Read results and data to build IPMs
Titambere.data <- readRDS("Titambere_data.rds")

vitalrates.Titambere.samples <- readRDS("sims_list_Titambere.rds")

vitalrates.Titambere.df <- read.csv("results_vitalrates_Titambere_df_100000iters.csv")

pradel.Titambere.df <- read.csv("results_pradel_itambere_df_400000_noecophys.csv")

fec.Titambere.df <- read.csv("results_fecund_itambere_df.csv")

f.vitalrates <- vitalrates.Titambere.df[grep(pattern = "f", x = vitalrates.Titambere.df$X)[1:850],]
f.pradel <- pradel.Titambere.df[grep(pattern = "f", x = pradel.Titambere.df$X)[1:850],]

rho.vitalrates <- vitalrates.Titambere.df[grep(pattern = "rho", x = vitalrates.Titambere.df$X)[1:850],]
rho.pradel <- pradel.Titambere.df[grep(pattern = "rho", x = pradel.Titambere.df$X)[1:850],]

phi.pradel <- pradel.Titambere.df[grep(pattern = "phi", x = pradel.Titambere.df$X)[1:850],]
phi.vitalrates <- vitalrates.Titambere.df[grep(pattern = "phi", x = vitalrates.Titambere.df$X)[1:850],]



f.vitalrates$plot <- rep(1:5,170)
f.pradel$plot <- rep(1:5,170)

f.vitalrates$time <- rep(1:170, each = 5)
f.pradel$time <- rep(1:170, each = 5)

rho.vitalrates$plot <- rep(1:5,170)
rho.pradel$plot <- rep(1:5,170)

rho.vitalrates$time <- rep(1:170, each = 5)
rho.pradel$time <- rep(1:170, each = 5)

phi.pradel$plot <- rep(1:5,170)
phi.pradel$time <- rep(1:170, each = 5)

table(f.vitalrates$X == f.pradel$X)
table(rho.vitalrates$X == rho.pradel$X)

f.diffmean <- f.vitalrates$mean - f.pradel$mean

rho.diffmean <- rho.vitalrates$mean - rho.pradel$mean

phi.diffmean <- phi.vitalrates$mean - phi.pradel$mean

f.vitalrates[f.diffmean > 1,]
rho.vitalrates[rho.diffmean > 1,]


summary(f.diffmean)
boxplot(f.diffmean ~ rep(1:5,170))

summary(rho.diffmean)
boxplot(rho.diffmean ~ rep(1:5,170))

summary(phi.diffmean)
boxplot(phi.diffmean ~ rep(1:5,170))

tapply(rho.pradel$mean, rho.pradel$plot, FUN = function(x) exp(mean(log(x), na.rm = T)))
tapply(rho.vitalrates$mean, rho.vitalrates$plot, FUN = function(x) exp(mean(log(x), na.rm = T)))

f.vitalrates$diff <- f.vitalrates$mean - f.pradel$mean
f.pradel$plot <- as.factor(f.pradel$plot)
f.vitalrates$plot <- as.factor(f.vitalrates$plot)
phi.pradel$plot <- as.factor(phi.pradel$plot)

ggplot(f.vitalrates, aes(x = time, y = mean, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  #geom_path(data=f.pradel[f.vitalrates$plot==3,], aes(x = time, y= mean, colour = plot), linetype = "dashed" )+
  ylim(c(0,5))+
  scale_color_manual(values=turbo(5))

ggplot(f.pradel, aes(x = time, y = mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill=plot), alpha=0.4, colour = NA)+
  ylim(c(0,2))+
  scale_color_manual(values=turbo(5))

ggplot(phi.pradel, aes(x = time, y = mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill=plot), alpha=0.4, colour = NA)+
  ylim(c(0.65,0.95))+
  scale_color_manual(values=turbo(5))

#Use the param option
#One with continuous varying environments and another with pre-defined ones

##################################################################################
#Simple deterministic IPM constructed from discretely varying parameter estimates#
##################################################################################

# Define some fixed parameters

fixed_list <- list(
  s_mu_slope   = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='beta.phi'],    #survival slope
  #s_mu_slope2   = vitalrates.Titambere.df['beta.phi2','mean'],    #survival slope
  
  s_sd_slope   = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='beta.phi'],    #survival slope
  #s_sd_slope2  = vitalrates.Titambere.df['beta.phi2','sd'],    #survival slope
  
  
  #Environmental slopes for survival
  s_mu_tmed2m  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[1]'],
  s_mu_RHmax   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[2]'],
  s_mu_sol     = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[3]'],
  s_mu_tmed0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[4]'],
  s_mu_tmin0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[5]'],
  s_mu_precip  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[6]'],
  s_mu_perf    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[7]'],
  s_mu_ha_90   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[8]'],
  s_mu_fire    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[9]'],
  s_mu_TSLF    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[10]'],
  
  s_sd_tmed2m  = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[1]'],
  s_sd_RHmax   = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[2]'],
  s_sd_sol     = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[3]'],
  s_sd_tmed0cm = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[4]'],
  s_sd_tmin0cm = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[5]'],
  s_sd_precip  = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[6]'],
  s_sd_perf    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[7]'],
  s_sd_ha_90   = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[8]'],
  s_sd_fire    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[9]'],
  s_sd_TSLF    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaphiJS[10]'],
  
  #Environmental slopes for reproduction
  r_f_mu_tmed2m  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[1]'],
  r_f_mu_RHmax   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[2]'],
  r_f_mu_sol     = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[3]'],
  r_f_mu_tmed0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[4]'],
  r_f_mu_tmin0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[5]'],
  r_f_mu_precip  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[6]'],
  r_f_mu_perf    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[7]'],
  r_f_mu_ha_90   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[8]'],
  r_f_mu_fire    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[9]'],
  r_f_mu_TSLF    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[10]'],
  
  r_f_sd_tmed2m  = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[1]'],
  r_f_sd_RHmax   = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[2]'],
  r_f_sd_sol     = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[3]'],
  r_f_sd_tmed0cm = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[4]'],
  r_f_sd_tmin0cm = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[5]'],
  r_f_sd_precip  = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[6]'],
  r_f_sd_perf    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[7]'],
  r_f_sd_ha_90   = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[8]'],
  r_f_sd_fire    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[9]'],
  r_f_sd_TSLF    = pradel.Titambere.df$sd[pradel.Titambere.df$X=='betaf[10]'],
  
  
  #Probability of reproduction
  r_r_mu_int   = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='alpha.prep'],
  r_r_mu_slope = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='beta1.prep'],  
  # r_r_mu_slope2 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='beta2.prep'], 
  
  r_r_sd_int   = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='alpha.prep'],
  r_r_sd_slope = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='beta1.prep'],  
  # r_r_sd_slope2 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='beta2.prep'], 
  
  #Number of eggs/embryos
  r_n_mu_int   = fec.Titambere.df$Mean[fec.Titambere.df$X=='alpha.fec'],
  r_n_mu_slope = fec.Titambere.df$Mean[fec.Titambere.df$X=='beta1.fec'],
  r_n_mu_slope2 =fec.Titambere.df$Mean[fec.Titambere.df$X=='beta2.fec'],

  r_n_sd_int   = fec.Titambere.df$SD[fec.Titambere.df$X=='alpha.fec'],
  r_n_sd_slope = fec.Titambere.df$SD[fec.Titambere.df$X=='beta1.fec'],
  r_n_sd_slope2 =fec.Titambere.df$SD[fec.Titambere.df$X=='beta2.fec'],
  
  
  #Size of newborns
  mu_rd     = Titambere.data$mu.L0,   
  sd_rd     = sqrt(Titambere.data$tau.L0),
  mu_LI = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.LI']
  
)


# Now, simulate some random intercepts for growth (g_), survival (s_),
# and offspring production (r_s_). This part is for the purpose of the example.

# First, we create vector of values that each random component can take.
s_params  <- list(
  s_g_mu_int_1 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[1]'],
  s_g_mu_int_2 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[2]'],
  s_g_mu_int_3 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[3]'],
  s_g_mu_int_4 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[4]'],
  s_g_mu_int_5 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[5]'],
  
  s_g_sd_int_1 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.phiJS[1]'],
  s_g_sd_int_2 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.phiJS[2]'],
  s_g_sd_int_3 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.phiJS[3]'],
  s_g_sd_int_4 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.phiJS[4]'],
  s_g_sd_int_5 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.phiJS[5]']
)

g_params <- list(
  
  g_g_mu_K_1 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[1]'],
  g_g_mu_K_2 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[2]'],
  g_g_mu_K_3 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[3]'],
  g_g_mu_K_4 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[4]'],
  g_g_mu_K_5 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[5]'],
  
  g_g_sd_K_1 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='mu.K[1]'],
  g_g_sd_K_2 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='mu.K[2]'],
  g_g_sd_K_3 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='mu.K[3]'],
  g_g_sd_K_4 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='mu.K[4]'],
  g_g_sd_K_5 = vitalrates.Titambere.df$sd[vitalrates.Titambere.df$X=='mu.K[5]']
)

r_params <- list(
  r_f_mu_int_1 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[1]'],
  r_f_mu_int_2 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[2]'],
  r_f_mu_int_3 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[3]'],
  r_f_mu_int_4 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[4]'],
  r_f_mu_int_5 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[5]'],
  
  r_f_sd_int_1 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.f[1]'],
  r_f_sd_int_2 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.f[2]'],
  r_f_sd_int_3 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.f[3]'],
  r_f_sd_int_4 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.f[4]'],
  r_f_sd_int_5 = pradel.Titambere.df$sd[pradel.Titambere.df$X=='alpha.f[5]'],
  
  r_p_mu_int_1 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.pJS[1]'],
  r_p_mu_int_2 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.pJS[2]'],
  r_p_mu_int_3 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.pJS[3]'],
  r_p_mu_int_4 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.pJS[4]'],
  r_p_mu_int_5 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.pJS[5]']
)

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
# add them all together using c()



inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}


pois_r <- function(x) {
  return(
    exp(x)
  )
}

#Function to estimate size from age
age_to_size <- function(x,mu.L0,mu.LI,K) mu.L0 + (mu.LI-mu.L0)*(1-inv_logit(K)^x)

#Function to estimate age from size
size_to_age <- function(x,mu.L0,mu.LI,K) log(1-((x - mu.L0)/(mu.LI - mu.L0)))/log(inv_logit(K))

#Function to estimate size in t1 from size in t0
sizet0_t1 <- function(x,mu.L0,mu.LI,K) age_to_size(size_to_age(x,mu.L0,mu.LI,K)+1,mu.L0,mu.LI,K)

#Variance in growth#
####################

sd_growth <- function(x,mu.L0,mu_LI,site){
  mean.values <- sizet0_t1(x,mu.L0,
                           mu_LI,
                           unlist(.Titambere.samples[,paste0('mu.K[',site,']')]))
  return(sd(mean.values,na.rm=T))
}


my_funs <- list(inv_logit   = inv_logit,
                pois_r      = pois_r,
                sizet0_t1 = sizet0_t1,
                sd_growth = sd_growth)

surv.pradel <- matrix(phi.pradel$mean,nrow = 5, ncol = 170)
recr.pradel <- matrix(f.pradel$mean,nrow = 5, ncol = 170)
recr.pradel[,170] <- 0.0001
env.states <- array(c(Titambere.data$amb,surv.pradel,recr.pradel), dim = c(5,170,12))

env.states[,,11]

all_params_list <- c(fixed_list, g_params, s_params, r_params)


################################################################################
my_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,s_g_sd_int_site) + 
                          rnorm(1, s_mu_tmed2m , s_sd_tmed2m ) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , s_sd_RHmax  ) * RHmax_site +
                          rnorm(1, s_mu_sol    , s_sd_sol    ) * sol_site +
                          rnorm(1, s_mu_tmed0cm, s_sd_tmed0cm) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, s_sd_tmin0cm) * tmin0cm_site+
                          rnorm(1, s_mu_precip , s_sd_precip ) * precip_site+
                          rnorm(1, s_mu_perf   , s_sd_perf   ) * perf_site+
                          rnorm(1, s_mu_ha_90  , s_sd_ha_90  ) * ha_90_site+
                          rnorm(1, s_mu_fire   , s_sd_fire   ) * fire_site+
                          rnorm(1, s_mu_TSLF   , s_sd_TSLF   ) * TSLF_site+
                          rnorm(1, s_mu_slope, s_sd_slope) * ht_1) ,
    s_sigma_site = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,r_r_sd_int) + 
                          rnorm(1, r_r_mu_slope, r_r_sd_slope) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, r_n_sd_int) +
                               rnorm(1,r_n_mu_slope,r_n_sd_slope) * ht_1 +
                               rnorm(1,r_n_mu_slope2,r_n_sd_slope2)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, r_f_sd_int_site) +
                               rnorm(1, r_f_mu_tmed2m , r_f_sd_tmed2m ) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  , r_f_sd_RHmax  ) * RHmax_site +
                               rnorm(1, r_f_mu_sol    , r_f_sd_sol    ) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm, r_f_sd_tmed0cm) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm, r_f_sd_tmin0cm) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip , r_f_sd_precip ) * precip_site+
                               rnorm(1, r_f_mu_perf   , r_f_sd_perf   ) * perf_site+
                               rnorm(1, r_f_mu_ha_90  , r_f_sd_ha_90  ) * ha_90_site+
                               rnorm(1, r_f_mu_fire   , r_f_sd_fire   ) * fire_site+
                               rnorm(1, r_f_mu_TSLF   , r_f_sd_TSLF   ) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site + r_f_sigma_site),
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"

sample_env <- function(env_states, site, iteration) {
  
  out <- as.list(env_states[site, iteration, ])
  names(out) <- c(paste0("tmed2m_",site), 
                  paste0("RHmax_",site),  
                  paste0("sol_",site),    
                  paste0("tmed0cm_",site),
                  paste0("tmin0cm_",site),
                  paste0("precip_",site) ,
                  paste0("perf_",site)   ,
                  paste0("ha_90_",site)  ,
                  paste0("fire_" ,site)  ,
                  paste0("TSLF_",site)   ,
                  paste0("surv_",site)   ,
                  paste0("f_",site))
                  
  
  return(out)
  
}


my_ipm <- my_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 1:5
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(1:5, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 1:5))
#Lambdas
lambda(my_ipm,type_lambda = 'all')
lambda(my_ipm,log = F)

#Kernels
library(fields)
quartz(8,8)
mean.kernel<-mean_kernel(my_ipm)

quartz(height=6,width=12)
par(mfrow=c(1,2))
plot(mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))


plot(mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))
par(mfrow=c(1,1))

######################
#Without sd estimates#
######################


my_ipm2 <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"

sample_env <- function(env_states, site, iteration) {
  
  out <- as.list(env_states[site, iteration, ])
  names(out) <- c(paste0("tmed2m_",site), 
                  paste0("RHmax_",site),  
                  paste0("sol_",site),    
                  paste0("tmed0cm_",site),
                  paste0("tmin0cm_",site),
                  paste0("precip_",site) ,
                  paste0("perf_",site)   ,
                  paste0("ha_90_",site)  ,
                  paste0("fire_" ,site)  ,
                  paste0("TSLF_",site)   ,
                  paste0("surv_",site)   ,
                  paste0("f_",site))
  
  
  return(out)
  
}


my_ipm2 <- my_ipm2 %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 1:5
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(1:5, each=170),
           iterations = 170,
           return_sub_kernels = TRUE)

lambda(my_ipm2,type_lambda = 'all')
lambda(my_ipm2,log = F)


mean.kernel<-mean_kernel(my_ipm2)

quartz(height = 6, width = 12)
par(mfrow=c(1,2))
plot(mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))


###############
#For each plot#
###############
#########
#Control#
#########

C_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 


# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"


C_ipm <- C_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 1
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(1, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 1),)

#Lambdas
lambda(C_ipm,type_lambda = 'all')
lambda(C_ipm,log = F)

#Kernels
C.mean.kernel<-mean_kernel(C_ipm)

quartz(8,8)
plot(C.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

quartz(8,8)
plot(C.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))



plot_ipm_sd_rv <- function(ipm_sd, xmin, xmax, ylab){
  plot(seq(xmin, xmax,length.out=100),ipm_sd[,1],type="l",ylim=c(0,max(ipm_sd[,-ncol(ipm_sd)])),
       xlab = "SVL (mm)", ylab = ylab, bty="n", col=rgb(0,0,0,.1))
for(i in 2:ncol(ipm_sd)-1){
  lines(seq(xmin, xmax,length.out=100),ipm_sd[,i],col=rgb(0,0,0,.1))
}
}

C_ipm_sd <- right_ev(C_ipm)
plot_ipm_sd_rv(C_ipm_sd$ht_w, 20, all_params_list$mu_LI, "Stable distribution")

C_ipm_rv <- left_ev(C_ipm, iterations = 170)
plot_ipm_sd_rv(C_ipm_rv$ht_v, 20, all_params_list$mu_LI, "Reproductive value")

par(mfrow=c(1,1))


#Lambdas IPM x PJS
quartz(10,10)
plot(lambda(C_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.9,1.8), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.pradel$mean[rho.pradel$plot==1][-c(1, 170)])

quartz(8,8)
plot(rho.pradel$mean[rho.pradel$plot==1][-c(1, 170)], lambda(C_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.pradel$mean[rho.pradel$plot==1][-c(1, 170)], lambda(C_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(C_ipm,type_lambda = 'all')[-c(1,170)], rho.pradel$mean[rho.pradel$plot==1][-c(1, 170)])

summary(rho.pradel$mean[rho.pradel$plot==1][-c(1, 170)] - lambda(C_ipm,type_lambda = 'all')[-c(1, 170)])

#############
#Quadrennial#
##############

Q_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 2),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 2),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 


# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"


Q_ipm <- Q_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 2
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(2, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 2),)

#Lambdas
lambda(Q_ipm,type_lambda = 'all')
lambda(Q_ipm,log = F)

#Kernels
Q.mean.kernel<-mean_kernel(Q_ipm)
par(mfrow = c(1,2))
quartz(8,8)
plot(Q.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

quartz(8,8)
plot(Q.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

#Stable size distributions and reproductive values
Q_ipm_sd <- right_ev(Q_ipm)
plot_ipm_sd_rv(Q_ipm_sd$ht_w, 
               Q_ipm$proto_ipm$domain[[1]][[1]][1], 
               Q_ipm$proto_ipm$domain[[1]][[1]][2], "Stable distribution")

Q_ipm_rv <- left_ev(Q_ipm, iterations = 170)
plot_ipm_sd_rv(Q_ipm_rv$ht_v, 
               Q_ipm$proto_ipm$domain[[1]][[1]][1], 
               Q_ipm$proto_ipm$domain[[1]][[1]][2], "Reproductive value")

par(mfrow = c(1,1))

#Lambdas IPM x PJS
quartz(h = 6, w = 8)
plot(lambda(Q_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.9,1.5), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.pradel$mean[rho.pradel$plot==2][-c(1, 170)])

quartz(h = 8, w = 8)
plot(rho.pradel$mean[rho.pradel$plot==2][-c(1, 170)], lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.pradel$mean[rho.pradel$plot==2][-c(1, 170)], lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(Q_ipm,type_lambda = 'all')[-c(1,170)], rho.pradel$mean[rho.pradel$plot==2][-c(1, 170)])

summary(rho.pradel$mean[rho.pradel$plot==2][-c(1, 170)] - lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)])

################
#Early biennial#
################
EB_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 3),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 3),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"


EB_ipm <- EB_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 3
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(3, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 3),)

#Lambdas
lambda(EB_ipm,type_lambda = 'all')
lambda(EB_ipm,log = F)

#Kernels
EB.mean.kernel<-mean_kernel(EB_ipm)

quartz(8,8)
plot(EB.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

quartz(8,8)
plot(EB.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

#Stable size distributions and reproductive values
EB_ipm_sd <- right_ev(EB_ipm)
plot_ipm_sd_rv(EB_ipm_sd$ht_w, 
               EB_ipm$proto_ipm$domain[[1]][[1]][1], 
               EB_ipm$proto_ipm$domain[[1]][[1]][2], "Stable distribution")

EB_ipm_rv <- left_ev(EB_ipm, iterations = 170)
plot_ipm_sd_rv(EB_ipm_rv$ht_v, 
               EB_ipm$proto_ipm$domain[[1]][[1]][1], 
               EB_ipm$proto_ipm$domain[[1]][[1]][2], "Reproductive value")

#Lambdas IPM x PJS
quartz(h = 6, w = 8)
plot(lambda(EB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.8,1.8), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.pradel$mean[rho.pradel$plot==3][-c(1, 170)])

quartz(h = 8, w = 8)
plot(rho.pradel$mean[rho.pradel$plot==3][-c(1, 170)], lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.pradel$mean[rho.pradel$plot==3][-c(1, 170)], lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(EB_ipm,type_lambda = 'all')[-c(1,170)], rho.pradel$mean[rho.pradel$plot==3][-c(1, 170)])

summary(rho.pradel$mean[rho.pradel$plot==3][-c(1, 170)] - lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)])

##############
#Mid biennial#
##############
MB_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 4),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 4),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"


MB_ipm <- MB_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 4
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(4, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 4),)

#Lambdas
lambda(MB_ipm,type_lambda = 'all')
lambda(MB_ipm,log = F)

#Kernels
MB.mean.kernel<-mean_kernel(MB_ipm)

par(mfrow = c(1,2))
quartz(8,8)
plot(MB.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

quartz(8,8)
plot(MB.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

MB_ipm_sd <- right_ev(MB_ipm)
plot_ipm_sd_rv(MB_ipm_sd$ht_w, 
               MB_ipm$proto_ipm$domain[[1]][[1]][1], 
               MB_ipm$proto_ipm$domain[[1]][[1]][2], "Stable distribution")

MB_ipm_rv <- left_ev(MB_ipm, iterations = 170)
plot_ipm_sd_rv(MB_ipm_rv$ht_v, 
               MB_ipm$proto_ipm$domain[[1]][[1]][1], 
               MB_ipm$proto_ipm$domain[[1]][[1]][2], "Reproductive value")

par(mfrow=c(1,1))

#Lambdas IPM x PJS
quartz(h = 6, w = 8)
plot(lambda(MB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.8,2.6), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.pradel$mean[rho.pradel$plot==4][-c(1, 170)])

quartz(h = 8, w = 8)
plot(rho.pradel$mean[rho.pradel$plot==4][-c(1, 170)], lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.pradel$mean[rho.pradel$plot==4][-c(1, 170)], lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(MB_ipm,type_lambda = 'all')[-c(1,170)], rho.pradel$mean[rho.pradel$plot==4][-c(1, 170)])

summary(rho.pradel$mean[rho.pradel$plot==4][-c(1, 170)] - lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)])

###############
#Late biennial#
###############
LB_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
  define_kernel(
    
    # Our P kernels will vary from site to site, so we index it with "_site"
    
    name             = 'P_site',
    
    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed
    
    formula          = s_site * g_site,
    family           = "CC",
    
    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.
    
    s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                          rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                          rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                          rnorm(1, s_mu_sol    , 0) * sol_site +
                          rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                          rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                          rnorm(1, s_mu_precip , 0) * precip_site+
                          rnorm(1, s_mu_perf   , 0) * perf_site+
                          rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                          rnorm(1, s_mu_fire   , 0) * fire_site+
                          rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                          rnorm(1, s_mu_slope, 0) * ht_1),
    s_sigma_site     = surv_site - inv_logit(s_lin_site), 
    s_site           =  inv_logit(s_lin_site) + s_sigma_site,
    
    # Again, we modify the vital rate expression to include "_site".
    
    g_site           = dnorm(ht_2, 
                             mean = sizet0_t1(ht_1, 
                                              mu_rd,
                                              mu_LI,
                                              g_g_mu_K_site), 
                             sd = sd_growth(ht_1,
                                            mu_rd,
                                            mu_LI,
                                            site)),
    
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 5),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                               r_f_lin_site +
                               rnorm(1,r_n_mu_slope,0) * ht_1 +
                               rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    r_n_site         = pois_r(r_n_lin_site),
    r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                               rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                               rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                               rnorm(1, r_f_mu_sol    ,0) * sol_site +
                               rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                               rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                               rnorm(1, r_f_mu_precip ,0) * precip_site+
                               rnorm(1, r_f_mu_perf   ,0) * perf_site+
                               rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                               rnorm(1, r_f_mu_fire   ,0) * fire_site+
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
    r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
    r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    
    # As in the P kernel, we specify the values the index can have.
    
    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      
      # The impl_args are also modified with the index
      
      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(20, all_params_list$mu_LI, 100)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"


LB_ipm <- LB_ipm %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 5
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(100))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(5, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 5),)

#Lambdas
lambda(LB_ipm,type_lambda = 'all')
lambda(LB_ipm,log = F)

#Kernels
LB.mean.kernel<-mean_kernel(LB_ipm)

par(mfrow=c(1,2))
plot(LB.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(LB.mean.kernel$mean_F_site, do_contour=T, col=turbo(1000))

#Stable size distributions and reproductive values
LB_ipm_sd <- right_ev(LB_ipm)
plot_ipm_sd_rv(LB_ipm_sd$ht_w, 
               LB_ipm$proto_ipm$domain[[1]][[1]][1], 
               LB_ipm$proto_ipm$domain[[1]][[1]][2], "Stable distribution")

LB_ipm_rv <- left_ev(LB_ipm, iterations = 170)
plot_ipm_sd_rv(LB_ipm_rv$ht_v, 
               LB_ipm$proto_ipm$domain[[1]][[1]][1], 
               LB_ipm$proto_ipm$domain[[1]][[1]][2], "Reproductive value")

par(mfrow=c(1,1))
#Lambdas IPM x PJS
quartz(h = 6, w = 8)
plot(lambda(LB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.8,2.6), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.pradel$mean[rho.pradel$plot==5][-c(1, 170)])

quartz(h = 8, w = 8)
plot(rho.pradel$mean[rho.pradel$plot==5][-c(1, 170)], lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.pradel$mean[rho.pradel$plot==5][-c(1, 170)], lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(LB_ipm,type_lambda = 'all')[-c(1,170)], rho.pradel$mean[rho.pradel$plot==5][-c(1, 170)], method = "spearman")

summary(rho.pradel$mean[rho.pradel$plot==5][-c(1, 170)] - lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)])

#######################
#Perturbation analyses#
#######################
mean.kernel$mean_K_site <-mean.kernel$mean_P_site + mean.kernel$mean_F_site
C.mean.kernel$mean_K_site <-C.mean.kernel$mean_P_site + C.mean.kernel$mean_F_site
Q.mean.kernel$mean_K_site <-Q.mean.kernel$mean_P_site + Q.mean.kernel$mean_F_site
EB.mean.kernel$mean_K_site <-EB.mean.kernel$mean_P_site + EB.mean.kernel$mean_F_site
MB.mean.kernel$mean_K_site <-MB.mean.kernel$mean_P_site + MB.mean.kernel$mean_F_site
LB.mean.kernel$mean_K_site <-LB.mean.kernel$mean_P_site + LB.mean.kernel$mean_F_site

#Sensitivity
sens.mean.K<-perturb_matrix(mean.kernel$mean_K_site)
plot(sens.mean.K)

sens.mean.K.C<-perturb_matrix(C.mean.kernel$mean_K_site)
plot(sens.mean.K.C)

sens.mean.K.Q<-perturb_matrix(Q.mean.kernel$mean_K_site)
plot(sens.mean.K.Q)

sens.mean.K.EB<-perturb_matrix(EB.mean.kernel$mean_K_site)
plot(sens.mean.K.EB)

sens.mean.K.MB<-perturb_matrix(MB.mean.kernel$mean_K_site)
plot(sens.mean.K.MB)

sens.mean.K.LB<-perturb_matrix(LB.mean.kernel$mean_K_site)
plot(sens.mean.K.LB)

#Elasticity
elas.mean.K<-perturb_matrix(mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K)

elas.mean.K.C<-perturb_matrix(C.mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K.C)

elas.mean.K.Q<-perturb_matrix(Q.mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K.Q)

elas.mean.K.EB<-perturb_matrix(EB.mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K.EB)

elas.mean.K.MB<-perturb_matrix(MB.mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K.MB)

elas.mean.K.LB<-perturb_matrix(LB.mean.kernel$mean_K_site,type = "elasticity")
plot(elas.mean.K.LB)


#Stochastic Elasticity

stoch_K <- function(ipm,sites,time){
  stoch.K <- array(0,dim = c(100,100,sites,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(100,100,sites*time*2))
  #stoch.K.list <- list()
  
  for(j in 1:time){
    
    for(i in 1:sites){
      stoch.K[,,i,j] <- ipm_array[,,i+seq.t[j]] + ipm_array[,,(i)+(sites+seq.t[j])]
      
    }
  }
  list.K <- apply(stoch.K,MARGIN = c(1,2),FUN=c)
  list.K <- aperm(list.K,c(2,3,1))
  list.K2 <- list()
  for(i in 1:(sites*time)){
    list.K2[[i]] <- list.K[,,i] 
  }
  return(list.K2)
}

stoch.K <- stoch_K(my_ipm2, sites = 5,time = 170)
stoch.K.C <- stoch_K(C_ipm, sites = 1,time = 170)
stoch.K.Q <- stoch_K(Q_ipm, sites = 1,time = 170)
stoch.K.EB <- stoch_K(EB_ipm, sites = 1,time = 170)
stoch.K.MB <- stoch_K(MB_ipm, sites = 1,time = 170)
stoch.K.LB <- stoch_K(LB_ipm, sites = 1,time = 170)


stoch.elas.mean.K<-perturb_stochastic(stoch.K,pop_vectors(stoch.K))
stoch.elas.mean.K
imagePlot(t(stoch.elas.mean.K$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K$E_sigma),ylim=c(1,0),col=turbo(100))


stoch.elas.mean.K.C<-perturb_stochastic(stoch.K.C,pop_vectors(stoch.K.C))
stoch.elas.mean.K.C
imagePlot(t(stoch.elas.mean.K.C$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.C$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.C$E_sigma),ylim=c(1,0),col=turbo(100))

stoch.elas.mean.K.Q<-perturb_stochastic(stoch.K.Q,pop_vectors(stoch.K.Q))
stoch.elas.mean.K.Q
imagePlot(t(stoch.elas.mean.K.Q$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.Q$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.Q$E_sigma),ylim=c(1,0),col=turbo(100))

stoch.elas.mean.K.EB<-perturb_stochastic(stoch.K.EB,pop_vectors(stoch.K.EB))
stoch.elas.mean.K.EB
imagePlot(t(stoch.elas.mean.K.EB$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.EB$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.EB$E_sigma),ylim=c(1,0),col=turbo(100))

stoch.elas.mean.K.MB<-perturb_stochastic(stoch.K.MB,pop_vectors(stoch.K.MB))
stoch.elas.mean.K.MB
imagePlot(t(stoch.elas.mean.K.MB$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.MB$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.MB$E_sigma),ylim=c(1,0),col=turbo(100))

stoch.elas.mean.K.LB<-perturb_stochastic(stoch.K.LB,pop_vectors(stoch.K.LB))
stoch.elas.mean.K.LB
imagePlot(t(stoch.elas.mean.K.LB$E),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.LB$E_mu),ylim=c(1,0),col=turbo(100))
imagePlot(t(stoch.elas.mean.K.LB$E_sigma),ylim=c(1,0),col=turbo(100))



######################################
#Life-history and Resilience measures#
######################################

#Survival and lifespan traits
life_expect_mean(matU = mean.kernel$mean_P_site, start = 1)  # mean life expectancy
life_expect_mean(matU = C.mean.kernel$mean_P_site, start = 1)  # mean life expectancy 
life_expect_mean(matU = Q.mean.kernel$mean_P_site, start = 1)  # mean life expectancy 
life_expect_mean(matU = EB.mean.kernel$mean_P_site, start = 1)  # mean life expectancy 
life_expect_mean(matU = MB.mean.kernel$mean_P_site, start = 1)  # mean life expectancy 
life_expect_mean(matU = LB.mean.kernel$mean_P_site, start = 1)  # mean life expectancy

longevity(matU = mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)   
longevity(matU = C.mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)   
longevity(matU = Q.mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)   
longevity(matU = EB.mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)   
longevity(matU = MB.mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)   
longevity(matU = LB.mean.kernel$mean_P_site, start = 1, lx_crit = 0.01)

gen_time(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site)
gen_time(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site)
gen_time(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site)
gen_time(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site)
gen_time(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site)
gen_time(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site)


#Reproduction and maturation traits
net_repro_rate(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site)   
net_repro_rate(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site)   
net_repro_rate(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site)   
net_repro_rate(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site)   
net_repro_rate(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site)   
net_repro_rate(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site) 

net_repro_rate(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site, method = "start")   
net_repro_rate(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site, method = "start")   
net_repro_rate(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site, method = "start")   
net_repro_rate(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site, method = "start")   
net_repro_rate(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site, method = "start")   
net_repro_rate(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site, method = "start") 

mature_age(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site, start = 1)
mature_age(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site, start = 1)
mature_age(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site, start = 1)
mature_age(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site, start = 1)
mature_age(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site, start = 1)
mature_age(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site, start = 1)

mature_prob(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site, start = 1)
mature_prob(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site, start = 1)
mature_prob(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site, start = 1)
mature_prob(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site, start = 1)
mature_prob(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site, start = 1)
mature_prob(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site, start = 1)

#Life-table component traits
lx <- mpm_to_lx(matU = mean.kernel$mean_P_site, start = 1)
px <- mpm_to_px(matU = mean.kernel$mean_P_site, start = 1)
hx <- mpm_to_hx(matU = mean.kernel$mean_P_site, start = 1)
mx <- mpm_to_mx(matU = mean.kernel$mean_P_site, matR = mean.kernel$mean_F_site, start = 1)

lx.C <- mpm_to_lx(matU = C.mean.kernel$mean_P_site, start = 1)
px.C <- mpm_to_px(matU = C.mean.kernel$mean_P_site, start = 1)
hx.C <- mpm_to_hx(matU = C.mean.kernel$mean_P_site, start = 1)
mx.C <- mpm_to_mx(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site, start = 1)

lx.C <- mpm_to_lx(matU = C.mean.kernel$mean_P_site, start = 1)
px.C <- mpm_to_px(matU = C.mean.kernel$mean_P_site, start = 1)
hx.C <- mpm_to_hx(matU = C.mean.kernel$mean_P_site, start = 1)
mx.C <- mpm_to_mx(matU = C.mean.kernel$mean_P_site, matR = C.mean.kernel$mean_F_site, start = 1)

lx.Q <- mpm_to_lx(matU = Q.mean.kernel$mean_P_site, start = 1)
px.Q <- mpm_to_px(matU = Q.mean.kernel$mean_P_site, start = 1)
hx.Q <- mpm_to_hx(matU = Q.mean.kernel$mean_P_site, start = 1)
mx.Q <- mpm_to_mx(matU = Q.mean.kernel$mean_P_site, matR = Q.mean.kernel$mean_F_site, start = 1)

lx.EB <- mpm_to_lx(matU = EB.mean.kernel$mean_P_site, start = 1)
px.EB <- mpm_to_px(matU = EB.mean.kernel$mean_P_site, start = 1)
hx.EB <- mpm_to_hx(matU = EB.mean.kernel$mean_P_site, start = 1)
mx.EB <- mpm_to_mx(matU = EB.mean.kernel$mean_P_site, matR = EB.mean.kernel$mean_F_site, start = 1)

lx.MB <- mpm_to_lx(matU = MB.mean.kernel$mean_P_site, start = 1)
px.MB <- mpm_to_px(matU = MB.mean.kernel$mean_P_site, start = 1)
hx.MB <- mpm_to_hx(matU = MB.mean.kernel$mean_P_site, start = 1)
mx.MB <- mpm_to_mx(matU = MB.mean.kernel$mean_P_site, matR = MB.mean.kernel$mean_F_site, start = 1)

lx.LB <- mpm_to_lx(matU = LB.mean.kernel$mean_P_site, start = 1)
px.LB <- mpm_to_px(matU = LB.mean.kernel$mean_P_site, start = 1)
hx.LB <- mpm_to_hx(matU = LB.mean.kernel$mean_P_site, start = 1)
mx.LB <- mpm_to_mx(matU = LB.mean.kernel$mean_P_site, matR = LB.mean.kernel$mean_F_site, start = 1)

quartz(8,8)
plot(lx, xlab="Survival time (years)", ylab="Survivorship", type="s", col="black")
lines(lx.C, type="s", col="forestgreen")
lines(lx.Q, type="s", col="darkcyan")
lines(lx.EB, type="s", col="gold")
lines(lx.MB, type="s", col="orange")
lines(lx.LB, type="s", col="red")

#A value of 0 indicates only one reproduction event in the entire life time (also known as semelparity), 
#and larger values indicate a more even distribution of reproductive events over the life cycle
entropy_d(lx, mx)   # Demetrius' entropy
entropy_d(lx.C, mx.C)   # Demetrius' entropy
entropy_d(lx.Q, mx.Q)   # Demetrius' entropy
entropy_d(lx.EB, mx.EB)   # Demetrius' entropy
entropy_d(lx.MB, mx.MB)   # Demetrius' entropy
entropy_d(lx.LB, mx.LB)   # Demetrius' entropy

#Values greater than 1 correspond to survival curves type I (high survival in early part of life and decreasing over time), 
#a value of 1 to survival curve type II (constant survival), 
#and values less than 1 to survival curve type III (low survival early and increasing over time) (Salguero-Gómez et al. 2016)
entropy_k(lx)       # Keyfitz' entropy
entropy_k(lx.C)       # Keyfitz' entropy
entropy_k(lx.Q)       # Keyfitz' entropy
entropy_k(lx.EB)       # Keyfitz' entropy
entropy_k(lx.MB)       # Keyfitz' entropy
entropy_k(lx.LB)       # Keyfitz' entropy

#Shape pf mortality and fecundity
shape_surv(lx)      # shape of survival/mortality trajectory
shape_surv(lx.C)      # shape of survival/mortality trajectory
shape_surv(lx.Q)      # shape of survival/mortality trajectory
shape_surv(lx.EB)      # shape of survival/mortality trajectory
shape_surv(lx.MB)      # shape of survival/mortality trajectory
shape_surv(lx.LB)      # shape of survival/mortality trajectory


shape_rep(mx)       # shape of fecundity trajectory
shape_rep(mx.C)       # shape of fecundity trajectory
shape_rep(mx.Q)       # shape of fecundity trajectory
shape_rep(mx.EB)       # shape of fecundity trajectory
shape_rep(mx.MB)       # shape of fecundity trajectory
shape_rep(mx.LB)       # shape of fecundity trajectory

#####################
#For monthly kernels#
#####################
P_F_array <- function(ipm,sites,time){
  stoch.P <- array(0,dim = c(100,100,time))
  stoch.F <- array(0,dim = c(100,100,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(100,100,sites*time*2))
  
  for(j in 1:time){
    
    for(i in 1:sites){
      stoch.P[,,j] <- ipm_array[,,i+seq.t[j]] 
      stoch.F[,,j] <- ipm_array[,,(i)+(sites+seq.t[j])]
    }
  }
  
  return(list(stoch.P=stoch.P, stoch.F=stoch.F))
}
P_F_C_ipm <- P_F_array(C_ipm,1,170)
P_F_Q_ipm <- P_F_array(Q_ipm,1,170)
P_F_EB_ipm <- P_F_array(EB_ipm,1,170)
P_F_MB_ipm <- P_F_array(MB_ipm,1,170)
P_F_LB_ipm <- P_F_array(LB_ipm,1,170)

#Survival and lifespan traits
(life.expect.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),life_expect_mean,  start = 1)))  # mean life expectancy
(life.expect.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),life_expect_mean,  start = 1)))  # mean life expectancy
(life.expect.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),life_expect_mean,  start = 1)))  # mean life expectancy
(life.expect.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),life_expect_mean,  start = 1)))  # mean life expectancy
(life.expect.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),life_expect_mean,  start = 1)))  # mean life expectancy

(life.expect.var.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),life_expect_var,  start = 1)))  # mean life expectancy
(life.expect.var.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),life_expect_var,  start = 1)))  # mean life expectancy
(life.expect.var.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),life_expect_var,  start = 1)))  # mean life expectancy
(life.expect.var.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),life_expect_var,  start = 1)))  # mean life expectancy
(life.expect.var.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),life_expect_var,  start = 1)))  # mean life expectancy

(longev.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),longevity,  start = 1)))  # longevity
(longev.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),longevity,  start = 1)))  # longevity
(longev.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),longevity,  start = 1)))  # longevity
(longev.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),longevity,  start = 1)))  # longevity
(longev.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),longevity,  start = 1)))  # longevity

(gen.time.C <- sapply(1:170,function(i) gen_time(P_F_C_ipm$stoch.P[,,i],P_F_C_ipm$stoch.F[,,i])))  # Generation time
(gen.time.Q <- sapply(1:170,function(i) gen_time(P_F_Q_ipm$stoch.P[,,i],P_F_Q_ipm$stoch.F[,,i])))  # Generation time
(gen.time.EB <- sapply(1:170,function(i) gen_time(P_F_EB_ipm$stoch.P[,,i],P_F_EB_ipm$stoch.F[,,i])))  # Generation time
(gen.time.MB <- sapply(1:170,function(i) gen_time(P_F_MB_ipm$stoch.P[,,i],P_F_MB_ipm$stoch.F[,,i])))  # Generation time
(gen.time.LB <- sapply(1:170,function(i) gen_time(P_F_LB_ipm$stoch.P[,,i],P_F_LB_ipm$stoch.F[,,i])))  # Generation time


#Reproduction and maturation traits
(repro.value.C <- sapply(1:170,function(i) net_repro_rate(P_F_C_ipm$stoch.P[,,i],P_F_C_ipm$stoch.F[,,i])))
(repro.value.Q <- sapply(1:170,function(i) net_repro_rate(P_F_Q_ipm$stoch.P[,,i],P_F_Q_ipm$stoch.F[,,i])))
(repro.value.EB <- sapply(1:170,function(i) net_repro_rate(P_F_EB_ipm$stoch.P[,,i],P_F_EB_ipm$stoch.F[,,i])))
(repro.value.MB <- sapply(1:170,function(i) net_repro_rate(P_F_MB_ipm$stoch.P[,,i],P_F_MB_ipm$stoch.F[,,i])))
(repro.value.LB <- sapply(1:170,function(i) net_repro_rate(P_F_LB_ipm$stoch.P[,,i],P_F_LB_ipm$stoch.F[,,i])))

(net.repro.rate.C <- sapply(1:170,function(i) net_repro_rate(P_F_C_ipm$stoch.P[,,i],P_F_C_ipm$stoch.F[,,i], method = "start")))
(net.repro.rate.Q <- sapply(1:170,function(i) net_repro_rate(P_F_Q_ipm$stoch.P[,,i],P_F_Q_ipm$stoch.F[,,i], method = "start")))
(net.repro.rate.EB <- sapply(1:170,function(i) net_repro_rate(P_F_EB_ipm$stoch.P[,,i],P_F_EB_ipm$stoch.F[,,i], method = "start")))
(net.repro.rate.MB <- sapply(1:170,function(i) net_repro_rate(P_F_MB_ipm$stoch.P[,,i],P_F_MB_ipm$stoch.F[,,i], method = "start")))
(net.repro.rate.LB <- sapply(1:170,function(i) net_repro_rate(P_F_LB_ipm$stoch.P[,,i],P_F_LB_ipm$stoch.F[,,i], method = "start")))

#Life-table component traits
(lx.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),mpm_to_lx,  start = 1)))
(px.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),mpm_to_px,  start = 1)))
(hx.C <- (apply(P_F_C_ipm$stoch.P, MARGIN = c(3),mpm_to_hx,  start = 1)))
(mx.C <- sapply(1:170,function(i) mpm_to_mx(P_F_C_ipm$stoch.P[,,i],P_F_C_ipm$stoch.F[,,i])))

(lx.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),mpm_to_lx,  start = 1)))
(px.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),mpm_to_px,  start = 1)))
(hx.Q <- (apply(P_F_Q_ipm$stoch.P, MARGIN = c(3),mpm_to_hx,  start = 1)))
(mx.Q <- sapply(1:170,function(i) mpm_to_mx(P_F_Q_ipm$stoch.P[,,i],P_F_Q_ipm$stoch.F[,,i])))

(lx.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),mpm_to_lx,  start = 1)))
(px.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),mpm_to_px,  start = 1)))
(hx.EB <- (apply(P_F_EB_ipm$stoch.P, MARGIN = c(3),mpm_to_hx,  start = 1)))
(mx.EB <- sapply(1:170,function(i) mpm_to_mx(P_F_EB_ipm$stoch.P[,,i],P_F_EB_ipm$stoch.F[,,i])))

(lx.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),mpm_to_lx,  start = 1)))
(px.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),mpm_to_px,  start = 1)))
(hx.MB <- (apply(P_F_MB_ipm$stoch.P, MARGIN = c(3),mpm_to_hx,  start = 1)))
(mx.MB <- sapply(1:170,function(i) mpm_to_mx(P_F_MB_ipm$stoch.P[,,i],P_F_MB_ipm$stoch.F[,,i])))

(lx.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),mpm_to_lx,  start = 1)))
(px.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),mpm_to_px,  start = 1)))
(hx.LB <- (apply(P_F_LB_ipm$stoch.P, MARGIN = c(3),mpm_to_hx,  start = 1)))
(mx.LB <- sapply(1:170,function(i) mpm_to_mx(P_F_LB_ipm$stoch.P[,,i],P_F_LB_ipm$stoch.F[,,i])))

# Demetrius' entropy
#A value of 0 indicates only one reproduction event in the entire life time (also known as semelparity), 
#and larger values indicate a more even distribution of reproductive events over the life cycle
(semel.C <- sapply(1:170,function(i) entropy_d(lx.C[[i]],mx.C[[i]])))
(semel.Q <- sapply(1:170,function(i) entropy_d(lx.Q[[i]],mx.Q[[i]])))
(semel.EB <- sapply(1:170,function(i) entropy_d(lx.EB[[i]],mx.EB[[i]])))
(semel.MB <- sapply(1:170,function(i) entropy_d(lx.MB[[i]],mx.MB[[i]])))
(semel.LB <- sapply(1:170,function(i) entropy_d(lx.LB[[i]],mx.LB[[i]])))

# Keyfitz' entropy
#Values greater than 1 correspond to survival curves type I (high survival in early part of life and decreasing over time), 
#a value of 1 to survival curve type II (constant survival), 
#and values less than 1 to survival curve type III (low survival early and increasing over time) (Salguero-Gómez et al. 2016)
(surv.curv.C <- sapply(1:170,function(i) entropy_k(lx.C[[i]])))
(surv.curv.Q <- sapply(1:170,function(i) entropy_k(lx.Q[[i]])))
(surv.curv.EB <- sapply(1:170,function(i) entropy_k(lx.EB[[i]])))
(surv.curv.MB <- sapply(1:170,function(i) entropy_k(lx.MB[[i]])))
(surv.curv.LB <- sapply(1:170,function(i) entropy_k(lx.LB[[i]])))

# shape of survival/mortality trajectory
(shape.surv.C <- (sapply(1:170,function(i) if(length(lx.C[[i]])>2){
  shape_surv(lx.C[[i]])}
  else{NA}
)))
(shape.surv.Q <- (sapply(1:170,function(i) if(length(lx.Q[[i]])>2){
  shape_surv(lx.Q[[i]])}
  else{NA}
)))
(shape.surv.EB <- (sapply(1:170,function(i) if(length(lx.EB[[i]])>2){
  shape_surv(lx.EB[[i]])}
  else{NA}
)))
(shape.surv.MB <- (sapply(1:170,function(i) if(length(lx.MB[[i]])>2){
  shape_surv(lx.MB[[i]])}
  else{NA}
)))

(shape.surv.LB <- (sapply(1:170,function(i) if(length(lx.LB[[i]])>2){
  shape_surv(lx.LB[[i]])}
  else{NA}
)))

# shape of fecundity trajectory
(shape.rep.C <- sapply(1:170,function(i) if(length(mx.C[[i]])>2){
  shape_rep(mx.C[[i]])}
  else{NA}
))
(shape.rep.Q <- sapply(1:170,function(i) if(length(mx.Q[[i]])>2){
  shape_rep(mx.Q[[i]])}
  else{NA}
))
(shape.rep.EB <- sapply(1:170,function(i) if(length(mx.EB[[i]])>2){
  shape_rep(mx.EB[[i]])}
  else{NA}
))
(shape.rep.MB <- sapply(1:170,function(i) if(length(mx.MB[[i]])>2){
  shape_rep(mx.MB[[i]])}
  else{NA}
))
(shape.rep.LB <- sapply(1:170,function(i) if(length(mx.LB[[i]])>2){
  shape_rep(mx.LB[[i]])}
  else{NA}
))


#Resilience parameters#
#######################

#Testing assumptions
isErgodic(mean.kernel$mean_K_site)
isIrreducible(mean.kernel$mean_K_site)
isPrimitive(matrix(unlist(mean.kernel$mean_K_site),100,100))

isErgodic(C.mean.kernel$mean_K_site)
isIrreducible(C.mean.kernel$mean_K_site)
isPrimitive(C.mean.kernel$mean_K_site)

isErgodic(Q.mean.kernel$mean_K_site)
isIrreducible(Q.mean.kernel$mean_K_site)
isPrimitive(Q.mean.kernel$mean_K_site)

isErgodic(EB.mean.kernel$mean_K_site)
isIrreducible(EB.mean.kernel$mean_K_site)
isPrimitive(EB.mean.kernel$mean_K_site)

isErgodic(MB.mean.kernel$mean_K_site)
isIrreducible(MB.mean.kernel$mean_K_site)
isPrimitive(MB.mean.kernel$mean_K_site)

isErgodic(LB.mean.kernel$mean_K_site)
isIrreducible(LB.mean.kernel$mean_K_site)
isPrimitive(as.matrix(unlist(LB.mean.kernel$mean_K_site),100,100))

#One is imprimitive!!!

#Reactivity (first-timestep amplification) and first-time step attenuation
(r.up.K <- reac(mean.kernel$mean_K_site, bound = "upper"))
(r.low.K <- reac(mean.kernel$mean_K_site, bound = "lower"))

(r.up.C.K <- reac(C.mean.kernel$mean_K_site, bound = "upper"))
(r.low.C.K <- reac(C.mean.kernel$mean_K_site, bound = "lower"))

(r.up.Q.K <- reac(Q.mean.kernel$mean_K_site, bound = "upper"))
(r.low.Q.K <- reac(Q.mean.kernel$mean_K_site, bound = "lower"))

(r.up.EB.K <- reac(EB.mean.kernel$mean_K_site, bound = "upper"))
(r.low.EB.K <- reac(EB.mean.kernel$mean_K_site, bound = "lower"))

(r.up.MB.K <- reac(MB.mean.kernel$mean_K_site, bound = "upper"))
(r.low.MB.K <- reac(MB.mean.kernel$mean_K_site, bound = "lower"))

(r.up.LB.K <- reac(LB.mean.kernel$mean_K_site, bound = "upper"))
(r.low.LB.K <- reac(LB.mean.kernel$mean_K_site, bound = "lower"))

#Upper and Lower Inertia (does not run for imprimitive matrices)
(in.up.K <- inertia(mean.kernel$mean_K_site, bound = "upper"))
(in.low.K <- inertia(mean.kernel$mean_K_site, bound = "lower"))

(in.up.C.K <- inertia(C.mean.kernel$mean_K_site, bound = "upper"))
(in.low.C.K <- inertia(C.mean.kernel$mean_K_site, bound = "lower"))

(in.up.Q.K <- inertia(Q.mean.kernel$mean_K_site, bound = "upper"))
(in.low.Q.K <- inertia(Q.mean.kernel$mean_K_site, bound = "lower"))

(in.up.EB.K <- inertia(EB.mean.kernel$mean_K_site, bound = "upper"))
(in.low.EB.K <- inertia(EB.mean.kernel$mean_K_site, bound = "lower"))

(in.up.MB.K <- inertia(MB.mean.kernel$mean_K_site, bound = "upper"))
(in.low.MB.K <- inertia(MB.mean.kernel$mean_K_site, bound = "lower"))

(in.up.LB.K <- inertia(LB.mean.kernel$mean_K_site, bound = "upper"))
(in.low.LB.K <- inertia(LB.mean.kernel$mean_K_site, bound = "lower"))

#Maximum amplification and attenuation
(max.up.K <- maxamp(matrix(unlist(mean.kernel$mean_K_site),100,100)))
(max.low.K <- maxatt(matrix(unlist(mean.kernel$mean_K_site),100,100)))

(max.up.C.K <- maxamp(matrix(unlist(C.mean.kernel$mean_K_site),100,100)))
(max.low.C.K <- maxatt(matrix(unlist(C.mean.kernel$mean_K_site),100,100)))

(max.up.Q.K <- maxamp(matrix(unlist(Q.mean.kernel$mean_K_site),100,100)))
(max.low.Q.K <- maxatt(matrix(unlist(Q.mean.kernel$mean_K_site),100,100)))

(max.up.EB.K <- maxamp(matrix(unlist(EB.mean.kernel$mean_K_site),100,100)))
(max.low.EB.K <- maxatt(matrix(unlist(EB.mean.kernel$mean_K_site),100,100)))

(max.up.MB.K <- maxamp(matrix(unlist(MB.mean.kernel$mean_K_site),100,100)))
(max.low.MB.K <- maxatt(matrix(unlist(MB.mean.kernel$mean_K_site),100,100)))

(max.up.LB.K <- maxamp(matrix(unlist(LB.mean.kernel$mean_K_site),100,100)))
(max.low.LB.K <- maxatt(matrix(unlist(LB.mean.kernel$mean_K_site),100,100)))


#Damping ratio and convergence time
(dr.K <- dr(mean.kernel$mean_K_site, return.time = T))

(dr.K.C <- dr(C.mean.kernel$mean_K_site, return.time = T))
(dr.K.Q <- dr(Q.mean.kernel$mean_K_site, return.time = T))
(dr.K.EB <- dr(EB.mean.kernel$mean_K_site, return.time = T))
(dr.K.MB <- dr(MB.mean.kernel$mean_K_site, return.time = T))
(dr.K.LB <- dr(LB.mean.kernel$mean_K_site, return.time = T))

#Kreiss bounds
(kreiss.up.K<- Kreiss(mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K<- Kreiss(mean.kernel$mean_K_site, bound = "lower"))

(kreiss.up.K.C<- Kreiss(C.mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K.C<- Kreiss(C.mean.kernel$mean_K_site, bound = "lower"))

(kreiss.up.K.Q<- Kreiss(Q.mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K.Q<- Kreiss(Q.mean.kernel$mean_K_site, bound = "lower"))

(kreiss.up.K.EB<- Kreiss(EB.mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K.EB<- Kreiss(EB.mean.kernel$mean_K_site, bound = "lower"))

(kreiss.up.K.MB<- Kreiss(MB.mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K.MB<- Kreiss(MB.mean.kernel$mean_K_site, bound = "lower"))

(kreiss.up.K.LB<- Kreiss(LB.mean.kernel$mean_K_site, bound = "upper"))
(kreiss.low.K.LB<- Kreiss(LB.mean.kernel$mean_K_site, bound = "lower"))

#Using monthly stochastic kernels#
##################################

#Reactivity (first-timestep amplification) and first-time step attenuation
(r.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))


# #Maximum amplification and attenuation
# (max.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))


#Damping ratio and convergence time
dr.stoch <- function(K, n){
  dr.K.stoch <- apply(array(unlist(K),dim = c(100,100,n)), MARGIN = c(3), FUN =  dr , return.time = T)
  dr.K.stoch.dr <- dr.K.stoch.t<- rep(NA, n)
  for(i in 1:n){
    dr.K.stoch.dr[i] <- dr.K.stoch[[i]]$dr
    dr.K.stoch.t[i] <- dr.K.stoch[[i]]$t
    dr.K.stoch.list <- list(dr = dr.K.stoch.dr, t = dr.K.stoch.t)
  }
  return(dr.K.stoch.list)
}

(dr.K.stoch <- dr.stoch(stoch.K, 850))

(dr.K.C.stoch <- dr.stoch(stoch.K.C, 170))
(dr.K.Q.stoch <- dr.stoch(stoch.K.Q, 170))
(dr.K.EB.stoch <- dr.stoch(stoch.K.EB, 170))
(dr.K.MB.stoch <- dr.stoch(stoch.K.MB, 170))
(dr.K.LB.stoch <- dr.stoch(stoch.K.LB, 170))

# #Kreiss bounds
# (kreiss.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(100,100,850)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.C.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.C.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.Q.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.Q.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.EB.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.EB.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.MB.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.MB.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.LB.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.LB.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(100,100,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))


par(mfrow=c(1,1))
Lizpd <- project(mean.kernel$mean_K_site, "diri", time = 12,
                 standard.A = TRUE)
plot(Lizpd, plottype = "shady", bounds = T, log = "y", bty="n")

Lizpd.C <- project(C.mean.kernel$mean_K_site, "diri", time = 12,
                   standard.A = TRUE)
plot(Lizpd.C, plottype = "shady", bounds = T, log = "y", bty="n")

Lizpd.Q <- project(Q.mean.kernel$mean_K_site, "diri", time = 12,
                   standard.A = TRUE)
plot(Lizpd.Q, plottype = "shady", bounds = T, log = "y", bty="n")

Lizpd.EB <- project(EB.mean.kernel$mean_K_site, "diri", time = 12,
                    standard.A = TRUE)
plot(Lizpd.EB, plottype = "shady", bounds = T, log = "y", bty="n")

Lizpd.MB <- project(MB.mean.kernel$mean_K_site, "diri", time = 12,
                    standard.A = TRUE)
plot(Lizpd.MB, plottype = "shady", bounds = T, log = "y", bty="n")

Lizpd.LB <- project(LB.mean.kernel$mean_K_site, "diri", time = 12,
                    standard.A = TRUE)
plot(Lizpd.LB, plottype = "shady", bounds = T, log = "y", bty="n")

bounds.ipms <- as.data.frame(rbind(bounds(Lizpd.C),bounds(Lizpd.Q),bounds(Lizpd.EB),
                                   bounds(Lizpd.MB),bounds(Lizpd.LB)))

bounds.ipms$plot <- rep(c(1:5), each = 13)
bounds.ipms$time <- rep(c(1:13), 5)

ggplot(bounds.ipms, aes(x = time, y = V1, color = as.factor(plot))) +
  geom_line()+
  geom_line(aes(y = V2))+
  scale_color_viridis_d(option = "H") +
  labs(x = "Time (months)", y = "Population growth")

(res.lh.param.Ti <- data.frame(species = c(rep("T_itambere", 850)),
                               plot = rep(c("C", "Q", "EB", "MB", "LB"),each = 170),
                               life.expect = c(life.expect.C, life.expect.Q, life.expect.EB, life.expect.MB, life.expect.LB),
                               life.expect.var = c(life.expect.var.C, life.expect.var.Q, life.expect.var.EB, life.expect.var.MB, life.expect.var.LB),
                               longev = c(longev.C, longev.Q, longev.EB, longev.MB, longev.LB),
                               gen.time = c(gen.time.C, gen.time.Q, gen.time.EB, gen.time.MB, gen.time.LB),
                               net.repro = c(net.repro.rate.C, net.repro.rate.Q, net.repro.rate.EB, net.repro.rate.MB, net.repro.rate.LB),
                               repro.value = c(repro.value.C, repro.value.Q, repro.value.EB, repro.value.MB, repro.value.LB),
                               semel = c(semel.C, semel.Q, semel.EB, semel.MB, semel.LB),
                               surv.curv = c(surv.curv.C, surv.curv.Q, surv.curv.EB, surv.curv.MB, surv.curv.LB),
                               shape.surv = c(shape.surv.C, shape.surv.Q, shape.surv.EB, shape.surv.MB, shape.surv.LB),
                               shape.rep = c(shape.rep.C, shape.rep.Q, shape.rep.EB, shape.rep.MB, shape.rep.LB),
                               fst.amp = c(r.up.C.K.stoch, r.up.Q.K.stoch, r.up.EB.K.stoch, r.up.MB.K.stoch, r.up.LB.K.stoch),
                               fst.att = c(r.low.C.K.stoch, r.low.Q.K.stoch, r.low.EB.K.stoch, r.low.MB.K.stoch, r.low.LB.K.stoch),
                               # max.amp = c(max.up.C.K.stoch, max.up.Q.K.stoch, max.up.EB.K.stoch, max.up.MB.K.stoch, max.up.LB.K.stoch),
                               # max.att = c(max.low.C.K.stoch, max.low.Q.K.stoch, max.low.EB.K.stoch, max.low.MB.K.stoch, max.low.LB.K.stoch),
                               #kreiss.up = c(kreiss.up.K.C.stoch, kreiss.up.K.Q.stoch, kreiss.up.K.EB.stoch, kreiss.up.K.MB.stoch, kreiss.up.K.LB.stoch),
                               #kreiss.low = c(kreiss.low.K.C.stoch, kreiss.low.K.Q.stoch, kreiss.low.K.EB.stoch, kreiss.low.K.MB.stoch, kreiss.low.K.LB.stoch),
                               damping.ratio = c(dr.K.C.stoch$dr, dr.K.Q.stoch$dr, dr.K.EB.stoch$dr, dr.K.MB.stoch$dr, dr.K.LB.stoch$dr),
                               recovery.time = c(dr.K.C.stoch$t, dr.K.Q.stoch$t, dr.K.EB.stoch$t, dr.K.MB.stoch$t, dr.K.LB.stoch$t)
))

cor(res.lh.param.Ti[,c(3:12)],use="na.or.complete", method = "spearman") 
cor(res.lh.param.Ti[,c(13:16)],use="na.or.complete", method = "spearman")     
saveRDS(res.lh.param.Ti, "res.lh.param.Ti.rds")

#########################################################################
#Bring more perturbation analyses from IPM Book - parameter perturbation#
#########################################################################
rm(list = ls())

Titambere.data <- readRDS("Titambere.data.rds")

vitalrates.Titambere.samples <- readRDS("sims_list_Titambere.rds")

vitalrates.Titambere.df <- read.csv("results.vitalrates.Titambere.df_100000iters.csv")

pradel.Titambere.df <- read.csv("results.pradel.itambere.df_400000_noecophys.csv")

fec.Titambere.df <- read.csv("results_fecund_itambere_df.csv")


inv_logit <- function(x) {
  return(
    1/(1 + exp(-(x)))
  )
}


pois_r <- function(x) {
  return(
    exp(x)
  )
}

#Function to estimate size from age
age_to_size <- function(x,mu.L0,mu.LI,K) mu.L0 + (mu.LI-mu.L0)*(1-inv_logit(K)^x)

#Function to estimate age from size
size_to_age <- function(x,mu.L0,mu.LI,K) log(1-((x - mu.L0)/(mu.LI - mu.L0)))/log(inv_logit(K))

#Function to estimate size in t1 from size in t0
sizet0_t1 <- function(x,mu.L0,mu.LI,K) age_to_size(size_to_age(x,mu.L0,mu.LI,K)+1,mu.L0,mu.LI,K)

#Variance in growth#
####################

sd_growth <- function(x,mu.L0,mu_LI,site){
  mean.values <- sizet0_t1(x,mu.L0,
                           mu_LI,
                           unlist(vitalrates.Titambere.samples[,paste0('mu.K[',site,']')]))
  return(sd(mean.values,na.rm=T))
}


my_funs <- list(inv_logit   = inv_logit,
                pois_r      = pois_r,
                sizet0_t1 = sizet0_t1,
                sd_growth = sd_growth)

#Environmental variation
f.pradel <- pradel.Titambere.df[grep(pattern = "f", x = pradel.Titambere.df$X)[1:850],]
phi.pradel <- pradel.Titambere.df[grep(pattern = "phi", x = pradel.Titambere.df$X)[1:850],]

surv.pradel <- matrix(phi.pradel$mean,nrow = 5, ncol = 170)
recr.pradel <- matrix(f.pradel$mean,nrow = 5, ncol = 170)
recr.pradel[,170] <- 0.0001
env.states <- array(c(Titambere.data$amb,surv.pradel,recr.pradel), dim = c(5,170,12))

env.states[,,11]

sample_env <- function(env_states, site, iteration) {
  
  out <- as.list(env_states[site, iteration, ])
  names(out) <- c(paste0("tmed2m_",site), 
                  paste0("RHmax_",site),  
                  paste0("sol_",site),    
                  paste0("tmed0cm_",site),
                  paste0("tmin0cm_",site),
                  paste0("precip_",site) ,
                  paste0("perf_",site)   ,
                  paste0("ha_90_",site)  ,
                  paste0("fire_" ,site)  ,
                  paste0("TSLF_",site)   ,
                  paste0("surv_",site)   ,
                  paste0("f_",site))
  
  
  return(out)
  
}

stoch_K <- function(ipm,sites,time){
  stoch.K <- array(0,dim = c(100,100,sites,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(100,100,sites*time*2))
  #stoch.K.list <- list()
  
  for(j in 1:time){
    
    for(i in 1:sites){
      stoch.K[,,i,j] <- ipm_array[,,i+seq.t[j]] + ipm_array[,,(i)+(sites+seq.t[j])]
      
    }
  }
  list.K <- apply(stoch.K,MARGIN = c(1,2),FUN=c)
  list.K <- aperm(list.K,c(2,3,1))
  list.K2 <- list()
  for(i in 1:(sites*time)){
    list.K2[[i]] <- list.K[,,i] 
  }
  return(list.K2)
}

dr.stoch <- function(K, n){
  dr.K.stoch <- apply(array(unlist(K),dim = c(100,100,n)), MARGIN = c(3), FUN =  dr , return.time = T)
  dr.K.stoch.dr <- dr.K.stoch.t<- rep(NA, n)
  for(i in 1:n){
    dr.K.stoch.dr[i] <- dr.K.stoch[[i]]$dr
    dr.K.stoch.t[i] <- dr.K.stoch[[i]]$t
    dr.K.stoch.list <- list(dr = dr.K.stoch.dr, t = dr.K.stoch.t)
  }
  return(dr.K.stoch.list)
}

#Function to perform Parameter perturbation
###########################################
res_param_perturb <- function(plot,nkernel){
  add.s <- seq(0.0,0.01,0.001)
  ord <- array(c(rep(1,11),rep(0,29*11),
                 rep(0,11),rep(1,11),rep(0,28*11),
                 rep(0,11*2),rep(1,11),rep(0,27*11),
                 rep(0,11*3),rep(1,11),rep(0,26*11),
                 rep(0,11*4),rep(1,11),rep(0,25*11),
                 rep(0,11*5),rep(1,11),rep(0,24*11),
                 rep(0,11*6),rep(1,11),rep(0,23*11),
                 rep(0,11*7),rep(1,11),rep(0,22*11),
                 rep(0,11*8),rep(1,11),rep(0,21*11),
                 rep(0,11*9),rep(1,11),rep(0,20*11),
                 rep(0,11*10),rep(1,11),rep(0,19*11),
                 rep(0,11*11),rep(1,11),rep(0,18*11),
                 rep(0,11*12),rep(1,11),rep(0,17*11),
                 rep(0,11*13),rep(1,11),rep(0,16*11),
                 rep(0,11*14),rep(1,11),rep(0,15*11),
                 rep(0,11*15),rep(1,11),rep(0,14*11),
                 rep(0,11*16),rep(1,11),rep(0,13*11),
                 rep(0,11*17),rep(1,11),rep(0,12*11),
                 rep(0,11*18),rep(1,11),rep(0,11*11),
                 rep(0,11*19),rep(1,11),rep(0,10*11),
                 rep(0,11*20),rep(1,11),rep(0,9*11),
                 rep(0,11*21),rep(1,11),rep(0,8*11),
                 rep(0,11*22),rep(1,11),rep(0,7*11),
                 rep(0,11*23),rep(1,11),rep(0,6*11),
                 rep(0,11*24),rep(1,11),rep(0,5*11),
                 rep(0,11*25),rep(1,11),rep(0,4*11),
                 rep(0,11*26),rep(1,11),rep(0,3*11),
                 rep(0,11*27),rep(1,11),rep(0,2*11),
                 rep(0,11*28),rep(1,11),rep(0,1*11),
                 rep(0,11*29),rep(1,11)),
               dim=c(11,30,30))
  
  res.sens <- list(fst.amp = array(NA,dim=c(170,11,30)), 
                   fst.att = array(NA,dim=c(170,11,30)), 
                   recov.t = array(NA,dim=c(170,11,30)))
  
  for(i in 1:30){
    for(j in 1:11){
      
      
      fixed_list <- list(
        s_mu_slope   = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='beta.phi'] + add.s[j]*ord[j,12,i],    #survival slope
        # s_mu_slope2   = vitalrates.itambere.df['beta2.phi','mean']+ add.s[j]*ord[j,13,i],    #survival slope
        
        #Environmental slopes for survival
        s_mu_tmed2m  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[1]'] + add.s[j]*ord[j,2,i],
        s_mu_RHmax   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[2]']+ add.s[j] *ord[j,3,i],
        s_mu_sol     = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[3]']+ add.s[j] *ord[j,4,i],
        s_mu_tmed0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[4]']+ add.s[j]*ord[j,5,i],
        s_mu_tmin0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[5]']+ add.s[j]*ord[j,6,i],
        s_mu_precip  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[6]']+ add.s[j]*ord[j,7,i],
        s_mu_perf    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[7]']+ add.s[j]*ord[j,8,i],
        s_mu_ha_90   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[8]']+ add.s[j]*ord[j,9,i],
        s_mu_fire    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[9]']+ add.s[j]*ord[j,10,i],
        s_mu_TSLF    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaphiJS[10]']+ add.s[j]*ord[j,11,i],
        
        #Environmental slopes for reproduction
        r_f_mu_tmed2m  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[1]']+ add.s[j]*ord[j,13,i],
        r_f_mu_RHmax   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[2]']+ add.s[j]*ord[j,14,i],
        r_f_mu_sol     = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[3]']+ add.s[j]*ord[j,15,i],
        r_f_mu_tmed0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[4]']+ add.s[j]*ord[j,16,i],
        r_f_mu_tmin0cm = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[5]']+ add.s[j]*ord[j,17,i],
        r_f_mu_precip  = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[6]']+ add.s[j]*ord[j,18,i],
        r_f_mu_perf    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[7]']+ add.s[j]*ord[j,19,i],
        r_f_mu_ha_90   = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[8]']+ add.s[j]*ord[j,20,i],
        r_f_mu_fire    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[9]']+ add.s[j]*ord[j,21,i],
        r_f_mu_TSLF    = pradel.Titambere.df$mean[pradel.Titambere.df$X=='betaf[10]']+ add.s[j]*ord[j,22,i],
        
        #Probability of reproduction
        r_r_mu_int   = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='alpha.prep'] + add.s[j]*ord[j,23,i],
        r_r_mu_slope = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='beta1.prep'] + add.s[j]*ord[j,24,i],  
        # r_r_mu_slope2 = vitalrates.itambere.df['beta2.prep','mean']+ add.s[j]*ord[j,26,i], 
        
        #Number of eggs/embryos
        r_n_mu_int   = fec.Titambere.df$Mean[fec.Titambere.df$X=='alpha.fec']+ add.s[j]*ord[j,25,i], 
        r_n_mu_slope = fec.Titambere.df$Mean[fec.Titambere.df$X=='beta1.fec']+ add.s[j]*ord[j,26,i],   
        r_n_mu_slope2 =fec.Titambere.df$Mean[fec.Titambere.df$X=='beta2.fec']+ add.s[j]*ord[j,27,i],
        
        
        #Size of newborns
        mu_rd     = Titambere.data$mu.L0,   
        sd_rd     = sqrt(Titambere.data$tau.L0),
        mu_LI = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.LI']
        
      )
      
      
      # Now, simulate some random intercepts for growth (g_), survival (s_),
      # and offspring production (r_s_). This part is for the purpose of the example.
      
      # First, we create vector of values that each random component can take.
      s_params  <- list(
        s_g_mu_int_1 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[1]'] + add.s[j]*ord[j,28,i],
        s_g_mu_int_2 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[2]']+ add.s[j]*ord[j,28,i],
        s_g_mu_int_3 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[3]']+ add.s[j]*ord[j,28,i],
        s_g_mu_int_4 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[4]']+ add.s[j]*ord[j,28,i],
        s_g_mu_int_5 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.phiJS[5]']+ add.s[j]*ord[j,28,i]
      )
      
      g_params <- list(

        g_g_mu_K_1 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[1]'] + add.s[j]*ord[j,29,i],
        g_g_mu_K_2 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[2]'] + add.s[j]*ord[j,29,i],
        g_g_mu_K_3 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[3]'] + add.s[j]*ord[j,29,i],
        g_g_mu_K_4 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[4]'] + add.s[j]*ord[j,29,i],
        g_g_mu_K_5 = vitalrates.Titambere.df$mean[vitalrates.Titambere.df$X=='mu.K[5]'] + add.s[j]*ord[j,29,i]
      )
      
      r_params <- list(
        r_f_mu_int_1 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[1]'] + add.s[j]*ord[j,30,i],
        r_f_mu_int_2 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[2]'] + add.s[j]*ord[j,30,i],
        r_f_mu_int_3 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[3]'] + add.s[j]*ord[j,30,i],
        r_f_mu_int_4 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[4]'] + add.s[j]*ord[j,30,i],
        r_f_mu_int_5 = pradel.Titambere.df$mean[pradel.Titambere.df$X=='alpha.f[5]'] + add.s[j]*ord[j,30,i]
      )
      
      all_params_list <- c(fixed_list, g_params, s_params, r_params)
      
      ipm <-  init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "stoch",kern_param = "param") %>%
        define_kernel(
          
          # Our P kernels will vary from site to site, so we index it with "_site"
          
          name             = 'P_site',
          
          # Similarly, our survival and growth functions will vary from site to site
          # so these are also indexed
          
          formula          = s_site * g_site,
          family           = "CC",
          
          # The linear predictor for the survival function can be split out
          # into its own expression as well. This might help keep track of things.
          # Survival is indexed by site as well.
          
          s_lin_site       = (rnorm(1,s_g_mu_int_site,0) + 
                                rnorm(1, s_mu_tmed2m , 0) * tmed2m_site +
                                rnorm(1, s_mu_RHmax  , 0) * RHmax_site +
                                rnorm(1, s_mu_sol    , 0) * sol_site +
                                rnorm(1, s_mu_tmed0cm, 0) * tmed0cm_site +
                                rnorm(1, s_mu_tmin0cm, 0) * tmin0cm_site+
                                rnorm(1, s_mu_precip , 0) * precip_site+
                                rnorm(1, s_mu_perf   , 0) * perf_site+
                                rnorm(1, s_mu_ha_90  , 0) * ha_90_site+
                                rnorm(1, s_mu_fire   , 0) * fire_site+
                                rnorm(1, s_mu_TSLF   , 0) * TSLF_site+
                                rnorm(1, s_mu_slope, 0) * ht_1),
          s_sigma_site     = surv_site - inv_logit(s_lin_site), 
          s_site           =  inv_logit(s_lin_site) + s_sigma_site,
          
          # Again, we modify the vital rate expression to include "_site".
          
          g_site           = dnorm(ht_2, 
                                   mean = sizet0_t1(ht_1, 
                                                    mu_rd,
                                                    mu_LI,
                                                    g_g_mu_K_site), 
                                   sd = sd_growth(ht_1,
                                                  mu_rd,
                                                  mu_LI,
                                                  site)),
          
          data_list        = all_params_list,
          states           = list(c('ht')),
          
          # Here, we tell ipmr that the model has some parameter sets, and
          # provide a list describing the values the index can take. The values in
          # par_set_indices are substituted for "site" everywhere in the model, except
          # for the data list. This is why we had to make sure that the names there
          # matched the levels we supply here.
          
          uses_par_sets    = TRUE,
          par_set_indices  = list(site = plot),
          
          # We must also index the variables in the eviction function
          
          evict_cor        = TRUE,
          evict_fun        = truncated_distributions("norm", "g_site")
          
        ) %>%
        define_kernel(
          
          # The F kernel also varies from site to site
          
          name             = "F_site",
          formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * r_n_site * r_d,
          family           = "CC",
          
          # We didn't include a site level effect for probability
          # of reproduction. Thus, this expression is NOT indexed.
          
          r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                                rnorm(1, r_r_mu_slope, 0) * ht_1),
          r_r_site              = inv_logit(r_r_lin),
          
          # We index the seed production expression with the site effect
          
          r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
                                     rnorm(1,r_n_mu_slope,0) * ht_1 +
                                     rnorm(1,r_n_mu_slope2,0)* ht_1^2),
          r_n_site         = pois_r(r_n_lin_site),
          r_f_lin_site         = ( rnorm(1, r_f_mu_int_site, 0) +
                                     rnorm(1, r_f_mu_tmed2m ,0) * tmed2m_site +
                                     rnorm(1, r_f_mu_RHmax  ,0) * RHmax_site +
                                     rnorm(1, r_f_mu_sol    ,0) * sol_site +
                                     rnorm(1, r_f_mu_tmed0cm,0) * tmed0cm_site +
                                     rnorm(1, r_f_mu_tmin0cm,0) * tmin0cm_site+
                                     rnorm(1, r_f_mu_precip ,0) * precip_site+
                                     rnorm(1, r_f_mu_perf   ,0) * perf_site+
                                     rnorm(1, r_f_mu_ha_90  ,0) * ha_90_site+
                                     rnorm(1, r_f_mu_fire   ,0) * fire_site+
                                     rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site),
          r_f_sigma_site = f_site - pois_r(r_f_lin_site), 
          r_f_site = pois_r(r_f_lin_site) + r_f_sigma_site,
          r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
          data_list        = all_params_list,
          states           = list(c('ht')),
          
          # As in the P kernel, we specify the values the index can have.
          
          uses_par_sets    = TRUE,
          par_set_indices  = list(site = plot),
          evict_cor        = TRUE,
          evict_fun        = truncated_distributions("norm", "r_d")
        ) %>%
        define_impl(
          make_impl_args_list(
            
            # The impl_args are also modified with the index
            
            kernel_names = c("P_site", "F_site"),
            int_rule     = rep("midpoint", 2),
            state_start    = rep("ht", 2),
            state_end      = rep("ht", 2)
          )
        ) %>%
        define_domains(ht = c(20, all_params_list$mu_LI, nkernel)) 
      
      # We also append the suffix in define_pop_state(). THis will create a deterministic
      # simulation for every "site"
      
      
      ipm <- ipm %>%
        define_env_state(env_params = sample_env(env.states, site=site,
                                                 iteration = t), # "t" indexes the current model iteration
                         
                         
                         data_list = list(
                           env.states = env.states,
                           sample_env = sample_env,
                           site = plot
                         )) %>%
        
        define_pop_state(pop_vectors = list(n_ht = runif(nkernel))) %>%
        make_ipm(usr_funs = my_funs,
                 iterate  = TRUE,
                 kernel_seq = rep(plot, each=170),
                 iterations = 170,
                 return_sub_kernels = TRUE,
                 uses_par_sets    = TRUE,
                 par_set_indices  = list(site = plot),)
      
      stoch.K <- stoch_K(ipm, sites = 1,time = 170)
      
      res.sens$fst.amp[,j,i] <- unlist(apply(array(unlist(stoch.K),
                                                   dim = c(nkernel,nkernel,170)),
                                             MARGIN = c(3),
                                             FUN =  reac, bound = "upper", simplify = F))
      res.sens$fst.att[,j,i] <- unlist(apply(array(unlist(stoch.K),
                                                   dim = c(nkernel,nkernel,170)), 
                                             MARGIN = c(3),
                                             FUN =  reac, bound = "lower", simplify = F))
      res.sens$recov.t[,j,i] <- dr.stoch(stoch.K, 170)$t
      
    }
  }
  return(res.sens)
}

res.Titambere.param.sens.C  <- res_param_perturb(1, 100)
res.Titambere.param.sens.Q  <- res_param_perturb(2, 100)
res.Titambere.param.sens.EB <- res_param_perturb(3, 100)
res.Titambere.param.sens.MB <- res_param_perturb(4, 100)
res.Titambere.param.sens.LB <- res_param_perturb(5, 100)

saveRDS(res.Titambere.param.sens.C,"res.Titambere.param.sens.C.rds")
saveRDS(res.Titambere.param.sens.Q,"res.Titambere.param.sens.Q.rds")
saveRDS(res.Titambere.param.sens.EB,"res.Titambere.param.sens.EB.rds")
saveRDS(res.Titambere.param.sens.MB,"res.Titambere.param.sens.MB.rds")
saveRDS(res.Titambere.param.sens.LB,"res.Titambere.param.sens.LB.rds")

res.Titambere.param.sens.C  <- readRDS("res.Titambere.param.sens.C.rds")
res.Titambere.param.sens.Q  <- readRDS("res.Titambere.param.sens.Q.rds")
res.Titambere.param.sens.EB <- readRDS("res.Titambere.param.sens.EB.rds")
res.Titambere.param.sens.MB <- readRDS("res.Titambere.param.sens.MB.rds")
res.Titambere.param.sens.LB <- readRDS("res.Titambere.param.sens.LB.rds")


sens.res <- function(res.array){
  add.s <- seq(0,0.01,0.001)
  res.sens <- array(NA, dim = c(169,10,30))
  for(i in 1:169){
    for(j in 2:11){
      for(k in 1:30){
        res.sens[i,j-1,k] <- (res.array[i,j,k] - res.array[i,1,k]) / add.s[j]
      }
    }
  }
  return(data.frame(X1 = c(res.sens[,,1]),
                    X2 = c(res.sens[,,2]),
                    X3 = c(res.sens[,,3]),
                    X4 = c(res.sens[,,4]),
                    X5 = c(res.sens[,,5]),
                    X6 = c(res.sens[,,6]),
                    X7 = c(res.sens[,,7]),
                    X8 = c(res.sens[,,8]),
                    X9 = c(res.sens[,,9]),
                    X10 = c(res.sens[,,10]),
                    X11 = c(res.sens[,,11]),
                    X12 = c(res.sens[,,12]),
                    X13 = c(res.sens[,,13]),
                    X14 = c(res.sens[,,14]),
                    X15 = c(res.sens[,,15]),
                    X16 = c(res.sens[,,16]),
                    X17 = c(res.sens[,,17]),
                    X18 = c(res.sens[,,18]),
                    X19 = c(res.sens[,,19]),
                    X20 = c(res.sens[,,20]),
                    X21 = c(res.sens[,,21]),
                    X22 = c(res.sens[,,22]),
                    X23 = c(res.sens[,,23]),
                    X24 = c(res.sens[,,24]),
                    X25 = c(res.sens[,,25]),
                    X26 = c(res.sens[,,26]),
                    X27 = c(res.sens[,,27]),
                    X28 = c(res.sens[,,28]),
                    X29 = c(res.sens[,,29]),
                    X30 = c(res.sens[,,30])
                    )
         )
}


(sens.Titambere.fst.amp.C <- sens.res(log10(res.Titambere.param.sens.C$fst.amp)))
(sens.Titambere.fst.amp.Q <- sens.res(log10(res.Titambere.param.sens.Q$fst.amp)))
(sens.Titambere.fst.amp.EB <- sens.res(log10(res.Titambere.param.sens.EB$fst.amp)))
(sens.Titambere.fst.amp.MB <- sens.res(log10(res.Titambere.param.sens.MB$fst.amp)))
(sens.Titambere.fst.amp.LB <- sens.res(log10(res.Titambere.param.sens.LB$fst.amp)))

(sens.Titambere.fst.att.C <- sens.res(res.Titambere.param.sens.C$fst.att))
(sens.Titambere.fst.att.Q <- sens.res(res.Titambere.param.sens.Q$fst.att))
(sens.Titambere.fst.att.EB <- sens.res(res.Titambere.param.sens.EB$fst.att))
(sens.Titambere.fst.att.MB <- sens.res(res.Titambere.param.sens.MB$fst.att))
(sens.Titambere.fst.att.LB <- sens.res(res.Titambere.param.sens.LB$fst.att))

(sens.Titambere.recov.t.C <- sens.res(res.Titambere.param.sens.C$recov.t))
(sens.Titambere.recov.t.Q <- sens.res(res.Titambere.param.sens.Q$recov.t))
(sens.Titambere.recov.t.EB <- sens.res(res.Titambere.param.sens.EB$recov.t))
(sens.Titambere.recov.t.MB <- sens.res(res.Titambere.param.sens.MB$recov.t))
(sens.Titambere.recov.t.LB <- sens.res(res.Titambere.param.sens.LB$recov.t))

#Summary statistics
library(psych)
sens.Titambere.fst.amp.summary <- print(describe(rbind(sens.Titambere.fst.amp.C, 
                                                    sens.Titambere.fst.amp.Q, 
                                                    sens.Titambere.fst.amp.EB,
                                                    sens.Titambere.fst.amp.MB,
                                                    sens.Titambere.fst.amp.LB),
                                              quant = c(0.25, 0.75)),3)


write.csv(sens.Titambere.fst.amp.summary , "sens_Titambere_fst_amp_summary.csv")

sens.Titambere.fst.att.summary <- print(describe(rbind(sens.Titambere.fst.att.C, 
                                                    sens.Titambere.fst.att.Q, 
                                                    sens.Titambere.fst.att.EB,
                                                    sens.Titambere.fst.att.MB,
                                                    sens.Titambere.fst.att.LB),
                                              quant = c(0.25, 0.75)),3)

write.csv(sens.Titambere.fst.att.summary , "sens_Titambere_fst_att_summary.csv")

sens.Titambere.recov.t.summary <- print(describe(rbind(sens.Titambere.recov.t.C, 
                                                    sens.Titambere.recov.t.Q, 
                                                    sens.Titambere.recov.t.EB,
                                                    sens.Titambere.recov.t.MB,
                                                    sens.Titambere.recov.t.LB),
                                              quant = c(0.25, 0.75)),3)

write.csv(sens.Titambere.recov.t.summary , "sens_Titambere_recov_t_summary.csv")

#Plots
#Resilience x Perturbation magnitude
quartz(height = 8, width = 12)
fst.amp.Titambere.param.sens.mean <- (res.Titambere.param.sens.C$fst.amp+
                                     res.Titambere.param.sens.Q$fst.amp +
                                     res.Titambere.param.sens.EB$fst.amp+
                                     res.Titambere.param.sens.MB$fst.amp+
                                     res.Titambere.param.sens.LB$fst.amp)/5


(fst.amp.Titambere.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              log10(apply(fst.amp.Titambere.param.sens.mean[,,c(which(sens.Titambere.fst.amp.summary$mean!=0))], 
                                                          MARGIN = c(2,3), mean))))

fst.amp.Titambere.param.sens.mean <- pivot_longer(data = fst.amp.Titambere.param.sens.mean, 
                                               cols = `X1`:`X6`,
                                               names_to = "parameter",
                                               values_to = "compensation")
quartz(height = 6, width = 8)
ggplot(fst.amp.Titambere.param.sens.mean,
       aes(x = perturb, y = compensation, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(6),name = 'Parameter', 
                      labels = c(expression(beta[phi]*"tmed2m"),
                                 expression(beta[phi]*"RHmax"),
                                 expression(beta[phi]*"sol"),
                                 expression(beta[phi]*"tmed0cm"),
                                 expression(beta[phi]*"tmin0cm"),
                                 expression(beta[phi]*"precip"),
                                 expression(beta[phi]*"perf"),
                                 expression(beta[phi]*"ha"),
                                 expression(beta[phi]*"fire"),
                                 expression(beta[phi]*"TSLF"),
                                 expression(beta[phi]*"SVL"),
                                 #expression(beta[phi]*"SVL"^2),
                                 expression(beta["f"]*"tmed2m"),
                                 expression(beta["f"]*"RHmax"),
                                 expression(beta["f"]*"sol"),
                                 expression(beta["f"]*"tmed0cm"),
                                 expression(beta["f"]*"tmin0cm"),
                                 expression(beta["f"]*"precip"),
                                 expression(beta["f"]*"perf"),
                                 expression(beta["f"]*"ha"),
                                 expression(beta["f"]*"fire"),
                                 expression(beta["f"]*"TSLF"),
                                 expression(alpha["prep"]),
                                 expression(beta["prep"]*"SVL"),
                                 #expression(beta["prep"]*"SVL"^2),
                                 expression(alpha["nb"]),
                                 expression(beta["nb"]*"SVL"),
                                 expression(beta["nb"]*"SVL"^2),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Titambere.fst.amp.summary$mean!=0))-1])+
  scale_y_log10()+
  labs(x = "Perturbation magnitude", y = "Compensation")

fst.att.Titambere.param.sens.mean <- (res.Titambere.param.sens.C$fst.att+
                                     res.Titambere.param.sens.Q$fst.att +
                                     res.Titambere.param.sens.EB$fst.att+
                                     res.Titambere.param.sens.MB$fst.att+
                                     res.Titambere.param.sens.LB$fst.att)/5

(fst.att.Titambere.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              apply(fst.att.Titambere.param.sens.mean[,,c(which(sens.Titambere.fst.att.summary$mean!=0))], 
                                                    MARGIN = c(2,3), mean)))

fst.att.Titambere.param.sens.mean <- pivot_longer(data = fst.att.Titambere.param.sens.mean, 
                                               cols = `X1`:`X6`,
                                               names_to = "parameter",
                                               values_to = "resistance")
quartz(height = 6, width = 8)
ggplot(fst.att.Titambere.param.sens.mean,
       aes(x = perturb, y = resistance, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(6),name = 'Parameter', 
                      labels = c(expression(beta[phi]*"tmed2m"),
                                 expression(beta[phi]*"RHmax"),
                                 expression(beta[phi]*"sol"),
                                 expression(beta[phi]*"tmed0cm"),
                                 expression(beta[phi]*"tmin0cm"),
                                 expression(beta[phi]*"precip"),
                                 expression(beta[phi]*"perf"),
                                 expression(beta[phi]*"ha"),
                                 expression(beta[phi]*"fire"),
                                 expression(beta[phi]*"TSLF"),
                                 expression(beta[phi]*"SVL"),
                                 #expression(beta[phi]*"SVL"^2),
                                 expression(beta["f"]*"tmed2m"),
                                 expression(beta["f"]*"RHmax"),
                                 expression(beta["f"]*"sol"),
                                 expression(beta["f"]*"tmed0cm"),
                                 expression(beta["f"]*"tmin0cm"),
                                 expression(beta["f"]*"precip"),
                                 expression(beta["f"]*"perf"),
                                 expression(beta["f"]*"ha"),
                                 expression(beta["f"]*"fire"),
                                 expression(beta["f"]*"TSLF"),
                                 expression(alpha["prep"]),
                                 expression(beta["prep"]*"SVL"),
                                 #expression(beta["prep"]*"SVL"^2),
                                 expression(alpha["nb"]),
                                 expression(beta["nb"]*"SVL"),
                                 expression(beta["nb"]*"SVL"^2),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Titambere.fst.att.summary$mean!=0))-1])+
  labs(x = "Perturbation magnitude", y = "Resistance")

recov.t.Titambere.param.sens.mean <- (res.Titambere.param.sens.C$recov.t+
                                     res.Titambere.param.sens.Q$recov.t +
                                     res.Titambere.param.sens.EB$recov.t+
                                     res.Titambere.param.sens.MB$recov.t+
                                     res.Titambere.param.sens.LB$recov.t)/5


(recov.t.Titambere.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              apply(recov.t.Titambere.param.sens.mean[,,c(which(sens.Titambere.recov.t.summary$mean!=0))], 
                                                    MARGIN = c(2,3), mean)))

recov.t.Titambere.param.sens.mean <- pivot_longer(data = recov.t.Titambere.param.sens.mean, 
                                               cols = `X1`:`X6`,
                                               names_to = "parameter",
                                               values_to = "recov.t")
quartz(height = 6, width = 8)
ggplot(recov.t.Titambere.param.sens.mean,
       aes(x = perturb, y = recov.t, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(6),name = 'Parameter', 
                      labels = c(expression(beta[phi]*"tmed2m"),
                                 expression(beta[phi]*"RHmax"),
                                 expression(beta[phi]*"sol"),
                                 expression(beta[phi]*"tmed0cm"),
                                 expression(beta[phi]*"tmin0cm"),
                                 expression(beta[phi]*"precip"),
                                 expression(beta[phi]*"perf"),
                                 expression(beta[phi]*"ha"),
                                 expression(beta[phi]*"fire"),
                                 expression(beta[phi]*"TSLF"),
                                 expression(beta[phi]*"SVL"),
                                 #expression(beta[phi]*"SVL"^2),
                                 expression(beta["f"]*"tmed2m"),
                                 expression(beta["f"]*"RHmax"),
                                 expression(beta["f"]*"sol"),
                                 expression(beta["f"]*"tmed0cm"),
                                 expression(beta["f"]*"tmin0cm"),
                                 expression(beta["f"]*"precip"),
                                 expression(beta["f"]*"perf"),
                                 expression(beta["f"]*"ha"),
                                 expression(beta["f"]*"fire"),
                                 expression(beta["f"]*"TSLF"),
                                 expression(alpha["prep"]),
                                 expression(beta["prep"]*"SVL"),
                                 #expression(beta["prep"]*"SVL"^2),
                                 expression(alpha["nb"]),
                                 expression(beta["nb"]*"SVL"),
                                 expression(beta["nb"]*"SVL"^2),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Titambere.recov.t.summary$mean!=0))-1])+
  labs(x = "Perturbation magnitude", y = "Recovery time")

#Exlcuding beta.SVL2
quartz(height = 6, width = 8)
ggplot(recov.t.Titambere.param.sens.mean[recov.t.Titambere.param.sens.mean$parameter!="X5",],
       aes(x = perturb, y = recov.t, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(6)[-5],name = 'Parameter', 
                      labels = c(expression(beta[phi]*"tmed2m"),
                                 expression(beta[phi]*"RHmax"),
                                 expression(beta[phi]*"sol"),
                                 expression(beta[phi]*"tmed0cm"),
                                 expression(beta[phi]*"tmin0cm"),
                                 expression(beta[phi]*"precip"),
                                 expression(beta[phi]*"perf"),
                                 expression(beta[phi]*"ha"),
                                 expression(beta[phi]*"fire"),
                                 expression(beta[phi]*"TSLF"),
                                 expression(beta[phi]*"SVL"),
                                 #expression(beta[phi]*"SVL"^2),
                                 expression(beta["f"]*"tmed2m"),
                                 expression(beta["f"]*"RHmax"),
                                 expression(beta["f"]*"sol"),
                                 expression(beta["f"]*"tmed0cm"),
                                 expression(beta["f"]*"tmin0cm"),
                                 expression(beta["f"]*"precip"),
                                 expression(beta["f"]*"perf"),
                                 expression(beta["f"]*"ha"),
                                 expression(beta["f"]*"fire"),
                                 expression(beta["f"]*"TSLF"),
                                 expression(alpha["prep"]),
                                 expression(beta["prep"]*"SVL"),
                                 #expression(beta["prep"]*"SVL"^2),
                                 expression(alpha["nb"]),
                                 expression(beta["nb"]*"SVL"),
                                 expression(beta["nb"]*"SVL"^2),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Titambere.recov.t.summary$mean!=0)[-5])-1])+
  labs(x = "Perturbation magnitude", y = "Recovery time")


#Barplots
quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.C[,-1],las=2,
        ylab="Compensation",
        main=expression("C - "*italic("T. itambere")),
        names = c(expression(beta[phi]*"tmed2m"),
                      expression(beta[phi]*"RHmax"),
                      expression(beta[phi]*"sol"),
                      expression(beta[phi]*"tmed0cm"),
                      expression(beta[phi]*"tmin0cm"),
                      expression(beta[phi]*"precip"),
                      expression(beta[phi]*"perf"),
                      expression(beta[phi]*"ha"),
                      expression(beta[phi]*"fire"),
                      expression(beta[phi]*"TSLF"),
                      expression(beta[phi]*"SVL"),
                      #expression(beta[phi]*"SVL"^2),
                      expression(beta["f"]*"tmed2m"),
                      expression(beta["f"]*"RHmax"),
                      expression(beta["f"]*"sol"),
                      expression(beta["f"]*"tmed0cm"),
                      expression(beta["f"]*"tmin0cm"),
                      expression(beta["f"]*"precip"),
                      expression(beta["f"]*"perf"),
                      expression(beta["f"]*"ha"),
                      expression(beta["f"]*"fire"),
                      expression(beta["f"]*"TSLF"),
                      expression(alpha["prep"]),
                      expression(beta["prep"]*"SVL"),
                      #expression(beta["prep"]*"SVL"^2),
                      expression(alpha["nb"]),
                      expression(beta["nb"]*"SVL"),
                      expression(beta["nb"]*"SVL"^2),
                      expression(alpha[phi]),
                      expression(mu["K"]),
                      expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.C[,-1],las=2,
        ylab="Compensation", ylim = c(0,0.25),
        main=expression("C - "*italic("T. itambere")),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.Q[, -1],las=2,ylab="Compensation",main=expression("Q - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.Q[, -1],las=2,ylab="Compensation",main=expression("Q - "*italic("T. itambere")),
        ylim=c(0,0.25),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.EB[, -1],las=2,ylab="Compensation",main=expression("EB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.EB[, -1],las=2,ylab="Compensation",main=expression("EB - "*italic("T. itambere")),
        ylim = c(0,.25),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.MB[, -1],las=2,ylab="Compensation",main=expression("MB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.MB[, -1],las=2,ylab="Compensation",main=expression("MB - "*italic("T. itambere")),
        ylim = c(0,.3),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.LB[, -1],las=2,ylab="Compensation",main=expression("LB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.amp.LB[, -1],las=2,ylab="Compensation",main=expression("LB - "*italic("T. itambere")),
        ylim = c(0,.25),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.C[, -1],las=2,ylab="Resistance",main=expression("C - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.C[, -1],las=2,ylab="Resistance",main=expression("C - "*italic("T. itambere")),
        ylim = c(-4,0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.Q[, -1],las=2,ylab="Resistance",main=expression("Q - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.Q[, -1],las=2,ylab="Resistance",main=expression("Q - "*italic("T. itambere")),
        ylim = c(-4, 0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.EB[, -1],las=2,ylab="Resistance",main=expression("EB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.EB[, -1],las=2,ylab="Resistance",main=expression("EB - "*italic("T. itambere")),
        ylim = c(-4,0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.MB[, -1],las=2,ylab="Resistance",main=expression("MB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.MB[, -1],las=2,ylab="Resistance",main=expression("MB - "*italic("T. itambere")),
        ylim = c(-4,0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.LB[, -1],las=2,ylab="Resistance",main=expression("LB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.fst.att.LB[, -1],las=2,ylab="Resistance",main=expression("LB - "*italic("T. itambere")),
        ylim = c(-4,0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.C[, -1],las=2,ylab="Recovery time",main=expression("C - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.C[, -1],las=2,ylab="Recovery time",main=expression("C - "*italic("T. itambere")),
        ylim = c(-250,250),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.Q[, -1],las=2,ylab="Recovery time",main=expression("Q - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.Q[, -1],las=2,ylab="Recovery time",main=expression("Q - "*italic("T. itambere")),
        ylim = c(-250,250),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.EB[, -1],las=2,ylab="Recovery time",main=expression("EB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.EB[, -1],las=2,ylab="Recovery time",main=expression("EB - "*italic("T. itambere")),
        ylim = c(-250,250),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.MB[, -1],las=2,ylab="Recovery time",main=expression("MB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.MB[, -1],las=2,ylab="Recovery time",main=expression("MB - "*italic("T. itambere")),
        ylim = c(-250, 250),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.LB[, -1],las=2,ylab="Recovery time",main=expression("LB - "*italic("T. itambere")),
           names = c(expression(beta[phi]*"tmed2m"),
                         expression(beta[phi]*"RHmax"),
                         expression(beta[phi]*"sol"),
                         expression(beta[phi]*"tmed0cm"),
                         expression(beta[phi]*"tmin0cm"),
                         expression(beta[phi]*"precip"),
                         expression(beta[phi]*"perf"),
                         expression(beta[phi]*"ha"),
                         expression(beta[phi]*"fire"),
                         expression(beta[phi]*"TSLF"),
                         expression(beta[phi]*"SVL"),
                         #expression(beta[phi]*"SVL"^2),
                         expression(beta["f"]*"tmed2m"),
                         expression(beta["f"]*"RHmax"),
                         expression(beta["f"]*"sol"),
                         expression(beta["f"]*"tmed0cm"),
                         expression(beta["f"]*"tmin0cm"),
                         expression(beta["f"]*"precip"),
                         expression(beta["f"]*"perf"),
                         expression(beta["f"]*"ha"),
                         expression(beta["f"]*"fire"),
                         expression(beta["f"]*"TSLF"),
                         expression(alpha["prep"]),
                         expression(beta["prep"]*"SVL"),
                         #expression(beta["prep"]*"SVL"^2),
                         expression(alpha["nb"]),
                         expression(beta["nb"]*"SVL"),
                         expression(beta["nb"]*"SVL"^2),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Titambere.recov.t.LB[, -1],las=2,ylab="Recovery time",main=expression("LB - "*italic("T. itambere")),
        ylim = c(-250, 250),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

#Average among plots
quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Titambere.fst.amp.C[, -1],
              sens.Titambere.fst.amp.Q[, -1],
              sens.Titambere.fst.amp.EB[, -1],
              sens.Titambere.fst.amp.MB[, -1],
              sens.Titambere.fst.amp.LB[, -1]),
        las=2,ylab="log10(Compensation)",main=expression("Average - "*italic("T. itambere")),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Titambere.fst.amp.C[, -1],
              sens.Titambere.fst.amp.Q[, -1],
              sens.Titambere.fst.amp.EB[, -1],
              sens.Titambere.fst.amp.MB[, -1],
              sens.Titambere.fst.amp.LB[, -1]),
        las=2,ylab="log10(Compensation)",main=expression("Average - "*italic("T. itambere")),
        ylim = c(0,25),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Titambere.fst.amp.C[, -1],
              sens.Titambere.fst.amp.Q[, -1],
              sens.Titambere.fst.amp.EB[, -1],
              sens.Titambere.fst.amp.MB[, -1],
              sens.Titambere.fst.amp.LB[, -1]),
        las=2,ylab="log10(Compensation)",main=expression("Average - "*italic("T. itambere")),
        ylim = c(0,0.3),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Titambere.fst.att.C[, -1],
              sens.Titambere.fst.att.Q[, -1],
              sens.Titambere.fst.att.EB[, -1],
              sens.Titambere.fst.att.MB[, -1],
              sens.Titambere.fst.att.LB[, -1]),
        las=2,ylab="log10(Resistance)",main=expression("Average - "*italic("T. itambere")),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Titambere.fst.att.C[, -1],
              sens.Titambere.fst.att.Q[, -1],
              sens.Titambere.fst.att.EB[, -1],
              sens.Titambere.fst.att.MB[, -1],
              sens.Titambere.fst.att.LB[, -1]),
        las=2,ylab="log10(Resistance)",main=expression("Average - "*italic("T. itambere")),
        ylim = c(-4,0),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(log10(rbind(sens.Titambere.recov.t.C[, -1]+1,
              sens.Titambere.recov.t.Q[, -1]+1,
              sens.Titambere.recov.t.EB[, -1]+1,
              sens.Titambere.recov.t.MB[, -1]+1,
              sens.Titambere.recov.t.LB[, -1]+1)),
        las=2,ylab="log10(Recovery time +1)",main=expression("Average - "*italic("T. itambere")),
        names = c(expression(beta[phi]*"tmed2m"),
                  expression(beta[phi]*"RHmax"),
                  expression(beta[phi]*"sol"),
                  expression(beta[phi]*"tmed0cm"),
                  expression(beta[phi]*"tmin0cm"),
                  expression(beta[phi]*"precip"),
                  expression(beta[phi]*"perf"),
                  expression(beta[phi]*"ha"),
                  expression(beta[phi]*"fire"),
                  expression(beta[phi]*"TSLF"),
                  expression(beta[phi]*"SVL"),
                  #expression(beta[phi]*"SVL"^2),
                  expression(beta["f"]*"tmed2m"),
                  expression(beta["f"]*"RHmax"),
                  expression(beta["f"]*"sol"),
                  expression(beta["f"]*"tmed0cm"),
                  expression(beta["f"]*"tmin0cm"),
                  expression(beta["f"]*"precip"),
                  expression(beta["f"]*"perf"),
                  expression(beta["f"]*"ha"),
                  expression(beta["f"]*"fire"),
                  expression(beta["f"]*"TSLF"),
                  expression(alpha["prep"]),
                  expression(beta["prep"]*"SVL"),
                  #expression(beta["prep"]*"SVL"^2),
                  expression(alpha["nb"]),
                  expression(beta["nb"]*"SVL"),
                  expression(beta["nb"]*"SVL"^2),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",15),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",15),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))
