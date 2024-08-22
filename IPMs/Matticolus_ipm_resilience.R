# Build Integral Projection Models (IPMs) and estimate life-history traits and demographic resilience components -------------------------------------------------------------------------

rm(list = ls())

options(timeout = 10000) # increase the time to download large files

#Set working directory
setwd("~/Documents/GitHub/DemographicResilience_CerradoLizards/IPMs")

#Load packages
library(jagsUI)
library(rjags)
library(ipmr)
library(Rage)
library(popdemo)
library(tidyverse)
library(tidybayes)
library(BayesPostEst)
library(MCMCvis)
library(mcmcr)
library(viridis)
library(psych)

#Read results and data to build IPMs
Matticolus.data <- readRDS("Matticolus_data.rds")
vitalrates.Matticolus <- readRDS("results_vitalrates_Matticolus.rds")
#vitalrates.Matticolus.samples <- as.mcmc.list(vitalrates.Matticolus$samples)

vitalrates.Matticolus.samples <- vitalrates.Matticolus$samples

muK.Matticolus.samples <- vitalrates.Matticolus$sims.list$mu.K

rm(vitalrates.Matticolus)
gc()
# vitalrates.Matticolus.df <- MCMCsummary(vitalrates.Matticolus.samples)

# write.csv(vitalrates.Matticolus.df, "results.vitalrates.Matticolus.df_100000iters.csv")
vitalrates.Matticolus.df <- read.csv("results.vitalrates.Matticolus.df_100000iters.csv")

View(vitalrates.Matticolus.df)

f.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "f", x = vitalrates.Matticolus.df$X)[1:850],]

rho.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "rho", x = vitalrates.Matticolus.df$X)[1:850],]

phi.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "phi", x = vitalrates.Matticolus.df$X)[1:850],]

f.vitalrates$plot <- rep(1:5,170)
f.vitalrates$time <- rep(1:170, each = 5)

rho.vitalrates$plot <- rep(1:5,170)
rho.vitalrates$time <- rep(1:170, each = 5)

phi.vitalrates$plot <- rep(1:5,170)
phi.vitalrates$time <- rep(1:170, each = 5)

tapply(rho.vitalrates$mean, rho.vitalrates$plot, FUN = function(x) exp(mean(log(x), na.rm = T)))
table(Matticolus.data$plot)
plot(y = colSums(Matticolus.data$y), x = (1:170)/12, type = "l")

f.vitalrates$plot <- as.factor(f.vitalrates$plot)
phi.vitalrates$plot <- as.factor(phi.vitalrates$plot)

ggplot(f.vitalrates, aes(x = time, y = mean, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  #geom_path(data=f.pradel[f.vitalrates$plot==3,], aes(x = time, y= mean, colour = plot), linetype = "dashed" )+
  ylim(c(0,2.5))+
  scale_color_manual(values=turbo(5))

ggplot(phi.vitalrates, aes(x = time, y = mean,fill = plot, colour = plot))+
  geom_line(aes(x = time, y = mean,colour=plot), alpha=0.5, linewidth = 2) +
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill=plot), alpha=0.4, colour = NA)+
  ylim(c(0.3,1))+
  scale_color_manual(values=turbo(5))

#Use the param option with predefined environments

##################################################################################
#Simple deterministic IPM constructed from discretely varying parameter estimates#
##################################################################################

# Define some fixed parameters

fixed_list <- list(
  s_mu_slope   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='beta.phi'],    #survival slope
  #s_mu_slope2   = vitalrates.Matticolus.df['beta.phi2','mean'],    #survival slope
  
  s_sd_slope   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='beta.phi'],    #survival slope
  #s_sd_slope2  = vitalrates.Matticolus.df['beta.phi2','sd'],    #survival slope
  
  
  #Environmental slopes for survival
  s_mu_tmed2m  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[1]'],
  s_mu_RHmax   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[2]'],
  s_mu_sol     = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[3]'],
  s_mu_tmed0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[4]'],
  s_mu_tmin0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[5]'],
  s_mu_precip  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[6]'],
  s_mu_perf    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[7]'],
  s_mu_ha_90   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[8]'],
  s_mu_fire    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[9]'],
  s_mu_TSLF    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[10]'],
  
  s_sd_tmed2m  = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[1]'],
  s_sd_RHmax   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[2]'],
  s_sd_sol     = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[3]'],
  s_sd_tmed0cm = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[4]'],
  s_sd_tmin0cm = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[5]'],
  s_sd_precip  = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[6]'],
  s_sd_perf    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[7]'],
  s_sd_ha_90   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[8]'],
  s_sd_fire    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[9]'],
  s_sd_TSLF    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaphiJS[10]'],
  
  sigma.phiJS = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='sigma.phiJS'],
  
  #Environmental slopes for reproduction
  r_f_mu_tmed2m  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[1]'],
  r_f_mu_RHmax   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[2]'],
  r_f_mu_sol     = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[3]'],
  r_f_mu_tmed0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[4]'],
  r_f_mu_tmin0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[5]'],
  r_f_mu_precip  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[6]'],
  r_f_mu_perf    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[7]'],
  r_f_mu_ha_90   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[8]'],
  r_f_mu_fire    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[9]'],
  r_f_mu_TSLF    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[10]'],
  
  r_f_sd_tmed2m  = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[1]'],
  r_f_sd_RHmax   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[2]'],
  r_f_sd_sol     = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[3]'],
  r_f_sd_tmed0cm = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[4]'],
  r_f_sd_tmin0cm = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[5]'],
  r_f_sd_precip  = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[6]'],
  r_f_sd_perf    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[7]'],
  r_f_sd_ha_90   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[8]'],
  r_f_sd_fire    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[9]'],
  r_f_sd_TSLF    = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='betaf[10]'],
  
  sigma.f = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='sigma.f'],
  
  
  #Probability of reproduction
  r_r_mu_int   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.prep'],
  r_r_mu_slope = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='beta1.prep'],  
  #r_r_mu_slope2 = vitalrates.Matticolus.df['beta2.fec','mean'], 
  
  r_r_sd_int   = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.prep'],
  r_r_sd_slope = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='beta1.prep'],  
  #r_r_sd_slope2 = vitalrates.Matticolus.df['beta2.prep','sd'], 
  
  #Number of eggs/embryos
  #r_n_mu_int   = vitalrates.Matticolus.df['alpha.fec','mean'], 
  #r_n_mu_slope = vitalrates.Matticolus.df['beta1.fec','mean'],   
  #r_n_mu_slope2 =vitalrates.Matticolus.df['beta2.fec','mean'],  
  
  #r_n_sd_int   = vitalrates.Matticolus.df['alpha.fec','sd'], 
  #r_n_sd_slope = vitalrates.Matticolus.df['beta1.fec','sd'],   
  #r_n_sd_slope2 =vitalrates.Matticolus.df['beta2.fec','sd'],  
  
  
  #Size of newborns
  mu_rd     = Matticolus.data$mu.L0,   
  sd_rd     = sqrt(Matticolus.data$tau.L0),
  mu_LI = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.LI']
  
)


# Now, simulate some random intercepts for growth (g_), survival (s_),
# and offspring production (r_s_). This part is for the purpose of the example.

# First, we create vector of values that each random component can take.
s_params  <- list(
  s_g_mu_int_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[1]'],
  s_g_mu_int_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[2]'],
  s_g_mu_int_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[3]'],
  s_g_mu_int_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[4]'],
  s_g_mu_int_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[5]'],
  
  s_g_sd_int_1 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.phiJS[1]'],
  s_g_sd_int_2 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.phiJS[2]'],
  s_g_sd_int_3 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.phiJS[3]'],
  s_g_sd_int_4 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.phiJS[4]'],
  s_g_sd_int_5 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.phiJS[5]']
)

g_params <- list(
  
  g_g_mu_K_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[1]'],
  g_g_mu_K_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[2]'],
  g_g_mu_K_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[3]'],
  g_g_mu_K_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[4]'],
  g_g_mu_K_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[5]'],
  
  g_g_sd_K_1 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='mu.K[1]'],
  g_g_sd_K_2 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='mu.K[2]'],
  g_g_sd_K_3 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='mu.K[3]'],
  g_g_sd_K_4 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='mu.K[4]'],
  g_g_sd_K_5 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='mu.K[5]']
)

r_params <- list(
  r_f_mu_int_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[1]'],
  r_f_mu_int_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[2]'],
  r_f_mu_int_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[3]'],
  r_f_mu_int_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[4]'],
  r_f_mu_int_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[5]'],
  
  r_f_sd_int_1 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.f[1]'],
  r_f_sd_int_2 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.f[2]'],
  r_f_sd_int_3 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.f[3]'],
  r_f_sd_int_4 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.f[4]'],
  r_f_sd_int_5 = vitalrates.Matticolus.df$sd[vitalrates.Matticolus.df$X=='alpha.f[5]'],
  
  r_p_mu_int_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.pJS[1]'],
  r_p_mu_int_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.pJS[2]'],
  r_p_mu_int_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.pJS[3]'],
  r_p_mu_int_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.pJS[4]'],
  r_p_mu_int_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.pJS[5]']
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

sd_growth <- function(x,mu.L0,mu.LI,site){
  mean.values <- sizet0_t1(x,mu.L0,mu.LI, muK.Matticolus.samples[, site])
  return(sd(mean.values,na.rm=T))
}


my_funs <- list(inv_logit   = inv_logit,
                pois_r      = pois_r,
                sizet0_t1 = sizet0_t1,
                sd_growth = sd_growth)

f.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "f", x = vitalrates.Matticolus.df$X)[1:850],]

phi.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "phi", x = vitalrates.Matticolus.df$X)[1:850],]


f.vitalrates$plot <- rep(1:5,170)
f.vitalrates$time <- rep(1:170, each = 5)

phi.vitalrates$plot <- rep(1:5,170)
phi.vitalrates$time <- rep(1:170, each = 5)

surv.pradel <- matrix(phi.vitalrates$mean,nrow = 5, ncol = 170)
recr.pradel <- matrix(f.vitalrates$mean,nrow = 5, ncol = 170)
recr.pradel[,170] <- 0.0001
env.states <- array(c(Matticolus.data$amb,surv.pradel,recr.pradel), dim = c(5,170,12))

env.states[,,12]

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
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,r_r_sd_int) + 
                          rnorm(1, r_r_mu_slope, r_r_sd_slope) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, r_n_sd_int) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,r_n_sd_slope) * ht_1),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 

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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(1:5, each=170),
           iterations = 170,
           return_sub_kernels = TRUE,
           uses_par_sets    = TRUE,
           par_set_indices  = list(site = 1:5))

lambda(my_ipm,type_lambda = 'all')
lambda(my_ipm,log = F)

library(fields)
quartz(8,12)
mean.kernel<-mean_kernel(my_ipm)

par(mfrow=c(1,2))
plot(mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

par(mfrow=c(1,1))
plot(20:38,sizet0_t1(c(20:38),all_params_list$mu_rd,all_params_list$mu_LI,all_params_list$g_g_mu_K_5))

######################
#Without sd estimates#
######################


################################################################################
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
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 

# We also append the suffix in define_pop_state(). THis will create a deterministic
# simulation for every "site"

my_ipm2 <- my_ipm2 %>%
  define_env_state(env_params = sample_env(env.states, site=site,
                                           iteration = t), # "t" indexes the current model iteration
                   
                   
                   data_list = list(
                     env.states = env.states,
                     sample_env = sample_env,
                     site = 1:5
                   )) %>%
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
  make_ipm(usr_funs = my_funs,
           iterate  = TRUE,
           kernel_seq = rep(1:5, each=170),
           iterations = 170,
           return_sub_kernels = TRUE)

lambda(my_ipm2,type_lambda = 'all')
lambda(my_ipm2,log = F)
quartz(8,12)
par(mfrow=c(1,2))
mean.kernel<-mean_kernel(my_ipm2)
plot(mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))


###############
#For each plot#
###############

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
    par_set_indices  = list(site = 1),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 

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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
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

quartz(8,12)
par(mfrow=c(1,2))
plot(C.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(C.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

#Stable size distributions and reproductive values
plot_ipm_sd_rv <- function(ipm_sd, xmin, xmax, ylab){
  plot(seq(xmin, xmax,length.out=60),ipm_sd[,1],type="l",ylim=c(0,max(ipm_sd[,-ncol(ipm_sd)])),
       xlab = "SVL (mm)", ylab = ylab, bty="n")
  for(i in 2:ncol(ipm_sd)-1){
    lines(seq(xmin, xmax,length.out=60),ipm_sd[,i],col=i)
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
lines(rho.vitalrates$mean[rho.vitalrates$plot==1][-c(1, 170)])

quartz(8,8)
plot(rho.vitalrates$mean[rho.vitalrates$plot==1][-c(1, 170)], lambda(C_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.vitalrates$mean[rho.vitalrates$plot==1][-c(1, 170)], lambda(C_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(C_ipm,type_lambda = 'all')[-c(1,170)], rho.vitalrates$mean[rho.vitalrates$plot==1][-c(1, 170)])

summary(rho.vitalrates$mean[rho.vitalrates$plot==1][-c(1, 170)] - lambda(C_ipm,type_lambda = 'all')[-c(1, 170)])

##############
#Quadriennial#
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
    par_set_indices  = list(site = 2),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 


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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
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

quartz(8,12)
par(mfrow=c(1,2))
plot(Q.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

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

#Lambdas IPM x PJS
quartz(10,10)
plot(lambda(Q_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.6,2.2), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)])

quartz(8,8)
plot(rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)], lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)], lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(Q_ipm,type_lambda = 'all')[-c(1,170)], rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)])

summary(rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)] - lambda(Q_ipm,type_lambda = 'all')[-c(1, 170)])

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
    par_set_indices  = list(site = 3),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 


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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
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

quartz(8,12)
par(mfrow=c(1,2))
plot(EB.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

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
quartz(10,10)
plot(lambda(EB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.5,2.2), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.vitalrates$mean[rho.vitalrates$plot==3][-c(1, 170)])

quartz(8,8)
plot(rho.vitalrates$mean[rho.vitalrates$plot==3][-c(1, 170)], lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.vitalrates$mean[rho.vitalrates$plot==3][-c(1, 170)], lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(EB_ipm,type_lambda = 'all')[-c(1,170)], rho.vitalrates$mean[rho.vitalrates$plot==3][-c(1, 170)])

summary(rho.vitalrates$mean[rho.vitalrates$plot==3][-c(1, 170)] - lambda(EB_ipm,type_lambda = 'all')[-c(1, 170)])

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
    par_set_indices  = list(site = 4),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 


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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
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

quartz(8,12)
par(mfrow=c(1,2))
plot(MB.mean.kernel$mean_P_site, do_contour=T,col=turbo(1000))

plot(MB.mean.kernel$mean_F_site, do_contour=T,col=turbo(1000))

#Stable size distributions and reproductive values
MB_ipm_sd <- right_ev(MB_ipm)
plot_ipm_sd_rv(MB_ipm_sd$ht_w, 
               MB_ipm$proto_ipm$domain[[1]][[1]][1], 
               MB_ipm$proto_ipm$domain[[1]][[1]][2], "Stable distribution")

MB_ipm_rv <- left_ev(MB_ipm, iterations = 170)
plot_ipm_sd_rv(MB_ipm_rv$ht_v, 
               MB_ipm$proto_ipm$domain[[1]][[1]][1], 
               MB_ipm$proto_ipm$domain[[1]][[1]][2], "Reproductive value")
par(mfrow = c(1,1))

#Lambdas IPM x PJS
quartz(10,10)
plot(lambda(MB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.5,2.5), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.vitalrates$mean[rho.vitalrates$plot==4][-c(1, 170)])

quartz(8,8)
plot(rho.vitalrates$mean[rho.vitalrates$plot==4][-c(1, 170)], lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.vitalrates$mean[rho.vitalrates$plot==4][-c(1, 170)], lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(MB_ipm,type_lambda = 'all')[-c(1,170)], rho.vitalrates$mean[rho.vitalrates$plot==4][-c(1, 170)])

summary(rho.vitalrates$mean[rho.vitalrates$plot==4][-c(1, 170)] - lambda(MB_ipm,type_lambda = 'all')[-c(1, 170)])

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
    par_set_indices  = list(site = 5),
    
    # We must also index the variables in the eviction function
    
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
    
  ) %>%
  define_kernel(
    
    # The F kernel also varies from site to site
    
    name             = "F_site",
    formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
    family           = "CC",
    
    # We didn't include a site level effect for probability
    # of reproduction. Thus, this expression is NOT indexed.
    
    r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                          rnorm(1, r_r_mu_slope, 0) * ht_1),
    r_r_site              = inv_logit(r_r_lin),
    
    # We index the seed production expression with the site effect
    # 
    # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
    #                            r_f_lin_site +
    #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
    #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
    # 2         = pois_r(r_n_lin_site),
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
                               rnorm(1, r_f_mu_TSLF   ,0) * TSLF_site +
                               rnorm(1,0,sigma.f)),
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
  define_domains(ht = c(20, all_params_list$mu_LI, 60)) 


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
  
  define_pop_state(pop_vectors = list(n_ht = runif(60))) %>%
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

quartz(8,12)
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
quartz(10,10)
plot(lambda(LB_ipm,type_lambda = 'all')[-c(1,170)], type = "l", ylim = c(0.5,2.5), 
     bty = "n", col = "red", ylab = "Population growth", xlab = "Time (months)")
lines(rho.vitalrates$mean[rho.vitalrates$plot==5][-c(1, 170)])

quartz(8,8)
plot(rho.vitalrates$mean[rho.vitalrates$plot==5][-c(1, 170)], lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)], bty = "n",
     ylab = "Population growth (IPM)", xlab = "Population growth (PJS)", col = rgb(0,0,0,0.5), pch = 19)

ccf(rho.vitalrates$mean[rho.vitalrates$plot==5][-c(1, 170)], lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)])
cor.test(lambda(Q_ipm,type_lambda = 'all')[-c(1,170)], rho.vitalrates$mean[rho.vitalrates$plot==2][-c(1, 170)])

summary(rho.vitalrates$mean[rho.vitalrates$plot==5][-c(1, 170)] - lambda(LB_ipm,type_lambda = 'all')[-c(1, 170)])

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
  stoch.K <- array(0,dim = c(60,60,sites,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(60,60,sites*time*2))
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
  stoch.P <- array(0,dim = c(60,60,time))
  stoch.F <- array(0,dim = c(60,60,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(60,60,sites*time*2))
  
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


#######################
#Resilience parameters#
#######################

#Testing assumptions
isErgodic(mean.kernel$mean_K_site)
isIrreducible(mean.kernel$mean_K_site)
isPrimitive(matrix(unlist(mean.kernel$mean_K_site),60,60))

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
isPrimitive(as.matrix(unlist(LB.mean.kernel$mean_K_site),60,60))

#Perfect!!!

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
(max.up.K <- maxamp(matrix(unlist(mean.kernel$mean_K_site),60,60)))
(max.low.K <- maxatt(matrix(unlist(mean.kernel$mean_K_site),60,60)))

(max.up.C.K <- maxamp(matrix(unlist(C.mean.kernel$mean_K_site),60,60)))
(max.low.C.K <- maxatt(matrix(unlist(C.mean.kernel$mean_K_site),60,60)))

(max.up.Q.K <- maxamp(matrix(unlist(Q.mean.kernel$mean_K_site),60,60)))
(max.low.Q.K <- maxatt(matrix(unlist(Q.mean.kernel$mean_K_site),60,60)))

(max.up.EB.K <- maxamp(matrix(unlist(EB.mean.kernel$mean_K_site),60,60)))
(max.low.EB.K <- maxatt(matrix(unlist(EB.mean.kernel$mean_K_site),60,60)))

(max.up.MB.K <- maxamp(matrix(unlist(MB.mean.kernel$mean_K_site),60,60)))
(max.low.MB.K <- maxatt(matrix(unlist(MB.mean.kernel$mean_K_site),60,60)))

(max.up.LB.K <- maxamp(matrix(unlist(LB.mean.kernel$mean_K_site),60,60)))
(max.low.LB.K <- maxatt(matrix(unlist(LB.mean.kernel$mean_K_site),60,60)))


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
(r.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))

(r.up.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "upper", simplify = F)))
(r.low.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN=  reac, bound = "lower", simplify = F)))


# #Maximum amplification and attenuation
# (max.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.C.K.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.Q.K.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.EB.K.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.MB.K.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))
# 
# (max.up.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxamp, simplify = F)))
# (max.low.LB.K.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN =  maxatt, simplify = F)))


#Damping ratio and convergence time
dr.stoch <- function(K, n){
  dr.K.stoch <- apply(array(unlist(K),dim = c(60,60,n)), MARGIN = c(3), FUN =  dr , return.time = T)
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
# (kreiss.up.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.stoch <- unlist(apply(array(unlist(stoch.K),dim = c(60,60,850)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.C.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.C.stoch <- unlist(apply(array(unlist(stoch.K.C),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.Q.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.Q.stoch <- unlist(apply(array(unlist(stoch.K.Q),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.EB.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.EB.stoch <- unlist(apply(array(unlist(stoch.K.EB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.MB.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.MB.stoch <- unlist(apply(array(unlist(stoch.K.MB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))
# 
# (kreiss.up.K.LB.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "upper", simplify = F)))
# (kreiss.low.K.LB.stoch <- unlist(apply(array(unlist(stoch.K.LB),dim = c(60,60,170)), MARGIN = c(3),FUN=  Kreiss, bound = "lower", simplify = F)))


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

(res.lh.param.Ma <- data.frame(species = c(rep("M_atticolus", 850)),
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

cor(res.lh.param.Ma[,c(3:12)],use="na.or.complete") 
cor(res.lh.param.Ma[,c(13:16)],use="na.or.complete")  

saveRDS(res.lh.param.Ma, "res.lh.param.Ma.rds")

#########################################################################
#Bring more perturbation analyses from IPM Book - parameter perturbation#
#########################################################################
rm(list = ls())

# setwd("/Volumes/Extreme SSD/Heitor/Doutorado/Analises/Cap2_LizardsDemography_Cerrado/Analysis")
library(jagsUI)
library(rjags)
library(ipmr)
library(Rage)
library(popdemo)
library(tidyverse)
library(tidybayes)
library(BayesPostEst)
library(MCMCvis)
library(mcmcr)
library(viridis)

Matticolus.data <- readRDS("Matticolus.data.rds")
vitalrates.Matticolus <- readRDS("results_vitalrates_Matticolus.rds")
#vitalrates.Matticolus.samples <- as.mcmc.list(vitalrates.Matticolus$samples)

muK.Matticolus.samples <- vitalrates.Matticolus$sims.list$mu.K
rm(vitalrates.Matticolus)
gc()
# vitalrates.Matticolus.df <- MCMCsummary(vitalrates.Matticolus.samples)

# write.csv(vitalrates.Matticolus.df, "results.vitalrates.Matticolus.df_100000iters.csv")
vitalrates.Matticolus.df <- read.csv("results.vitalrates.Matticolus.df_100000iters.csv")

#Functions
##########
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
sd_growth <- function(x,mu.L0,mu.LI,site){
  mean.values <- sizet0_t1(x,mu.L0,mu.LI, muK.Matticolus.samples[, site])
  return(sd(mean.values,na.rm=T))
}


my_funs <- list(inv_logit   = inv_logit,
                pois_r      = pois_r,
                sizet0_t1 = sizet0_t1,
                sd_growth = sd_growth)

f.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "f", x = vitalrates.Matticolus.df$X)[1:850],]
phi.vitalrates <- vitalrates.Matticolus.df[grep(pattern = "phi", x = vitalrates.Matticolus.df$X)[1:850],]


f.vitalrates$plot <- rep(1:5,170)
f.vitalrates$time <- rep(1:170, each = 5)

phi.vitalrates$plot <- rep(1:5,170)
phi.vitalrates$time <- rep(1:170, each = 5)

surv.pradel <- matrix(phi.vitalrates$mean,nrow = 5, ncol = 170)
recr.pradel <- matrix(f.vitalrates$mean,nrow = 5, ncol = 170)
recr.pradel[,170] <- 0.0001
env.states <- array(c(Matticolus.data$amb,surv.pradel,recr.pradel), dim = c(5,170,12))

env.states[,,12]

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
  stoch.K <- array(0,dim = c(60,60,sites,time))
  seq.t <- seq(0,(sites*time*2)-(sites*2),sites*2)
  ipm_array <- array(unlist(ipm$sub_kernels),dim=c(60,60,sites*time*2))
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
  dr.K.stoch <- apply(array(unlist(K),dim = c(60,60,n)), MARGIN = c(3), FUN =  dr , return.time = T)
  dr.K.stoch.dr <- dr.K.stoch.t<- rep(NA, n)
  for(i in 1:n){
    dr.K.stoch.dr[i] <- dr.K.stoch[[i]]$dr
    dr.K.stoch.t[i] <- dr.K.stoch[[i]]$t
    dr.K.stoch.list <- list(dr = dr.K.stoch.dr, t = dr.K.stoch.t)
  }
  return(dr.K.stoch.list)
}

############################################
#Function to perform parameter perturbations
############################################
res_param_perturb <- function(plot,mu_LI_plot,nkernel){
  add.s <- seq(0.0,0.01,0.001)
  ord <- array(c(rep(1,11),rep(0,26*11),
                 rep(0,11),rep(1,11),rep(0,25*11),
                 rep(0,11*2),rep(1,11),rep(0,24*11),
                 rep(0,11*3),rep(1,11),rep(0,23*11),
                 rep(0,11*4),rep(1,11),rep(0,22*11),
                 rep(0,11*5),rep(1,11),rep(0,21*11),
                 rep(0,11*6),rep(1,11),rep(0,20*11),
                 rep(0,11*7),rep(1,11),rep(0,19*11),
                 rep(0,11*8),rep(1,11),rep(0,18*11),
                 rep(0,11*9),rep(1,11),rep(0,17*11),
                 rep(0,11*10),rep(1,11),rep(0,16*11),
                 rep(0,11*11),rep(1,11),rep(0,15*11),
                 rep(0,11*12),rep(1,11),rep(0,14*11),
                 rep(0,11*13),rep(1,11),rep(0,13*11),
                 rep(0,11*14),rep(1,11),rep(0,12*11),
                 rep(0,11*15),rep(1,11),rep(0,11*11),
                 rep(0,11*16),rep(1,11),rep(0,10*11),
                 rep(0,11*17),rep(1,11),rep(0,9*11),
                 rep(0,11*18),rep(1,11),rep(0,8*11),
                 rep(0,11*19),rep(1,11),rep(0,7*11),
                 rep(0,11*20),rep(1,11),rep(0,6*11),
                 rep(0,11*21),rep(1,11),rep(0,5*11),
                 rep(0,11*22),rep(1,11),rep(0,4*11),
                 rep(0,11*23),rep(1,11),rep(0,3*11),
                 rep(0,11*24),rep(1,11),rep(0,2*11),
                 rep(0,11*25),rep(1,11),rep(0,1*11),
                 rep(0,11*26),rep(1,11)),
               dim=c(11,27,27))
  
  res.sens <- list(fst.amp = array(NA,dim=c(170,11,27)), fst.att = array(NA,dim=c(170,11,27)), recov.t = array(NA,dim=c(170,11,27)))
  
  for(i in 1:27){
    for(j in 1:11){
      
      
      fixed_list <- list(
        s_mu_slope   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='beta.phi'] + add.s[j]*ord[j,12,i],    #survival slope
        
        #Environmental slopes for survival
        s_mu_tmed2m  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[1]'] + add.s[j]*ord[j,2,i],
        s_mu_RHmax   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[2]']+ add.s[j] *ord[j,3,i],
        s_mu_sol     = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[3]']+ add.s[j] *ord[j,4,i],
        s_mu_tmed0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[4]']+ add.s[j]*ord[j,5,i],
        s_mu_tmin0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[5]']+ add.s[j]*ord[j,6,i],
        s_mu_precip  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[6]']+ add.s[j]*ord[j,7,i],
        s_mu_perf    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[7]']+ add.s[j]*ord[j,8,i],
        s_mu_ha_90   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[8]']+ add.s[j]*ord[j,9,i],
        s_mu_fire    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[9]']+ add.s[j]*ord[j,10,i],
        s_mu_TSLF    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaphiJS[10]']+ add.s[j]*ord[j,11,i],

        
        #sigma.phiJS = vitalrates.Matticolus.df['sigma.phiJS','mean'],
        
        #Environmental slopes for reproduction
        r_f_mu_tmed2m  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[1]']+ add.s[j]*ord[j,13,i],
        r_f_mu_RHmax   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[2]']+ add.s[j]*ord[j,14,i],
        r_f_mu_sol     = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[3]']+ add.s[j]*ord[j,15,i],
        r_f_mu_tmed0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[4]']+ add.s[j]*ord[j,16,i],
        r_f_mu_tmin0cm = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[5]']+ add.s[j]*ord[j,17,i],
        r_f_mu_precip  = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[6]']+ add.s[j]*ord[j,18,i],
        r_f_mu_perf    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[7]']+ add.s[j]*ord[j,19,i],
        r_f_mu_ha_90   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[8]']+ add.s[j]*ord[j,20,i],
        r_f_mu_fire    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[9]']+ add.s[j]*ord[j,21,i],
        r_f_mu_TSLF    = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='betaf[10]']+ add.s[j]*ord[j,22,i],

        
        #sigma.f = vitalrates.Matticolus.df['sigma.f','mean'],
        
        
        #Probability of reproduction
        r_r_mu_int   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.prep']+ add.s[j]*ord[j,23,i],
        r_r_mu_slope = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='beta1.prep']+ add.s[j]*ord[j,24,i],  
        # r_r_mu_slope2 = vitalrates.Matticolus.df['beta2.prep','mean']+ add.s[j]*ord[j,26,i], 
        # r_r_sd_slope2 = vitalrates.Matticolus.df['beta2.prep','sd'], 
        
        #Number of eggs/embryos
        # r_n_mu_int   = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.fec']+ add.s[j]*ord[j,25,i], 
        # r_n_mu_slope = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='beta1.fec']+ add.s[j]*ord[j,26,i],   
        #r_n_mu_slope2 =vitalrates.Matticolus.df['beta2.fec','mean'],  
        
        #r_n_sd_slope2 =vitalrates.Matticolus.df['beta2.fec','sd'],  
        
        
        #Size of newborns
        mu_rd     = Matticolus.data$mu.L0,   
        sd_rd     = sqrt(Matticolus.data$tau.L0),
        mu_LI = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.LI']
        
      )
      
      
      # Now, simulate some random intercepts for growth (g_), survival (s_),
      # and offspring production (r_s_). This part is for the purpose of the example.
      
      # First, we create vector of values that each random component can take.
      s_params  <- list(
        s_g_mu_int_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[1]']+ add.s[j]*ord[j,25,i],
        s_g_mu_int_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[2]']+ add.s[j]*ord[j,25,i],
        s_g_mu_int_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[3]']+ add.s[j]*ord[j,25,i],
        s_g_mu_int_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[4]']+ add.s[j]*ord[j,25,i],
        s_g_mu_int_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.phiJS[5]']+ add.s[j]*ord[j,25,i]

      )
      
      g_params <- list(
        
        g_g_mu_K_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[1]']+ add.s[j]*ord[j,26,i],
        g_g_mu_K_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[2]']+ add.s[j]*ord[j,26,i],
        g_g_mu_K_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[3]']+ add.s[j]*ord[j,26,i],
        g_g_mu_K_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[4]']+ add.s[j]*ord[j,26,i],
        g_g_mu_K_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='mu.K[5]']+ add.s[j]*ord[j,26,i]
        
      )
      
      r_params <- list(
        r_f_mu_int_1 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[1]']+ add.s[j]*ord[j,27,i],
        r_f_mu_int_2 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[2]']+ add.s[j]*ord[j,27,i],
        r_f_mu_int_3 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[3]']+ add.s[j]*ord[j,27,i],
        r_f_mu_int_4 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[4]']+ add.s[j]*ord[j,27,i],
        r_f_mu_int_5 = vitalrates.Matticolus.df$mean[vitalrates.Matticolus.df$X=='alpha.f[5]']+ add.s[j]*ord[j,27,i]
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
          par_set_indices  = list(site = plot),
          
          # We must also index the variables in the eviction function
          
          evict_cor        = TRUE,
          evict_fun        = truncated_distributions("norm", "g_site")
          
        ) %>%
        define_kernel(
          
          # The F kernel also varies from site to site
          
          name             = "F_site",
          formula          = ((1-(surv_site/(surv_site+r_f_site)))+r_r_site) * 2 * r_d,
          family           = "CC",
          
          # We didn't include a site level effect for probability
          # of reproduction. Thus, this expression is NOT indexed.
          
          r_r_lin          = (rnorm(1,r_r_mu_int,0) + 
                                rnorm(1, r_r_mu_slope, 0) * ht_1),
          r_r_site              = inv_logit(r_r_lin),
          
          # We index the seed production expression with the site effect
          # 
          # r_n_lin_site          = (rnorm(1,r_n_mu_int, 0) +
          #                            r_f_lin_site +
          #                            rnorm(1,r_n_mu_slope,0) * ht_1 +
          #                            rnorm(1,r_n_mu_slope2,0)* ht_1^2),
          # 2         = pois_r(r_n_lin_site),
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
                                                   dim = c(60,60,170)), 
                                             MARGIN = c(3),FUN =  reac, 
                                             bound = "upper", 
                                             simplify = F))
      res.sens$fst.att[,j,i] <- unlist(apply(array(unlist(stoch.K),
                                                   dim = c(60,60,170)), 
                                             MARGIN = c(3),
                                             FUN =  reac, 
                                             bound = "lower", 
                                             simplify = F))
      res.sens$recov.t[,j,i] <- dr.stoch(stoch.K, 170)$t
      
    }
  }
  return(res.sens)
}

res.Matticolus.param.sens.C <- res_param_perturb(1, all_params_list$mu_LI, 60)
res.Matticolus.param.sens.Q <- res_param_perturb(2, all_params_list$mu_LI, 60)
res.Matticolus.param.sens.EB <- res_param_perturb(3, all_params_list$mu_LI, 60)
res.Matticolus.param.sens.MB <- res_param_perturb(4, all_params_list$mu_LI, 60)
res.Matticolus.param.sens.LB <- res_param_perturb(5, all_params_list$mu_LI, 60)

saveRDS(res.Matticolus.param.sens.C,"res.param.sens.C.Matticolus.rds")
saveRDS(res.Matticolus.param.sens.Q,"res.param.sens.Q.Matticolus.rds")
saveRDS(res.Matticolus.param.sens.EB,"res.param.sens.EB.Matticolus.rds")
saveRDS(res.Matticolus.param.sens.MB,"res.param.sens.MB.Matticolus.rds")
saveRDS(res.Matticolus.param.sens.LB,"res.param.sens.LB.Matticolus.rds")

res.Matticolus.param.sens.C  <- readRDS("res.param.sens.C.Matticolus.rds")
res.Matticolus.param.sens.Q  <- readRDS("res.param.sens.Q.Matticolus.rds")
res.Matticolus.param.sens.EB <- readRDS("res.param.sens.EB.Matticolus.rds")
res.Matticolus.param.sens.MB <- readRDS("res.param.sens.MB.Matticolus.rds")
res.Matticolus.param.sens.LB <- readRDS("res.param.sens.LB.Matticolus.rds")


sens.res <- function(res.array){
  add.s <- seq(0,0.01,0.001)
  res.sens <- array(NA, dim = c(169,10,27))
  for(i in 1:169){
    for(j in 2:11){
      for(k in 1:27){
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
                    X27 = c(res.sens[,,27])
                    )
         )
}

(sens.Matticolus.fst.amp.C <- sens.res(log10(res.Matticolus.param.sens.C$fst.amp)))
(sens.Matticolus.fst.amp.Q <- sens.res(log10(res.Matticolus.param.sens.Q$fst.amp)))
(sens.Matticolus.fst.amp.EB <- sens.res(log10(res.Matticolus.param.sens.EB$fst.amp)))
(sens.Matticolus.fst.amp.MB <- sens.res(log10(res.Matticolus.param.sens.MB$fst.amp)))
(sens.Matticolus.fst.amp.LB <- sens.res(log10(res.Matticolus.param.sens.LB$fst.amp)))

(sens.Matticolus.fst.att.C <- sens.res(res.Matticolus.param.sens.C$fst.att))
(sens.Matticolus.fst.att.Q <- sens.res(res.Matticolus.param.sens.Q$fst.att))
(sens.Matticolus.fst.att.EB <- sens.res(res.Matticolus.param.sens.EB$fst.att))
(sens.Matticolus.fst.att.MB <- sens.res(res.Matticolus.param.sens.MB$fst.att))
(sens.Matticolus.fst.att.LB <- sens.res(res.Matticolus.param.sens.LB$fst.att))

(sens.Matticolus.recov.t.C <- sens.res(res.Matticolus.param.sens.C$recov.t))
(sens.Matticolus.recov.t.Q <- sens.res(res.Matticolus.param.sens.Q$recov.t))
(sens.Matticolus.recov.t.EB <- sens.res(res.Matticolus.param.sens.EB$recov.t))
(sens.Matticolus.recov.t.MB <- sens.res(res.Matticolus.param.sens.MB$recov.t))
(sens.Matticolus.recov.t.LB <- sens.res(res.Matticolus.param.sens.LB$recov.t))

#Summary statistics
library(psych)
sens.Matticolus.fst.amp.summary <- print(describe(rbind(sens.Matticolus.fst.amp.C, 
                                             sens.Matticolus.fst.amp.Q, 
                                             sens.Matticolus.fst.amp.EB,
                                             sens.Matticolus.fst.amp.MB,
                                             sens.Matticolus.fst.amp.LB),
                                            quant = c(0.25, 0.75)),3)


write.csv(sens.Matticolus.fst.amp.summary , "sens_Matticolus_fst_amp_summary.csv")

sens.Matticolus.fst.att.summary <- print(describe(rbind(sens.Matticolus.fst.att.C, 
                                             sens.Matticolus.fst.att.Q, 
                                             sens.Matticolus.fst.att.EB,
                                             sens.Matticolus.fst.att.MB,
                                             sens.Matticolus.fst.att.LB),
                                            quant = c(0.25, 0.75)),3)


write.csv(sens.Matticolus.fst.att.summary , "sens_Matticolus_fst_att_summary.csv")

sens.Matticolus.recov.t.summary <- print(describe(rbind(sens.Matticolus.recov.t.C, 
                                             sens.Matticolus.recov.t.Q, 
                                             sens.Matticolus.recov.t.EB,
                                             sens.Matticolus.recov.t.MB,
                                             sens.Matticolus.recov.t.LB),
                                             quant = c(0.25, 0.75)),3)


write.csv(sens.Matticolus.recov.t.summary , "sens_Matticolus_recov_t_summary.csv")

#Plots
#Resilience x Perturbation magnitude
fst.amp.Matticolus.param.sens.mean <- (res.Matticolus.param.sens.C$fst.amp+
                                     res.Matticolus.param.sens.Q$fst.amp +
                                     res.Matticolus.param.sens.EB$fst.amp+
                                     res.Matticolus.param.sens.MB$fst.amp+
                                     res.Matticolus.param.sens.LB$fst.amp)/5


(fst.amp.Matticolus.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              log10(apply(fst.amp.Matticolus.param.sens.mean[,,c(which(sens.Matticolus.fst.amp.summary$mean!=0))], 
                                                          MARGIN = c(2,3), mean))))

fst.amp.Matticolus.param.sens.mean <- pivot_longer(data = fst.amp.Matticolus.param.sens.mean, 
                                               cols = `X1`:`X3`,
                                               names_to = "parameter",
                                               values_to = "compensation")
quartz(height = 6, width = 8)
ggplot(fst.amp.Matticolus.param.sens.mean,
       aes(x = perturb, y = compensation, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(3),name = 'Parameter', 
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
                                 # expression(alpha["nb"]),
                                 # expression(beta["nb"]*"SVL"),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Matticolus.fst.amp.summary$mean!=0))-1])+
  labs(x = "Perturbation magnitude", y = "Compensation")

fst.att.Matticolus.param.sens.mean <- (res.Matticolus.param.sens.C$fst.att+
                                     res.Matticolus.param.sens.Q$fst.att +
                                     res.Matticolus.param.sens.EB$fst.att+
                                     res.Matticolus.param.sens.MB$fst.att+
                                     res.Matticolus.param.sens.LB$fst.att)/5

(fst.att.Matticolus.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              apply(fst.att.Matticolus.param.sens.mean[,,c(which(sens.Matticolus.fst.att.summary$mean!=0))], 
                                                    MARGIN = c(2,3), mean)))

fst.att.Matticolus.param.sens.mean <- pivot_longer(data = fst.att.Matticolus.param.sens.mean, 
                                               cols = `X1`:`X3`,
                                               names_to = "parameter",
                                               values_to = "resistance")
quartz(height = 6, width = 8)
ggplot(fst.att.Matticolus.param.sens.mean,
       aes(x = perturb, y = resistance, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(3),name = 'Parameter', 
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
                                 # expression(alpha["nb"]),
                                 # expression(beta["nb"]*"SVL"),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Matticolus.fst.att.summary$mean!=0))-1])+
  labs(x = "Perturbation magnitude", y = "Resistance")

recov.t.Matticolus.param.sens.mean <- (res.Matticolus.param.sens.C$recov.t+
                                     res.Matticolus.param.sens.Q$recov.t +
                                     res.Matticolus.param.sens.EB$recov.t+
                                     res.Matticolus.param.sens.MB$recov.t+
                                     res.Matticolus.param.sens.LB$recov.t)/5


(recov.t.Matticolus.param.sens.mean <- data.frame(perturb = seq(0,0.01,0.001),
                                              apply(recov.t.Matticolus.param.sens.mean[,,c(which(sens.Matticolus.recov.t.summary$mean!=0))], 
                                                    MARGIN = c(2,3), mean)))

recov.t.Matticolus.param.sens.mean <- pivot_longer(data = recov.t.Matticolus.param.sens.mean, 
                                               cols = `X1`:`X3`,
                                               names_to = "parameter",
                                               values_to = "recov.t")
quartz(height = 6, width = 8)
ggplot(recov.t.Matticolus.param.sens.mean,
       aes(x = perturb, y = recov.t, colour = parameter))+
  geom_line(linewidth = 1.5, alpha = 0.5)+
  scale_colour_manual(values = turbo(3),name = 'Parameter', 
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
                                 # expression(alpha["nb"]),
                                 # expression(beta["nb"]*"SVL"),
                                 expression(alpha[phi]),
                                 expression(mu["K"]),
                                 expression(alpha["f"]))[c(which(sens.Matticolus.recov.t.summary$mean!=0))-1])+
  labs(x = "Perturbation magnitude", y = "Recovery time")


#Barplots
quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.amp.C[, -1],las=2,ylab="Compensation",main=expression("C - "*italic("M. atticolus")),
        ylim = c(0,2),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
           col = c(rep("brown",11),
                   rep("darkcyan",12),
                   rep("brown",2),
                   "darkcyan"),
           border = c(rep("brown",11),
                     rep("darkcyan",12),
                     rep("brown",2),
                     "darkcyan"),
                     notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.amp.Q[, -1],las=2,ylab="Compensation",main=expression("Q - "*italic("M. atticolus")),
        ylim = c(0,1.5),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.amp.EB[, -1],las=2,ylab="Compensation",main=expression("EB - "*italic("M. atticolus")),
        ylim = c(0,1.2),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.amp.MB[, -1],las=2,ylab="Compensation",main=expression("MB - "*italic("M. atticolus")),
        ylim = c(0,1.2),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.amp.LB[, -1],las=2,ylab="Compensation",main=expression("LB - "*italic("M. atticolus")),
        ylim = c(0,2),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.att.C[, -1],las=2,ylab="Resistance",main=expression("C - "*italic("M. atticolus")),
        ylim = c(-.7,0),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.att.Q[, -1],las=2,ylab="Resistance",main=expression("Q - "*italic("M. atticolus")),
        ylim = c(-.8,0),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.att.EB[, -1],las=2,ylab="Resistance",main=expression("EB - "*italic("M. atticolus")),
        ylim = c(-1,0),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.att.MB[, -1],las=2,ylab="Resistance",main=expression("MB - "*italic("M. atticolus")),
        ylim = c(-1.2,0),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.fst.att.LB[, -1],las=2,ylab="Resistance",main=expression("LB - "*italic("M. atticolus")),
        ylim = c(-1.1,0),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.recov.t.C[, -1],las=2,ylab="Recovery time",main=expression("C - "*italic("M. atticolus")),
        ylim = c(-25,2),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.recov.t.Q[, -1],las=2,ylab="Recovery time",main=expression("Q - "*italic("M. atticolus")),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.recov.t.EB[, -1],las=2,ylab="Recovery time",main=expression("EB - "*italic("M. atticolus")),
        ylim = c(-60,10),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.recov.t.MB[, -1],las=2,ylab="Recovery time",main=expression("MB - "*italic("M. atticolus")),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(sens.Matticolus.recov.t.LB[, -1],las=2,ylab="Recovery time",main=expression("LB - "*italic("M. atticolus")),
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
                         # expression(alpha["nb"]),
                         # expression(beta["nb"]*"SVL"),
                         expression(alpha[phi]),
                         expression(mu["K"]),
                         expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))


#Average among plots
quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Matticolus.fst.amp.C[, -1],
              sens.Matticolus.fst.amp.Q[, -1],
              sens.Matticolus.fst.amp.EB[, -1],
              sens.Matticolus.fst.amp.MB[, -1],
              sens.Matticolus.fst.amp.LB[, -1]),
        las=2,ylab="Compensation",main=expression("Average - "*italic("M. atticolus")),
        ylim = c(0,1.5),
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
                  # expression(alpha["nb"]),
                  # expression(beta["nb"]*"SVL"),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("topleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Matticolus.fst.att.C[, -1],
              sens.Matticolus.fst.att.Q[, -1],
              sens.Matticolus.fst.att.EB[, -1],
              sens.Matticolus.fst.att.MB[, -1],
              sens.Matticolus.fst.att.LB[, -1]),
        las=2,ylab="Resistance",main=expression("Average - "*italic("M. atticolus")),
        ylim = c(-1.2,0),
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
                  # expression(alpha["nb"]),
                  # expression(beta["nb"]*"SVL"),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))

quartz(8,12)
par(mar=c(7,4,1,.5))
boxplot(rbind(sens.Matticolus.recov.t.C[, -1],
              sens.Matticolus.recov.t.Q[, -1],
              sens.Matticolus.recov.t.EB[, -1],
              sens.Matticolus.recov.t.MB[, -1],
              sens.Matticolus.recov.t.LB[, -1]),las=2,ylab="Recovery time",main=expression("Average - "*italic("M. atticolus")),
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
                  # expression(alpha["nb"]),
                  # expression(beta["nb"]*"SVL"),
                  expression(alpha[phi]),
                  expression(mu["K"]),
                  expression(alpha["f"])),
        col = c(rep("brown",11),
                rep("darkcyan",12),
                rep("brown",2),
                "darkcyan"),
                border = c(rep("brown",11),
                           rep("darkcyan",12),
                           rep("brown",2),
                           "darkcyan"),
                           notch = T)

legend("bottomleft",legend = c("P kernel","F kernel"), fill = c("brown","darkcyan"))


