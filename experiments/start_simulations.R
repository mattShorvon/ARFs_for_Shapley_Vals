# init:
library(data.table)
library(MASS)
library(ggplot2)
# library(devtools)
# install_github("bips-hb/arf@tweaks") # doesn't work on cluster w/o sudo access
# library(arf) just installs the main branch which we don't want for this
library(party)
library(partykit)
library(doParallel)

print(getwd())

# load the functions from the shapr package:
source("shapr_original/R/explanation.R")
source("shapr_original/R/features.R")
source("shapr_original/R/models.R")
source("shapr_original/R/observations.R")
source("shapr_original/R/plot.R")
source("shapr_original/R/predictions.R")
source("shapr_original/R/preprocess_data.R")
source("shapr_original/R/RcppExports.R")
source("shapr_original/R/sampling.R")
source("shapr_original/R/shapley.R")
source("shapr_original/R/transformation.R")
source("shapr_original/R/utils.R")

# load the new ARF functions:
# library(ranger) # required if you're using source() to install arf
source("VectorisedMAP7.R")
source("arf_shap_sampling.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/adversarial_rf.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/forde.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/forge.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/utils.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/map.R")
# source("ARFs_for_Shapley_Vals/arf@tweaks/lik.R")
# install.packages("~/ARFs_for_Shapley_Vals/arf-tweaks", repos = NULL, type="source") run if arf not installed
library(arf)

# load the c++ functions
library(Rcpp)
library(RcppArmadillo)
sourceCpp("shapr_original/src/AICc.cpp")
sourceCpp("shapr_original/src/distance.cpp")
sourceCpp("shapr_original/src/features.cpp")
sourceCpp("shapr_original/src/impute_data.cpp")
sourceCpp("shapr_original/src/weighted_matrix.cpp")
# may throw warning messages, but they don't seem to affect things

source("experiments/3-calculate_true_shapley_withdatatable.R")

# doing continuous and categorical features:
source("experiments/source_mixed_cat2cont4.R")

directory <- "experiments/experiment_data/nbcat2cont4_n1000/"
dir.create(directory)
no_experiments <- 1

no_cores <- 6 # need to set to correct value
cl <- makeCluster(no_cores)
registerDoParallel(cl)

for (i in 1:no_experiments) {
  tm0 <- proc.time()
  tod_date0 <- format(Sys.Date(), "%d_%m_%y")
  
  dim <- 4
  no_categories <- 4
  
  clock_seed_0 <- round(as.numeric(Sys.time()) * 1000)
  clock_seed <- signif(clock_seed_0) - clock_seed_0
  set.seed(clock_seed)
  rand_string <- stringi::stri_rand_strings(1,5, pattern = "[A-HJKM-Za-hjkm-z]") # avoiding i and l
  folder <- paste0(tod_date0, "_", rand_string, "_nbcat2cont2_n1000_sampling")
  
  dir.create(paste(directory, folder, sep = ""))
  
  parameters_list <- list()
  
  corr <- c(0, 0.1, 0.3, 0.5, 0.8, 0.9)
  k <- 1
  for(j in corr){
    parameters_list[[k]] <- list(methods = c("kernelSHAP","ctree","arf","arf_sampling"),
                                 No_sample_gaussian = c(10),
                                 No_cont_var = 4,
                                 No_cat_var = 2,
                                 No_levels = 4,
                                 Sigma_diag = 1,
                                 corr = j,
                                 No_train_obs = 1000,
                                 No_test_obs = 500,
                                 cat_cutoff = c(-200, -0.5, 0, 1, 200),
                                 noise = TRUE,
                                 name = 'testing',
                                 seed = clock_seed,
                                 beta_0 = 1,
                                 beta_cont = c(1, -1, 1, -1),
                                 beta_cat = c(1, 0, -1, 0.5,
                                              2, 3, -1, -0.5),
                                 arf_alpha = 0,
                                 No_cores = no_cores)
    k <- k + 1
  }
  all_methods <- list()
  for(i in 1:length(parameters_list)){
    all_methods[[i]] <- compute_shapley_mixed_data(parameters_list[[i]])
    
    nm = paste(folder, '_rho_', parameters_list[[i]]$corr, ".rds", sep = "")
    
    saveRDS(all_methods, file = paste(directory, folder, nm, sep = "/"))
  }
  tm1 <- proc.time()
  print(tm1 - tm0)
}


