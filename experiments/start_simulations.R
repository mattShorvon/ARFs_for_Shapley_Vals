
# make sure environment is initialised properly with functions from init.R
library(data.table)
library(MASS)
library(ggplot2)
library(arf)
library(party)
library(partykit)
library(doParallel)

source("experiments/3-calculate_true_shapley_withdatatable.R")
source("experiments/source_mixed_data.R")

## Start simulation study
tod_date0 <- format(Sys.Date(), "%d_%m_%y")

dim <- 4
no_categories <- 4

clock_seed_0 <- round(as.numeric(Sys.time()) * 1000)
clock_seed <- signif(clock_seed_0) - clock_seed_0
set.seed(clock_seed)
rand_string <- stringi::stri_rand_strings(1,5)
folder <- paste0(tod_date0, "_", rand_string, "_dim", dim, "_nbcat", no_categories)

dir.create(paste("experiments/experiment_data", folder, sep = ""))

##

parameters_list <- list()

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
registerDoParallel(cl)

corr <- c(0, 0.1, 0.3, 0.5, 0.8, 0.9)
k <- 1
for(j in corr){
  parameters_list[[k]] <- list(methods = c("kernelSHAP", "ctree", "arf"),
                          No_sample_gaussian = c(10),
                          No_cont_var = 2,
                          No_cat_var = 2,
                          No_levels = 8,
                          Sigma_diag = 1,
                          corr = j,
                          No_train_obs = 1000,
                          No_test_obs = 500,
                          cat_cutoff = c(-200, -100, -50, -0.5, 0, 1, 50, 100, 200),
                          noise = FALSE,
                          name = 'testing',
                          seed = 123,
                          No_mc_cores = 1,
                          beta_0 = 1,
                          beta_cont = c(1, -1),
                          beta_cat = c(1, 0, -1, 0.5, 2, 3, -1, 1.5,
                                        2, 3, -1, -0.5, 2, 3, -1, 1.5))
  k <- k + 1
}
# parameters_list <- parameters_list[[1]] # this is for working out what the code does
all_methods <- list()
for(i in 1:length(parameters_list)){
  all_methods[[i]] <- compute_shapley_mixed_data(parameters_list[[i]])

  nm = paste(folder, '_rho_', parameters_list[[i]]$corr, ".rds", sep = "")

  saveRDS(all_methods, file = paste("C:/Users/mshor/OneDrive - King's College London/PhD/ShapleyValExptData", folder, nm, sep = "/"))
}

for(i in 1:length(parameters_list)){
  all_methods[[i]] <- compute_shapley_mixed_data(parameters_list[[i]])

  nm = paste(folder, '_rho_', parameters_list[[i]]$corr, ".rds", sep = "")

  saveRDS(all_methods, file = paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData", folder, nm, sep = "/"))
}

# running several experiments to get error bars:

# doing continuous and categorical features:
source("C:/Users/w21113599/OneDrive - King's College London/PhD/shapr-ctree_paper/inst/paper_simulations/source_mixed_data.R")

no_experiments <- 10
for (i in 1:no_experiments) {
  tod_date0 <- format(Sys.Date(), "%d_%m_%y")

  dim <- 4
  no_categories <- 4

  clock_seed_0 <- round(as.numeric(Sys.time()) * 1000)
  clock_seed <- signif(clock_seed_0) - clock_seed_0
  set.seed(clock_seed)
  rand_string <- stringi::stri_rand_strings(1,5)
  folder <- paste0(tod_date0, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories)

  dir.create(paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData/dim4_nbcat2cont2_n1000/", folder, sep = ""))

  parameters_list <- list()

  corr <- c(0, 0.1, 0.3, 0.5, 0.8, 0.9)
  k <- 1
  for(j in corr){
    parameters_list[[k]] <- list(methods = c("kernelSHAP", "ctree", "arf"),
                                 No_sample_gaussian = c(10),
                                 No_cont_var = 2,
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
                                 No_mc_cores = 1)
    k <- k + 1
  }
  all_methods <- list()
  for(i in 1:length(parameters_list)){
    all_methods[[i]] <- compute_shapley_mixed_data(parameters_list[[i]])

    nm = paste(folder, '_rho_', parameters_list[[i]]$corr, ".rds", sep = "")

    saveRDS(all_methods, file = paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData/dim4_nbcat2cont2_n1000", folder, nm, sep = "/"))
  }
}


# doing continuous features only:
source("C:/Users/w21113599/OneDrive - King's College London/PhD/shapr-ctree_paper/inst/paper_simulations/source_mixedorcont_data.R")

no_experiments <- 10
for (i in 1:no_experiments) {
  tod_date0 <- format(Sys.Date(), "%d_%m_%y")

  dim <- 3
  no_categories <- 0

  clock_seed_0 <- round(as.numeric(Sys.time()) * 1000)
  clock_seed <- signif(clock_seed_0) - clock_seed_0
  set.seed(clock_seed)
  rand_string <- stringi::stri_rand_strings(1,5)
  folder <- paste0(tod_date0, "_", rand_string, "_dim", dim, "_nbcont3", no_categories)

  dir.create(paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData/dim3_nbcont3_n1000/", folder, sep = ""))

  parameters_list <- list()

  corr <- c(0, 0.1, 0.3, 0.5, 0.8, 0.9)
  k <- 1
  for(j in corr){
    parameters_list[[k]] <- list(methods = c("gaussian","empirical", "kernelSHAP", "ctree", "arf"),
                                 No_sample_gaussian = c(10),
                                 No_cont_var = 3,
                                 No_cat_var = 0,
                                 No_levels = 4,
                                 Sigma_diag = 1,
                                 corr = j,
                                 No_train_obs = 1000,
                                 No_test_obs = 500,
                                 cat_cutoff = c(-200, -0.5, 0, 1, 200),
                                 noise = TRUE,
                                 name = 'testing',
                                 seed = clock_seed,
                                 No_mc_cores = 1)
    k <- k + 1
  }
  all_methods <- list()
  for(i in 1:length(parameters_list)){
    all_methods[[i]] <- compute_shapley_mixedorcont_data(parameters_list[[i]])

    nm = paste(folder, '_rho_', parameters_list[[i]]$corr, ".rds", sep = "")

    saveRDS(all_methods, file = paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData/dim3_nbcont3_n1000/", folder, nm, sep = "/"))
  }
}
