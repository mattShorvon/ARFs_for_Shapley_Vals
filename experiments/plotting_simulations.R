# library(shapr)
library(data.table)
library(MASS)
library(lqmm) ## to check if Sigma is positive definite
library(rapportools) # for testing booleans
library(ggplot2)
library(stringr)

source("C:/Users/mshor/OneDrive - King's College London/PhD/ARFs_for_Shapley_Vals/experiments/3-calculate_true_shapley_withdatatable.R")
source("C:/Users/w21113599/OneDrive - King's College London/PhD/shapr-ctree_paper/inst/paper_simulations/3-calculate_true_shapley_withdatatable.R")


tod_date <- "04_07_23"
rand_string <- "5U05Z"
dim <- 4
no_categories <- 4

folder <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories)

##

## load data
nm <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories, "_rho_0.9.rds")
all_methods <- readRDS(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ARFs_for_Shapley_Vals/experiments/experiment_data/dim4_nbcat2cont2_n1000", folder, nm, sep = "/"))


## for gaussian stuff
if(any(grepl("gaussian", names(all_methods[[1]]$methods)))){
  gauss_mat <- list()
  for(i in 1:length(all_methods)){
    gauss_mat[[i]] <- all_methods[[i]]$methods$gaussian_nsamples100$dt_sum[, id := rep(1:dim^no_categories, each = 50)]
    gauss_mat[[i]] <- gauss_mat[[i]][, .(none = mean(none), feat_1_ = mean(feat_1_), feat_2_ = mean(feat_2_), feat_3_ = mean(feat_3_)), by = .(id)]
    gauss_mat[[i]] = gauss_mat[[i]][, id := NULL]
  }

  gauss_mat2 <- list()
  for(i in 1:length(all_methods)){
    gauss_mat2[[i]] <- all_methods[[i]]$methods$gaussian_nsamples1000$dt_sum[, id := rep(1:dim^no_categories, each = 50)]
    gauss_mat2[[i]] <- gauss_mat2[[i]][, .(none = mean(none), feat_1_ = mean(feat_1_), feat_2_ = mean(feat_2_), feat_3_ = mean(feat_3_)), by = .(id)]
    gauss_mat2[[i]] = gauss_mat2[[i]][, id := NULL]
  }
}


MAE_truth <- NULL
MAE_methods <- NULL
MAE_methods_names <- NULL
MAE_parameters <- NULL
MAE_seed <- NULL
N_test <- 500 # don't hard code?


MAE <- function(true_shapley, shapley_method, weights = rep(1 / nrow(true_shapley), nrow(true_shapley))){
  mean(apply(abs((true_shapley - shapley_method) * weights), 2, sum)[-1])
}

# if you need to remove gaussian10 trials
for (i in 1:length(all_methods)) {
  all_methods[[i]]$methods <- all_methods[[i]]$methods[-c(2,3,4,5,6,7,8,9,10)]
}


for(i in 1:length(all_methods)){
  if(!is.null(all_methods[[i]][['true_linear']])){
    MAE_truth <- c(MAE_truth, MAE(all_methods[[i]][['true_shapley']], all_methods[[i]][['true_linear']], N_test = N_test))
  }
  for(m in names(all_methods[[1]]$methods)){

    if(m == 'gaussian_nsamples100'){
      MAE_methods <- c(MAE_methods, MAE(all_methods[[i]][['true_shapley']], gauss_mat[[i]], N_test = N_test))
      MAE_methods_names <- c(MAE_methods_names, m)
      MAE_parameters <- c(MAE_parameters, all_methods[[i]]$parameters$name)
      MAE_seed <- c(MAE_seed, all_methods[[i]]$seed)
    } else if (m == 'gaussian_nsamples1000'){
      MAE_methods <- c(MAE_methods, MAE(all_methods[[i]][['true_shapley']], gauss_mat2[[i]], N_test = N_test))
      MAE_methods_names <- c(MAE_methods_names, m)
      MAE_parameters <- c(MAE_parameters, all_methods[[i]]$parameters$name)
      MAE_seed <- c(MAE_seed, all_methods[[i]]$seed)
    }
    else if(m == "empirical"){
      MAE_methods <- c(MAE_methods, MAE(all_methods[[i]][['true_shapley']], all_methods[[i]][['methods']][[m]]$dt_sum))
      MAE_methods_names <- c(MAE_methods_names, m)
      MAE_parameters <- c(MAE_parameters, all_methods[[i]]$parameters$corr)
      MAE_seed <- c(MAE_seed, all_methods[[i]]$parameters$seed)
    }
    else{
      # Mean_abserror <- mean(colSums(abs(all_methods[[i]][['true_shapley']] - all_methods[[i]][['methods']][[m]]$dt))[-1])/N_test
      MAE_methods <- c(MAE_methods, MAE(all_methods[[i]][['true_shapley']], all_methods[[i]][['methods']][[m]]$dt))
      MAE_methods_names <- c(MAE_methods_names, m)
      MAE_parameters <- c(MAE_parameters, all_methods[[i]]$parameters$corr)
      MAE_seed <- c(MAE_seed, all_methods[[i]]$parameters$seed)
    }
  }
}


results <- data.table(MAE_methods, MAE_methods_names, MAE_parameters, MAE_seed)
colnames(results)[3] <- "correlation"

nm = paste(tod_date, '_results', '.rds', sep = "")
folder <- paste0(tod_date, "_WbMlD", "_dim", dim, "_nbcat", no_categories)
saveRDS(results0, file = paste("/nr/project/stat/BigInsight/Projects/Fraud/Subprojects/NAV/Annabelle/results/paper_simulations_Annabelle", folder, nm, sep = "/"))


# p1 <- ggplot(data = results0, aes(y = MAE_methods, x = MAE_parameters, col = as.factor(MAE_methods_names ))) +
#   geom_point(size = 4, stroke = 1.5) +
#   scale_x_discrete(labels = c("corr0" = "0", "corr0.05" = "0.05", "corr0.1" = "0.1", "corr0.3" = "0.3", "corr0.5" = "0.5", "corr0.8" = "0.8", "corr0.9" = "0.9")) +
#   theme_bw(base_size = 22) + xlab("correlation") +
#   ylab("Mean absolute error (MAE)") +
#   scale_color_discrete(name = "Method" ) +
#   ggtitle("")
# #  labels = c("Ctree", "Ctree one-hot", "Empirical", "Empirical independence", "Gaussian_100", "Gaussian_1000", "kernelSHAP")

p1 <- ggplot(data = results, aes(y = MAE_methods, x = correlation, col = as.factor(MAE_methods_names))) +
  geom_point(size = 4, stroke = 1.5) +
  geom_line() +
  theme_bw(base_size = 22) +
  xlab("correlation") + ylab("Mean Absolute Error (MAE)") +
  scale_color_discrete(name = "Method" ) + ggtitle("")

plot(p1)

img = paste(tod_date, '_MAE', '.png', sep = "")
ggsave(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ShapleyValExptData", folder, img, sep = "/"), plot = p1, device = NULL, path = NULL,
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 300, limitsize = TRUE)

ggsave(paste("C:/Users/w21113599/OneDrive - King's College London/PhD/ShapleyValExptData", folder, img, sep = "/"), plot = p1, device = NULL, path = NULL,
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 300, limitsize = TRUE)

ggsave(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ARFs_for_Shapley_Vals/experiments/experiment_data/dim4_nbcat2cont2_n1000", folder, img, sep = "/"), plot = p1, device = NULL, path = NULL,
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 300, limitsize = TRUE)

# getting error bars:
rand_string_list <- c("7pl43","aoeYf","esgA8","flEXJ","GsKYu","iCk3y","UfyJX","ulX58","VjvuQ")
results_1 <- list()
results_2 <- list()
results_3 <- list()
results_4 <- list()
results_5 <- list()
results_6 <- list()
results_7 <- list()
results_8 <- list()
results_9 <- list()
results_10 <- list()
for (i in 1:length(rand_string_list)) {
  tod_date <- "06_07_23"
  rand_string <- rand_string_list[i]
  dim <- 4
  no_categories <- 4
  print(rand_string)
  folder <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories)
  nm <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories, "_rho_0.9.rds")
  all_methods <- readRDS(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ARFs_for_Shapley_Vals/experiments/experiment_data/dim4_nbcat2cont2_n1000", folder, nm, sep = "/"))
  MAE_truth <- NULL
  MAE_methods <- NULL
  MAE_methods_names <- NULL
  MAE_parameters <- NULL
  MAE_seed <- NULL
  for (j in 1:length(all_methods)) {
    all_methods[[j]]$methods <- all_methods[[j]]$methods
    if(!is.null(all_methods[[j]][['true_linear']])){
      MAE_truth <- c(MAE_truth, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['true_linear']], N_test = N_test))
    }
    for(m in names(all_methods[[1]]$methods)){

      if(m == 'gaussian_nsamples100'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      } else if (m == 'gaussian_nsamples1000'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat2[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      }
      else if(m == "empirical"){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt_sum))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$corr)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
      else{
        # Mean_abserror <- mean(colSums(abs(all_methods[[i]][['true_shapley']] - all_methods[[i]][['methods']][[m]]$dt))[-1])/N_test
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$corr)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
    }
  }
  results <- data.table(MAE_methods, MAE_methods_names, MAE_parameters, MAE_seed)
  colnames(results)[3] <- "correlation"
  assign(paste0("results_",i), results)
}
big_results <- rbind(results_1,results_2,results_3,results_4, results_5, results_6, results_7, results_8, results_9)
big_results <- na.omit(big_results)
big_results[, std_dev := apply(.SD,2,sd), by = c("correlation", "MAE_methods_names"), .SDcols = "MAE_methods"]
big_results[, avg := apply(.SD,2,mean), by = c("correlation", "MAE_methods_names"), .SDcols = "MAE_methods"]
final_results <- big_results[, !c("MAE_methods","MAE_seed")]
final_results <- unique(final_results)
plot <- ggplot(final_results,
               aes(x = correlation,y = avg,group = MAE_methods_names)) +
  geom_line(aes(color = MAE_methods_names)) +
  geom_point(aes(color = MAE_methods_names)) +
  geom_errorbar(aes(ymin=avg-std_dev, ymax=avg+std_dev, width = 0.03)) +
  ylab("Average MAE") +
  scale_colour_discrete(name = "Approach")

# img save doesn't work for error bars for some reason
img = paste(tod_date, 'dim4nbcat2cont2_toeplitz', '.png', sep = "")
ggsave(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ShapleyValExptData", folder, img, sep = "/"), plot = p1, device = NULL, path = NULL,
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 300, limitsize = TRUE)

# adding results from before leaf_posteriors were fixed

rand_string_list <- c("5o7yv","D3Wmd","dzQjw","EkFQQ","g4Rg5","iMUcv","ipgIj","kpXTZ","KTA0B","V2kxn")
results_1 <- list()
results_2 <- list()
results_3 <- list()
results_4 <- list()
results_5 <- list()
results_6 <- list()
results_7 <- list()
results_8 <- list()
results_9 <- list()
results_10 <- list()
for (i in 1:length(rand_string_list)) {
  tod_date <- "17_06_23"
  rand_string <- rand_string_list[i]
  dim <- 4
  no_categories <- 4
  print(rand_string)
  folder <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories)
  nm <- paste0(tod_date, "_", rand_string, "_dim", dim, "_nbcat2cont2_n1000_", no_categories, "_rho_0.9.rds")
  all_methods <- readRDS(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ShapleyValExptData/dim4_nbcat2cont2_n1000", folder, nm, sep = "/"))
  MAE_truth <- NULL
  MAE_methods <- NULL
  MAE_methods_names <- NULL
  MAE_parameters <- NULL
  MAE_seed <- NULL
  for (j in 1:length(all_methods)) {
    all_methods[[j]]$methods <- all_methods[[j]]$methods
    if(!is.null(all_methods[[j]][['true_linear']])){
      MAE_truth <- c(MAE_truth, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['true_linear']], N_test = N_test))
    }
    for(m in names(all_methods[[1]]$methods)){
      
      if(m == 'gaussian_nsamples100'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      } else if (m == 'gaussian_nsamples1000'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat2[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      }
      else if(m == "empirical"){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt_sum))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$corr)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
      else{
        # Mean_abserror <- mean(colSums(abs(all_methods[[i]][['true_shapley']] - all_methods[[i]][['methods']][[m]]$dt))[-1])/N_test
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$corr)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
    }
  }
  results <- data.table(MAE_methods, MAE_methods_names, MAE_parameters, MAE_seed)
  colnames(results)[3] <- "correlation"
  assign(paste0("results_",i), results)
}

old_results <- rbind(results_1,results_2,results_3,results_4, results_5, results_6, results_7, results_8, results_9)
old_results <- na.omit(old_results)
old_results <- old_results[MAE_methods_names == "arf"]
old_results[, MAE_methods_names := "arf (old)"]
old_results[, std_dev := apply(.SD,2,sd), by = c("correlation", "MAE_methods_names"), .SDcols = "MAE_methods"]
old_results[, avg := apply(.SD,2,mean), by = c("correlation", "MAE_methods_names"), .SDcols = "MAE_methods"]
new_n_old_results <- rbind(big_results, old_results)
final_results <- new_n_old_results[, !c("MAE_methods","MAE_seed")]
final_results <- unique(final_results)
plot <- ggplot(final_results,
               aes(x = correlation,y = avg,group = MAE_methods_names)) +
  geom_line(aes(color = MAE_methods_names)) +
  geom_point(aes(color = MAE_methods_names)) +
  geom_errorbar(aes(ymin=avg-std_dev, ymax=avg+std_dev, width = 0.03, color = MAE_methods_names)) +
  ylab("Average MAE") +
  scale_colour_discrete(name = "Approach")


# convergence results:
rand_string_list <- c("5BGuO","cPfZW","rTNVA","vhsEQ")
results_1 <- list()
results_2 <- list()
results_3 <- list()
results_4 <- list()
for (i in 1:length(rand_string_list)) {
  tod_date <- "20_06_23"
  rand_string <- rand_string_list[i]
  dim <- 4
  no_categories <- 4
  print(rand_string)
  folder <- paste0(tod_date, "_", rand_string, "_dim", dim, "nbcat2cont2_convergence", no_categories)
  nm <- paste0(tod_date, "_", rand_string, "_dim", dim, "nbcat2cont2_convergence", no_categories, "_rho_0.5.rds")
  all_methods <- readRDS(paste("C:/Users/mshor/OneDrive - King's College London/PhD/ShapleyValExptData/dim4_nbcat2cont2_convergence", folder, nm, sep = "/"))
  MAE_truth <- NULL
  MAE_methods <- NULL
  MAE_methods_names <- NULL
  MAE_parameters <- NULL
  MAE_seed <- NULL
  for (j in 1:length(all_methods)) {
    all_methods[[j]]$methods <- all_methods[[j]]$methods
    if(!is.null(all_methods[[j]][['true_linear']])){
      MAE_truth <- c(MAE_truth, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['true_linear']], N_test = N_test))
    }
    for(m in names(all_methods[[1]]$methods)){

      if(m == 'gaussian_nsamples100'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      } else if (m == 'gaussian_nsamples1000'){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], gauss_mat2[[j]], N_test = N_test))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$name)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$seed)
      }
      else if(m == "empirical"){
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt_sum))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$corr)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
      else{
        # Mean_abserror <- mean(colSums(abs(all_methods[[i]][['true_shapley']] - all_methods[[i]][['methods']][[m]]$dt))[-1])/N_test
        MAE_methods <- c(MAE_methods, MAE(all_methods[[j]][['true_shapley']], all_methods[[j]][['methods']][[m]]$dt))
        MAE_methods_names <- c(MAE_methods_names, m)
        MAE_parameters <- c(MAE_parameters, all_methods[[j]]$parameters$No_train_obs)
        MAE_seed <- c(MAE_seed, all_methods[[j]]$parameters$seed)
      }
    }
  }
  results <- data.table(MAE_methods, MAE_methods_names, MAE_parameters, MAE_seed)
  colnames(results)[3] <- "Number_of_training_datapoints"
  assign(paste0("results_",i), results)
}
big_results <- rbind(results_1,results_2,results_3,results_4)
big_results <- na.omit(big_results)
big_results[, std_dev := apply(.SD,2,sd), by = c("Number_of_training_datapoints", "MAE_methods_names"), .SDcols = "MAE_methods"]
big_results[, avg := apply(.SD,2,mean), by = c("Number_of_training_datapoints", "MAE_methods_names"), .SDcols = "MAE_methods"]
final_results <- big_results[, !c("MAE_methods","MAE_seed")]
final_results <- unique(final_results)
plot <- ggplot(final_results,
               aes(x = Number_of_training_datapoints,y = avg,group = MAE_methods_names)) +
  geom_line(aes(color = MAE_methods_names)) +
  geom_point(aes(color = MAE_methods_names)) +
  geom_errorbar(aes(ymin=avg-std_dev, ymax=avg+std_dev, width = 0.03, color = MAE_methods_names)) +
  ylab("Average MAE") +
  scale_colour_discrete(name = "Approach")
plot

