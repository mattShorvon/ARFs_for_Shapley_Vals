library(MASS)
library(data.table)
library(arf)
library(doParallel)

# make sure your working directory is set to wherever you have saved the
# ARFs_for_Shapley_Vals repository folder, and change your working directory
# with setwd() if necessary
getwd()

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
source("VectorisedMAP5.R")

# load the c++ functions
library(Rcpp)
library(RcppArmadillo)
sourceCpp("shapr_original/src/AICc.cpp")
sourceCpp("shapr_original/src/distance.cpp")
sourceCpp("shapr_original/src/features.cpp")
sourceCpp("shapr_original/src/impute_data.cpp")
sourceCpp("shapr_original/src/weighted_matrix.cpp")

# creating synthetic train and test datasets:
No_cont_var <- 2
No_cat_var <- 2
No_levels <- 4
No_tot_var <- No_cont_var + No_cat_var
Sigma_diag <- 1
corr <- 0.5
mu <- rep(0, No_tot_var)
No_train_obs <- 100
No_test_obs <- 50
Sigma <- matrix(corr, No_tot_var, No_tot_var)
diag(Sigma) <- Sigma_diag
seed <- 123
set.seed(seed)
x_train <- mvrnorm(n =  No_train_obs, mu = mu, Sigma = Sigma)
dt_train <- data.table(x_train[, 1:No_cont_var])
cat_cutoff <- c(-200, -0.5, 0, 1, 200)
for (i in (No_cont_var + 1):No_tot_var){
  dt_train <- cbind(dt_train, cut(x_train[, i], cat_cutoff, labels = 1:No_levels))
}
names(dt_train) <- c(paste0("cont_", 1:No_cont_var, "_"), paste0("cat_", 1:No_cat_var, "_"))
dt_train[, epsilon := rnorm(No_train_obs, 0, 0.1^2)]
x_test <- mvrnorm(n =  No_test_obs, mu = mu, Sigma = Sigma)
dt_test <- data.table(x_test[, 1:No_cont_var])
for (i in (No_cont_var + 1):No_tot_var){
  dt_test <- cbind(dt_test, cut(x_test[, i], cat_cutoff, labels = 1:No_levels))
}
names(dt_test) <- c(paste0("cont_", 1:No_cont_var, "_"), paste0("cat_", 1:No_cat_var, "_"))
dt_test[, epsilon := rnorm(No_test_obs, 0, 0.1^2)]
dt <- rbind(dt_train, dt_test)
cont_cols <- names(dt)[grep("cont", names(dt))]
cat_cols <- names(dt)[grep("cat", names(dt))]
feat_names <- c(cont_cols, cat_cols)
mod_matrix <- model.matrix(~.-1, data = dt[, ..feat_names], contrasts.arg = lapply(dt[, ..cat_cols], contrasts, contrasts = FALSE))
dt <- cbind(dt, data.table(mod_matrix[, -(1:No_cont_var)]))

# creating a simple linear model to be explained:
beta_0 <- 1
beta_cont <- c(1, -1)
beta_cat <- c(1, 0, -1, 0.5,
              2, 3, -1, -0.5)
beta <- c(beta_0, beta_cont, beta_cat)
train_obs <- 1:No_train_obs
test_obs <- (No_train_obs + 1):(No_train_obs + No_test_obs)
dt[, response := as.vector(cbind(1, mod_matrix) %*% beta) + epsilon]
fmla <- as.formula(paste0("response~", paste(feat_names, collapse = "+")))
model <- lm(formula = fmla, data = dt[train_obs, ])

# create the explainer object using shapr and the linear model
x_train <- dt[train_obs, ..feat_names]
x_test <- dt[test_obs, ..feat_names]
y_train <- dt[train_obs, .(response)]
explainer <- shapr(x_train, model)

# train an ARF to get the psi object
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
registerDoParallel(cl)
arf <- adversarial_rf(x_train, num_trees = 100)
psi <- forde(arf, x_train)
arf_shaps <- explain(x_test,explainer = explainer, approach = "arf", prediction_zero = mean(y_train), psi = psi)











