###### GOING THROUGH HOW GAUSSIAN WORKS TO SEE THE SAMPLING IN ACTION ######

library(xgboost)

data("airquality")
data <- data.table::as.data.table(airquality)
data <- data[complete.cases(data), ]

x_var <- c("Solar.R", "Wind", "Temp", "Month")
y_var <- "Ozone"

ind_x_explain <- 1:6
x_train <- data[-ind_x_explain, ..x_var]
y_train <- data[-ind_x_explain, get(y_var)]
x_test <- data[ind_x_explain, ..x_var]

model <- xgboost(
  data = as.matrix(x_train),
  label = y_train,
  nround = 20,
  verbose = FALSE
)

explainer <- shapr(x_train, model)

# setting inputs:
x <- x_test
approach <- "gaussian"
n_samples <- 1e3
n_batches <- 1
seed <- 1
prediction_zero <- mean(y_train)

################################# inside explain() in explanation.R


# Check input for x
if (!is.matrix(x) & !is.data.frame(x)) {
  stop("x should be a matrix or a data.frame/data.table.")
}

if (n_batches < 1 || n_batches > nrow(explainer$S)) {
  stop("`n_batches` is smaller than 1 or greater than the number of rows in explainer$S.")
}
# Check input for approach
if (!(is.vector(approach) &&
      is.atomic(approach) &&
      (length(approach) == 1 | length(approach) == length(explainer$feature_list$labels)) &&
      all(is.element(approach, c("empirical", "gaussian", "copula", "ctree", "independence", "arf"))))
) {
  stop(
    paste(
      "It seems that you passed a non-valid value for approach.",
      "It should be either 'empirical', 'gaussian', 'copula', 'ctree', 'independence' or",
      "a vector of length=ncol(x) with only the above characters."
    )
  )
}

# adding this here, cos i think its been removed from current cran version of shapr, which will cause
# it not to work with the github version functions
explainer$is_groupwise <- FALSE
print(!explainer$is_groupwise)

################################# inside explain.gaussian in explanation.R

if (!is.null(seed)) set.seed(seed)

# Add arguments to explainer object
explainer$x_test <- as.matrix(preprocess_data(x, explainer$feature_list)$x_dt)
explainer$approach <- approach
explainer$n_samples <- n_samples


# If mu is not provided directly, use mean of training data
explainer$mu <- unname(colMeans(explainer$x_train))

# If cov_mat is not provided directly, use sample covariance of training data
cov_mat = NULL
if (is.null(cov_mat)) {
  cov_mat <- stats::cov(explainer$x_train)
}

# Make sure that covariance matrix is positive-definite
eigen_values <- eigen(cov_mat)$values
if (any(eigen_values <= 1e-06)) {
  explainer$cov_mat <- as.matrix(Matrix::nearPD(cov_mat)$mat)
} else {
  explainer$cov_mat <- cov_mat
}

################################# inside prepare_and_predict:

# For R CMD check
row_id <- NULL

index_S <- NULL
only_return_contrib_dt <- NULL
if(is.null(only_return_contrib_dt)) only_return_contrib_dt <- FALSE

S_batch <- create_S_batch(explainer, n_batches, index_S)
pred_batch <- list()
r_batch <- list()
p <- NA

################################# inside prepare_data.gaussian:

id <- id_combination <- w <- NULL # due to NSE notes in R CMD check

n_xtest <- nrow(explainer$x_test)
dt_l <- list()
index_features = NULL

if (is.null(index_features)) {
  features <- explainer$X$features
} else {
  features <- explainer$X$features[index_features]
}

# THIS LOOP IS WHERE THE SAMPLES ARE TAKEN: 
for (i in seq(n_xtest)) {
  l <- lapply(
    X = features,
    FUN = sample_gaussian,
    n_samples = explainer$n_samples,
    mu = explainer$mu,
    cov_mat = explainer$cov_mat,
    m = ncol(explainer$x_test),
    x_test = explainer$x_test[i, , drop = FALSE]
  )
  
  dt_l[[i]] <- data.table::rbindlist(l, idcol = "id_combination")
  dt_l[[i]][, w := 1 / explainer$n_samples]
  dt_l[[i]][, id := i]
  if (!is.null(index_features)) dt_l[[i]][, id_combination := index_features[id_combination]]
}

dt_airqual <- data.table::rbindlist(dt_l, use.names = TRUE, fill = TRUE)

################################# inside prediction(): 

id <- w <- id_combination <- p_hat <- NULL # due to NSE notes in R CMD check
stopifnot(
  data.table::is.data.table(dt),
  !is.null(dt[["id"]]),
  !is.null(dt[["id_combination"]]),
  !is.null(dt[["w"]])
)

# Setup
feature_names <- colnames(explainer$x_test)
data.table::setkeyv(dt, c("id", "id_combination"))

# Check that the number of test observations equals max(id)
stopifnot(nrow(explainer$x_test) == dt[, max(id)])

# Reducing the prediction data.table
max_id_combination <- nrow(explainer$S)
V1 <- keep <- NULL # due to NSE notes in R CMD check
dt[, keep := TRUE] # adds keep column, all values set to TRUE first
first_element <- dt[, tail(.I, 1), .(id, id_combination)][id_combination %in% c(1, max_id_combination), V1]
dt[id_combination %in% c(1, max_id_combination), keep := FALSE]
# full and empty coalitions set to false (the rows with id_comb = 1 or 16)
dt[first_element, c("keep", "w") := list(TRUE, 1.0)] # why are they set back to true here?
dt <- dt[keep == TRUE][, keep := NULL] 

# Predictions
if (!all(dt[, unique(id_combination)] == 1)) { # Avoid warnings when predicting with empty newdata
  dt[id_combination != 1, p_hat := predict_model(explainer$model, newdata = .SD), .SDcols = feature_names]
  # code should go in here.
}
dt[id_combination == 1, p_hat := prediction_zero]

if (dt[, max(id_combination)] < max_id_combination) {
  p_all <- NULL
} else {
  p_all <- dt[id_combination == max_id_combination, p_hat] # takes all the empty coalition p_hat values
  names(p_all) <- 1:nrow(explainer$x_test)
}

# Calculate contributions, THIS IS WHERE THE INTEGRAL APPROXIMATION IS DONE
dt_res <- dt[, .(k = sum((p_hat * w) / sum(w))), .(id, id_combination)]
# dt res basically ends up being the predictions col of dt_mat, w doesn't seem to affect anything when set to 1 or 1/n, 
# just makes the integral approximation into an average.
dt_mat <- data.table::dcast(dt_res, id_combination ~ id, value.var = "k")
dt_mat[, id_combination := NULL]
# these lines have made dt -  a table of coalitions with corresponding subset values - 
# into the matrix of subset values dt_mat.

r <- list(p = p_all, dt_mat = dt_mat)

#################################  getting samples from forge: 

# run init.R to get the dummy dataset. 

library(arf)
arf <- adversarial_rf(x_train, num_trees = 100)
psi <- forde(arf, x_train)

coalition <- dt[row_id == 7]

variable <- operator <- value <- prob <- f_idx <- cvg <- wt <- . <- NULL

# this is all we need: 
evi = data.frame(variable = c(as.character(coalition[!is.na(value), variable])), 
                 relation = "==", 
                 value = c(coalition[!is.na(value), value]))
forge(psi, 1000, evidence = evi)

# as a loop
dt <- data.table(x_test)
dt[,id := .I]
dt <- dt[rep(seq_len(nrow(x_test)), each = 2^ncol(x_test))]
dt <- as.matrix(dt)
S <- matrix(rep(t(explainer$S),nrow(x_test)),ncol=ncol(explainer$S),byrow=TRUE)
coal <- matrix(NA, nrow = nrow(S), ncol = ncol(S))
coal[,][S == 1] <- dt[,1:(ncol(dt)-1)][S == 1]
dt <- as.data.table(coal)
colnames(dt) <- explainer$feature_list$labels
dt[, id := rep(seq_along(1:nrow(x_test)), each = 2^ncol(x_test))]
dt[, id_combination := rep(seq_along(1:2^ncol(x_test)),nrow(x_test))]
dt <- melt(dt, id.vars = c("id", "id_combination"))
dt <- dt[order(id, id_combination)]
dt[, row_id := .GRP, by = c("id", "id_combination")]
to_eval <- unique(dt[id_combination != 1 & id_combination != 2^ncol(x_test), row_id])
forge_wrapper <- function(params = psi, n_synth = 1000, datapoint) {
  if (datapoint[, unique(id_combination)] == 1 | datapoint[, unique(id_combination == 16)]) {
    out <- dcast(datapoint, row_id + id + id_combination ~ variable, value.var = "value")
    return(out)
  }
  evi = data.frame(variable = c(as.character(datapoint[!is.na(value), variable])), 
                   relation = "==", 
                   value = c(datapoint[!is.na(value), value]))
  out <- forge(params, n_synth, evi)
  out[, id := datapoint[, unique(id)]]
  out[, id_combination := datapoint[, unique(id_combination)]]
  out[, row_id := datapoint[, unique(row_id)]]
  return(out)
}
dt_out <- dt
dt_out[row_id %in% to_eval,] <- foreach(i = seq_len(nrow(x_test)*2^ncol(x_test)), .combine = rbind) %do% forge(psi, 100, evi = data.frame(variable = c(as.character(dt[row_id == i,.(variable, value)][!is.na(value), variable])),
                                                                                  relation = "==",
                                                                                  value = c(dt[row_id == i,.(variable, value)][!is.na(value), value])))
foreach(i = seq_len(nrow(x_test)*2^ncol(x_test)), .combine = rbind) %do% forge_wrapper(psi, 100, dt[row_id == i, ])
