#' Calculate Shapley weight (More specifically the shapley kernel weights! See top of pg 5 of Aas et al.)
#'
#' @param m Positive integer. Total number of features/feature groups.
#' @param n_components Positive integer. Represents the number of features/feature groups you want to sample from
#' a feature space consisting of \code{m} unique features/feature groups. Note that \code{ 0 < = n_components <= m}.
#' @param N Positive integer. The number of unique combinations when sampling \code{n_components} features/feature
#' groups, without replacement, from a sample space consisting of \code{m} different features/feature groups.
#' @param weight_zero_m Positive integer. Represents the Shapley weight for two special
#' cases, i.e. the case where you have either \code{0} or \code{m} features/feature groups.
#'
#' @return Numeric
#' @keywords internal
#'
#' @author Nikolai Sellereite
shapley_weights <- function(m, N, n_components, weight_zero_m = 10^6) {
  x <- (m - 1) / (N * n_components * (m - n_components))
  x[!is.finite(x)] <- weight_zero_m
  x
}

#' Calculate weighted matrix
#'
#' @param X data.table
#' @param normalize_W_weights Logical. Whether to normalize the weights for the combinations to sum to 1 for
#' increased numerical stability before solving the WLS (weighted least squares). Applies to all combinations
#' except combination \code{1} and \code{2^m}.
#' @param is_groupwise Logical. Indicating whether group wise Shapley values are to be computed.
#'
#' @return Numeric matrix. See \code{\link{weight_matrix_cpp}} for more information.
#' @keywords internal
#'
#' @author Nikolai Sellereite, Martin Jullum
weight_matrix <- function(X, normalize_W_weights = TRUE, is_groupwise = FALSE) {

  # Fetch weights
  w <- X[["shapley_weight"]]

  if (normalize_W_weights) {
    w[-c(1, length(w))] <- w[-c(1, length(w))] / sum(w[-c(1, length(w))])
  }

  if (!is_groupwise) {
    W <- weight_matrix_cpp(
      subsets = X[["features"]],
      m = X[.N][["n_features"]],
      n = X[, .N],
      w = w
    )
  } else {
    W <- weight_matrix_cpp(
      subsets = X[["groups"]],
      m = X[.N][["n_groups"]],
      n = X[, .N],
      w = w
    )
  }

  return(W)
}

#' Create an explainer object with Shapley weights for test data.
#'
#' @param x Numeric matrix or data.frame/data.table. Contains the data used to estimate the (conditional)
#' distributions for the features needed to properly estimate the conditional expectations in the Shapley formula.
#'
#' @param model The model whose predictions we want to explain. Run
#' \code{\link[shapr:get_supported_models]{shapr:::get_supported_models()}}
#' for a table of which models \code{shapr} supports natively.
#'
#' @param n_combinations Integer. The number of feature combinations to sample. If \code{NULL},
#' the exact method is used and all combinations are considered. The maximum number of
#' combinations equals \code{2^ncol(x)}.
#'
#' @param group List. If \code{NULL} regular feature wise Shapley values are computed.
#' If provided, group wise Shapley values are computed. \code{group} then has length equal to
#' the number of groups. The list element contains character vectors with the features included
#' in each of the different groups.
#'
#' @return Named list that contains the following items:
#' \describe{
#'   \item{exact}{Boolean. Equals \code{TRUE} if \code{n_combinations = NULL} or
#'   \code{n_combinations < 2^ncol(x)}, otherwise \code{FALSE}.}
#'   \item{n_features}{Positive integer. The number of columns in \code{x}}
#'   \item{S}{Binary matrix. The number of rows equals the number of unique combinations, and
#'   the number of columns equals the total number of features. I.e. let's say we have a case with
#'   three features. In that case we have \code{2^3 = 8} unique combinations. If the j-th
#'   observation for the i-th row equals \code{1} it indicates that the j-th feature is present in
#'   the i-th combination. Otherwise it equals \code{0}.}
#'   \item{W}{Matrix. This matrix is equal to the matrix \code{R_D} in Equation 7 in the reference
#'   of \code{link{explain}}. The Shapley value for a test observation will be equal to the matrix-vector product
#'   of \code{W} and the contribution vector.}
#'   \item{X}{data.table. Returned object from \code{\link{feature_combinations}}}
#'   \item{x_train}{data.table. Transformed \code{x} into a data.table.}
#'   \item{feature_list}{List. The \code{updated_feature_list} output from
#'   \code{\link[shapr:preprocess_data]{preprocess_data}}}
#' }
#'
#' In addition to the items above, \code{model} and \code{n_combinations} are also present in the returned object.
#'
#' @export
#'
#' @author Nikolai Sellereite
#'
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   # Load example data
#'   data("Boston", package = "MASS")
#'   df <- Boston
#'
#'   # Example using the exact method
#'   x_var <- c("lstat", "rm", "dis", "indus")
#'   y_var <- "medv"
#'   df0 <- df[, x_var]
#'   model <- lm(medv ~ lstat + rm + dis + indus, data = df)
#'   explainer <- shapr(df0, model)
#'
#'   print(nrow(explainer$X))
#'   # 16 (which equals 2^4)
#'
#'   # Example using approximation
#'   y_var <- "medv"
#'   model <- lm(medv ~ ., data = df)
#'   explainer <- shapr(df, model, n_combinations = 1e3)
#'
#'   print(nrow(explainer$X))
#'
#'   # Example using approximation where n_combinations > 2^m
#'   x_var <- c("lstat", "rm", "dis", "indus")
#'   y_var <- "medv"
#'   model <- lm(medv ~ lstat + rm + dis + indus, data = df)
#'   explainer <- shapr(df0, model, n_combinations = 1e3)
#'
#'   print(nrow(explainer$X))
#'   # 16 (which equals 2^4)
#'
#'   # Example using groups
#'   group <- list(A=x_var[1:2], B=x_var[3:4])
#'
#'   explainer_group <- shapr(df0, model, group = group)
#'   print(nrow(explainer_group$X))
#'   # 4 (which equals 2^(#groups))
#' }
shapr <- function(x,
                  model,
                  n_combinations = NULL,
                  group = NULL) {

  # Checks input argument
  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x should be a matrix or a dataframe.")
  }

  # Setup
  explainer <- as.list(environment())
  explainer$exact <- ifelse(is.null(n_combinations), TRUE, FALSE)


  # Check features of training data against model specification
  feature_list_model <- get_model_specs(model)

  processed_list <- preprocess_data(
    x = x,
    feature_list = feature_list_model
  )



  x_train <- processed_list$x_dt
  updated_feature_list <- processed_list$updated_feature_list

  explainer$n_features <- ncol(x_train)


  # Processes groups if specified. Otherwise do nothing
  is_groupwise <- !is.null(group)
  if (is_groupwise) {
    group_list <- process_groups(
      group = group,
      feature_labels = updated_feature_list$labels
    )
    group <- group_list$group
    group_num <- group_list$group_num
  } else {
    group_num <- NULL
  }

  # Checking that the prediction function works
  tmp <- predict_model(model, head(x_train, 2))
  if (!(all(is.numeric(tmp)) & length(tmp) == 2)) {
    stop(
      paste0(
        "The predict_model function of class ", class(model), " is invalid.\n",
        "See the 'Advanced usage' section of the vignette:\n",
        "vignette('understanding_shapr', package = 'shapr')\n",
        "for more information on running shapr with custom models.\n"
      )
    )
  }

  # Get all combinations ----------------
  dt_combinations <- feature_combinations(
    m = explainer$n_features,
    exact = explainer$exact,
    n_combinations = n_combinations,
    weight_zero_m = 10^6,
    group_num = group_num
  )
  # View(dt_combinations, "dt_combinations") # the shapley weights col contains the 
  #                                          # shapley kernel weight values (aas et al. top of pg5 has deets on this)
  # Get weighted matrix ----------------
  weighted_mat <- weight_matrix(
    X = dt_combinations,
    normalize_W_weights = TRUE,
    is_groupwise = is_groupwise
  )
  # View(weighted_mat, "weight matrix") # this is the matrix R, or Rd, described in eq5 Aas et al. pg5. 
  ## Get feature matrix ---------
  feature_matrix <- feature_matrix_cpp(
    features = dt_combinations[["features"]],
    m = explainer$n_features
  )
  # View(feature_matrix, "feature matrix") # this is just the S matrix
  # Updating explainer$exact as done in feature_combinations
  if (!explainer$exact && n_combinations > (2^explainer$n_features - 2)) {
    explainer$exact <- TRUE
  }

  explainer$S <- feature_matrix
  explainer$W <- weighted_mat
  explainer$X <- dt_combinations
  explainer$x_train <- x_train
  explainer$x <- NULL
  explainer$feature_list <- updated_feature_list
  explainer$group <- group
  explainer$is_groupwise <- is_groupwise
  explainer$n_combinations <- nrow(feature_matrix)

  attr(explainer, "class") <- c("explainer", "list")

  return(explainer)
}

#' @keywords internal
distance_matrix <- function(x_train, x_test = NULL, list_features, mcov) {
  if (is.null(x_test)) {
    return(NULL)
  }

  if (is.null(dim(x_test))) {
    x_test <- t(as.matrix(x_test))
  }
  # Note that D equals D_S(,)^2 in the paper
  print(class(as.matrix(x_train)))
  D <- mahalanobis_distance_cpp(
    featureList = list_features,
    Xtrain_mat = as.matrix(x_train),
    Xtest_mat = as.matrix(x_test),
    mcov = mcov,
    S_scale_dist = TRUE
  )

  # Normalize distance rows to ensure numerical stability in later operations
  colmin <- apply(X = D, MARGIN = c(2, 3), FUN = min)
  for (i in 1:dim(D)[3]) {
    D[, , i] <- t(t(D[, , i]) - colmin[, i])
  }

  return(D)
}
