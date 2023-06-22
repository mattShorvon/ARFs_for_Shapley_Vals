#' Calculate Shapley weights for test data
#'
#' @description This function should only be called internally, and not be used as
#' a stand-alone function.
#'
#' @param dt data.table
#' @param prediction_zero Numeric. The value to use for \code{phi_0}.
#' @param explainer An object of class \code{explainer}. See \code{\link{shapr}}.
#'
#' @details If \code{dt} does not contain three columns called \code{id}, \code{id_combination} and \code{w}
#' the function will fail. \code{id} represents a unique key for a given test observation,
#' and \code{id_combination} is a unique key for which feature combination the row represents. \code{w}
#' represents the Shapley value of feature combination given by \code{id_combination}. In addition
#' to these three columns, \code{dt} should also have columns which matches the variables used
#' when training the model.
#'
#' I.e. you have fitted a linear model using the features \code{x1},
#' \code{x2} and \code{x3}, and you want to explain 5 test observations using the exact method, i.e.
#' setting \code{exact = TRUE} in \code{\link{shapr}}, the following properties should be satisfied
#' \enumerate{
#' \item \code{colnames(dt)} equals \code{c("x1", "x2", "x3", "id", "id_combination", ""w)}
#' \item \code{dt[, max(id)]} equals the number of test observations
#' \item \code{dt[, min(id)]} equals 1L.
#' \item \code{dt[, max(id_combination)]} equals \code{2^m} where m equals the number of features.
#' \item \code{dt[, min(id_combination)]} equals 1L.
#' \item \code{dt[, type(w)]} equals \code{double}.
#' }
#'
#'
#' @return An object of class \code{c("shapr", "list")}. For more details see \code{\link{explain}}.
#'
#' @keywords internal
#'
#' @author Nikolai Sellereite
prediction <- function(dt, prediction_zero, explainer) {
  # View(dt, "dt line 36 prediction.R")
  # Checks on input data
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

  # View(dt, "dt line 53 prediction.R")
  write.csv(dt, file = "shaprcoalitions.csv")
  # Reducing the prediction data.table
  max_id_combination <- nrow(explainer$S)
  V1 <- keep <- NULL # due to NSE notes in R CMD check
  dt[, keep := TRUE] # adds keep column, all values set to TRUE first
  first_element <- dt[, tail(.I, 1), .(id, id_combination)][id_combination %in% c(1, max_id_combination), V1]
  dt[id_combination %in% c(1, max_id_combination), keep := FALSE]
  # full and empty coalitions set to false (the rows with id_comb = 1 or 16)
  dt[first_element, c("keep", "w") := list(TRUE, 1.0)]
  dt <- dt[keep == TRUE][, keep := NULL] # then the full and empty ones are set back to true in arf! is this correct? def nottt
                                         # for empirical where there are obv a lot more points the last
                                         # point of the id_comb = 1 group is set to true with w = 1 (rest are still false) then
                                         # all others are discarded basically. Same with id_comb = 16. Also weirdly,
                                         # id_comb = 16 is the full set, whilst id_comb = 1 is the empty. THIS IS IMPORTANT THO
  # View(dt, "dt line 65 prediction.R")

  # Predictions
  if (!all(dt[, unique(id_combination)] == 1)) { # Avoid warnings when predicting with empty newdata
    dt[id_combination != 1, p_hat := predict_model(explainer$model, newdata = .SD), .SDcols = feature_names]
    # code should go in here.
  }
  dt[id_combination == 1, p_hat := prediction_zero]
  # View(dt, "dt after predict_model()")

  if (dt[, max(id_combination)] < max_id_combination) {
    p_all <- NULL
  } else {
    p_all <- dt[id_combination == max_id_combination, p_hat] # takes all the empty coalition p_hat values
    names(p_all) <- 1:nrow(explainer$x_test)
  }

  # Calculate contributions
  dt_res <- dt[, .(k = sum((p_hat * w) / sum(w))), .(id, id_combination)]
  # dt res basically ends up being the predictions col of dt_mat, w doesn't seem to affect anything when set to 1.
  dt_mat <- data.table::dcast(dt_res, id_combination ~ id, value.var = "k")
  dt_mat[, id_combination := NULL]
  # these lines have made dt, a table of coalitions with corresponding subset values,
  # into the matrix of subset values dt_mat.

  r <- list(p = p_all, dt_mat = dt_mat)

  return(r)
}


#' Compute shapley values
#' @param explainer An \code{explain} object.
#' @param contribution_mat The contribution matrix.
#' @return A \code{data.table} with shapley values for each test observation.
#' @export
#' @keywords internal
compute_shapley <- function(explainer, contribution_mat) {

  feature_names <- colnames(explainer$x_test)
  if (!explainer$is_groupwise) {
    print('success')
    shap_names <- feature_names
  } else {
    shap_names <- names(explainer$group)
  }

  # View(contribution_mat, "contrib mat in computeshapley()")
  kshap <- t(explainer$W %*% contribution_mat)
  dt_kshap <- data.table::as.data.table(kshap)
  colnames(dt_kshap) <- c("none", shap_names)

  return(dt_kshap)

}
