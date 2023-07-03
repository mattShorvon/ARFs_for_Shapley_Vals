#' Generates coalition data used to calculate coalition values. Called from
#' prepare_and_predict() in explanation.R, in turn called from explain.arf(),
#' also in explanation.R.
#'
#' Takes the test data, uses a binary matrix mask to create the coalitions containing
#' different subsets of features, imputes the missing values for each coalition using
#' prep() (where the imputed values are MAP estimates of the missing features, given
#' the values of the observed features) and returns a data table of coalitions with
#' all features filled, ready to be used in coalition value calculations.
#'
#' N.B. the current implemented method of finding the MAP estimates does this using
#' weighted means. This works if the distributions in the ARF leaves are Gaussian,
#' which may not always be true.
#'
#' @param x explainer object created from a call to shapr(), containing the test
#' data.
#'
#' @param index_features Positive integer vector. Specifies the indices of
#' combinations to apply to the present method.\code{NULL} means all combinations.
#' Only used internally.
#'
#' @param psi object containing parameters of the ARF's leaves, created from a call
#' to adversarial_rf()
#'
#' @return A `data.table` containing filled coalition data, passed to
#' \code{\link{prediction}} in the prepare_and_predict() function.
prepare_data.arf <- function(x = explainer, index_features = NULL, psi = psi , ...) {
  print("Using vectorised version")
  dt <- data.table(x$x_test)
  dt[,id := .I]
  dt <- dt[rep(seq_len(nrow(x$x_test)), each = 2^ncol(x$x_test))]
  dt <- as.matrix(dt)
  S <- matrix(rep(t(x$S),nrow(x$x_test)),ncol=ncol(x$S),byrow=TRUE)
  coal <- matrix(NA, nrow = nrow(S), ncol = ncol(S))
  coal[,][S == 1] <- dt[,1:(ncol(dt)-1)][S == 1]
  dt <- as.data.table(coal)
  colnames(dt) <- x$feature_list$labels
  dt[, id := rep(seq_along(1:nrow(x$x_test)), each = 2^ncol(x$x_test))]
  dt[,id_combination := rep(seq_along(1:2^ncol(x$x_test)),nrow(x$x_test))]
  dt <- melt(dt, id.vars = c("id", "id_combination"))
  dt <- dt[order(id, id_combination)]
  var_names_cnt = psi[["cnt"]][, unique(variable)]
  is_cnt <- x$feature_list$labels %in% var_names_cnt
  dt[, type := fifelse(variable %in% var_names_cnt, "cnt", "cat")]
  dt[, row_id := .GRP, by = c("id", "id_combination")]
  dt[, mu_cvg_sum := as.numeric(NA)]
  dt[, cvg_sum := as.numeric(NA)]
  to_eval <- unique(dt[id_combination != 1 & id_combination != 2^ncol(x$x_test), row_id])
  dt[row_id %in% to_eval,] <- foreach(i = 1:length(to_eval), .combine = rbind,
                                      .packages = "data.table", .export = "prep") %dopar%
    {prep(dt[row_id == to_eval[i],], psi)
    }
  dt[is.na(value) & type == "cnt" & !is.na(cvg_sum), value := mu_cvg_sum/cvg_sum]
  dt[, type := NULL]
  dt <- dcast(dt, row_id + id + id_combination ~ variable, value.var = "value")
  for (column in x$feature_list$labels) {
    col_class <- x$feature_list$classes[column]
    dt[, (column) := switch(col_class,
                            numeric = as.numeric(get(column)),
                            character = as.character(get(column)),
                            factor = as.factor(get(column)))]
  }
  dt[, w := 1]
  dt[, row_id := NULL]
  return(dt)
}

#' Takes a coalition containing missing features, finds leaves in the ARF that
#' contain data that adhere to constraints created by the observed features
#' (specifically, the observed feature values are between min and max of the leaf)
#' and imputes values for the missing features by taking the MAP estimates of
#' the distributions formed by combining the distributions of these leaves.
#'
#' For continuous features, the mu and cvg of all leaves adhering to the
#' continuous constraint are stored in leaves_info_cnt and the sum of mu's and
#' cvg's are then calculated and added to the coalition's data, to be used to
#' calculate the imputed value. For categorical features, the prob of each possible
#' value in each leaf and the cvg of each leaf is stored in leaves_info_cat, and
#' the prob times cvg sum values and cvg sum values are calculated for each
#' feature value, and the value with the highest prob_cvg_sum / cvg_sum is taken
#' as the imputed value.
#'
#' @param coalition a row of data created from a test point that contains some
#' subset of feature values missing, inputted in melted form with feature names
#' and values as columns.
#'
#' @param psi object containing parameters of the ARF's leaves, created from a call
#' to adversarial_rf()
#'
#' @return coalition in melted form containing data for the cvg_sum and mu times
#' cvg sum columns for continuous features, and imputed values for the categorical
#' features
#'
prep <- function(coalition, psi) {
  # Find the ids of leaves that satisfy constraints from the input coalition:
  coal <- coalition[!is.na(value)]
  leafIDs = as.list(rep(NA, nrow(coal)))
  leaves = numeric()
  for (i in sequence(nrow(coal))) {
    ifelse(coal[i, type] == "cnt", {
      leafIDs[[i]] <- psi$cnt[variable == coal[i, variable]
                              & min < as.numeric(coal[i, value])
                              & max > as.numeric(coal[i, value]), f_idx]},
      {leafIDs[[i]] <- psi$cat[variable == coal[i, variable]
                               & val == coal[i, value]
                               & prob > 0, f_idx]}
    )
  }
  leaves <- Reduce(intersect, leafIDs)
  if (length(leaves) == 0) {
    leaves <- unique(psi$cnt[,f_idx])
  }

  # Find the continuous data info for leaves that satisfy these constraints, add it to coalition
  if (nrow(coalition[is.na(value) & type == "cnt"]) > 0) {
    leaves_info_cnt <- merge(psi[["cnt"]][f_idx %in% leaves
                                          & variable %in% coalition[, unique(variable)],
                                          .(variable, mu, f_idx)],
                             psi$forest[f_idx %in% leaves,
                                        .(cvg, f_idx)], by = 'f_idx', sort = F)

    leaves_info_cnt[, mu_x_cvg := mu*cvg]
    coalition[is.na(value) & type == "cnt", cvg_sum := leaves_info_cnt[variable == .SD[, variable],
                                                                       sum(cvg),
                                                                       by = variable]$V1,
              by = variable, .SDcols = "variable"]
    coalition[is.na(value) & type == "cnt", mu_cvg_sum := leaves_info_cnt[variable == .SD[, variable],
                                                                          sum(mu_x_cvg),
                                                                          by = variable]$V1,
              by = variable, .SDcols = "variable"]
  }

  # Find the categorical data info for leaves that satisfy these constraints, add it to the coalition
  if (nrow(coalition[is.na(value) & type == "cat"]) > 0) {
    leaves_info_cat <- merge(psi[["cat"]][f_idx %in% leaves & variable %in% coalition[, unique(variable)],
                                          .(variable, val, prob, f_idx)],
                             psi$forest[f_idx %in% leaves, .(cvg, f_idx)],
                             by = 'f_idx', sort = F)
    leaves_info_cat[, prob_x_cvg := prob*cvg]
    leaves_info_cat[, cvg_sum := sum(cvg), by = .(val, variable)]
    leaves_info_cat[, prob_cvg_sum := sum(prob_x_cvg ), by = .(val, variable)]
    leaves_info_cat[, division := prob_cvg_sum/cvg_sum]
    coalition[is.na(value) & type == "cat", value := unique(leaves_info_cat[variable == .SD[, variable]]
                                                            [division == max(division), val])[1],
              by = variable, .SDcols = "variable"]
  }
  return(coalition)
}
