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
  to_eval <- unique(dt[id_combination != 1 & id_combination != 2^ncol(x$x_test), row_id])
  dt[row_id %in% to_eval,] <- foreach(i = 1:length(to_eval), .combine = rbind,
                                      .packages = "data.table", .export = "prep") %dopar%
    {prep(dt[row_id == to_eval[i],], psi)
    }
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
  # get table of p(th | x) values where x is the set of observed features in the coalition
  # start_t <- Sys.time()
  leaves_info <- leaf_posterior_shapley(params = psi, coalition = coalition)
  leaves_info <- unique(leaves_info[, .(f_idx, p_th_given_x)])
  # print(Sys.time() - start_t)

  # loop through unobserved features, add column of mu's to p(th | x) table,
  # do an element-wise multiplication between these two columns, sum the resulting
  # column to get the MAP.
  # start_t <- Sys.time()
  unobserved <- coalition[is.na(value)]
  for (i in 1:nrow(coalition[is.na(value)])) {
  coal <- unobserved[i]
    var <- coal[, variable]
    if (coal$type == "cnt") {
      leaves_info_cnt <- merge(psi$cnt[variable == var, .(mu, f_idx)],
                              leaves_info[, .(f_idx, p_th_given_x)],
                              by = 'f_idx', sort = F)
      leaves_info_cnt[, mu_x_p := mu * p_th_given_x]
      coalition[variable == var, value := sum(leaves_info_cnt[, mu_x_p])]
    } else {
      leaves_info_cat <- psi$cat[variable == var,]
      leaves_info_cat <- merge(leaves_info_cat, leaves_info, by = 'f_idx', sort = F)
      leaves_info_cat[, p_y_th_times_p_th_x := prob * p_th_given_x]
      leaves_info_cat[, total := sum(.SD[, p_y_th_times_p_th_x]), by = val]
      coalition[variable == var, value := leaves_info_cat[which.max(total), val]]
    }
  }
  # print(Sys.time() - start_t)
  return(coalition)
}

leaf_posterior_shapley <- function(params, coalition) {

  # To avoid data.table check issues
  variable <- operator <- value <- prob <- f_idx <- cvg <- wt <- . <- NULL

  # Likelihood per leaf-event
  coal <- coalition[!is.na(value),]
  leaf_list <- lapply(seq_len(nrow(coal)), function(k) {
    j <- coal$variable[k]
    value <- coal$value[k]
    if (params$meta[variable == j, class == 'numeric']) {
      value <- as.numeric(value)
      leaves <- params$cnt[variable == j]
      leaves[, prob := truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    } else {
      leaves <- params$cat[variable == j]
      grd <- expand.grid(f_idx = params$forest$f_idx, val = leaves[, unique(val)])
      leaves <- merge(leaves, grd, by = c('f_idx', 'val'), all.y = TRUE, sort = FALSE)
      leaves[is.na(prob), prob := 0][is.na(variable), variable := j]
      leaves <- leaves[val == value]
    }
    out <- leaves[, .(f_idx, variable, prob)]
    return(out)
  })

  # Weight is proportional to coverage times product of likelihoods
  leaves <- rbindlist(leaf_list)
  leaves <- merge(leaves, params$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  leaves[, p_th_given_x := cvg * prod(prob), by = f_idx] # Worth doing in log space?

  # Normalize, export
  out <- unique(leaves)
  out[, p_th_given_x := p_th_given_x / sum(p_th_given_x)]
  return(out)
}
