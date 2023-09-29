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
#' weighted means. This works if the distributions in the ARF leaves are trunc norms,
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
#' \code{\link{prediction}} in the prepare_and_predict() function. The prediction
#' function then calculated values for each of these coalitions. The vector of 
#' these coalition values is used in the matrix multiplication to find the shapley
#' values. 
prepare_data.arf <- function(x = explainer, index_features = NULL, psi = psi , ...) {
  print("Using vectorised version with new MAP method")
  dt <- data.table(x$x_test)
  dt[,id := .I]
  dt_out <- foreach(i = seq_len(nrow(x$x_test)),.combine = rbind) %dopar% map_wrap1(point_to_explain = dt[id == i],
                                                                               psi = psi)
  dt_out[, id := rep(seq_along(1:nrow(x$x_test)), each = 2^ncol(x$x_test))]
  dt_out[, id_combination := rep(seq_along(1:2^ncol(x$x_test)),nrow(x$x_test))]
  dt_out[, w := 1]
  for (column in x$feature_list$labels) {
    col_class <- x$feature_list$classes[column]
    dt_out[, (column) := switch(col_class,
                            numeric = as.numeric(get(column)),
                            character = as.character(get(column)),
                            factor = as.factor(get(column)))]
  }
  return(dt_out)
}

map_wrap1 <- function(point_to_explain, psi = psi) {
  # Precompute all marginal likelihoods
  psi_cnt <- psi_cat <- NULL
  point_to_explain[, id := NULL]
  factor_cols <- psi$meta[, family == 'multinom']
  cont_cols <- psi$meta[, family != 'multinom']
  # Continuous case
  if (any(!factor_cols)) {
    fam <- psi$meta[class == 'numeric', unique(family)]
    x_long <- melt(
      point_to_explain[, ..cont_cols, drop = FALSE], 
      measure.vars = seq_len(sum(!factor_cols)), variable.factor = FALSE
    )
    psi_cnt <- merge(psi$cnt, x_long, by = 'variable', sort = FALSE)
    if (fam == 'truncnorm') {
      psi_cnt[, lik := truncnorm::dtruncnorm(value, a = min, b = max, 
                                             mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi_cnt[, lik := stats::dunif(value, min = min, max = max)]
    }
    psi_cnt[value == min, lik := 0]
    psi_cnt[, c('min', 'max', 'mu', 'sigma', 'value') := NULL]
  }
  
  # Categorical case
  if (any(factor_cols)) {
    suppressWarnings(
      x_long <- melt(
        point_to_explain[, ..factor_cols, drop = FALSE], 
        measure.vars = seq_len(sum(factor_cols)), value.name = 'val', 
        variable.factor = FALSE
      )
    )
    grd <- rbindlist(lapply(which(factor_cols), function(j) {
      expand.grid('f_idx' = psi$forest$f_idx, 
                  'variable' = colnames(point_to_explain)[j],
                  'val' = x_long[variable == colnames(point_to_explain)[j], val],
                  stringsAsFactors = FALSE)
    }))
    psi_cat <- merge(psi$cat, grd, by = c('f_idx', 'variable', 'val'), 
                     sort = FALSE, all.y = TRUE)
    psi_cat[is.na(prob), prob := 0]
    psi_cat[, lik := prob]
    psi_cat[, c('val', 'prob') := NULL]
  }
  
  # Put it all together
  psi_x <- rbind(psi_cnt, psi_cat)
  psi_x <- merge(psi_x, psi$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  
  # Create all 2^d bit vectors of length d
  d <- length(psi$meta$variable)
  o <- expand.grid(rep(list(c(0, 1)), times = d))
  colnames(o) <- psi$meta$variable
  
  # Call map_wrap2 to obtain table of coalitions with filled (imputed) values
  out1 <- foreach(idx = seq_len(2^d), .combine = rbind) %do% map_wrap2(i = idx,
                                                                       o = o,
                                                                       psi = psi,
                                                                       point_to_explain = point_to_explain,
                                                                       psi_x = psi_x)
  # 
  # for (idx in 1:seq_len(2^d)) {
  #   print(idx)
  #   out1 <- map_wrap2(i = idx, o = o, psi = psi, point_to_explain = point_to_explain, psi_x = psi_x)
  # }
  return(out1)
}

map_wrap2 <- function(i, o, psi = psi, point_to_explain, psi_x) {
  bitv <- as.logical(o[i, ])
  if (all(bitv)) {
    # Impute nothing
    out2 <- point_to_explain # just return the original sample/coalition
  } else {
    if (all(!bitv)) {
      # Impute everything
      q <- colnames(point_to_explain)
      e <- NULL
    } else {
      # Impute some but not all things
      q <- colnames(point_to_explain)[!bitv]
      vars <- colnames(point_to_explain)[bitv]
      # Compute posterior over leaves from the marginals over features
      e <- psi_x[variable %in% vars]
      e[, wt := cvg * prod(lik), by = f_idx]
      e <- unique(e[wt > 0, .(f_idx, wt)])
      e[, wt := wt / sum(wt)] 
    }
    imputed <- map(psi, q, e)
    out2 <- cbind(point_to_explain[, ..bitv, drop = FALSE], imputed)
    setcolorder(out2, psi$meta$variable)
  }
  return(out2)
}
