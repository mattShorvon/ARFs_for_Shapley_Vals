#' Compute leaf posterior
#'
#' This function returns a posterior distribution on leaves, conditional on some
#' evidence.
#'
#' @param params Parameters learned via \code{\link{forde}}.
#' @param evidence Data frame of conditioning event(s).
#'
#' @import data.table
#' @importFrom truncnorm dtruncnorm ptruncnorm
#'
leaf_posterior <- function(params, evidence) {

  # To avoid data.table check issues
  variable <- operator <- value <- prob <- f_idx <- cvg <- wt <- . <- NULL

  # Likelihood per leaf-event
  leaf_list <- lapply(seq_len(nrow(evidence)), function(k) {
    j <- evidence$variable[k]
    op <- evidence$operator[k]
    value <- evidence$value[k]
    if (params$meta[variable == j, class == 'numeric']) {
      value <- as.numeric(value)
      psi <- params$cnt[variable == j]
      if (op == '==') {
        psi[, prob := truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      } else if (op %in% c('<', '<=')) {
        psi[, prob := truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      } else if (op %in% c('>', '>=')) {
        psi[, prob := 1 - truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      }
    } else {
      psi <- params$cat[variable == j]
      grd <- expand.grid(f_idx = params$forest$f_idx, val = psi[, unique(val)])
      psi <- merge(psi, grd, by = c('f_idx', 'val'), all.y = TRUE, sort = FALSE)
      psi[is.na(prob), prob := 0][is.na(variable), variable := j]
      if (op == '==') {
        psi <- psi[val == value]
      } else if (op == '!=') {
        psi <- psi[val != value]
        psi[, prob := sum(prob), by = f_idx]
        psi <- unique(psi[, .(f_idx, variable, prob)])
      }
    }
    out <- psi[, .(f_idx, variable, prob)]
    return(out)
  })

  # Weight is proportional to coverage times product of likelihoods
  psi <- rbindlist(leaf_list)
  psi <- merge(psi, params$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  psi[, wt := cvg * prod(prob), by = f_idx] # Worth doing in log space?

  # Normalize, export
  out <- unique(psi[, .(f_idx, wt)])
  out[, wt := wt / sum(wt)]
  return(out)
}

# example query:
evidence <- data.frame(variable = "cont_1_", operator = "==", value = 1)
evidence <- data.frame(variable = "cat_1_", operator = "==", value = 1)


# plugging in coalition instead:
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

