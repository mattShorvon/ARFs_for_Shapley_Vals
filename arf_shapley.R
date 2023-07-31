### Shapley Applications ###

# Load libraries
library(data.table)
library(arf)
library(palmerpenguins)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123)

# Prep dataset
df <- as.data.table(na.omit(penguins))
df[, flipper_length_mm := as.numeric(flipper_length_mm)]
df[, body_mass_g := as.numeric(body_mass_g)]
d <- ncol(df)

# Fit ARF, learn circuit
arf <- adversarial_rf(df, num_trees = 50)
psi <- forde(arf, df, alpha = 0.1)
factor_cols <- psi$meta[, family == 'multinom']

# Just taking the first point as an example
sample_i <- setDF(df[1, ])

# Precompute all marginal likelihoods
psi_cnt <- psi_cat <- NULL

# Continuous case
if (any(!factor_cols)) {
  fam <- psi$meta[class == 'numeric', unique(family)]
  x_long <- melt(
    setDT(sample_i[, !factor_cols, drop = FALSE]), 
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
      setDT(sample_i[, factor_cols, drop = FALSE]), 
      measure.vars = seq_len(sum(factor_cols)), value.name = 'val', 
      variable.factor = FALSE
    )
  )
  grd <- rbindlist(lapply(which(factor_cols), function(j) {
    expand.grid('f_idx' = psi$forest$f_idx, 
                'variable' = colnames(sample_i)[j],
                'val' = x_long[variable == colnames(sample_i)[j], val],
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
o <- expand.grid(rep(list(c(0, 1)), times = d))
colnames(o) <- colnames(df)

# Wrapper function for MAP calls
wrap_fn <- function(i) {
  bitv <- as.logical(o[i, ])
  if (all(bitv)) {
    # Impute nothing
    out <- sample_i
  } else {
    if (all(!bitv)) {
      # Impute everything
      q <- colnames(df)
      e <- NULL
    } else {
      # Impute some but not all things
      q <- colnames(sample_i)[!bitv]
      vars <- colnames(sample_i)[bitv]
      # Compute posterior over leaves from the marginals over features
      e <- psi_x[variable %in% vars]
      e[, wt := cvg * prod(lik), by = f_idx]
      e <- unique(e[wt > 0, .(f_idx, wt)])
      e[, wt := wt / sum(wt)] 
    }
    imputed <- map(psi, q, e)
    out <- setDT(cbind(sample_i[, bitv, drop = FALSE], imputed))
    setcolorder(out, colnames(df))
  }
  return(out)
}
out <- foreach(idx = seq_len(2^d), .combine = rbind) %do% wrap_fn(idx)

# This executes in about 5.6 seconds on my machine (no parallelization)
library(microbenchmark)
microbenchmark(
  'tst' = foreach(idx = seq_len(2^d), .combine = rbind) %do% wrap_fn(idx),
  times = 10L
)

# Ramblings...
# This suggests an alternative approach, more in line with the GeF method of 
# Correia et al., where we actually marginalize instead of just imputing
# missing values. That is, our model for Y is secretly 2^d models, one for
# each missingness pattern. We call whichever we need for each query. Is this
# faster than imputing? Maybe. It's not model agnostic, of course. 







