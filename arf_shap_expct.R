# requires arf from https://github.com/bips-hb/arf/tree/conditional_sampling

explain.arf_expct <- function(x, explainer, approach, psi, prediction_zero, n_batches = 1, n_samples = 1e3) {
  
  # Add arguments to explainer object
  explainer$x_test <- as.matrix(preprocess_data(x, explainer$feature_list)$x_dt)
  explainer$approach <- approach
  explainer$n_samples <- n_samples
  
  r <- prepare_and_predict(explainer, n_batches, prediction_zero, psi)
}

prepare_data.arf_expct <- function(x = explainer, index_features = NULL, psi = psi , ...) {
  print("Using arf_expct method")
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
  dt[, id_combination := rep(seq_along(1:2^ncol(x$x_test)),nrow(x$x_test))]
  dt <- melt(dt, id.vars = c("id", "id_combination"))
  dt <- dt[order(id, id_combination)]
  dt[, row_id := .GRP, by = c("id", "id_combination")]
  full_coals <- dt[id_combination == 2^ncol(x$x_test),]
  dt_out <- foreach(i = seq_len(nrow(x$x_test)*2^ncol(x$x_test)), .combine = rbind) %dopar%
    expct_wrapper(params = psi, datapoint = dt[row_id == i, ])
  dt_out[, w := 1]
  
  for (column in x$feature_list$labels) {
    col_class <- x$feature_list$classes[column]
    dt_out[, (column) := switch(col_class,
                                numeric = as.numeric(get(column)),
                                character = as.character(get(column)),
                                factor = as.factor(get(column)))]
  }
  dt_out[, row_id := NULL]
  setcolorder(dt_out, neworder =  c("id_combination", psi$meta$variable,"w", "id"))
  print(dt_out[1:50,])
  print("any NAs in arf expct dt_out:")
  print(anyNA(dt_out[id_combination != 64,]))
  return(dt_out)
}

expct_wrapper <- function(params, datapoint) {
  # if datapoint is the empty coalition, compute the expectation for all variables:
  if (datapoint[, unique(id_combination)] == 1) {
    out <- expct(params = params)
    out[, id := datapoint[, unique(id)]]
    out[, id_combination := datapoint[,unique(id_combination)]]
    out[, row_id := datapoint[,unique(row_id)]]
    return(out)
  }
  # if the datapoint is the full coalition, return the full coalition in wide form
  last_coal <- 2^(nrow(datapoint))
  if (datapoint[, unique(id_combination == last_coal)]) {
    out <- dcast(datapoint, row_id + id + id_combination ~ variable, value.var = "value")
    return(out)
  }
  # otherwise, construct evidence dataframe and compute expected values for missing variables
  evi = data.frame(variable = c(as.character(datapoint[!is.na(value), variable])),
                   relation = "==",
                   value = c(datapoint[!is.na(value), value]))
  out <- expct(params = params, evidence = evi)
  out[, id := datapoint[, unique(id)]]
  out[, id_combination := datapoint[, unique(id_combination)]]
  out[, row_id := datapoint[, unique(row_id)]]
  out_melt <- melt(out, id.vars = c("id", "id_combination", "row_id"))
  out <- rbind(datapoint[!is.na(value)], out_melt)
  out <- dcast(out, row_id + id + id_combination ~ variable, value.var = "value")
  return(out)
}