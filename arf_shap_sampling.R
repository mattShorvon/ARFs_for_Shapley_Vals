explain.arf_sampling <- function(x, explainer, approach, psi, prediction_zero, n_batches = 1) {
  
  # Don't think I'm doing anything that needs a seed to be set?
  
  # Add arguments to explainer object
  explainer$x_test <- as.matrix(preprocess_data(x, explainer$feature_list)$x_dt)
  explainer$approach <- approach
  explainer$n_samples <- n_samples
  
  r <- prepare_and_predict(explainer, n_batches, prediction_zero, psi)
}

prepare_data.arf_sampling <- function(x = explainer, index_features = NULL, psi = psi , ...) {
  print("Using arf sampling method")
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
  dt_out <- foreach(i = seq_len(nrow(x$x_test)*2^ncol(x$x_test)), .combine = rbind) %do% 
    forge_wrapper(params = psi, n_synth = x$n_samples, datapoint = dt[row_id == i, ],full_coals = full_coals)
  dt_out[, w := 1/x$n_samples]
  for (column in x$feature_list$labels) {
    col_class <- x$feature_list$classes[column]
    dt_out[, (column) := switch(col_class,
                                numeric = as.numeric(get(column)),
                                character = as.character(get(column)),
                                factor = as.factor(get(column)))]
  }
  dt_out[, row_id := NULL]
  setcolorder(dt_out, neworder =  c("id_combination", psi$meta$variable,"w", "id"))
  return(dt_out)
}

forge_wrapper <- function(params = psi, n_synth, datapoint, full_coals) {
  if (datapoint[, unique(id_combination)] == 1) {
    out <- full_coals[id == datapoint[,unique(id)]]
    out[, id_combination := datapoint[,unique(id_combination)]]
    out[, row_id := datapoint[,unique(row_id)]]
    out <- dcast(out, row_id + id + id_combination ~ variable, value.var = "value")
    return(out)
  }
  if (datapoint[, unique(id_combination == 16)]) {
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
