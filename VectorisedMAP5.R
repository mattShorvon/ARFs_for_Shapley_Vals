
# compared to 4: eliminated another for loop and used ifelse once more, uses fifelse line 21.
prepare_data.arf <- function(x = explainer, index_features = NULL, psi = psi , ...) {
  # start_t = Sys.time()
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
  # print(Sys.time() - start_t) # around 0.035 seconds on old laptop
  dt[, row_id := .GRP, by = c("id", "id_combination")]
  dt[, mu_cvg_sum := as.numeric(NA)]
  dt[, cvg_sum := as.numeric(NA)]

  to_eval <- unique(dt[id_combination != 1 & id_combination != 2^ncol(x$x_test), row_id])
  dt[row_id %in% to_eval,] <- foreach(i = 1:length(to_eval), .combine = rbind, .packages = "data.table", .export = "prep") %dopar% {
    prep(dt[row_id == to_eval[i],], psi)
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

prep <- function(coalition, psi) {
  # Find leaves that satisfy constraints from the input coalition:
  coal <- coalition[!is.na(value)]
  leafIDs = as.list(rep(NA, nrow(coal)))
  leaves = numeric()
  for (i in sequence(nrow(coal))) {
    ifelse(coal[i, type] == "cnt", {
      leafIDs[[i]] <- psi$cnt[variable == coal[i, variable] & min < as.numeric(coal[i, value]) & max > as.numeric(coal[i, value]), f_idx]},
      {leafIDs[[i]] <- psi$cat[variable == coal[i, variable] & val == coal[i, value] & prob > 0, f_idx]}
      )
  } # identical came out true when compared to finalMAP version.
  leaves <- Reduce(intersect, leafIDs)
  if (length(leaves) == 0) {
    leaves <- unique(psi$cnt[,f_idx])
  }
  # Find the continuous data info for leaves that satisfy these constraints, add it to coalition
  # start_t <- Sys.time()
  if (nrow(coalition[is.na(value) & type == "cnt"]) > 0) {
    leaves_info_cnt <- merge(psi[["cnt"]][f_idx %in% leaves & variable %in% coalition[, unique(variable)], .(variable, mu, f_idx)],
                             psi$forest[f_idx %in% leaves, .(cvg, f_idx)],by = 'f_idx', sort = F) # set sort to false, you don't need to sort by f_idx, by call in the merge!

    leaves_info_cnt[, mu_x_cvg := mu*cvg]
    # got rid of the for loop here, replacing with a data table operation.
    coalition[is.na(value) & type == "cnt", cvg_sum := leaves_info_cnt[variable == .SD[, variable], sum(cvg),
                                                                       by = variable]$V1, by = variable, .SDcols = "variable"]
    coalition[is.na(value) & type == "cnt", mu_cvg_sum := leaves_info_cnt[variable == .SD[, variable], sum(mu_x_cvg),
                                                                          by = variable]$V1, by = variable, .SDcols = "variable"]
  }
  # Find the categorical data info for leaves that satisfy these constraints, store in the cat_prep table.
  if (nrow(coalition[is.na(value) & type == "cat"]) > 0) {
    leaves_info_cat <- merge(psi[["cat"]][f_idx %in% leaves & variable %in% coalition[, unique(variable)], .(variable, val, prob, f_idx)],
                             psi$forest[f_idx %in% leaves, .(cvg, f_idx)],
                             by = 'f_idx', sort = F)
    leaves_info_cat[, prob_x_cvg := prob*cvg]
    leaves_info_cat[, cvg_sum := sum(cvg), by = .(val, variable)]
    leaves_info_cat[, prob_cvg_sum := sum(prob_x_cvg ), by = .(val, variable)]
    leaves_info_cat[, division := prob_cvg_sum/cvg_sum]
    # for (var in coalition[is.na(value) & type == "cat", variable]) {
    #   # print(var)
    #   coalition[variable == var, value := unique(leaves_info_cat[variable == var][division == max(division), val])]
    # }
    coalition[is.na(value) & type == "cat", value := unique(leaves_info_cat[variable == .SD[, variable]]
                                                            [division == max(division), val])[1], by = variable, .SDcols = "variable"]
  }
  # print(Sys.time() - start_t)
  # store the coalition and cat_prep
  return(coalition)
}
