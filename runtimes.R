library(microbenchmark)
library(doParallel)

# arf methods and no parallelisation:
registerDoSEQ()
microbenchmark(
  'MAP' = prepare_data.arf(x = explainer, psi = psi),
  'Sampling' = prepare_data.arf_sampling(x = explainer, psi = psi),
  times = 3L
)
# ~80s for MAP, ~90s for Sampling. 

arf_shaps <- explain(x_test,
                     explainer = explainer, 
                     approach = "arf_sampling", 
                     prediction_zero = mean(y_train$response), 
                     n_samples = 1000,
                     psi = psi)
