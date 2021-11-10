# DPMNIG
Dirichlet process mixture of MNIGs

"dpmnig_fun.R" contains the code for functions: dpMNIG(data,true_lab,maxiter,dp.alpha,verbs=FALSE,plots=FALSE,outputs) is the core function.
Input:
  data: data matrix
  true_lab: true class label
  maxiter: stop iteration if reached the maximum of iteration
  dp.alpha: numeric indicating a number as constant dp alpha parameter; or a vector(a,b) as updating along the algorithm with a gamma(a,b) prior.
  verbs: logic, print the chain numbers, iterations and current table of classifications or not
  plots: logic, trace plot or not
  outputs: "vector" or "lists", format of output estimated parameter
Output: 
  Parameter estimation for Mu, Beta, Gamma, Sigma
  final.G: final model numer of component selected 
  ARI 
  BIC 
  runtime
  
"dpmnig_run.R" contains an example of fitting Dirichlet process mixture of MNIGs to a dataset.