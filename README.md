This is an R package to estimate revealed preferences based on observed matchings in large bipartite graphs (large two-sided markets). The estimation for one-to-one matching uses a logit model with latent opportunity sets, where the Independence of Irrelevant Alternatives (IIA) assumption is satisfied. The estimation method is a re-parameterization of Menzel's model (Menzel 2015).

Dependencies: abind, nloptr, numDeriv, MASS

Exported functions: fitrpm_R_CP, Gale_Shapley, check_CP

Example: one_to_one_matching_CP.R is an example (a test driver) that sets the ground truth for the utility coefficients, along with various user-specified parameters (e.g. sample size, gender proportion, number of runs for each sample size). It then creates a simulated matching using the Gale-Shapley stable matching algorithm. Finally, it invokes the algorithm entry function in RPM, fitrpm_R_CP(), to estimate the utility coefficients. At the end, it prints out the median of the estimated values of utility coefficients, beta_med, and then reconstructs the joint density of the matching using the simulated data and the estimated utility coefficients from the last run. The estimated and the observed joint density matrices are printed at the end. The matrix "beta_sim" stores the following from each run: the estimated coefficients, the expected utility of the opportunity set for each type (the Gammas from Menzel's paper), the values returned from the equality constraints, and the status from the optimization procedure. 

Note: 
1) The estimates of the utility coefficients tend to be more accurate with roughly even proportion of men and women (up to around 20%/80%, and large sample size (total sample size >= 200).
2) TODO: add hypothesis testing using likelihood ratio. Currently, both the log-likelihood from the estimated and null models are returned from fitrpm_R_CP(), but they are not being utilized in the test driver yet.
3) TODO: all code is in R currently, but the likelihood and equality constraint functions will be implemented in C++.

Reference:
K Menzel. Large Matching Markets as Two-Sided Demand Systems. Econometrica, 83(3):987-941, 2015. URL https://bfi.uchicago.edu/sites/default/files/research/large_matching.pdf 

