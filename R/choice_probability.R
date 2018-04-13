logitProb_latent_opp_set = function(beta, GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric){
  
  Xdim = dim(Xd)
  nXu = Xdim[1];
  nZu = Xdim[2];
  W = Xdim[3];
  
  CPW = matrix(0, nrow=nXu, ncol=nZu+1) # women decision makers indexed by row
  CPM = matrix(0, nrow=nZu, ncol=nXu+1) # men decision makers indexed by row
  
  # single option
  for (i in 1:nXu) {
    CPW[i,(1+nZu)] = 1/(1+GammaW[i])
  }
  for (i in 1:nZu) {
    CPM[i,(1+nXu)] = 1/(1+GammaM[i])
  }
  
  if (symmetric) {
    beta_w = beta
    beta_m = beta
  } else {
    beta_w = beta[1:W]
    beta_m = beta[(W+1):length(beta)]
  }
  # married women
  for (j in 1:nXu) {
    for (k in 1:nZu) {
      Ustar = 0.0
      for (i in 1:length(beta_w)) {
        Ustar = Ustar + beta_w[i]*Xd[j,k,i]
      }
      CPW[j,k] = 1/sqrt(n) * exp(Ustar) / (1.0 + GammaW[j])
    }
  }
  
  # married men
  for (j in 1:nZu) {
    for (k in 1:nXu) {
      Vstar = 0.0
      for (i in 1:length(beta_m)) {
        Vstar = Vstar + beta_m[i]*Zd[j,k,i]
      }
      CPM[j,k] = 1/sqrt(n) * exp(Vstar) / (1.0 + GammaM[j])
    }
  }
  
  return(list(CPW = CPW, CPM = CPM))

}


