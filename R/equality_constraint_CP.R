equality_constraint_CP = function(beta, GammaW, GammaM, Xd, Zd, pmfW, pmfM, gw, gm, n, symmetric, sampling) {
  
  NumGammaW = length(GammaW)
  NumGammaM = length(GammaM)
  NumGamma =  NumGammaW + NumGammaM
  
  # # omp_set_num_threads(4);
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  sumVec = numeric(NumGamma)
  diffVec = numeric(NumGamma)
  
  # compute f(x_j,diamond) for all j women, where f(x_j,z_k) = w(x_j)*Px_jk * m(z_k)*Pz_kj
  for (j in 1:NumGammaW) {
    sum=0.0
    for (k in 1:NumGammaM) {
      sum = sum + CPW[j,k]*pmfM[k]*exp(gm) *  CPM[k,j]*pmfW[j]*exp(gw) * n
    }
    sumVec[j]=sum
  }

  # compute f(diamond, z_j) for all j men, where f(z_j,x_k) = m(z_j)*Pz_jk * w(x_k)*Px_kj
  for (j in 1:NumGammaM) {
    sum=0.0
    for (k in 1:NumGammaW) {
      sum = sum + CPM[j,k]*pmfW[k]*exp(gw) * CPW[k,j]*pmfM[j]*exp(gm) * n
    }
    sumVec[j+NumGammaW]=sum
  }

  # GammaW = f(x_j,diamond)/f(x_j,*), where f(x_j,*) = w(x_j)*Px_j0
  for (k in 1:NumGammaW) {
    diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - (1-CPW[k,1+NumGammaM])/CPW[k,1+NumGammaM] # works
    # diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - GammaW[k] # works
  }
  # GammaM = f(diamond, z_j)/f(*,z_j), where f(*,z_j) = m(z_j)*Pz_j0
  for (k in 1:NumGammaM) {
    diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - (1-CPM[k,1+NumGammaW])/CPM[k,1+NumGammaW] # works
    # diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - GammaM[k] # works
  }
  
  return(diffVec)
  
}

equality_constraint_null_CP = function(GammaW, GammaM, Xd, Zd, pmfW, pmfM, gw, gm, n, symmetric, sampling) {
  
  NumGammaW = length(GammaW)
  NumGammaM = length(GammaM)
  NumGamma =  NumGammaW + NumGammaM
  if (symmetric) {
    numBeta = dim(Xd)[3] # assume same dimension for both Xd and Zd (same)
  } else {
    numBeta = dim(Xd)[3]+dim(Zd)[3]
  }
  
  # # omp_set_num_threads(4);
  
  probs = logitProb_latent_opp_set(rep(0,numBeta), GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  sumVec = numeric(NumGamma)
  diffVec = numeric(NumGamma)
  
  # compute f(x_j,diamond) for all j women, where f(x_j,z_k) = w(x_j)*Px_jk * m(z_k)*Pz_kj
  for (j in 1:NumGammaW) {
    sum=0.0
    for (k in 1:NumGammaM) {
      sum = sum + CPW[j,1+NumGammaM]*pmfM[k]*exp(gm) *  CPM[k,1+NumGammaW]*pmfW[j]*exp(gw)
    }
    sumVec[j]=sum
  }
  
  # compute f(diamond, z_j) for all j men, where f(z_j,x_k) = m(z_j)*Pz_jk * w(x_k)*Px_kj
  for (j in 1:NumGammaM) {
    sum=0.0
    for (k in 1:NumGammaW) {
      sum = sum + CPM[j,1+NumGammaW]*pmfW[k]*exp(gw) * CPW[k,1+NumGammaM]*pmfM[j]*exp(gm)
    }
    sumVec[j+NumGammaW]=sum
  }
  
  # GammaW = f(x_j,diamond)/f(x_j,*), where f(x_j,*) = w(x_j)*Px_j0
  for (k in 1:NumGammaW) {
    diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - (1-CPW[k,1+NumGammaM])/CPW[k,1+NumGammaM] # works
    # diffVec[k] = sumVec[k]/(pmfW[k]*exp(gw)*CPW[k,1+NumGammaM]) - GammaW[k] # works
  }
  # GammaM = f(diamond, z_j)/f(*,z_j), where f(*,z_j) = m(z_j)*Pz_j0
  for (k in 1:NumGammaM) {
    diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - (1-CPM[k,1+NumGammaW])/CPM[k,1+NumGammaW] # works
    # diffVec[k+NumGammaW] = sumVec[k+NumGammaW]/(pmfM[k]*exp(gm)*CPM[k,1+NumGammaW]) - GammaM[k] # works
  }
  
  return(diffVec)
}

