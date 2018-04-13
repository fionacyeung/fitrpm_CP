loglikelihood_CP = function(beta, GammaW, GammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
  
  NumGammaW=length(GammaW)
  NumGammaM=length(GammaM)
  
  probs = logitProb_latent_opp_set(beta, GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  llik = 0.0
  
  
  if (sampling == "HOUSEHOLD") {
    llik=0.0
  
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        llik = llik + pmfj[j,k] * log(CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n) 
      }
    }
    
    pfps=0.0
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        pfps = pfps + CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n
      }
    }
    llik = llik - gw - gw + log(pfps)
    
    
  } else if (sampling == "COUPLE") {
    
    llik=0.0
    
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        llik = llik + pmfj[j,k] * log(CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n) 
      }
    }
    
    pfps=0.0
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        pfps = pfps + CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n
      }
    }
    llik = llik - log(pfps)
    
    
  } else { # assume sampling == "INDIV"
    
      llik=0.0
      
      # couples
      for (j in 1:NumGammaW) {
        for (k in 1:NumGammaM) {
          # llik = llik + pmfj[j,k] * log(CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n / (exp(gw)+exp(gm)))
          llik = llik + pmfj[j,k] * log(2*CPW[j,k]*pmfM[k]*exp(gm) * CPM[k,j]*pmfW[j]*exp(gw) * n / (exp(gw)+exp(gm)))
        }
      }
      
      # singles women (X each row)
      for (j in 1:NumGammaW) {
        # llik = llik + pmfj[j,(1+NumGammaM)] * log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw)/2 / (exp(gw)+exp(gm)))
        llik = llik + pmfj[j,(1+NumGammaM)] * log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw) / (exp(gw)+exp(gm)))
      }
      # single men (Z each col)
      for (k in 1:NumGammaM) {
        # llik = llik + pmfj[(1+NumGammaW),k] * log(CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm)/2 / (exp(gw)+exp(gm))) 
        llik = llik + pmfj[(1+NumGammaW),k] * log(CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm) / (exp(gw)+exp(gm))) 
      }
  }
  
  return (llik)
 
}

loglikelihood_null_CP = function(GammaW, GammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
  
  NumGammaW=length(GammaW)
  NumGammaM=length(GammaM)
  if (symmetric) {
    numBeta = dim(Xd)[3] # assume same dimension for both Xd and Zd (same)
  } else {
    numBeta = dim(Xd)[3]+dim(Zd)[3]
  }
  
  pfps=0.0
  llik=0.0
  
  probs = logitProb_latent_opp_set(rep(0,numBeta), GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric)
  CPW = probs$CPW
  CPM = probs$CPM
  
  if (sampling == "HOUSEHOLD") {
    llik=0.0
    
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        llik = llik + pmfj[j,k] * log(CPW[j,(1+NumGammaM)]*pmfM[k]*exp(gm) * CPM[k,(1+NumGammaW)]*pmfW[j]*exp(gw)) 
      }
    }
    
    pfps=0.0
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        pfps = pfps + CPW[j,(1+NumGammaM)]*pmfM[k]*exp(gm) * CPM[k,(1+NumGammaW)]*pmfW[j]*exp(gw)
      }
    }
    llik = llik - gw - gw + log(pfps)
    
    
  } else if (sampling == "COUPLE") {
    
    llik=0.0
    
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        llik = llik + pmfj[j,k] * log(CPW[j,(1+NumGammaM)]*pmfM[k]*exp(gm) * CPM[k,(1+NumGammaW)]*pmfW[j]*exp(gw)) 
      }
    }
    
    pfps=0.0
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        pfps = pfps + CPW[j,(1+NumGammaM)]*pmfM[k]*exp(gm) * CPM[k,(1+NumGammaW)]*pmfW[j]*exp(gw)
      }
    }
    llik = llik - log(pfps)
    
    
  } else { # assume sampling == "INDIV"
    
    # couples
    for (j in 1:NumGammaW) {
      for (k in 1:NumGammaM) {
        # llik = llik +  pmfj[j,k]*(log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw)*CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm) / (exp(gw)+exp(gm)) ));
        llik = llik +  pmfj[j,k]*(log(2*CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw)*CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm) / (exp(gw)+exp(gm)) ));
      }
    }
    
    # singles women (X each row)
    for (j in 1:NumGammaW) {
      # llik = llik + pmfj[j,(1+NumGammaM)] * log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw)/2 / (exp(gw)+exp(gm)) ) 
      llik = llik + pmfj[j,(1+NumGammaM)] * log(CPW[j,(1+NumGammaM)]*pmfW[j]*exp(gw) / (exp(gw)+exp(gm)) )
    }
    # single men (Z each col)
    for (k in 1:NumGammaM) {
      # llik = llik + pmfj[(1+NumGammaW),k] * log(CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm)/2 / (exp(gw)+exp(gm)) )
      llik = llik + pmfj[(1+NumGammaW),k] * log(CPM[k,(1+NumGammaW)]*pmfM[k]*exp(gm) / (exp(gw)+exp(gm)) ) 
    }
    
  }
  
  return (llik)
  
}



