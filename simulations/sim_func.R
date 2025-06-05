#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances,
# which are used to support numerical simulation studies.
#####################################################################################

#####################################################################################
# Functions for generation of simulated data
#####################################################################################
generate.tensor = function(n, p.vec, Omega.list, mu=array(0,p.vec)){
  M = length(p.vec)
  Sigma.list = list()
  for ( m in 1:M) {
      Sigma.list[[m]] = solve(Omega.list[[m]])
  }

  kro.sigma = 1
  for (m in M:1){
    kro.sigma = kronecker(kro.sigma, Sigma.list[[m]])
  }
  
  ncols = ncol(kro.sigma)
  vecdata = matrix(rnorm(n * ncols), ncol = ncols) %*% chol(kro.sigma)
  Data = array(rep(mu,n), c(p.vec,n))+array(t(vecdata), c(p.vec,n))
  
  return(Data)
  
}

generate.Omega.A = function(t.Omega.true.list, A.orac, n, p.vec, K, 
                            c0=1, prob0=0.1, 
                            c1=10, prob1=0.5){
  M = length(t.Omega.true.list)
  Sigma.list = list()
  for ( m in 1:M) {
    Sigma.list[[m]] = solve(t.Omega.true.list[[m]])
  }
  A.Omega.true.list = list()
  
  for (k in 1:K) {
    A.Omega.true.list.k = list()
    if(k %in% A.orac){
      for (m in 1:M) {
        Sigma0.m = Sigma.list[[m]]
        pm = p.vec[m]
        s0 = max(apply(t.Omega.true.list[[m]], 2, function(x) sum(x!=0)))
        uni.up = c0 * sqrt( pm*log(pm) / prod(p.vec) / n )
        delta.m.k = matrix(rbinom(pm^2,size=1,prob=prob0)*runif(pm^2,-uni.up,uni.up), ncol=pm)
        Sig.m.k = Sigma0.m %*% (delta.m.k + diag(1,pm))
        Sig.m.k = (Sig.m.k+t(Sig.m.k))/2
        if(min(eigen(Sig.m.k)$values)<0.05){
            Sig.m.k = Sig.m.k + diag(0.1-min(eigen(Sig.m.k)$values),pm)
        }
        A.Omega.true.list.k[[m]] = solve(Sig.m.k)
      }
      A.Omega.true.list[[k]] = A.Omega.true.list.k
    } else {
      for (m in 1:M) {
        Sigma0.m = Sigma.list[[m]]
        pm = p.vec[m]
        s0 = max(apply(t.Omega.true.list[[m]], 2, function(x) sum(x!=0)))
        uni.up = c1 * s0 * sqrt( pm*log(pm) / prod(p.vec) / n )
        delta.m.k = matrix(rbinom(pm^2,size=1,prob=prob1)*runif(pm^2,-uni.up,uni.up), ncol=pm)
        Sig.m.k = Sigma0.m %*% (delta.m.k + diag(1,pm))
        Sig.m.k = (Sig.m.k+t(Sig.m.k))/2
        if(min(eigen(Sig.m.k)$values)<0.05){
          Sig.m.k = Sig.m.k + diag(0.1-min(eigen(Sig.m.k)$values),pm)
        }
        A.Omega.true.list.k[[m]] = solve(Sig.m.k)
      }
      A.Omega.true.list[[k]] = A.Omega.true.list.k
    }
    
    }
  
  return(A.Omega.true.list)
}

generate.data = function(M, t.Omega.true.list, A.Omega.true.list, nA.vec, K){
  if (M == 2){
    # tensor data of the target domain
    sigmaS = solve(t.Omega.true.list[[1]])
    sigmaT = solve(t.Omega.true.list[[2]])
    t.data = rmatnorm(n, sigmaS=sigmaS,sigmaT=sigmaT,method="chol")
    # tensor data of auxiliary domains
    A.data = list()
    for (k in 1:K) {
      sigmaS = solve(A.Omega.true.list[[k]][[1]])
      sigmaT = solve(A.Omega.true.list[[k]][[2]])
      A.data[[k]] = rmatnorm(nA.vec[k], sigmaS=sigmaS,sigmaT=sigmaT,method="chol")
    }
  }
  if (M != 2){
    # tensor data of the target domain
    t.data = generate.tensor(n, p.vec, t.Omega.true.list) 
    # tensor data of auxiliary domains
    A.data = list()
    for (k in 1:K) {
      A.data[[k]] = generate.tensor(nA.vec[k], p.vec, A.Omega.true.list[[k]]) 
    }
  }
  
  data.list = list(t.data=t.data, A.data=A.data)
  return(data.list)
}


################################################################################
# Functions for evaluating performances of proposed methods
################################################################################
evaluation = function(res.final, t.Omega.true.list, symmetric=F, test = F){
  
  # after model selection & auxiliary covariance matrices is weighted by sample sizes
  Omega.list = res.final$Omega.list
  # if(symmetric){Omega.list = res.final$Omega.sym.list}
  i.Omega = est.analysis(Omega.list, t.Omega.true.list)
  i.Omega = as.data.frame(t(unlist(i.Omega)))

  # after model selection & auxiliary covariance matrices is weighted by the differences with the target domain
  Omega.list.diff = res.final$Omega.list.diff
  # if(symmetric){Omega.list.diff = res.final$Omega.sym.list.diff}
  i.Omega.diff = est.analysis(Omega.list.diff, t.Omega.true.list)
  i.Omega.diff = as.data.frame(t(unlist(i.Omega.diff)))
  
  # before model selection & using only informative auxiliary covariance matrices
  Theta.list.o = res.final$res.trans.list$Theta.hat.list.o
  if(is.null(Theta.list.o)){
    i.Theta.o = NULL
  }else{
    i.Theta.o = est.analysis(Theta.list.o, t.Omega.true.list)
    i.Theta.o = as.data.frame(t(unlist(i.Theta.o)))
  }
  
  
  if(test){
    # before model selection & auxiliary covariance matrices is weighted by sample sizes
    Theta.list = res.final$res.trans.list$Theta.hat.list
    i.Theta = est.analysis(Theta.list, t.Omega.true.list)
    i.Theta = as.data.frame(t(unlist(i.Theta)))
    
    # before model selection & auxiliary covariance matrices is weighted by the differences with the target domain
    Theta.list.diff = res.final$res.trans.list$Theta.hat.diff.list
    i.Theta.diff = est.analysis(Theta.list.diff, t.Omega.true.list)
    i.Theta.diff = as.data.frame(t(unlist(i.Theta.diff)))
    
    # only based on the target domain using "sepa"
    sepa.Omega.list0 = res.final$res.trans.list$t.Omega.hat.list
    i.sepa.Omega0 = est.analysis(sepa.Omega.list0, t.Omega.true.list)
    i.sepa.Omega0 = as.data.frame(t(unlist(i.sepa.Omega0)))
    
    index.list = list(i.Omega=i.Omega, i.Omega.diff=i.Omega.diff, i.Theta.o=i.Theta.o,
                      i.Theta=i.Theta, i.Theta.diff=i.Theta.diff, i.sepa.Omega0=i.sepa.Omega0)
  } else {
    index.list = list(i.Omega=i.Omega, i.Omega.diff=i.Omega.diff, i.Theta.o=i.Theta.o)
  }

  
  return(index.list)
}

chg.eva = function(index.list){
  L = length(index.list)
  index.names = names(index.list[[1]])
  res.all = as.data.frame(matrix(0, nrow = L, ncol = length(index.names)))
  names(res.all) = index.names
  res.all.mean = res.all
  res.all.mean.sd = res.all
  for (l in 1:L) {
    index = index.list[[l]]
    res.all.mean[l,] = round(apply(index, 2, mean),3)
    res.all.mean.sd[l,] = paste0(round(apply(index, 2, mean),3), "(", round(apply(index, 2, sd),3), ")")
  }
  return(list(res.all.mean=res.all.mean, res.all.mean.sd=res.all.mean.sd))
}


eva.cha = function(index){
  M = length(index$error.f)
  index.form = as.data.frame(t(unlist(index)))
  return(index.form)
}


## Function: rmatnorm
## paper: Brain connectivity alteration detection via matrix‐variate differential network model, Biometrics, 2020
## R package: MVDN
## DOI: 10.1111/biom.13359
## Author: Jiadong Ji, Yong He, Lei Liu & Lei Xie.
## Download in https://github.com/jijiadong/MVDN
rmatnorm <- function (n, mean = matrix(0, nrow=nrow(sigmaS),ncol=ncol(sigmaT)),
                      sigmaS, sigmaT, method = c("eigen", "svd", "chol")){
  
  # vec(X) ~ N(vec(M), kron(sigmaT, sigmaS))
  # sigmaS Covariance matrix between rows
  # sigmaT Covariance matrix between columns
  
  nr <- nrow(sigmaS)
  nc <- ncol(sigmaT)
  
  if (!isSymmetric(sigmaS, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaS must be a symmetric matrix")
  }
  if (nrow(mean) != nr)
    stop("nrow of mean and nrow of sigmaS have non-conforming size")
  
  if (!isSymmetric(sigmaT, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaT must be a symmetric matrix")
  }
  if (ncol(mean) != nc)
    stop("ncol of mean and ncol of sigmaT have non-conforming size")
  
  mat = array( dim = c(nr, nc, n) )
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigmaT, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaT)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    R <- chol(sigmaT, pivot = TRUE)
    R <- R[, order(attr(R, "pivot"))]
  }
  
  ##
  if (method == "eigen") {
    ev <- eigen(sigmaS, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaS)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    Q <- chol(sigmaS, pivot = TRUE)
    Q <- Q[, order(attr(Q, "pivot"))]
  }
  
  for(i in 1:n){
    mat[,,i] = mean + t(Q) %*% matrix(rnorm(nr*nc), nrow=nr) %*% R
  }
  
  if(n==1){
    mat <- mat[,,1]
  }
  mat
}


