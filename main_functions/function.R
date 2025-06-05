################################################################################
# This document includes main functions of proposed methods,
# which are used to support numerical simulation studies and real data analysis
# in the paper
# ""
################################################################################


############################# Functions for main algorithms ####################
tensor.GGM.trans = function(t.data, A.data, A.lambda, A.orac = NULL, c=0.6,
                            t.lambda.int.trans=NULL, t.lambda.int.aggr=NULL,
                            theta.algm="cd", cov.select="inverse",
                            cov.select.agg.size = "inverse",
                            cov.select.agg.diff = "tensor.prod",
                            symmetric = T, init.method="sepa",
                            init.method.aux="sepa", mode.set = NULL,
                            init.iter.Tlasso=2, cn.lam2=seq(0.1,2,length.out =10),
                            c.lam.Tlasso=20, c.lam.sepa=20, adjust.BIC=F,
                            normalize = T, inti.the=T, sel.ind="fit"){

  p.vec = dim(t.data)
  M = length(p.vec) - 1
  if(M < 2){
    print("Error: M is less than 2!")
    # break
  }
  n = p.vec[M+1]
  p.vec = p.vec[1:M]

  if(is.null(mode.set)){
    mode.set = 1:M
  } else {
    mode.set = mode.set
  }

  # split the target data
  nc = n*c
  t.data.tran=0
  eval(parse(text=paste('t.data.tran=t.data[',paste(rep(',',M),collapse=''),'1:floor(nc)]')))
  t.data.aggr=0
  eval(parse(text=paste('t.data.aggr=t.data[',paste(rep(',',M),collapse=''),'setdiff(1:n,1:floor(nc))]')))


  res.trans = trans.estimate(t.data.tran, A.data, A.lambda, A.orac = A.orac,
                             t.lambda.int=t.lambda.int.trans, adjust.BIC=adjust.BIC,
                             mode.set = mode.set, init.iter=init.iter.Tlasso,
                             init.method=init.method, init.method.aux=init.method.aux,
                             theta.algm=theta.algm, cov.select=cov.select,
                             c.lam.sepa=c.lam.sepa, c.lam.Tlasso=c.lam.Tlasso,
                             cn.lam2=cn.lam2, normalize = normalize,
                             inti.the=inti.the)
  # initial precision matrices of all modes in the target domain
  t.Omega.hat.list = res.trans$t.Omega.hat.list

  # transfer learning-based precision matrices of all modes
  # using the weight (in Sigma.A.m) determined by the sample sizes
  Theta.hat.list.trans = res.trans$Theta.hat.list
  # using the weight (in Sigma.A.m) determined by the differences
  Theta.hat.list.trans.diff = res.trans$Theta.hat.diff.list

  if(sel.ind == "predict"){
    n.a = dim(t.data.aggr)[M+1]
    X = matrix(0, nrow = n.a, ncol = prod(p.vec))
    for (i in 1:n.a) {
      t.data.aggr.i=0
      eval(parse(text=paste('t.data.aggr.i=t.data.aggr[',paste(rep(',',M),collapse=''),'i]')))
      t.data.i.vec = as.vector(t.data.aggr.i)
      X[i,] = t.data.i.vec
    }
    Sigma.X = cov(X)

    KOmega.hat=1; KOmega.hat.diff=1; KOmega.hat0 = 1;
    for (m in M:1){
      KOmega.hat0 = kronecker(KOmega.hat0, t.Omega.hat.list[[m]])
      # Theta.hat using the weight (in Sigma.A.m) determined by the sample sizes
      KOmega.hat = kronecker(KOmega.hat, Theta.hat.list.trans[[m]])
      # Theta.hat using the weight (in Sigma.A.m) determined by the differences
      KOmega.hat.diff = kronecker(KOmega.hat.diff, Theta.hat.list.trans.diff[[m]])
    }

    pre0 = sum(diag(Sigma.X %*% KOmega.hat0)) - log(det(KOmega.hat0))
    pre1 = sum(diag(Sigma.X %*% KOmega.hat)) - log(det(KOmega.hat))
    pre1.diff = sum(diag(Sigma.X %*% KOmega.hat.diff)) - log(det(KOmega.hat.diff))

    # Theta.hat using the weight (in Sigma.A.m) determined by the sample sizes
    W.list = which.min(c(pre0, pre1))
    if(W.list == 2){
      Omega.hat.final.list = Theta.hat.list.trans
    } else {
      Omega.hat.final.list = t.Omega.hat.list
    }
    Omega.hat.final.sym.list = Omega.hat.final.list
    for (m in 1:M) {
      Omega.hat.final.sym.list[[m]] = Omega.hat.final.list[[m]]
    }

    # Theta.hat using the weight (in Sigma.A.m) determined by the differences
    W.list.diff = which.min(c(pre0, pre1.diff))
    if(W.list == 2){
      Omega.hat.final.list.diff = Theta.hat.list.trans.diff
    } else {
      Omega.hat.final.list.diff = t.Omega.hat.list
    }
    Omega.hat.final.sym.list.diff = Omega.hat.final.list.diff
    for (m in 1:M) {
      Omega.hat.final.sym.list.diff[[m]] = Omega.hat.final.list.diff[[m]]
    }

  }

  if(sel.ind == "fit"){
    # covariance matrices for aggregation & Theta.hat using the weight (in Sigma.A.m) determined by the sample sizes
    init.aggr = Initial.aggr(t.data.aggr, t.lambda.int=t.lambda.int.aggr,
                             method = init.method, cov.select=cov.select.agg.size,
                             c.lam.sepa=c.lam.sepa, c.lam.Tlasso=c.lam.Tlasso,
                             normalize = normalize)
    t.sigma.tilde.list = init.aggr$t.S.hat.list

    Omega.res.list = select.1(t.sigma.tilde.list, t.Omega.hat.list, Theta.hat.list.trans, mode.set)
    Omega.hat.final.list = Omega.res.list$Omega.hat.final.list
    Omega.hat.final.sym.list = Omega.res.list$Omega.hat.final.sym.list
    W.list = Omega.res.list$W.list

    # covariance matrices for aggregation & Theta.hat using the weight (in Sigma.A.m) determined by the differences
    init.aggr = Initial.aggr(t.data.aggr, t.lambda.int=t.lambda.int.aggr,
                             method = init.method, cov.select=cov.select.agg.diff,
                             c.lam.sepa=c.lam.sepa, c.lam.Tlasso=c.lam.Tlasso,
                             normalize = normalize)
    t.sigma.tilde.list = init.aggr$t.S.hat.list

    Omega.res.list.diff = select.1(t.sigma.tilde.list, t.Omega.hat.list, Theta.hat.list.trans.diff, mode.set)
    Omega.hat.final.list.diff = Omega.res.list.diff$Omega.hat.final.list
    Omega.hat.final.sym.list.diff = Omega.res.list.diff$Omega.hat.final.sym.list
    W.list.diff = Omega.res.list.diff$W.list
  }

  tensor.GGM.trans.res = list(Omega.list = Omega.hat.final.list,
                              Omega.sym.list = Omega.hat.final.sym.list,
                              Omega.list.diff = Omega.hat.final.list.diff,
                              Omega.sym.list.diff = Omega.hat.final.sym.list.diff,
                              res.trans.list = res.trans,
                              W.list=W.list, W.list.diff=W.list.diff)
  return(tensor.GGM.trans.res)
}

trans.estimate = function(t.data.tran, A.data, A.lambda, A.orac = NULL,
                          t.lambda.int=NULL, adjust.BIC=F, mode.set = NULL,
                          init.method="sepa", init.method.aux="sepa",
                          init.iter=3, normalize = T,
                          theta.algm="cd", cov.select="tensor.prod",
                          c.lam.sepa=20, c.lam.Tlasso=20, cn.lam2=1,
                          inti.the=T){

  p.vec = dim(t.data.tran)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  if(is.null(mode.set)){
    mode.set = 1:M
  } else {
    mode.set = mode.set
  }

  K = length(A.data)
  nA.vec = rep(0, K)
  for (k in 1:K) {
    p.vec.A = dim(A.data[[k]])
    nA.vec[k] = p.vec.A[length(p.vec.A)]
  }

  # Initialization
  if(is.null(t.lambda.int)){
    if(init.method == "sepa"){
      t.lambda.tran = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(init.method == "Tlasso"){
      t.lambda.tran = c.lam.Tlasso*sqrt( log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(M==2){
      t.lambda.tran = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
  }else{
    t.lambda.tran = t.lambda.int
  }

  init.res = Initial(t.data.tran, t.lambda.tran, A.data, A.lambda,
                     A.orac = A.orac, method = init.method, method.aux = init.method.aux,
                     TT=init.iter, normalize = normalize, mode.set = mode.set)
  init.time = init.res$time
  t.Omega.hat.list = init.res$t.Omega.hat.list


  if(cov.select=="inverse"){
    # 0 Covariance matrix: directly inverting precision matrix
    S.hat.A.list = init.res$S.hat.A.M0
    S.hat.A.diff.list = init.res$S.hat.A.M.diff0
  }
  if(cov.select=="tensor.prod"){
    # 1 Covariance matrix: multiplication by tensors and precision matrices
    S.hat.A.list = init.res$S.hat.A.M1
    S.hat.A.diff.list = init.res$S.hat.A.M.diff1
  }


  t1 = proc.time()
  # using the weight (in Sigma.A.m) determined by the differences
  Theta.hat.diff.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A = S.hat.A.diff.list[[mi]]
    Omega.hat0 = t.Omega.hat.list[[m]]
    lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
    delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)

    h.hat = max(apply(delta.hat, 2, function(x) sum(abs(x))))
    lam2.del = min(h.hat*sqrt(p.vec[m]*log(p.vec[m]) / n.da / prod(p.vec)), h.hat^2)
    lam2.N = sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec) / prod(p.vec))
    lambda2 = cn.lam2*min(max(lam2.del, lam2.N), n.da*lam2.N)

    Omega.hat00 = Omega.hat0 * inti.the + 0 * (!inti.the)
    Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat00,
                                    n.A=sum(nA.vec), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
    Theta.hat.m = Theta.tuning.res$Theta.hat.m
    Theta.hat.diff.list[[mi]] = Theta.hat.m
  }

  ## using the weight (in Sigma.A.m) determined by the sample sizes
  Theta.hat.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A = S.hat.A.list[[mi]]
    Omega.hat0 = t.Omega.hat.list[[m]]
    lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
    delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)

    h.hat = max(apply(delta.hat, 2, function(x) sum(abs(x))))
    lam2.del = min(h.hat*sqrt(p.vec[m]*log(p.vec[m]) / n.da / prod(p.vec)), h.hat^2)
    lam2.N = sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec) / prod(p.vec))
    lambda2 = cn.lam2*min(max(lam2.del, lam2.N), n.da*lam2.N)

    Omega.hat00 = Omega.hat0 * inti.the + 0 * (!inti.the)
    Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat00,
                                    n.A=sum(nA.vec), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
    Theta.hat.m = Theta.tuning.res$Theta.hat.m
    Theta.hat.list[[mi]] = Theta.hat.m
  }
  t.theta = proc.time() - t1

  if(sum(A.orac) > 0){
    ########### oracle auxiliary domains
    if(cov.select=="inverse"){
      S.hat.A.list.o = init.res$S.hat.A.M0.o   # 0 Covariance matrix: directly inverting precision matrix
    }
    if(cov.select=="tensor.prod"){
      S.hat.A.list.o = init.res$S.hat.A.M1.o   # 1 Covariance matrix: multiplication by tensors and precision matrices
    }

    Theta.hat.list.o = list()
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A = S.hat.A.list.o[[mi]]
      Omega.hat0 = t.Omega.hat.list[[m]]
      lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
      delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)

      lambda2 = cn.lam2*sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec[A.orac]) / prod(p.vec))
      Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat0,
                                      n.A=sum(nA.vec[A.orac]), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
      Theta.hat.m = Theta.tuning.res$Theta.hat.m
      Theta.hat.list.o[[mi]] = Theta.hat.m
    }
    res.trans = list(Theta.hat.list=Theta.hat.list, Theta.hat.diff.list=Theta.hat.diff.list,
                     S.hat.A.list=S.hat.A.list, S.hat.A.diff.list=S.hat.A.diff.list,
                     t.Omega.hat.list=t.Omega.hat.list, init.res=init.res,
                     init.time=init.time, theta.time=t.theta,
                     Theta.hat.list.o=Theta.hat.list.o, S.hat.A.list.o=S.hat.A.list.o)
  } else {
    res.trans = list(Theta.hat.list=Theta.hat.list, Theta.hat.diff.list=Theta.hat.diff.list,
                     S.hat.A.list=S.hat.A.list, S.hat.A.diff.list=S.hat.A.diff.list,
                     t.Omega.hat.list=t.Omega.hat.list, init.res=init.res,
                     init.time=init.time, theta.time=t.theta)
  }

  return(res.trans)

}

Initial = function(t.data, t.lambda, A.data, A.lambda,
                   A.orac = NULL, method = "sepa", method.aux = "sepa",
                   TT=2, normalize = T, mode.set = NULL){
  # Initial: the function calculating initial precision matrices of the target domain
  #          and covariance matrices of the auxiliary domain,
  #          via two alternative methods:
  #          "Tlasso" (PAMI, 2020) & "sepa" (JCGS, 2022)

  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  K = length(A.data)
  nA.vec = rep(0, K)
  for (k in 1:K) {
    p.vec.A = dim(A.data[[k]])
    nA.vec[k] = p.vec.A[length(p.vec.A)]
  }
  if(is.null(mode.set)){
    mode.set = 1:M
  } else {
    mode.set = mode.set
  }

  t0 <- proc.time()
  if(method == "sepa"){
    # Initialization in target domain
    t.Omega.hat.list = Separate.fit(t.data, lambda.vec=t.lambda, normalize = normalize)$Omegahat
  }

  if(method == "Tlasso"){
    # Initialization in target domain
    t.Omega.hat.list = Tlasso.fit(t.data, T=TT, lambda.vec = t.lambda, norm.type = 1+as.numeric(normalize))
  }

  if(method.aux == "sepa"){
    # Initialization in auxiliary domains
    A.Omega.hat.list = list()
    for (k in 1:K) {
      A.Omega.hat.list[[k]] = Separate.fit(A.data[[k]], lambda.vec=A.lambda[[k]], normalize = normalize)$Omegahat
    }
  }
  if(method.aux == "Tlasso"){
    # Initialization in auxiliary domains
    A.Omega.hat.list = list()
    for (k in 1:K) {
      A.Omega.hat.list[[k]] = Tlasso.fit(A.data[[k]], lambda.vec = A.lambda[[k]], norm.type = 1+as.numeric(normalize))
    }
  }

  A.S.hat.list0 = list()
  A.S.hat.list1 = list()
  for (k in 1:K) {
    A.S.hat.list = S.est(A.data[[k]], A.Omega.hat.list[[k]])
    A.S.hat.list0[[k]] = A.S.hat.list$sig0
    A.S.hat.list1[[k]] = A.S.hat.list$sig1
  }

  S.hat.A.M = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A.M[[mi]] = diag(p.vec[m]) - diag(p.vec[m])
  }
  # weight determined by the differences
  S.hat.A.M.diff0 = S.hat.A.M
  S.hat.A.M.diff1 = S.hat.A.M
  weight.KM0 = matrix(0, ncol = K, nrow = length(mode.set))
  weight.KM1 = matrix(0, ncol = K, nrow = length(mode.set))
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      weight.KM0[mi,k] = 1/sum((A.S.hat.list0[[k]][[m]] %*% t.Omega.hat.list[[m]] - diag(p.vec[m]))^2)
      weight.KM1[mi,k] = 1/sum((A.S.hat.list1[[k]][[m]] %*% t.Omega.hat.list[[m]] - diag(p.vec[m]))^2)
    }
  }

  wed0 = t(t(weight.KM0)*nA.vec)
  alpha.k.diff0 = wed0 / apply(wed0, 1, sum)
  wed1 = t(t(weight.KM1)*nA.vec)
  alpha.k.diff1 = wed1 / apply(wed1, 1, sum)
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A.M.diff0[[mi]] = S.hat.A.M.diff0[[mi]] + A.S.hat.list0[[k]][[m]] * alpha.k.diff0[mi,k]
      S.hat.A.M.diff1[[mi]] = S.hat.A.M.diff1[[mi]] + A.S.hat.list1[[k]][[m]] * alpha.k.diff1[mi,k]
    }
  }


  # weight determined by the sample sizes
  S.hat.A.M0 = S.hat.A.M
  S.hat.A.M1 = S.hat.A.M
  alpha.k = nA.vec/sum(nA.vec)
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A.M0[[mi]] = S.hat.A.M0[[mi]] + A.S.hat.list0[[k]][[m]] * alpha.k[k]
      S.hat.A.M1[[mi]] = S.hat.A.M1[[mi]] + A.S.hat.list1[[k]][[m]] * alpha.k[k]
    }
  }
  t00 = proc.time() - t0


  if(sum(A.orac) > 0){
    S.hat.A.M0.o = S.hat.A.M
    S.hat.A.M1.o = S.hat.A.M
    alpha.k.o = nA.vec[A.orac]/sum(nA.vec[A.orac])
    A.S.hat.list0.o = list()
    A.S.hat.list1.o = list()
    for (k in A.orac) {
      ki = match(k, A.orac)
      A.S.hat.list.o = S.est(A.data[[k]], A.Omega.hat.list[[k]])
      A.S.hat.list0.o[[ki]] = A.S.hat.list.o$sig0
      A.S.hat.list1.o[[ki]] = A.S.hat.list.o$sig1
      for (m in mode.set) {
        mi = match(m, mode.set)
        S.hat.A.M0.o[[mi]] = S.hat.A.M0.o[[mi]] + A.S.hat.list0.o[[ki]][[m]] * alpha.k.o[ki]
        S.hat.A.M1.o[[mi]] = S.hat.A.M1.o[[mi]] + A.S.hat.list1.o[[ki]][[m]] * alpha.k.o[ki]
      }
    }
    Init.res = list(t.Omega.hat.list=t.Omega.hat.list, A.Omega.hat.list=A.Omega.hat.list,
                    S.hat.A.M0=S.hat.A.M0, S.hat.A.M1=S.hat.A.M1,
                    S.hat.A.M.diff0=S.hat.A.M.diff0, S.hat.A.M.diff1=S.hat.A.M.diff1,
                    A.S.hat.list0=A.S.hat.list0, A.S.hat.list1=A.S.hat.list1,
                    S.hat.A.M0.o=S.hat.A.M0.o, S.hat.A.M1.o=S.hat.A.M1.o,
                    A.S.hat.list0.o=A.S.hat.list0.o, A.S.hat.list1.o=A.S.hat.list1.o,
                    time = t00)
    return(Init.res)
  } else {
    Init.res = list(t.Omega.hat.list=t.Omega.hat.list, A.Omega.hat.list=A.Omega.hat.list,
                    S.hat.A.M0=S.hat.A.M0, S.hat.A.M1=S.hat.A.M1,
                    S.hat.A.M.diff0=S.hat.A.M.diff0, S.hat.A.M.diff1=S.hat.A.M.diff1,
                    A.S.hat.list0=A.S.hat.list0, A.S.hat.list1=A.S.hat.list1,
                    time = t00)
    return(Init.res)
  }

}

Initial.aggr = function(t.data, t.lambda.int=NULL, method = "sepa",
                        cov.select= "tensor.prod", TT=2,
                        c.lam.sepa=20, c.lam.Tlasso=20, normalize = T){
  # Initial.aggr: the function calculating initial covariance matrices of
  #               the target domain for the aggregation step,
  #               via two alternative methods:
  #               "Tlasso" (PAMI, 2020) & "sepa" (JCGS, 2022)

  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]

  if(is.null(t.lambda.int)){
    if(method == "sepa"){
      t.lambda = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(method == "Tlasso"){
      t.lambda = c.lam.Tlasso*sqrt( log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(M==2){
      t.lambda = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
  }else{
    t.lambda = t.lambda.int
  }


  if(method == "sepa"){
    # Initialization in target domain
    t.Omega.hat.list = Separate.fit(t.data, lambda.vec=t.lambda, normalize=normalize)$Omegahat
  }

  if(method == "Tlasso"){
    # Initialization in target domain
    t.Omega.hat.list = Tlasso.fit(t.data, T=TT, t.lambda, norm.type = 1+as.numeric(normalize))
  }

  t.S.hat.list = S.est(t.data, t.Omega.hat.list)

  t.S.hat.list0 = t.S.hat.list$sig0  # 0 Covariance matrix: directly inverting precision matrix
  t.S.hat.list1 = t.S.hat.list$sig1  # 1 Covariance matrix: multiplication by tensors and precision matrices

  if(cov.select=="inverse"){
    t.S.hat.list = t.S.hat.list0
  }
  if(cov.select=="tensor.prod"){
    t.S.hat.list = t.S.hat.list1
  }

  Init.res = list(t.Omega.hat.list=t.Omega.hat.list,
                  t.S.hat.list=t.S.hat.list)



}

Theta.tuning = function(lambda2, S.hat.A, delta.hat, Omega.hat0, n.A,
                        theta.algm="cd", adjust.BIC=F){

  L                  = length(lambda2)
  aBIC               = rep(0,L)
  Theta.list         = list()
  BIC.summary        = as.data.frame(matrix(0,ncol = 5, nrow = L))
  names(BIC.summary) = c("lambda", "BIC", "fitness", "BIC.penalty", "degree")

  for (l in 1:L) {
    lam2 = lambda2[l]
    Theta.hat.m = Theta.est(S.hat.A, delta.hat, lam2=lam2,
                            Omega.hat0=Omega.hat0, method = theta.algm)
    BIC.list = BIC.value(S.hat.A, delta.hat, Theta.hat.m, n.A, adjust.BIC)
    BIC.summary[l,] = c(lam2, unlist(BIC.list))
    Theta.list[[l]] = Theta.hat.m
  }
  Theta.hat.m.opt = Theta.list[[which.min(BIC.summary$BIC)]]
  Theta.tuning.res = list(Theta.hat.m=Theta.hat.m.opt, BIC.summary=BIC.summary,
                          Theta.hat.list.m=Theta.list)
  return(Theta.tuning.res)
}

S.est = function(data, Omega.hat.list){
  # S.est: the function calculating the covariance matrix of each mode

  p.vec = dim(data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]

  Omega.hat.list.sqrt = list()
  S.hat.list0 = list()
  for (m in 1:M) {
    Omega.hat.list.sqrt[[m]] = expm::sqrtm(Omega.hat.list[[m]])
    S.hat.list0[[m]] = solve(Omega.hat.list[[m]])
  }

  S.hat.list1 = list()
  for(m in 1:M){
    S.array = array(0,c(p.vec[m],p.vec[m],n.da))
    Omega.hat.list.sqrt.m = Omega.hat.list.sqrt
    Omega.hat.list.sqrt.m[[m]] = diag(p.vec[m])
    for(i in 1:n.da){
      d=0
      eval(parse(text=paste('d=data[',paste(rep(',',M),collapse=''),'i]')))
      Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.hat.list.sqrt.m , ms=1:M)@data) ,m=m)@data
      S.array[,,i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array,c(1,2),mean) * p.vec[m] / prod(p.vec)
    S.hat.list1[[m]] = S.mat
  }
  S.hat.list = list(sig0=S.hat.list0, sig1=S.hat.list1)
  return(S.hat.list)
}

delta.est = function(S.hat.A, Omega.hat0, lam1){
  # delta.est: the function estimating divergence matrix (Delta_m) of the mode
  #            corresponding to S.hat.A.
  pm = dim(S.hat.A)[1]
  B.hat = Omega.hat0 %*% S.hat.A - diag(pm)
  Z = array(rep(0, 2*pm^2), dim=c(pm,pm,2))
  Z[,,2] = abs(B.hat) - lam1
  delta.hat = sign(B.hat) * apply(Z, 1:2, max)
  return(delta.hat)
}

Theta.est = function(S.hat.A, delta.hat, lam2=0.1, Omega.hat0=NULL,
                     n=100, max_iter=10, eps=1e-3, method = "admm"){
  # Theta.est: the function estimating transfer learning-based estimator of
  #            precision matrix of the mode corresponding to S.hat.A, via
  #            two algorithms:
  #            ADMM (method = "admm") & Coordinate descent (method = "cd").
  pm = dim(S.hat.A)[1]
  deltaI = delta.hat + diag(pm)

  if(sum(abs(Omega.hat0)) == 0){
    if(abs(det(S.hat.A)) < 1e-3){
      Omega.hat0 = solve(S.hat.A + diag(pm)/n) %*% deltaI
    } else {
      Omega.hat0 = solve(S.hat.A) %*% deltaI
    }
  }

  if(method == "cd"){
    Theta.hat = Thetaest.cd(S.hat.A, deltaI, lam2, Omega.hat0, max_iter, eps)
  }
  if(method == "admm"){
    Theta.hat = Thetaest.admm(S.hat.A, deltaI, lam2, Omega.hat0, max_iter, eps)
  }

  return(Theta.hat)
}

Thetaest.cd = function(S.hat.A, deltaI, lam2, Omega.hat0, max_iter=10, eps=0.001){
  # Thetaest.cd: the function estimating transfer learning-based estimator of
  #              precision matrix of the mode corresponding to S.hat.A, via
  #              coordinate descent algorithm.
  p = dim(S.hat.A)[1]
  Theta_hat = Omega.hat0
  for (j in 1:p){
    thetaj = Omega.hat0[,j]
    iter = 0
    diff = 10
    while(iter < max_iter && diff > eps){
      thetaj0 = thetaj
      for (i in 1:p){
        Sj = S.hat.A[i,]
        thetaji = deltaI[i,j] - Sj %*% thetaj + Sj[i] * thetaj[i]
        if(i == j){
          thetaj[i] = S_soft(thetaji, 0) / Sj[i]
        }else{
          thetaj[i] = S_soft(thetaji, lam2) / Sj[i]
        }
      }
      diff = sqrt( sum((thetaj - thetaj0)^2) / p )
      iter = iter + 1
      # print(c(iter,diff))
    }
    Theta_hat[,j] = thetaj
  }
  return(Theta_hat)
}

Thetaest.admm = function(S.hat.A, deltaI, lam2, Omega.hat0,
                         max_iter=10, eps=1e-3, kappa = 1){
  # Thetaest.cd: the function estimating transfer learning-based estimator of
  #              precision matrix of the mode corresponding to S.hat.A, via
  #              ADMM algorithm.
  p = dim(S.hat.A)[1]
  SI = S.hat.A+kappa*diag(p)
  Theta_hat = Omega.hat0

  for (j in 1:p){
    thetaj = Omega.hat0[,j]
    v = thetaj
    ej = rep(1,p)
    ej[j] = 0
    gamma = rep(0,p)
    iter = 0
    diff = 10
    while(iter < max_iter && diff > eps){
      thetaj0 = thetaj
      thetaj = as.numeric(solve(SI) %*% ( deltaI[,j] + gamma + kappa*v ))
      v = S_soft.vec(thetaj - gamma/kappa, lam2, ej)
      gamma = gamma + kappa * ( v - thetaj )

      diff = sqrt( sum((thetaj - thetaj0)^2) / p )
      iter = iter + 1
    }
    Theta_hat[,j] = v
  }
  return(Theta_hat)
}

select.1 = function(t.sigma.tilde.list, t.Omega.hat.list, Theta.hat.list, mode.set, symmetric=T){
  Omega.hat.final.list = list()
  Omega.hat.final.sym.list = list()
  W.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    pm = p.vec[m]
    t.sigma.tilde.m = t.sigma.tilde.list[[m]]
    Theta.trans.m = Theta.hat.list[[mi]]
    t.Omega.hat.m = t.Omega.hat.list[[m]]

    resid.Omega0 = apply((t.sigma.tilde.m %*% t.Omega.hat.m - diag(pm))^2, 2, sum)
    resid.trans = apply((t.sigma.tilde.m %*% Theta.trans.m - diag(pm))^2, 2, sum)
    w.m = apply(rbind(resid.Omega0, resid.trans),2,which.min)
    Omega.final.m = t( t(t.Omega.hat.m) * c(w.m == 1) + t(Theta.trans.m) * c(w.m == 2) )

    Omega.hat.final.list[[mi]] = Omega.final.m
    W.list[[mi]] = w.m # 1: only target; 2: transfer

    if(symmetric){
      Omega.hat.final.sym.list[[mi]] = symmetric.mat(Omega.final.m)}
  }
  selec = list(Omega.hat.final.list = Omega.hat.final.list,
               Omega.hat.final.sym.list = Omega.hat.final.sym.list,
               W.list = W.list)
  return(selec)

}

############################# Some fundamental supporting functions ############################
symmetric.mat = function(Omega){
  pm = dim(Omega)[1]
  Z = array(rep(0, 2*pm^2), dim=c(pm,pm,2))
  Z[,,1] = Omega
  Z[,,2] = t(Omega)
  Omega.sym = apply(Z, 1:2, function(x) x[which.min(abs(x))])
  return(Omega.sym)
}

S_soft = function(z,lambda){
  # S_soft: single lasso shrinkage estimate
  norm.z = sqrt(sum(z^2))
  if(norm.z!=0){
    n.x = 1 - lambda/norm.z
    rho = n.x*(n.x > 0)*z
  } else{
    rho = z
  }
  return(rho)
}

S_soft.vec = function(z,lambda,ej=rep(1,length(z))){
  # S_soft.vec: single lasso shrinkage estimate for a vector
  n.z = abs(z) - lambda*ej
  return(sign(z) * (n.z > 0) * n.z)
}

BIC.value = function(S.hat.A, delta.hat, Theta.hat, n=100, adjust=F){
  pm = dim(S.hat.A)[1]
  deltaI = delta.hat + diag(pm)
  fitness = 0.5*sum(diag(t(Theta.hat) %*% S.hat.A %*% Theta.hat)) - sum(diag( t(deltaI) %*%  Theta.hat))
  degree = sum(Theta.hat != 0) - pm
  if(adjust){Cn = log(n*p)} else {Cn = 1}
  BIC.penalty = Cn*degree*log(n)  / n
  BICvalue = fitness + BIC.penalty

  return(list(BIC=BICvalue, fitness=fitness, BIC.penalty=BIC.penalty,degree=degree))
}










## Function: Separate.fit
## paper: Fast and Separable Estimation in High-Dimensional Tensor Gaussian Graphical Models, JCGS, 2022
## DOI: 10.1080/10618600.2021.1938086
## Author: Keqian Min, Qing Mai & Xin Zhang
## Download in https://www.tandfonline.com/doi/suppl/10.1080/10618600.2021.1938086?scroll=top
#### This function estimates the precision matrices in the tensor graphical model using the proposed parallel scheme.
# x: p1*p2*...*pM*n
# val: (Optional) validation set. If supplied, lambda.list should be provided
# est.mode: index set of precision matrices to be estimated. If not specified, all precision matrices will be estimated
# lambda.vec: the sequence of regularization parameters for each mode in est.mode. It is used when val is missing.
# lambda.list: A list of regularization parameters that provides a lambda sequence for each mode in `est.mode`. Must be supplied with validation set. When a validation set is supplied, the optimal
# tuning parameters will be chosen from lambda.list based on the log-likelihood calculated using validation set.
# Omegatilde.list: (Optional) a list of M matrices
# scale.vec: constants to scale the log-likelihood to avoid infinite value. Default is 1 for all modes
# normalize: indicates whether $\widetilde{\boldsymbol{\Omega}}_m$ and $\widehat{\boldsymbol{\Omega}}_m$ should be normalized to have unit Frobenius norm. Default is TRUE.
# thres: threshold for convergence. Default value is 1e-4
# maxit: maximum number of iterations for fitting glasso. Default value is 10,000
# njobs: number of nodes used to do parallel computing
Separate.fit = function(x, val = NULL, est.mode = NULL, lambda.vec = NULL, lambda.list = NULL, Omegatilde.list = NULL, scale.vec = NULL, normalize = TRUE, thres = 1.0e-4, maxit = 1e4, njobs = 4) {

  dimen = dim(x) # dimension of dataset
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n_val = dim(val)[K + 1] # sample size of validation set
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation

  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  if (is.null(val)) {
    if (is.null(lambda.vec) | length(lambda.vec) != length(est.mode)) {
      stop("lambda.vec is missing or does not have the correct length")
    }
  } else {
    if (is.null(lambda.list) | length(lambda.list) != length(est.mode)) {
      stop("lambda.list is missing or does not have the correct length")
    }
  }



  if (!(is.null(Omegatilde.list) | length(Omegatilde.list) == K)) {
    stop("argument Omegatilde.list should be a list of M matrices")
  }

  if (is.null(scale.vec)) {
    scale.vec = rep(1, length(est.mode))
  }

  ##### Calculate \tilde\Omega #####
  if (is.null(Omegatilde.list) == FALSE) {
    fit1 = Omegatilde.list  # user-specified value for \tilde\Omega
  } else {
    # Calculate \tilde\Omega by the definition in the paper
    c1 = makeCluster(njobs)
    registerDoParallel(c1)
    fit1 = foreach(k = 1:K, .export = c("x"), .combine = list, .multicombine = TRUE) %dopar% {
      # when sample size is small, use the identity matrix
      if (n * nvars < ((dimen[k]**2) * (dimen[k] - 1) / 2)) {
        Omega_tilde = diag(dimen[k])
      }
      else {
        # when sample size is large, calculate the sample estimator of the precision matrices
        S.array = array(0, c(m.vec[k], m.vec[k], n))
        for (i in 1:n) {
          d = 0
          eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
          Vi = rTensor::k_unfold(rTensor::as.tensor(d), m = k)@data  # unfold tensor
          S.array[, , i] = Vi %*% t(Vi)
        }
        S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # sample estimation of \Sigma_k
        Omega_tilde = solve(S.mat)
        # normalization
        if (normalize) {
          Omega_tilde = Omega_tilde / norm(Omega_tilde, type = "F")}
      }
      Omega_tilde
    }
    stopCluster(c1)
  }

  K1 = length(est.mode) # number of precision matrices to be estimated
  lam.best = rep(0, K1)
  loglik = list()

  ###### Tuning process ######
  # When validation set is supplied, the lambdas with the maximum log-likelihood will be chosen
  if (!(is.null(val))) {
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = fit1[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt

    for (mode_index in 1:K1) {
      k = est.mode[mode_index]

      # Calculate \tilde S_k using the training set
      S.array = array(0, c(m.vec[k], m.vec[k], n))
      Omega.list.sqrt[[k]] = diag(m.vec[k]) # set \tilde\Omega_k to identity matrix
      for (i in 1:n) {
        d = 0
        eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        S.array[, , i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k

      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]

      # fit model with a sequence of lambdas
      lamk = lambda.list[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      for (i in 1:lam.length) {
        Out1 = glasso(S.mat, rho = lamk[i], penalize.diagonal = FALSE, maxit = 1e4, thr = 1.0e-4)
        hat_Omega = Out1$wi
        loglik2[i] = -tr(testS.mat %*% hat_Omega) + log(det(hat_Omega * scale.vec[mode_index]))
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      ind = which.max(loglik2)
      lam.best[mode_index] = lamk[ind] # get the optimal lambda that maximizes the log-likelihood
      loglik[[mode_index]] = loglik2
    }
  } else {
    # if validation set is not provided, directly use lambda.vec to fit model
    lam.best = lambda.vec
  }


  ##### Model fitting using parallel computing #####
  # register cluster for parallel computing
  c1 = makeCluster(njobs)
  registerDoParallel(c1)
  K1 = length(est.mode)
  fit_result = foreach(mode_ind = 1:K1, .packages = c("glasso", "rTensor", "expm"), .export = c("x"), .combine = list, .multicombine = TRUE) %dopar% {
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(fit1[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    Out1 = glasso(S.mat, rho = lam.best[mode_ind], penalize.diagonal = FALSE, maxit = maxit, thr = thres)
    hat_Omega = as.matrix(Out1$wi)
    # normalization
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    hat_Omega
  }
  stopCluster(c1)

  result = list()
  result$Omegahat = fit_result
  result$lambda = lam.best
  result$loglik = loglik
  return(result)
}

