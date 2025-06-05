all.names = c("KKI","NeuroIMAGE","Peking","Pittsburgh","NYU","OHSU","WashU")
methods = c("Tlasso", "proposed", "proposed.v")
##### summary of sample size #####
A.data = fMRI.data.list$Typical
n.vec.Typ = rep(0, 7)
for (k in 1:length(A.data)) {
  p.vec.A = dim(A.data[[k]])
  n.vec.Typ[k] = p.vec.A[length(p.vec.A)]
}
A.data = fMRI.data.list$ADHD
n.vec.ADHD = rep(0, 7)
for (k in 1:length(A.data)) {
  p.vec.A = dim(A.data[[k]])
  n.vec.ADHD[k] = p.vec.A[length(p.vec.A)]
}

sample.size = as.data.frame(rbind(n.vec.Typ, n.vec.ADHD))
names(sample.size) = all.names
row.names(sample.size) = c("Typical","ADHD")
# write.csv(sample.size, "./results/sample.size.csv")



# compare prediction errors
mode.n = 2
c.coe = 20
c.test = 0.1
c = 0.6
cn = seq(0.1,2, length.out = 20)
L = 5

# Typical
type = "Typ"
all.data = fMRI.data.list$Typical
pred.error = as.data.frame(matrix(0, ncol=length(all.data), nrow=3))
names(pred.error) = all.names
row.names(pred.error) = methods
if(type == "ADHD"){
  all.data$data.Pittsburgh = NULL
  all.data$data.WashU = NULL
  all.data.na = c("KKI","NeuroIMAGE","Peking","NYU","OHSU")
  pred.error = as.data.frame(matrix(0, ncol=length(all.data.na), nrow=3))
  names(pred.error) = all.names
  row.names(pred.error) = methods
}
for (k in 1:length(all.data)) {
  t.data = all.data[[k]]
  A.data = all.data
  A.data[[k]] = NULL
  K = length(A.data)
  
  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.t = p.vec[M+1]
  p.vec = p.vec[1:M]
  K = length(A.data)

  n.n.k = 1:n.t
  for (l in 1:L) {
    n.test = n.n.k[(l-1)*floor(n.t/L) + 1:floor(n.t/L)]
    data.test = t.data[,,n.test]
    data.train = t.data[,,setdiff(n.n.k, n.test)]
    
    test.list = Initial.aggr(data.test, method = "Tlasso", cov.select= "tensor.prod", 
                             c.lam.sepa=c.coe, c.lam.Tlasso=c.coe, normalize = F)
    S.test.list = test.list$t.S.hat.list
    S.test = S.test.list[[mode.n]]
    
    
    t.lambda = c.coe*sqrt( p.vec*log(p.vec) / ( length(n.test) * prod(p.vec) ))
    t.Omega.hat.list = Tlasso.fit(data.train, lambda.vec = t.lambda, norm.type = 1)
    
    A.lambda = list()
    for (kk in 1:K) {
      A.lambda[[kk]] = c.coe*sqrt( log(p.vec) / ( dim(A.data[[kk]])[M+1] * prod(p.vec) ))
    }
    res.final = tensor.GGM.trans(data.train, A.data, A.lambda, c=c, cn.lam2=cn,
                                 init.method="Tlasso", init.method.aux="Tlasso",
                                 mode.set=c(mode.n), normalize = F)
    Omega.hat = res.final$Omega.list[[1]]
    Omega.hat.diff = res.final$Omega.list.diff[[1]]
    Omega.hat0 = t.Omega.hat.list[[mode.n]]
    
    
    pred.0 = sum(diag(S.test %*% Omega.hat0))/p.vec[mode.n] - log(det(Omega.hat0))/p.vec[mode.n]
    pred.pro = sum(diag(S.test %*% Omega.hat))/p.vec[mode.n] - log(det(Omega.hat))/p.vec[mode.n]
    pred.pro.diff = sum(diag(S.test %*% Omega.hat.diff))/p.vec[mode.n] - log(det(Omega.hat.diff))/p.vec[mode.n]
    
    pred.error[,k] = c(pred.0, pred.pro, pred.pro.diff)/L
    
    print(c(all.names[k], pred.error[,k]))
  }
  
}
relative.error = t( t(pred.error) / as.numeric(pred.error[1,]) )
size = as.character(sample.size[1,])
index.Typ = rbind(size, round(pred.error,6), round(relative.error[-1,],6))
pred.error.Typ = pred.error


# ADHD
type = "ADHD"
all.data = fMRI.data.list$ADHD
pred.error = as.data.frame(matrix(0, ncol=length(all.data), nrow=3))
names(pred.error) = all.names
row.names(pred.error) = methods
if(type == "ADHD"){
  all.data$data.Pittsburgh = NULL
  all.data$data.WashU = NULL
  all.names.na = c("KKI","NeuroIMAGE","Peking","NYU","OHSU")
  pred.error = as.data.frame(matrix(0, ncol=length(all.names.na), nrow=3))
  names(pred.error) = all.names.na
  row.names(pred.error) = methods
}
for (k in 1:length(all.data)) {
  t.data = all.data[[k]]
  A.data = all.data
  A.data[[k]] = NULL
  K = length(A.data)
  
  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.t = p.vec[M+1]
  p.vec = p.vec[1:M]
  K = length(A.data)

  n.n.k = 1:n.t
  for (l in 1:L) {
    n.test = n.n.k[(l-1)*floor(n.t/L) + 1:floor(n.t/L)]
    data.test = t.data[,,n.test]
    data.train = t.data[,,setdiff(n.n.k, n.test)]
    
    test.list = Initial.aggr(data.test, method = "Tlasso", cov.select= "tensor.prod", 
                             c.lam.sepa=c.coe, c.lam.Tlasso=c.coe, normalize = F)
    S.test.list = test.list$t.S.hat.list
    S.test = S.test.list[[mode.n]]
    
    
    t.lambda = c.coe*sqrt( p.vec*log(p.vec) / ( length(n.test) * prod(p.vec) ))
    t.Omega.hat.list = Tlasso.fit(data.train, lambda.vec = t.lambda, norm.type = 1)
    
    A.lambda = list()
    for (kk in 1:K) {
      A.lambda[[kk]] = c.coe*sqrt( log(p.vec) / ( dim(A.data[[kk]])[M+1] * prod(p.vec) ))
    }
    res.final = tensor.GGM.trans(data.train, A.data, A.lambda, c=c, cn.lam2=cn,
                                 init.method="Tlasso", init.method.aux="Tlasso",
                                 mode.set=c(mode.n), normalize = F)
    Omega.hat = res.final$Omega.list[[1]]
    Omega.hat.diff = res.final$Omega.list.diff[[1]]
    Omega.hat0 = t.Omega.hat.list[[mode.n]]
    
    
    pred.0 = sum(diag(S.test %*% Omega.hat0))/p.vec[mode.n] - log(det(Omega.hat0))/p.vec[mode.n]
    pred.pro = sum(diag(S.test %*% Omega.hat))/p.vec[mode.n] - log(det(Omega.hat))/p.vec[mode.n]
    pred.pro.diff = sum(diag(S.test %*% Omega.hat.diff))/p.vec[mode.n] - log(det(Omega.hat.diff))/p.vec[mode.n]
    
    pred.error[,k] = c(pred.0, pred.pro, pred.pro.diff)/L
    
    print(c(k, pred.error[,k]))
  }
  
}
relative.error = t( t(pred.error) / as.numeric(pred.error[1,]) )
size = as.character(sample.size[2,])
index.ADHD0 = rbind(size, round(pred.error,6), round(relative.error[-1,],6))
pred.error.ADHD = pred.error

index.ADHD = index.Typ
index.ADHD[,match(c("Pittsburgh","WashU"),all.names)] = "-"
index.ADHD[,-match(c("Pittsburgh","WashU"),all.names)] = index.ADHD0
index.ADHD[1,] = as.character(sample.size[2,])

name1 = c("Typical", rep("",5), "ADHD", rep("",5))
name2 = rep(c("sample size", "absolute error", "", "","relative error", ""), 2)
name3 = rep(c("",methods, methods[-1]), 2)
index.all = as.data.frame(cbind(name1, name2,name3,rbind(index.Typ,index.ADHD)))
names(index.all)[1:3] = c("","","Methods")


save.image("./results/prediction_error.RData")
write.csv(index.all, "./results/index.all.csv")



