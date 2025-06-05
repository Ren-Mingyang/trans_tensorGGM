################################################################################
# parameters supporting data generation
################################################################################
# precision matrices of the target domain
if(graph.type == "triangle"){
  t.Omega.true.list =  vector("list", length = M)
  for (m in 1:M) {
    t.Omega.true.list[[m]] = ChainOmega(p.vec[m], sd = m, norm.type = 1+as.numeric(normalize))
  }
}
if(graph.type == "nearest"){
  t.Omega.true.list =  vector("list", length = M)
  for (m in 1:M) {
    t.Omega.true.list[[m]] = NeighborOmega(p.vec[m], sd = m, norm.type = 1+as.numeric(normalize))
  }
}
# precision matrices of auxiliary domains
c0    = 1     # similarity between target and informative auxiliary domains
c1    = 10    # similarity between target and informative non-auxiliary domains
prob0 = 0.1
prob1 = 0.25
A.orac = c(1:K)
A.Omega.true.list = generate.Omega.A(t.Omega.true.list, A.orac, n, 
                                     p.vec, K, c0, prob0, c1, prob1)

################################################################################
# parameters supporting the proposed method
################################################################################
init.method="Tlasso"
c=0.6
# tuning parameters for initialization
cn = seq(0.1,2,length.out =10)
tla.lambda = 20*sqrt( p.vec*log(p.vec) / ( n * prod(p.vec) ))
t.lambda.int.trans = 20*sqrt( p.vec*log(p.vec) / ( n * prod(p.vec) ))
t.lambda.int.aggr = 20*sqrt( p.vec*log(p.vec) / ( n*(1-c) * prod(p.vec) ))
A.lambda = list()
for (k in 1:K) {
  A.lambda[[k]] = 20*sqrt( log(p.vec) / ( nA.vec[k] * prod(p.vec) ))
}


