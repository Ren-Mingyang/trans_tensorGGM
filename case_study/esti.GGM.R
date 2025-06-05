##### Specify the target domain and corresponding parameters #####
target.name = "OHSU"
all.names = c("KKI","NeuroIMAGE","Peking","Pittsburgh","NYU","OHSU","WashU")
t.num = match(target.name, all.names)
c = 0.6
cn = seq(0.1,2, length.out = 20)
mode.n = 2

### Typical
A.data = fMRI.data.list$Typical
t.data = A.data[[t.num]]
A.data[[t.num]] = NULL

p.vec = dim(t.data)
M = length(p.vec) - 1
n = p.vec[M+1]
p.vec = p.vec[1:M]
K = length(A.data)


A.lambda = list()
for (k in 1:K) {
  c.coe = 20
  A.lambda[[k]] = c.coe*sqrt( log(p.vec) / ( dim(A.data[[k]])[M+1] * prod(p.vec) ))
}

t1 = proc.time()
res.final = tensor.GGM.trans(t.data, A.data, A.lambda, c=c, cn.lam2=cn,
                             init.method="Tlasso", init.method.aux="Tlasso",
                             mode.set=c(mode.n), normalize = F)
t2 = proc.time() - t1
res.Typ = res.final


### ADHD
A.data = fMRI.data.list$ADHD
t.data = A.data[[t.num]]
A.data[[t.num]] = NULL
A.data$data.Pittsburgh = NULL
A.data$data.WashU = NULL

p.vec = dim(t.data)
M = length(p.vec) - 1
n = p.vec[M+1]
p.vec = p.vec[1:M]
K = length(A.data)

A.lambda = list()
for (k in 1:K) {
  A.lambda[[k]] = c.coe*sqrt( log(p.vec) / ( dim(A.data[[k]])[M+1] * prod(p.vec) ))
}

t1 = proc.time()
res.final = tensor.GGM.trans(t.data, A.data, A.lambda, c=c, cn.lam2=cn,
                             init.method="Tlasso", init.method.aux="Tlasso",
                             cov.select="tensor.prod",
                             mode.set=c(mode.n), normalize = F)
t2.A = proc.time() - t1
res.ADHD = res.final


Omega.Typ = res.Typ$Omega.sym.list.diff[[1]]
Omega.ADHD = res.ADHD$Omega.sym.list.diff[[1]]
Theta.Typ = res.Typ$res.trans.list$Theta.hat.diff.list[[1]]
Theta.ADHD = res.ADHD$res.trans.list$Theta.hat.diff.list[[1]]

summary = function(Omega.Typ, Omega.ADHD){
  Omega.Typ[Omega.Typ!=0] = 1
  Omega.ADHD[Omega.ADHD!=0] = 1
  diff.Typ = Omega.Typ - Omega.ADHD
  diff.Typ[diff.Typ < 0] = 0
  diff.ADHD = Omega.ADHD - Omega.Typ
  diff.ADHD[diff.ADHD < 0] = 0
  
  edge.num.Typ = apply(Omega.Typ,1,sum) - 1
  edge.num.ADHD = apply(Omega.ADHD,1,sum) - 1
  
  edge.num.Typ.diff = apply(diff.Typ,1,sum) - 1
  edge.num.ADHD.diff = apply(diff.ADHD,1,sum) - 1
  
  return(list(edge.Typ=Omega.Typ, edge.ADHD=Omega.ADHD, 
              diff.Typ=diff.Typ, diff.ADHD=diff.ADHD,
              edge.num.Typ=edge.num.Typ, edge.num.ADHD=edge.num.ADHD,
              edge.num.Typ.diff=edge.num.Typ.diff,
              edge.num.ADHD.diff=edge.num.ADHD.diff))
}

res.Omega = summary(Omega.Typ, Omega.ADHD)
edge.num.Typ = res.Omega$edge.num.Typ
edge.num.ADHD = res.Omega$edge.num.ADHD
edge.Typ = res.Omega$edge.Typ
edge.ADHD = res.Omega$edge.ADHD

edge.num.Typ.diff = res.Omega$edge.num.Typ.diff
edge.num.ADHD.diff = res.Omega$edge.num.ADHD.diff
diff.Typ = res.Omega$diff.Typ
diff.ADHD = res.Omega$diff.ADHD

node.cor = node.info[,c(1:5,9)]
node.cor[,6] = as.numeric(node.cor[,6])

node.cor[,5] = 2.5
node.cor.Typ = node.cor

node.cor[,5] = 2.5
node.cor.ADHD = node.cor



node.cor[,5] = 2.5
node.diff.Typ = node.cor
diff.mat = diff.Typ
edge.num.diff = apply(abs(diff.mat), 1, sum)
hub.nodes0 = cbind(edge.num.diff, node.info[,c(9,6,7)])
node.diff.Typ[order(hub.nodes0[,1], decreasing = T)[1:12],4] = 1


node.cor[,5] = 2.5
node.diff.ADHD = node.cor
diff.mat = diff.ADHD
edge.num.diff = apply(abs(diff.mat), 1, sum)
hub.nodes0 = cbind(edge.num.diff, node.info[,c(9,6,7)])
node.diff.ADHD[order(hub.nodes0[,1], decreasing = T)[1:12],4] = 1






save.image("./results/OHSU.GGM.RData")

write.table(node.cor.Typ, "./results/aTyp.node", row.names = F, col.names = F)
write.table(edge.Typ, "./results/aTyp.edge", row.names = F, col.names = F)

write.table(node.cor.ADHD, "./results/aADHD.node", row.names = F, col.names = F)
write.table(edge.ADHD, "./results/aADHD.edge", row.names = F, col.names = F)

write.table(node.diff.Typ, "./results/aTyp.diff.node", row.names = F, col.names = F)
write.table(diff.Typ, "./results/aTyp.diff.edge", row.names = F, col.names = F)

write.table(node.diff.ADHD, "./results/aADHD.diff.node", row.names = F, col.names = F)
write.table(diff.ADHD, "./results/aADHD.diff.edge", row.names = F, col.names = F)



diff.mat = diff.Typ
edge.num.diff = apply(abs(diff.mat), 1, sum)
hub.nodes0 = cbind(edge.num.diff, node.info[,c(9,6,7)])
hub.Typ.diff = hub.nodes0[order(hub.nodes0[,1], decreasing = T),]

diff.mat = diff.ADHD
edge.num.diff = apply(abs(diff.mat), 1, sum)
hub.nodes0 = cbind(edge.num.diff, node.info[,c(9,6,7)])
hub.ADHD.diff = hub.nodes0[order(hub.nodes0[,1], decreasing = T),]

write.csv(hub.Typ.diff, "./results/hub.Typ.diff.csv")
write.csv(hub.ADHD.diff, "./results/hub.ADHD.diff.csv")














