################################################################################################
#Codes for conducting simulations
################################################################################################
rm(list = ls(all = TRUE))
ls()
# setwd("./simulations")
library(Tlasso)
library(rTensor)
library(doParallel)
library(glasso)
source("../main_functions/function.R")
source("../simulations/sim_func.R")

n = 50                   # sample size of target domain
nk = 100                 # sample size of an auxiliary domain
K = 5                    # number of auxiliary domains
nA.vec = rep(nk, K)      # sample sizes of all auxiliary domains
normalize = T            # normalization of precision matrix

p.vec = c(10,10,10)      # dimensions of tensor
M = length(p.vec)        # number of tensor modes
graph.type = "triangle"  # the type of the target graph: "triangle" and "nearest"


source("../simulations/sim_func_para.R")
data.list = generate.data(M, t.Omega.true.list, A.Omega.true.list, nA.vec, K)
t.data = data.list$t.data
A.data = data.list$A.data

# the proposed method
res.final = tensor.GGM.trans(t.data, A.data, A.lambda, A.orac = A.orac, c=c,
                             cn.lam2=cn, init.method=init.method,
                             t.lambda.int.trans=t.lambda.int.trans,
                             t.lambda.int.aggr=t.lambda.int.aggr,
                             normalize = normalize)
# Tlasso
Tlasso.Omega.list = Tlasso.fit(t.data, lambda.vec = tla.lambda, norm.type = 1+as.numeric(normalize))
i.Tlasso = as.data.frame(t(unlist(est.analysis(Tlasso.Omega.list, t.Omega.true.list))))

# summary
index.list = evaluation(res.final, t.Omega.true.list)
index.list$i.Omega           # proposed
index.list$i.Omega.diff      # proposed.v
i.Tlasso                     # Tlasso
