################################################################################################
#Codes for conducting real data analysis
################################################################################################
library(Tlasso)
library(rTensor)
library(doParallel)
library(glasso)

rm(list = ls(all = TRUE))
ls()
# setwd("./case_study")
#############################################################################
## data-process
## Running time: about 1 minute
source("data-process.R")
#############################################################################


#############################################################################
## Use each site individually as a target domain and calculate 
## the prediction errors of different methods based on 5-fold cross-validation.
## Running time: about 55 minutes
rm(list = ls(all = TRUE))
ls()
load("./results/fMRI_data_list.RData")
source("../main_functions/function.R")
source("prediction-test.R")
#############################################################################


#############################################################################
## Apply the proposed method and output the tensor GGM 
## of ADHD brain functional connectivity networks
## Running time: about 2 minutes
rm(list = ls(all = TRUE))
ls()
load("./results/fMRI_data_list.RData")
source("../main_functions/function.R")
node.info = read.csv("./data/Node_AAL116_info.csv")
source("esti.GGM.R")
#############################################################################




