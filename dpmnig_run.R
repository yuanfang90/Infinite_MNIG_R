## remove previous list and set working directory
rm(list=ls())
## load packages
library(statmod)
library(MASS)
library(coda)
library(GeneralizedHyperbolic)
library(MCMCpack)
library(mclust)
library(class)
library(MixGHD)
library(mvtnorm)
library(truncnorm)


## load datasets
load("./data/sim1.RData") ## this contains 100 datasets in simulation 1 of our paper
set.seed(20211110)

## change the directory here
source("dpmnig_fun.R")

#### infinite MNIG
# data: data matrix
# true_lab: true class label
# maxiter: stop iteration if reached the maximum of iteration
# dp.alpha: numeric indicating a number as constant dp alpha parameter; or a vector(a,b) as updating along the algorithm with a gamma(a,b) prior.
# verbs: logic, print the chain numbers, iterations and current table of classifications or not
# plots: logic, trace plot or not
# outputs: "vector" or "lists", format of output estimated parameter

##################
procid <- 1

data=X[[procid]][,-3]
true_lab=X[[procid]][,3]

test1result <- dpMNIG(data=data,true_lab=true_lab,maxiter=10000,dp.alpha=0.5,verbs=TRUE,plots=TRUE,outputs="vectors")

