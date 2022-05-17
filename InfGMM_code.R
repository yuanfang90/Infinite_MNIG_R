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
load("/usr2/postdoc/yuanf/test/DPMNIG/sim1/sim1.RData")
set.seed(20210311)

##### Functions
lik <- function(mu,Sigma_inv,x){
  den <- dmvnorm(x,mean = mu,sigma = solve(Sigma_inv))
  return(den)
}

hyper_fun<-function(z,x){
  a0<-colSums(z)
  a1<-t(t(x)%*%z)
  return(list(a0=a0,a1=a1))
}

update_V<-function(a0,a1,a2,a3,a4,x,z,d){
  V0<-list()
  for (g in 1:length(a0)){
    A <- matrix(rep(0, d^2), nrow = d, ncol = d)
    obs<-which(z[,g]==1)
    A_0 <- matrix(nrow = d, ncol = d)
    mean_mu<-a1[g,]/a0[g]
    for (i in obs){
      A_0 = (x[i,])%*%t(x[i,])
      A = A_0 + A
    }
    V2<-as.matrix(mean_mu)%*%t(as.matrix(mean_mu))*a0[g]
    V <- A-V2
    
    if (det(V)<0.00000001) {V0[[g]]<-solve(cov(x))}else{V0[[g]]<-solve(V)}
  }
  return(V0)
}

update_ISig<-function(a0,V0,g,d){
  if(a0[g]>(d+1)){
    new_ISig<-matrix(rWishart(1,a0[g],Sigma=V0[[g]]),d,d)
  }else{
    new_ISig<-matrix(rWishart(1,a0[g]+d+1,Sigma=V0[[g]]),d,d)
  }
  return(new_ISig)
}

update_Mu<-function(a0,a1,ISig,g){
  T_mu<-ISig[[g]]*a0[g]
  for_cov<-solve(T_mu,tol=(.Machine$double.neg.eps*1e-10))
  sym_for_cov<-matrix(NA,d,d) 
  sym_for_cov[lower.tri(sym_for_cov,diag=TRUE)]<-for_cov[lower.tri(for_cov,diag=TRUE)]
  sym_for_cov[upper.tri(sym_for_cov)]<-t(for_cov)[upper.tri(for_cov)]
  if(isSymmetric(for_cov)==TRUE){cov_mu<-for_cov}else{cov_mu<-sym_for_cov}
  mean_mu<-a1[g,]/a0[g]
  Mu<-rmvnorm(1,mean=mean_mu,sigma =cov_mu)
  
  return(Mu=Mu)
}

sample_common<-function(up_p,V0){
  up_p<-up_p
  V0<-V0
  for (g in 1:1){
    ISig[[g]]<-update_ISig(a0=up_p$a0,V0=V0,g=g,d=d)
    Mu[[g]]<-update_Mu(a0=up_p$a0,a1=up_p$a1,ISig,g)
  }
  return(list(ISig=ISig[[1]],Mu=Mu[[1]]))
}

#### find if there is any bad chain
find_diff <- function(x){
  a <- rep(1,3)
  # if(length(unique(x))==1){a=c(1,1,1)}
  if(x[1]!=x[2]&&x[2]==x[3]){a[1]=0}
  if(x[1]!=x[2]&&x[1]==x[3]){a[2]=0}
  if(x[1]==x[2]&&x[2]!=x[3]){a[3]=0}
  return(a)
}

#### update the concentration parameter
Update_d_alpha <- function(old_d_alpha,G,n,a,b){
  x = rbeta(1,old_d_alpha+1,n)
  pi_x = (a+G-1)/((a+G-1)+n*(b-log(x)))
  r = runif(1,0,1)
  if(r <= pi_x){
    d_alpha = rgamma(1,shape=(a+G),rate=(b-log(x)))
  }else{
    d_alpha = rgamma(1,shape=(a+G-1),rate=(b-log(x)))
  }
  return(d_alpha)
}

#######################
#### run the model ####
#######################

procid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# procid <- 1

data=X[[procid]][,-3]
true_lab=X[[procid]][,3]
maxiter=2000
dp.alpha=1
verbs=FALSE
plots=FALSE

#################
#### DPGMM  #####
#################
x <- data
require(MixGHD)
Sys.time()
fit_MGHD <- MGHD(data=x,G=1:5,gpar0=NULL,method="kmeans",scale=FALSE,modelSel="BIC")
GMGHD <- summary(fit_MGHD)$Cluster
MGHDchoose <- as.integer(length(GMGHD))
MixGHD_ARI <- adjustedRandIndex(fit_MGHD@map,true_lab)
MixGHD_BIC <- fit_MGHD@BIC

ghd_results <- c(MGHDchoose,MixGHD_ARI,MixGHD_BIC)
output.filename = paste0("/usr2/postdoc/yuanf/test/DPMNIG/sim1/gmm/ghdsim",procid,".csv",sep='')
# write.csv(t(ghd_results),file=output.filename)
Sys.time()

n<-nrow(x)
d<-ncol(x)
it<-1
max_it<-600
check<-0
G<-1
G_lab<-1 #This will make sure that when a group is removed the label doesn't go down like G would

z<-matrix(1,ncol=1,nrow=n)

####Initialize parameters
ISig<-list()
Mu<-list()

Mu[[1]]<-t(sapply(1:G,function(g){colSums(x*z[,g])/sum(z[,g])}))
ISig[[1]]<-solve(cov(x))

########## Get ready for extra layer common prior
b0 = 1/n
b1 = cov(x)
c1 = colSums(x)
nu0 = d+1
L0 = (cov(x))/nu0

###### Get three chains ready
all_loglik <- list()
all_Class_matrix <- list()
all_parameters <- list()
all_classification <- list()

###### Get three chains ready
all_loglik <- list()
all_Class_matrix <- list()
all_parameters <- list()
all_classification <- list()

for(chain in 1:3){
  all_loglik[[chain]] <- numeric(max_it)
  all_Class_matrix[[chain]] <- matrix(0,ncol=max_it,nrow=n)
  all_parameters[[chain]] <- list()
  all_classification[[chain]] <- list()
}

multi_starts <- c()
multi_starts[1] <- 1
multi_starts[2] <- n
multi_starts[3] <- sample(1:n,1,replace=TRUE)

for(chain in 1:3){
# chain=1

  it<-1
  if(length(dp.alpha)==1){d_alpha <- dp.alpha}
  if(length(dp.alpha)==2){d_alpha <- rgamma(1,shape=dp.alpha[1],rate=dp.alpha[2])}

  G<-multi_starts[chain]
  G_lab<-G #This will make sure that when a group is removed the label doesn't go down like G would
  
  Class<-all_Class_matrix[[chain]]
  
  if(chain == 1){
    Class[,1]<-rep("c1",n) #All observations in the same group
  }
  if(chain == 2){
    Class[,1]<-paste("c",c(1:n),sep="") #Each obs on its own grp
  }
  if(chain == 3){
    Class[,1]<-paste("c",kmeans(x,centers=G)$cluster,sep="") #k-means
  }
  
  z<-as.matrix(unmap(Class[,1]))
  pi_g<-colMeans(z)
  Ng<-table(Class[,1]) 
  labels<-unique(Class[,1])
  
  ####Initialize the hyper parameters
  a0 = mean(rexp(10000,rate = b0))
  a1 = t(colMeans(as.matrix(mvrnorm(10000, c1, b1))))
  forV0 = rWishart(10000,nu0,L0)
  V0_co = solve(Reduce("+",lapply(seq(dim(forV0)[3]),function(x) forV0[,,x]))/10000)
  
  up_p<-list(a0=a0,a1=a1)
  V0_co<-list(V0_co)
  
  #### get parameters for iteration 1
  ISig<-list()
  Mu<-list()

  for (g in 1:G){
    com_par<-sample_common(up_p=up_p,V0=V0_co)
    ISig[[g]]<-com_par$ISig
    Mu[[g]]<-com_par$Mu
  }
  names(ISig)<-labels
  names(Mu)<-labels

  it<-it+1
  i<-1
  
  loglik <- NULL
  all_loglik[[chain]][1] <- -1400
  
  ####################### 
  ##### Start Iterations
  while (it<(max_it+1)){
    Class_old<-Class[,it-1]
    Class[,it]<-Class[,it-1]
    Ng_old<-table(Class[,it-1])
    # print(Ng_old)
    
    if(length(dp.alpha)==2){
      old_d_alpha <- d_alpha
      d_alpha <- Update_d_alpha(old_d_alpha,G,n,a=dp.alpha[1],b=dp.alpha[2])
    }
    
    for (i in 1:n){
      # print(i)
      ### For new group  
      sam_com<-sample_common(up_p=up_p,V0=V0_co)
      ##### Prob for new groups
      new<-d_alpha*lik(mu=sam_com$Mu,Sigma_inv=sam_com$ISig,x=t(x[i,]))
      ### Prob for old groups
      old<-NULL
      for (j in labels){
        if(Ng[j] == 1){
          old[j]<-d_alpha*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
        }else{
          old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
        }
      }
      if(all(c(new,old)==0)){old <- Ng}
      prob_g<-c(new,old)/sum(c(new,old))
      
      ### Updating class labels
      Class[i,it]<-sample(c("new",labels),1,prob=prob_g)
      
      ####Since the group didn't change, do nothing####
      # if (Class[i,it]==Class[i,it-1]){
      # print("Do Nothing")
      # }
      
      if (Class[i,it]=="new"){
        ###Checking if it was the only observation in that group. In that case, it would give a length of 0 because there are no other observations with that group label
        only_obs<-(length(which(Class[,it]==Class_old[i]))==0)
        if (only_obs){
          Class[i,it]<-Class_old[i]
          Mu[[Class_old[i]]]<-sam_com$Mu
          ISig[[Class_old[i]]]<-sam_com$ISig
        }
        if (!only_obs){
          # print("combine to new")
          G<-G+1
          G_lab<-G_lab+1
          Mu[[G]]<-sam_com$Mu
          ISig[[G]]<-sam_com$ISig

          Class[i,it]<-paste("c",G_lab,sep="")
          names(Mu)<-c(labels,paste("c",G_lab,sep=""))
          names(ISig)<-c(labels,paste("c",G_lab,sep=""))
          Ng<-table(Class[,it]) 
          labels<-names(table(Class[,it]))
        }
      }
      
      if (Class[i,it]!=Class[i,it-1]&Class[i,it]!="new"){
        empty<-(length(which(Class[,it]==Class_old[i]))==0) 
        if (!empty){
          # print("Update group membership totals only")
          Ng<-table(Class[,it]) 
        }
        ##If empties the previous group
        if (empty){
          # print("remove one grp")
          G<-G-1
          Mu[[Class_old[i]]]<-NULL
          ISig[[Class_old[i]]]<-NULL
          # V0[[Class_old[i]]]<-NULL
          Ng<-table(Class[,it])
          labels<-names(table(Class[,it]))
        }
      }
    }
    z<-as.matrix(unmap(Class[,it]))
    pi_g <- colMeans(z)
    # print(pi_g)
    # pi_g <- rdirichlet(G, d_alpha/G)
    # plot(x,col=as.factor(Class[,it]))
    
    ###Updating the hyperparameters
    hyper<-hyper_fun(z,x)
    V0<-update_V(a0=hyper$a0,a1=hyper$a1,x=x,z=z,d=d)
    # print(hyper)
    # print(V0)
    
    ###Updating the parameters
    for (g in 1:G){
      obs<-which(z[,g]==1)
      if (length(obs)>1){
        ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
        Mu[[g]]<-update_Mu(hyper$a0,hyper$a1,ISig,g=g)
      }else{
        sam_com<-sample_common(up_p=up_p,V0=V0_co)
        ISig[[g]]<-sam_com$ISig
        Mu[[g]]<-sam_com$Mu
      }
    }
    # print(ISig)
    # print(Mu)
    
    ### Store parameters and the values for updating common hyperparameters
    mu_isig<-matrix(nrow=G,ncol=d)
    mu_isig_mu<-numeric(G)
    for (g in 1:G){
      # print(g)
      mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
      mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
    }
    
    ###Updating common hyperparameter using their posterior distributions
    b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)+sum(mu_isig_mu))
    if(b0_new <= 0){b0_new=b0}
    b1_new=b1
    c1_new=c1+(colSums(mu_isig))%*%b1
    # print(b0_new)
    # print(b1_new)
    # print(c1_new)

    a0<-rexp(1,rate=b0_new)
    a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
    # print(a0)
    # print(a1)
    
    nu0_new=nu0+G*a0
    L0_new=solve(solve(L0)+Reduce("+",ISig))
    # print(nu0_new*L0_new)
    
    up_p<-list(a0=a0,a1=a1)
    V0_co<-list(riwish(nu0_new,L0_new))
    # print(up_p)
    # print(V0_co)
    
    if(verbs){
      cat("Chain:",multi_starts[chain],"start","\n")
      cat("Iteration", it, "\n")
      print(table(Class[,it]))
    }
    
    ##Calucluate the log likelihood
    likelihood <- matrix(0,nrow=n,ncol=G)
    for (g in 1:G){
      likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],Sigma_inv=ISig[[g]],x=t(x[i,]))})
    }
    loglik <- sum(log(rowSums(likelihood)))
    
    all_loglik[[chain]][it] <- loglik
    
    it<-it+1
  }
  all_Class_matrix[[chain]] <- Class
  all_parameters[[chain]] <- list(up_p,V0_co,Mu,ISig)
  all_classification[[chain]] <- list(z,labels,Ng,G,G_lab,d_alpha)
}

### check convergence
burnin <- 1:(max_it-500)
# burnin <- 1:(0.9*max_it)
object <- list()
negInfLik <- 0
for(chain in 1:3){
  if(any(all_loglik[[chain]] == -Inf)){negInfLik <- 1}
  object[[chain]]=mcmc(all_loglik[[chain]][-burnin])
}
obj <- mcmc.list(object)
conv <- gelman.diag(obj)
if(max_it == maxiter|negInfLik == 1){
  check.conv <- 1
}else{
  check.conv <- conv[[1]][1]
}
check.thresh <- 1.1

if(plots){
  # plot likelihood for all chains
  plot(all_loglik[[1]],type = 'l',#ylim=c(-1160,-1060),
       main=paste("check.conv is",check.conv))
  lines(all_loglik[[2]],col='red')
  lines(all_loglik[[3]],col='blue')
}

while(check.conv > check.thresh){
  it_old <- it
  max_it <- max_it+100
  
  for(chain in 1:3){
    up_p <- all_parameters[[chain]][[1]]
    V0_co <- all_parameters[[chain]][[2]]
    Mu <- all_parameters[[chain]][[3]]
    ISig <- all_parameters[[chain]][[4]]
    
    z <- all_classification[[chain]][[1]]
    labels <- all_classification[[chain]][[2]]
    Ng <- all_classification[[chain]][[3]]
    G <- all_classification[[chain]][[4]]
    G_lab <- all_classification[[chain]][[5]]
    d_alpha <- all_classification[[chain]][[6]]
    
    add_class_mat <- matrix(0,ncol=100,nrow=n)
    Class<-cbind(all_Class_matrix[[chain]],add_class_mat)
    
    it <- it_old
    while (it<(max_it+1)){
      # for(it in 2:102){
      Class_old<-Class[,it-1]
      Class[,it]<-Class[,it-1]
      Ng_old<-table(Class[,it-1])
      
      if(length(dp.alpha)==2){
        old_d_alpha <- d_alpha
        d_alpha <- Update_d_alpha(old_d_alpha,G,n,a=dp.alpha[1],b=dp.alpha[2])
      }
      
      for (i in 1:n){
        # print(i)
        ### For new group  
        sam_com<-sample_common(up_p=up_p,V0=V0_co)
        ##### Prob for new groups
        new<-d_alpha*lik(mu=sam_com$Mu,Sigma_inv=sam_com$ISig,x=t(x[i,]))
        ### Prob for old groups
        old<-NULL
        for (j in labels){
          if(Ng[j] == 1){
            old[j]<-d_alpha*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
          }else{
            old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
          }
        }
        if(all(c(new,old)==0)){old <- Ng}
        prob_g<-c(new,old)/sum(c(new,old))
        
        ### Updating class labels
        Class[i,it]<-sample(c("new",labels),1,prob=prob_g)
        
        ####Since the group didn't change, do nothing####
        # if (Class[i,it]==Class[i,it-1]){
        # print("Do Nothing")
        # }
        
        if (Class[i,it]=="new"){
          ###Checking if it was the only observation in that group. In that case, it would give a length of 0 because there are no other observations with that group label
          only_obs<-(length(which(Class[,it]==Class_old[i]))==0)
          if (only_obs){
            Class[i,it]<-Class_old[i]
            Mu[[Class_old[i]]]<-sam_com$Mu
            ISig[[Class_old[i]]]<-sam_com$ISig
          }
          if (!only_obs){
            # print("combine to new")
            G<-G+1
            G_lab<-G_lab+1
            Mu[[G]]<-sam_com$Mu
            ISig[[G]]<-sam_com$ISig
            
            Class[i,it]<-paste("c",G_lab,sep="")
            names(Mu)<-c(labels,paste("c",G_lab,sep=""))
            names(ISig)<-c(labels,paste("c",G_lab,sep=""))
            Ng<-table(Class[,it]) 
            labels<-names(table(Class[,it]))
          }
        }
        
        if (Class[i,it]!=Class[i,it-1]&Class[i,it]!="new"){
          empty<-(length(which(Class[,it]==Class_old[i]))==0) 
          if (!empty){
            # print("Update group membership totals only")
            Ng<-table(Class[,it]) 
          }
          ##If empties the previous group
          if (empty){
            # print("remove one grp")
            G<-G-1
            Mu[[Class_old[i]]]<-NULL
            ISig[[Class_old[i]]]<-NULL
            # V0[[Class_old[i]]]<-NULL
            Ng<-table(Class[,it])
            labels<-names(table(Class[,it]))
          }
        }
      }
      z<-as.matrix(unmap(Class[,it]))
      pi_g <- colMeans(z)
      # print(pi_g)
      # pi_g <- rdirichlet(G, d_alpha/G)
      # plot(x,col=as.factor(Class[,it]))
      
      ###Updating the hyperparameters
      hyper<-hyper_fun(z,x)
      V0<-update_V(a0=hyper$a0,a1=hyper$a1,x=x,z=z,d=d)
      # print(hyper)
      # print(V0)
      
      ###Updating the parameters
      for (g in 1:G){
        obs<-which(z[,g]==1)
        if (length(obs)>1){
          ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
          Mu[[g]]<-update_Mu(hyper$a0,hyper$a1,ISig,g=g)
        }else{
          sam_com<-sample_common(up_p=up_p,V0=V0_co)
          ISig[[g]]<-sam_com$ISig
          Mu[[g]]<-sam_com$Mu
        }
      }
      # print(ISig)
      # print(Mu)
      
      ### Store parameters and the values for updating common hyperparameters
      mu_isig<-matrix(nrow=G,ncol=d)
      mu_isig_mu<-numeric(G)
      for (g in 1:G){
        # print(g)
        mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
        mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
      }
      
      ###Updating common hyperparameter using their posterior distributions
      b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)+sum(mu_isig_mu))
      if(b0_new <= 0){b0_new=b0}
      b1_new=b1
      c1_new=c1+(colSums(mu_isig))%*%b1
      # print(b0_new)
      # print(b1_new)
      # print(c1_new)
      
      a0<-rexp(1,rate=b0_new)
      a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
      # print(a0)
      # print(a1)
      
      nu0_new=nu0+G*a0
      L0_new=solve(solve(L0)+Reduce("+",ISig))
      # print(nu0_new*L0_new)
      
      up_p<-list(a0=a0,a1=a1)
      V0_co<-list(riwish(nu0_new,L0_new))
      # print(up_p)
      # print(V0_co)
      
      if(verbs){
        cat("Chain:",multi_starts[chain],"start","\n")
        cat("Iteration", it, "\n")
        print(table(Class[,it]))
      }
      
      ##Calucluate the log likelihood
      likelihood <- matrix(0,nrow=n,ncol=G)
      for (g in 1:G){
        likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],Sigma_inv=ISig[[g]],x=t(x[i,]))})
      }
      loglik <- sum(log(rowSums(likelihood)))
      
      all_loglik[[chain]][it] <- loglik
      
      it<-it+1
    }
    all_Class_matrix[[chain]] <- Class
    all_parameters[[chain]] <- list(up_p,V0_co,Mu,ISig)
    all_classification[[chain]] <- list(z,labels,Ng,G,G_lab,d_alpha)
  }
  
  ### check convergence
  burnin <- 1:(max_it-500)
  # burnin <- 1:(0.9*max_it)
  object <- list()
  negInfLik <- 0
  for(chain in 1:3){
    if(any(all_loglik[[chain]] == -Inf)){negInfLik <- 1}
    object[[chain]]=mcmc(all_loglik[[chain]][-burnin])
  }
  obj <- mcmc.list(object)
  conv <- gelman.diag(obj)
  if(max_it == maxiter|negInfLik == 1){
    check.conv <- 1
  }else{
    check.conv <- conv[[1]][1]
  }
  check.thresh <- 1.1
  
  if(plots){
    # plot likelihood for all chains
    plot(all_loglik[[1]],type = 'l',#ylim=c(-1160,-1060),
         main=paste("check.conv is",check.conv))
    lines(all_loglik[[2]],col='red')
    lines(all_loglik[[3]],col='blue')
  }
}

####################################################################################
####### creat 500 more samples each chain after reaching convergence or maxiter ####
####################################################################################
it_old <- it
max_it <- max_it+500
store_all_pars <- list()
ord_Class_labels <- list()
# ord_Class_names <- list()

for(chain in 1:3){
  store_all_pars[[chain]] <- list()
  ord_Class_labels[[chain]] <- matrix(nrow=n,ncol=500)
  # ord_Class_names[[chain]] <- matrix(nrow=n,ncol=500)
}

for(chain in 1:3){
  up_p <- all_parameters[[chain]][[1]]
  V0_co <- all_parameters[[chain]][[2]]
  Mu <- all_parameters[[chain]][[3]]
  ISig <- all_parameters[[chain]][[4]]
  
  z <- all_classification[[chain]][[1]]
  labels <- all_classification[[chain]][[2]]
  Ng <- all_classification[[chain]][[3]]
  G <- all_classification[[chain]][[4]]
  G_lab <- all_classification[[chain]][[5]]
  d_alpha <- all_classification[[chain]][[6]]
  
  add_class_mat <- matrix(0,ncol=500,nrow=n)
  Class<-cbind(all_Class_matrix[[chain]],add_class_mat)
  
  it <- it_old
  while (it<(max_it+1)){
    # for(it in 2:102){
    Class_old<-Class[,it-1]
    Class[,it]<-Class[,it-1]
    Ng_old<-table(Class[,it-1])
    
    if(length(dp.alpha)==2){
      old_d_alpha <- d_alpha
      d_alpha <- Update_d_alpha(old_d_alpha,G,n,a=dp.alpha[1],b=dp.alpha[2])
    }
    
    for (i in 1:n){
      # print(i)
      ### For new group  
      sam_com<-sample_common(up_p=up_p,V0=V0_co)
      ##### Prob for new groups
      new<-d_alpha*lik(mu=sam_com$Mu,Sigma_inv=sam_com$ISig,x=t(x[i,]))
      ### Prob for old groups
      old<-NULL
      for (j in labels){
        if(Ng[j] == 1){
          old[j]<-d_alpha*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
        }else{
          old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],Sigma_inv=ISig[[j]],x=t(x[i,]))
        }
      }
      if(all(c(new,old)==0)){old <- Ng}
      prob_g<-c(new,old)/sum(c(new,old))
      
      ### Updating class labels
      Class[i,it]<-sample(c("new",labels),1,prob=prob_g)
      
      ####Since the group didn't change, do nothing####
      # if (Class[i,it]==Class[i,it-1]){
      # print("Do Nothing")
      # }
      
      if (Class[i,it]=="new"){
        ###Checking if it was the only observation in that group. In that case, it would give a length of 0 because there are no other observations with that group label
        only_obs<-(length(which(Class[,it]==Class_old[i]))==0)
        if (only_obs){
          Class[i,it]<-Class_old[i]
          Mu[[Class_old[i]]]<-sam_com$Mu
          ISig[[Class_old[i]]]<-sam_com$ISig
        }
        if (!only_obs){
          # print("combine to new")
          G<-G+1
          G_lab<-G_lab+1
          Mu[[G]]<-sam_com$Mu
          ISig[[G]]<-sam_com$ISig
          
          Class[i,it]<-paste("c",G_lab,sep="")
          names(Mu)<-c(labels,paste("c",G_lab,sep=""))
          names(ISig)<-c(labels,paste("c",G_lab,sep=""))
          Ng<-table(Class[,it]) 
          labels<-names(table(Class[,it]))
        }
      }
      
      if (Class[i,it]!=Class[i,it-1]&Class[i,it]!="new"){
        empty<-(length(which(Class[,it]==Class_old[i]))==0) 
        if (!empty){
          # print("Update group membership totals only")
          Ng<-table(Class[,it]) 
        }
        ##If empties the previous group
        if (empty){
          # print("remove one grp")
          G<-G-1
          Mu[[Class_old[i]]]<-NULL
          ISig[[Class_old[i]]]<-NULL
          # V0[[Class_old[i]]]<-NULL
          Ng<-table(Class[,it])
          labels<-names(table(Class[,it]))
        }
      }
    }
    z<-as.matrix(unmap(Class[,it]))
    pi_g <- colMeans(z)
    # print(pi_g)
    # pi_g <- rdirichlet(G, d_alpha/G)
    # plot(x,col=as.factor(Class[,it]))
    
    ###Updating the hyperparameters
    hyper<-hyper_fun(z,x)
    V0<-update_V(a0=hyper$a0,a1=hyper$a1,x=x,z=z,d=d)
    # print(hyper)
    # print(V0)
    
    ###Updating the parameters
    for (g in 1:G){
      obs<-which(z[,g]==1)
      if (length(obs)>1){
        ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
        Mu[[g]]<-update_Mu(hyper$a0,hyper$a1,ISig,g=g)
      }else{
        sam_com<-sample_common(up_p=up_p,V0=V0_co)
        ISig[[g]]<-sam_com$ISig
        Mu[[g]]<-sam_com$Mu
      }
    }
    # print(ISig)
    # print(Mu)
    
    ### Store parameters and the values for updating common hyperparameters
    mu_isig<-matrix(nrow=G,ncol=d)
    mu_isig_mu<-numeric(G)
    for (g in 1:G){
      # print(g)
      mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
      mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
    }
    
    ###Updating common hyperparameter using their posterior distributions
    b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)+sum(mu_isig_mu))
    if(b0_new <= 0){b0_new=b0}
    b1_new=b1
    c1_new=c1+(colSums(mu_isig))%*%b1
    # print(b0_new)
    # print(b1_new)
    # print(c1_new)
    
    a0<-rexp(1,rate=b0_new)
    a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
    # print(a0)
    # print(a1)
    
    nu0_new=nu0+G*a0
    L0_new=solve(solve(L0)+Reduce("+",ISig))
    # print(nu0_new*L0_new)
    
    up_p<-list(a0=a0,a1=a1)
    V0_co<-list(riwish(nu0_new,L0_new))
    # print(up_p)
    # print(V0_co)
    
    if(verbs){
      cat("Chain:",multi_starts[chain],"start","\n")
      cat("Iteration", it, "\n")
      print(table(Class[,it]))
    }
    
    ##Calucluate the log likelihood
    likelihood <- matrix(0,nrow=n,ncol=G)
    for (g in 1:G){
      likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],Sigma_inv=ISig[[g]],x=t(x[i,]))})
    }
    
    loglik <- sum(log(rowSums(likelihood)))
    all_loglik[[chain]] <- append(all_loglik[[chain]],loglik)
    
    ### First order the parameters according to the first element of Mu
    ### Then store the parameters for further calculation
    ord_grp_ind <- which(Ng > 0.03*n)
    res_ord_ind <- which(Ng <= 0.03*n)
    mu.temp <- Mu[names(ord_grp_ind)]
    
    ord.temp <- rank(-sapply(mu.temp,'[[',2))
    resid.ord <- c((length(ord.temp)+1):G)
    names(resid.ord) <- names(res_ord_ind)
    
    if(length(ord.temp)==G){orders<-ord.temp}else{orders<-c(ord.temp,resid.ord)}
    
    ord_Mu <- list()
    ord_ISig <- list()
    ord_class <- numeric(n)
    # ord_names <- numeric(n)
    for (g in 1:G){
      order_ind <- which(orders == g)
      ord_Mu[[g]] <- Mu[[order_ind]]
      ord_ISig[[g]] <- ISig[[order_ind]]
      ord_class[which(Class[,it] == names(which(orders == g)))] = g
      # ord_names[which(Class[,it] == names(which(orders == g)))] = names(order_ind)
    }
    store_all_pars[[chain]][[it]] <- list(ord_Mu,ord_ISig) 
    ord_Class_labels[[chain]][,(it+1-it_old)] <- ord_class
    # ord_Class_names[[chain]][,(it+1-it_old)] <- ord_names
    
    it<-it+1
  }
  all_Class_matrix[[chain]] <- Class
  all_parameters[[chain]] <- list(up_p,V0_co,Mu,ISig)
  all_classification[[chain]] <- list(z,labels,Ng,G,G_lab,d_alpha)
}

### check convergence
burnin <- 1:(max_it-500)
# burnin <- 1:(0.9*max_it)
object <- list()
negInfLik <- 0
for(chain in 1:3){
  if(any(all_loglik[[chain]] == -Inf)){negInfLik <- 1}
  object[[chain]]=mcmc(all_loglik[[chain]][-burnin])
}
obj <- mcmc.list(object)
conv <- gelman.diag(obj)
check.conv <- conv[[1]][1]

######################################################################
#### Summarize the clustering result and the parameter estimation ####
######################################################################
# finalize the clustering using last 400 class labels 
# for each chain separatetly
all_final_G <- numeric(3)
all_final_class <- matrix(nrow=n,ncol=3)
for(chain in 1:3){
  Class <- all_Class_matrix[[chain]]
  final_class<-NULL
  for (i in 1:n){
    final_class[i]<-names(which.max(table(Class[i,(max_it-399):max_it])))
  }
  all_final_class[,chain] <- final_class
  labels<-names(table(final_class))
  G = length(labels)
  all_final_G[chain] <- G
}

chain_conv_ind <- find_diff(all_final_G)
convchains <- c(1,2,3)[which(chain_conv_ind==1)]
# convchains <- c(1,2,3)[which(all_final_G==4)]

ord_class_mat <- Reduce(cbind,ord_Class_labels[convchains])
# ord_names_mat <- Reduce(cbind,ord_Class_names[convchains])
ord_final_class <- NULL
# ord_final_class_names <- NULL
for(i in 1:n){
  ord_final_class[i] <- names(which.max(table(ord_class_mat[i,])))
  # ord_final_class_names[i] <- names(which.max(table(ord_names_mat[i,(max_it-399):max_it])))
}
z<-as.matrix(unmap(ord_final_class))
pi_g <- colMeans(z)
ARI <- adjustedRandIndex(ord_final_class,true_lab)

# summarize parameter estimation
mu_exp <- list()
sig_exp <- list()
class_labs <- list()
loglik_approx_ls <- list()
  
for (chain in convchains){
  len = length(all_loglik[[chain]])
  allparas <- store_all_pars[[chain]][(len-399):len]
  G = all_final_G[chain]
  loglik_approx_ls[[chain]] <- all_loglik[[chain]][(len-399):len]
  # G=3
  
  mu_exp[[chain]]<-matrix(nrow=400,ncol=d*G)
  sig_exp[[chain]]<-matrix(nrow=400,ncol=d*d*G)
  
  for(i in 1:400){
    for(g in 1:G){
      mu.try <- try(allparas[[i]][[1]][[g]],silent=TRUE)
      if(is.numeric(mu.try)&!(class(mu.try) == "try-error")){
        mu_exp[[chain]][i,(1+d*(g-1)):(d*g)]<-mu.try
      }
      sig.try <- try(as.vector(solve(allparas[[i]][[2]][[g]])),silent=TRUE)
      if(is.numeric(sig.try)&!(class(sig.try) == "try-error")){
        sig_exp[[chain]][i,(1+d*d*(g-1)):(d*d*g)]<-sig.try
      }
    }
  }
}

loglik_approx <- mean(unlist(loglik_approx_ls))
mu <- apply(Reduce(rbind,mu_exp),2,function(x){mean(x,na.rm=TRUE)})
sig <- apply(Reduce(rbind,sig_exp),2,function(x){mean(x,na.rm=TRUE)})
G=unique(all_final_G[convchains])

Mu <- list()
Sigma <- list()

for(g in 1:G){
  Mu[[g]] <- mu[(1+d*(g-1)):(d*g)]
  Sigma[[g]] <- matrix(sig[(1+d*d*(g-1)):(d*d*g)],nrow=d,ncol=d)
}
ISig <- lapply(Sigma,solve)

npar <- (d)*G+0.5*(d+1)*(d)*G+(G-1)
BIC <- 2*loglik_approx - npar*log(n)
gmm_results <- c(G,ARI,BIC)
Sys.time()
###################################
##### Save all results ############
###################################
gmm_results <- c(gmm_results)
output.filename = paste0("/usr2/postdoc/yuanf/test/DPMNIG/sim1/gmm/gmmsim",procid,".csv",sep='')
write.csv(t(gmm_results),file=output.filename)

