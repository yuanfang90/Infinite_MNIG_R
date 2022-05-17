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
library(rrcov)
library(DAAG)


## load datasets
# load("/usr2/postdoc/yuanf/test/DPMNIG/sim1/sim1.RData")
set.seed(20210331)

set.tol <- .Machine$double.xmin

## change the directory here
# source("/usr2/postdoc/yuanf/test/DPMNIG/dpmnig_fun.r")

#####################
###           #######
### functions #######
###           #######
#####################
phi_fun<-function(dat,G,n,Mu,ISig){
  phi<-matrix(NA,ncol=G,nrow=n)
  for (g in 1:G){
    phi[,g]<-sqrt(1+mahalanobis(dat,center=Mu[[g]],cov=ISig[[g]],inverted=TRUE))	
  }
  return(phi)
}

alpha_fun<-function(G,Gamma,Beta,ISig){
  alpha<-sapply(1:G,function(g){sqrt(Gamma[[g]]^2+(Beta[[g]])%*%ISig[[g]]%*%t(Beta[[g]]))})
  return(alpha)
}

U_fun<-function(G,phi,alpha,lam){
  U<-sapply(1:G,function(g){phi[,g]*besselK(phi[,g]*alpha[g],lam+1,expon.scaled=TRUE)/(besselK(phi[,g]*alpha[g],lam,expon.scaled=TRUE)*alpha[g])})
  return(U)
}

Uinv_fun<-function(G,phi,alpha,lam){
  Uinv<-sapply(1:G,function(g){alpha[g]/phi[,g]*besselK(phi[,g]*alpha[g],lam-1,expon.scaled=TRUE)/besselK(phi[,g]*alpha[g],lam,expon.scaled=TRUE)})
  return(Uinv)
}

hyper_fun<-function(z,U,Uinv,x){
  a0<-colSums(z)
  a1<-t(t(x)%*%z)
  a2<-t(t(x)%*%(z*Uinv))
  a3<-colSums(z*U)
  a4<-colSums(z*Uinv)
  return(list(a0=a0,a1=a1,a2=a2,a3=a3,a4=a4))
}

update_V<-function(a0,a1,a2,a3,a4,x,z,U,Uinv,d){
  V0<-list()
  for_u_inv <- numeric(dim(x)[1])
  for (g in 1:length(a0)){
    for_u_inv[which(z[,g] == 1)]=Uinv[which(z[,g] == 1),g]
    A <- matrix(rep(0, d^2), nrow = d, ncol = d)
    obs<-which(z[,g]==1)
    A_0 <- matrix(nrow = d, ncol = d)
    mean_mu<-a2[g,]*a3[g]/(a3[g]*a4[g]-a0[g]^2)-a1[g,]*a0[g]/(a3[g]*a4[g]-a0[g]^2)
    mean_beta<-a1[g,]*a4[g]/(a3[g]*a4[g]-a0[g]^2)-a2[g,]*a0[g]/(a3[g]*a4[g]-a0[g]^2)
    for (i in obs){
      A_0 = (x[i,])%*%t(x[i,])*for_u_inv[i]
      A = A_0 + A
    }
    V2<-as.matrix(mean_mu)%*%t(as.matrix(mean_mu))*a4[g]
    V3<-as.matrix(mean_beta)%*%t(as.matrix(mean_beta))*a3[g]
    V4<-as.matrix(mean_mu)%*%t(as.matrix(mean_beta))*a0[g]
    V5<-as.matrix(mean_beta)%*%t(as.matrix(mean_mu))*a0[g]
    V <- A-V2-V3-V4-V5
    if (det(V)<1e-6) {V0[[g]]<-solve(cov(x))}else{V0[[g]]<-solve(V,tol=set.tol)}
  }
  return(V0)
}

update_Gamma<-function(a0,a3,g){
  mu_ga <- a0[g]/(a3[g])
  sigma_ga <- sqrt(1/(a3[g]))
  alpha_ga <- -mu_ga/sigma_ga
  Z_ga <- 1-pnorm(alpha_ga)
  mean_trunc_ga <- mu_ga+(dnorm(alpha_ga)-1)*sigma_ga/Z_ga
  sd_trunc_ga <- sigma_ga*sqrt(1+alpha_ga*dnorm(alpha_ga)/Z_ga-(dnorm(alpha_ga)/Z_ga)^2)
  
  trunc_norm_sample<-rtruncnorm(1, a=0,b=Inf,mean=mu_ga,sd=sigma_ga)
  if(trunc_norm_sample <= mean_trunc_ga+2*sd_trunc_ga){
    Gamma <- trunc_norm_sample
  }else{
    Gamma <- mu_ga
  }
  return(Gamma)
}

update_ISig<-function(a0,V0,g,d){
  if(a0[g]>(d+1)){
    new_ISig<-matrix(rWishart(1,a0[g],Sigma=V0[[g]]),d,d)
  }else{
    new_ISig<-matrix(rWishart(1,a0[g]+d+1,Sigma=V0[[g]]),d,d)
  }
}

update_MuBeta<-function(a0,a1,a2,a3,a4,ISig,g){
  T_mu<-ISig[[g]]*a4[g]
  T_mubeta<-ISig[[g]]*a0[g]
  T_beta<-ISig[[g]]*a3[g]
  precision_mubeta<-rbind(cbind(T_mu,T_mubeta),cbind(T_mubeta,T_beta))
  for_cov<-solve(precision_mubeta,tol=set.tol)
  sym_for_cov<-matrix(NA,2*d,2*d) 
  sym_for_cov[lower.tri(sym_for_cov,diag=TRUE)]<-for_cov[lower.tri(for_cov,diag=TRUE)]
  sym_for_cov[upper.tri(sym_for_cov)]<-t(for_cov)[upper.tri(for_cov)]
  if(isSymmetric(for_cov)==TRUE){cov_mubeta<-for_cov}else{cov_mubeta<-sym_for_cov}
  mean_mu<-a2[g,]*a3[g]/(a3[g]*a4[g]-a0[g]^2)-a1[g,]*a0[g]/(a3[g]*a4[g]-a0[g]^2)
  mean_beta<-a1[g,]*a4[g]/(a3[g]*a4[g]-a0[g]^2)-a2[g,]*a0[g]/(a3[g]*a4[g]-a0[g]^2)
  Mu_beta<-rmvnorm(1,mean=c(mean_mu,mean_beta),sigma =cov_mubeta )
  Mu<-Mu_beta[1:d]
  Beta<-matrix(Mu_beta[(d+1):(d*2)],ncol=d,nrow=1)
  return(list(Mu=Mu,Beta=Beta))
}

sample_common<-function(up_p,V0){
  up_p<-up_p
  V0<-V0
  for (g in 1:1){
    Gamma[[g]]<-update_Gamma(a0=up_p$a0,a3=up_p$a3,g=g)
    ISig[[g]]<-update_ISig(a0=up_p$a0,V0=V0,g=g,d=d)
    Mu_Beta<-update_MuBeta(a0=up_p$a0,a1=up_p$a1,a2=up_p$a2,a3=up_p$a3,a4=up_p$a4,ISig,g)
    Mu[[g]]<-Mu_Beta$Mu
    Beta[[g]]<-Mu_Beta$Beta
  }
  return(list(Gamma=Gamma[[1]],ISig=ISig[[1]],Mu=Mu[[1]],Beta=Beta[[1]]))
}

lik <- function(mu,beta,gamma,Sigma_inv,x){
  lam = (d+1)/2
  mu_o<-as.matrix(mu)
  del_o<-det(solve(Sigma_inv))^{1/(2*d)}
  gam_o<-gamma/del_o
  beta_o<-as.matrix(beta%*%Sigma_inv)
  Delta_o<-solve(Sigma_inv)/del_o^2
  alpha<-sqrt(gam_o^2+beta_o%*%Delta_o%*%t(beta_o))
  x_mu<-as.matrix(x-mu_o)
  q<-sqrt(del_o^2+t(x_mu)%*%solve(Delta_o)%*%(x_mu))
  
  first<-del_o/(2^((d-1)/2))*exp(del_o*gam_o+t(x-mu_o)%*%t(beta_o))
  if (first=="Inf") first<-1
  second<-c((alpha/(pi*q))^lam)
  third <- besselK(c(alpha*q),lam)
  if(third==0){product=0}else{
    product <- first*third*second}
  return(product)
}

#### find if there is any bad chain only valid for three chains
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

#### infinite MNIG
# data: data matrix
# true_lab: true class label
# maxiter: stop iteration if reached the maximum of iteration
# dp.alpha: numeric indicating a number as constant dp alpha parameter; or a vector(a,b) as updating along the algorithm with a gamma(a,b) prior.
# verbs: logic, print the chain numbers, iterations and current table of classifications or not
# plots: logic, trace plot or not
# outputs: "vector" or "lists", format of output estimated parameter
# CI.out: if output credible intervals or not
# dpMNIG <- function(data,true_lab,maxiter,dp.alpha,verbs=FALSE,plots=FALSE,outputs,CI.out=TRUE){
  ###########################################################################
  ###########################################################################

#######################
#### run the model ####
#######################

# procid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# procid <- 5

# data=X[[procid]][,-3]
# true_lab=X[[procid]][,3]
maxiter=1000
dp.alpha=c(2,4)
verbs=TRUE
plots=TRUE
# outputs="vectors"
# CI.out=TRUE

  #################
  #### DPMNIG #####
  #################
  # x<-data
  # true_lab<-true_lab

data(crabs)
x <- crabs[,c(4:8)]
x <- scale(x)

#   data(fish)
#   x <- fish[,c("Length2","Height","Width")]
#   x <- scale(x)
#   true_lab <- fish$Species
  
  # data(ais)
  # x <- as.matrix(ais[,c("bmi","pcBfat")])
#   x <- scale(x)
# true_lab <- ais$sex
  
  x <- as.matrix(x)
  n<-nrow(x)
  d<-ncol(x)
  it<-1
  max_it<-600
  check<-0
  G<-1
  G_lab<-1 #This will make sure that when a group is removed the label doesn't go down like G would
  lam<--(d+1)/2
  
  z<-matrix(1,ncol=1,nrow=n)
  
  ####Initialize parameters
  ISig<-list()
  Gamma<-list()
  Mu<-list()
  Beta<-list()
  
  Gamma[[1]]<-1
  Mu[[1]]<-t(sapply(1:G,function(g){colSums(x*z[,g])/sum(z[,g])}))
  Beta[[1]]<-matrix(0,ncol=d,nrow=G)
  ISig[[1]]<-solve(cov(x))
  
  ########## Get ready for extra layer common prior
  phi<-phi_fun(dat=x,G=G,Mu=Mu,ISig=ISig,n=n)
  alpha<-alpha_fun(G=G,Gamma=Gamma,Beta=Beta,ISig=ISig)
  U<-U_fun(G=G,phi=phi,alpha=alpha,lam=lam)
  Uinv<-Uinv_fun(G=G,phi=phi,alpha=alpha,lam=lam)
  
  b0 = 1/n
  b1 = cov(x)
  b2 = cov(x*as.vector(Uinv))
  b3 = 1/sum(U)
  b4 = 1/sum(Uinv)
  c1 = colSums(x)
  c2 = colSums(x*as.vector(Uinv))
  nu0 = d+1
  L0 = (cov(x))/nu0
  
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
  
  start_time <- Sys.time()
  for(chain in 1:3){
    it<-1
    if(length(dp.alpha)==1){d_alpha <- dp.alpha}
    if(length(dp.alpha)==2){d_alpha <- rgamma(1,shape=dp.alpha[1],rate=dp.alpha[2])}
    
    lam<--(d+1)/2
    
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
    a2 = t(colMeans(as.matrix(mvrnorm(10000, c2, b2))))
    a3 = mean(rexp(10000,rate = b3))
    a4 = mean(rexp(10000,rate = b4))
    forV0 = rWishart(10000,nu0,L0)
    V0_co = solve(Reduce("+",lapply(seq(dim(forV0)[3]),function(x) forV0[,,x]))/10000)
    
    up_p<-list(a0=a0,a1=a1,a2=a2,a3=a3,a4=a4)
    V0_co<-list(V0_co)
    
    #### get parameters for iteration 1
    ISig<-list()
    Gamma<-list()
    Mu<-list()
    Beta<-list()
    
    for (g in 1:G){
      com_par<-sample_common(up_p=up_p,V0=V0_co)
      Gamma[[g]]<-com_par$Gamma
      ISig[[g]]<-com_par$ISig
      Mu[[g]]<-com_par$Mu
      Beta[[g]]<-matrix(0.01,ncol=d,nrow=1)
    }
    names(Gamma)<-labels
    names(ISig)<-labels
    names(Mu)<-labels
    names(Beta)<-labels
    
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
      
      if(length(dp.alpha)==2){
        old_d_alpha <- d_alpha
        d_alpha <- Update_d_alpha(old_d_alpha,G,n,a=dp.alpha[1],b=dp.alpha[2])
      }
      
      for (i in 1:n){
        # print(i)
        ### For new group  
        sam_com<-sample_common(up_p=up_p,V0=V0_co)
        ##### Prob for new groups
        new<-d_alpha*lik(mu=sam_com$Mu,beta=sam_com$Beta,gamma=sam_com$Gamma,Sigma_inv=sam_com$ISig,x=as.matrix(x[i,]))
        ### Prob for old groups
        old<-NULL
        for (j in labels){
          if(Ng[j] == 1){
            old[j]<-d_alpha*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
          }else{
            old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
          }
        }
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
            Gamma[[Class_old[i]]]<-sam_com$Gamma
            Beta[[Class_old[i]]]<-matrix(sam_com$Beta,ncol=d)
          }
          if (!only_obs){
            # print("combine to new")
            G<-G+1
            G_lab<-G_lab+1
            Mu[[G]]<-sam_com$Mu
            ISig[[G]]<-sam_com$ISig
            Gamma[[G]]<-sam_com$Gamma
            Beta[[G]]<-matrix(sam_com$Beta,ncol=d)
            
            Class[i,it]<-paste("c",G_lab,sep="")
            names(Mu)<-c(labels,paste("c",G_lab,sep=""))
            names(ISig)<-c(labels,paste("c",G_lab,sep=""))
            names(Beta)<-c(labels,paste("c",G_lab,sep=""))
            names(Gamma)<-c(labels,paste("c",G_lab,sep=""))
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
            Beta[[Class_old[i]]]<-NULL
            Gamma[[Class_old[i]]]<-NULL
            # V0[[Class_old[i]]]<-NULL
            Ng<-table(Class[,it])
            labels<-names(table(Class[,it]))
          }
        }
      }
      z<-as.matrix(unmap(Class[,it]))
      pi_g <- colMeans(z)
      # pi_g <- rdirichlet(G, d_alpha/G)
      # plot(x,col=as.factor(Class[,it]))
      
      ###Updating U and U inverse 
      phi<-phi_fun(dat=x,G=length(Ng),n=n,Mu=Mu,ISig=ISig)
      alpha<-alpha_fun(G=G,Gamma=Gamma,Beta=Beta,ISig=ISig)
      
      U<-U_fun(G=G,phi=phi,alpha=alpha,lam=lam)
      Uinv<-Uinv_fun(G=G,phi=phi,alpha=alpha,lam=lam)
      
      ###Updating the hyperparameters
      hyper<-hyper_fun(z,U,Uinv,x)
      V0<-update_V(a0=hyper$a0,a1=hyper$a1,a2=hyper$a2,a3=hyper$a3,a4=hyper$a4,x=x,z=z,U=U,Uinv=Uinv,d=d)
      
      ###Updating the parameters
      for (g in 1:G){
        obs<-which(z[,g]==1)
        if (length(obs)>1){
          Gamma[[g]]<-update_Gamma(hyper$a0,hyper$a3,g=g)
          ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
          mu_beta<-update_MuBeta(hyper$a0,hyper$a1,hyper$a2,hyper$a3,hyper$a4,ISig,g=g)
          Mu[[g]]<-mu_beta$Mu
          Beta[[g]]<-matrix(mu_beta$Beta,ncol=d,nrow=1)
        } else {
          sam_com<-sample_common(up_p=up_p,V0=V0_co)
          Gamma[[g]]<-sam_com$Gamma
          ISig[[g]]<-sam_com$ISig
          Mu[[g]]<-sam_com$Mu
          Beta[[g]]<-sam_com$Beta
        }
      }
      
      ### Store parameters and the values for updating common hyperparameters
      gam_vec<-numeric(G)
      mu_isig_beta<-numeric(G)
      beta_isig<-matrix(nrow=G,ncol=d)
      mu_isig<-matrix(nrow=G,ncol=d)
      beta_isig_beta<-numeric(G)
      mu_isig_mu<-numeric(G)
      for (g in 1:G){
        # print(g)
        gam_vec[g]<-Gamma[[g]]
        mu_isig_beta[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
        beta_isig[g,]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]
        mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
        beta_isig_beta[g]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
        mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
      }
      
      ###Updating common hyperparameter using their posterior distributions
      b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)-sum(gam_vec)+sum(mu_isig_beta))
      if(b0_new <= 0){b0_new=b0}
      b1_new=b1
      b2_new=b2
      c1_new=c1+(colSums(beta_isig))%*%b1
      c2_new=c2+(colSums(mu_isig))%*%b2
      b3_new=1/(1/b3+0.5*(sum(beta_isig_beta)+sum(gam_vec^2)))
      b4_new=1/(1/b4+0.5*(sum(mu_isig_mu)+G))
      
      a0<-rexp(1,rate=b0_new)
      a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
      a2<-as.matrix(rmvnorm(1,mean=c2_new,sigma=b2_new),nrow=1)
      a3<-rexp(1,rate=b3_new)
      a4<-rexp(1,rate=b4_new)
      
      nu0_new=nu0+G*a0
      L0_new=solve(solve(L0)+Reduce("+",ISig))
      
      up_p<-list(a0=a0,a1=a1,a2=a2,a3=a3,a4=a4)
      V0_co<-list(riwish(nu0_new,L0_new))
      
      if(verbs){
        cat("Chain:",multi_starts[chain],"start","\n")
        cat("Iteration", it, "\n")
        print(table(Class[,it]))
      }
      
      ###Calucluate the log likelihood
      likelihood <- matrix(0,nrow=n,ncol=G)
      for (g in 1:G){
        likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],beta=Beta[[g]],gamma=Gamma[[g]],Sigma_inv=ISig[[g]],x=as.matrix(x[i,]))})
      }
      loglik <- sum(log(rowSums(likelihood)))
      all_loglik[[chain]][it] <- loglik
      
      it<-it+1
    }
    all_Class_matrix[[chain]] <- Class
    all_parameters[[chain]] <- list(phi,alpha,U,Uinv,up_p,V0_co,Gamma,Mu,Beta,ISig)
    all_classification[[chain]] <- list(z,labels,Ng,G,G_lab,d_alpha)
  }
  
  ########################
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
  if(negInfLik == 1){
    check.conv <- 1
  }else{
    check.conv <- conv[[1]][1]
  }
  check.thresh <- 1.1 
  
  if(plots){
    # plot likelihood for all chains
    plot(all_loglik[[1]],type = 'l',main=paste("check.conv is",check.conv))
    lines(all_loglik[[2]],col='red')
    lines(all_loglik[[3]],col='blue')
  }
  
  while(check.conv > check.thresh){
    it_old <- it
    max_it <- max_it+100
    
    for(chain in 1:3){
      phi <- all_parameters[[chain]][[1]]
      alpha <- all_parameters[[chain]][[2]]
      U <- all_parameters[[chain]][[3]]
      Uinv <- all_parameters[[chain]][[4]]
      up_p <- all_parameters[[chain]][[5]]
      V0_co <- all_parameters[[chain]][[6]]
      Gamma <- all_parameters[[chain]][[7]]
      Mu <- all_parameters[[chain]][[8]]
      Beta <- all_parameters[[chain]][[9]]
      ISig <- all_parameters[[chain]][[10]]
      
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
          new<-d_alpha*lik(mu=sam_com$Mu,beta=sam_com$Beta,gamma=sam_com$Gamma,Sigma_inv=sam_com$ISig,x=as.matrix(x[i,]))
          ### Prob for old groups
          old<-NULL
          # labels<-unique(Class[,it])
          for (j in labels){
            if(Ng[j] == 1){
              old[j]<-d_alpha*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
              
            }else{
              old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
            }
          }
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
              # print("generate new")
              Class[i,it]<-Class_old[i]
              Mu[[Class_old[i]]]<-sam_com$Mu
              ISig[[Class_old[i]]]<-sam_com$ISig
              Gamma[[Class_old[i]]]<-sam_com$Gamma
              Beta[[Class_old[i]]]<-matrix(sam_com$Beta,ncol=d)
            }
            if (!only_obs){
              # print("combine to new")
              G<-G+1
              G_lab<-G_lab+1
              Mu[[G]]<-sam_com$Mu
              ISig[[G]]<-sam_com$ISig
              Gamma[[G]]<-sam_com$Gamma
              Beta[[G]]<-matrix(sam_com$Beta,ncol=d)
              
              Class[i,it]<-paste("c",G_lab,sep="")
              names(Mu)<-c(labels,paste("c",G_lab,sep=""))
              names(ISig)<-c(labels,paste("c",G_lab,sep=""))
              names(Beta)<-c(labels,paste("c",G_lab,sep=""))
              names(Gamma)<-c(labels,paste("c",G_lab,sep=""))
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
              # G_lab<-G_lab-1
              Mu[[Class_old[i]]]<-NULL
              ISig[[Class_old[i]]]<-NULL
              Beta[[Class_old[i]]]<-NULL
              Gamma[[Class_old[i]]]<-NULL
              Ng<-table(Class[,it])
              labels<-names(table(Class[,it]))
            }
          }
        }
        
        z<-as.matrix(unmap(Class[,it]))
        pi_g <- colMeans(z)
        # plot(x,col=as.factor(Class[,it]))
        
        ###Updating U and U inverse 
        phi<-phi_fun(dat=x,G=length(Ng),n=n,Mu=Mu,ISig=ISig)
        alpha<-alpha_fun(G=G,Gamma=Gamma,Beta=Beta,ISig=ISig)
        
        U<-U_fun(G=G,phi=phi,alpha=alpha,lam=lam)
        Uinv<-Uinv_fun(G=G,phi=phi,alpha=alpha,lam=lam)
        
        ###Updating the hyperparameters
        hyper<-hyper_fun(z,U,Uinv,x)
        V0<-update_V(a0=hyper$a0,a1=hyper$a1,a2=hyper$a2,a3=hyper$a3,a4=hyper$a4,x=x,z=z,U=U,Uinv=Uinv,d=d)
        
        ###Updating the parameters
        for (g in 1:G){
          obs<-which(z[,g]==1)
          if (length(obs)>1){
            Gamma[[g]]<-update_Gamma(hyper$a0,hyper$a3,g=g)
            ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
            mu_beta<-update_MuBeta(hyper$a0,hyper$a1,hyper$a2,hyper$a3,hyper$a4,ISig,g=g)
            Mu[[g]]<-mu_beta$Mu
            Beta[[g]]<-matrix(mu_beta$Beta,ncol=d,nrow=1)
          } else {
            sam_com<-sample_common(up_p=up_p,V0=V0_co)
            Gamma[[g]]<-sam_com$Gamma
            ISig[[g]]<-sam_com$ISig
            Mu[[g]]<-sam_com$Mu
            Beta[[g]]<-sam_com$Beta
          }
        }
        
        ### Store parameters and the values for updating common hyperparameters
        gam_vec<-numeric(G)
        mu_isig_beta<-numeric(G)
        beta_isig<-matrix(nrow=G,ncol=d)
        mu_isig<-matrix(nrow=G,ncol=d)
        beta_isig_beta<-numeric(G)
        mu_isig_mu<-numeric(G)
        for (g in 1:G){
          # print(g)
          gam_vec[g]<-Gamma[[g]]
          mu_isig_beta[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
          beta_isig[g,]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]
          mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
          beta_isig_beta[g]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
          mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
        }
        
        ###Updating common hyperparameter using their posterior distributions
        b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)-sum(gam_vec)+sum(mu_isig_beta))
        if(b0_new <= 0){b0_new=b0}
        b1_new=b1
        b2_new=b2
        c1_new=c1+(colSums(beta_isig))%*%b1
        c2_new=c2+(colSums(mu_isig))%*%b2
        b3_new=1/(1/b3+0.5*(sum(beta_isig_beta)+sum(gam_vec^2)))
        b4_new=1/(1/b4+0.5*(sum(mu_isig_mu)+G))
        
        a0<-rexp(1,rate=b0_new)
        a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
        a2<-as.matrix(rmvnorm(1,mean=c2_new,sigma=b2_new),nrow=1)
        a3<-rexp(1,rate=b3_new)
        a4<-rexp(1,rate=b4_new)
        
        nu0_new=nu0+G*a0
        L0_new=solve(solve(L0)+Reduce("+",ISig))
        
        up_p<-list(a0=a0,a1=a1,a2=a2,a3=a3,a4=a4)
        V0_co<-list(riwish(nu0_new,L0_new))
        
        if(verbs){
          cat("Chain:",multi_starts[chain],"start","\n")
          cat("Iteration", it, "\n")
          print(table(Class[,it]))
        }
        ###Calucluate the log likelihood
        likelihood <- matrix(0,nrow=n,ncol=G)
        for (g in 1:G){
          likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],beta=Beta[[g]],gamma=Gamma[[g]],Sigma_inv=ISig[[g]],x=as.matrix(x[i,]))})
        }
        loglik <- sum(log(rowSums(likelihood)))
        all_loglik[[chain]] <- append(all_loglik[[chain]],loglik)
        
        it<-it+1
      }
      all_Class_matrix[[chain]] <- Class
      all_parameters[[chain]] <- list(phi,alpha,U,Uinv,up_p,V0_co,Gamma,Mu,Beta,ISig)
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
  for(chain in 1:3){
    store_all_pars[[chain]] <- list()
    ord_Class_labels[[chain]] <- matrix(nrow=n,ncol=500)
  }
  
  for(chain in 1:3){
    phi <- all_parameters[[chain]][[1]]
    alpha <- all_parameters[[chain]][[2]]
    U <- all_parameters[[chain]][[3]]
    Uinv <- all_parameters[[chain]][[4]]
    up_p <- all_parameters[[chain]][[5]]
    V0_co <- all_parameters[[chain]][[6]]
    Gamma <- all_parameters[[chain]][[7]]
    Mu <- all_parameters[[chain]][[8]]
    Beta <- all_parameters[[chain]][[9]]
    ISig <- all_parameters[[chain]][[10]]
    
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
        new<-d_alpha*lik(mu=sam_com$Mu,beta=sam_com$Beta,gamma=sam_com$Gamma,Sigma_inv=sam_com$ISig,x=as.matrix(x[i,]))
        ### Prob for old groups
        old<-NULL
        # labels<-unique(Class[,it])
        for (j in labels){
          if(Ng[j] == 1){
            old[j]<-d_alpha*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
            
          }else{
            old[j]<-(Ng[j]-1)*lik(mu=Mu[[j]],beta=Beta[[j]],gamma=Gamma[[j]],Sigma_inv=ISig[[j]],x=as.matrix(x[i,]))
          }
        }
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
            # print("generate new")
            Class[i,it]<-Class_old[i]
            Mu[[Class_old[i]]]<-sam_com$Mu
            ISig[[Class_old[i]]]<-sam_com$ISig
            Gamma[[Class_old[i]]]<-sam_com$Gamma
            Beta[[Class_old[i]]]<-matrix(sam_com$Beta,ncol=d)
          }
          if (!only_obs){
            # print("combine to new")
            G<-G+1
            G_lab<-G_lab+1
            Mu[[G]]<-sam_com$Mu
            ISig[[G]]<-sam_com$ISig
            Gamma[[G]]<-sam_com$Gamma
            Beta[[G]]<-matrix(sam_com$Beta,ncol=d)
            
            Class[i,it]<-paste("c",G_lab,sep="")
            names(Mu)<-c(labels,paste("c",G_lab,sep=""))
            names(ISig)<-c(labels,paste("c",G_lab,sep=""))
            names(Beta)<-c(labels,paste("c",G_lab,sep=""))
            names(Gamma)<-c(labels,paste("c",G_lab,sep=""))
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
            Beta[[Class_old[i]]]<-NULL
            Gamma[[Class_old[i]]]<-NULL
            Ng<-table(Class[,it])
            labels<-names(table(Class[,it]))
          }
        }
      }
      
      z<-as.matrix(unmap(Class[,it]))
      pi_g <- colMeans(z)
      
      ###Updating U and U inverse 
      phi<-phi_fun(dat=x,G=length(Ng),n=n,Mu=Mu,ISig=ISig)
      alpha<-alpha_fun(G=G,Gamma=Gamma,Beta=Beta,ISig=ISig)
      
      U<-U_fun(G=G,phi=phi,alpha=alpha,lam=lam)
      Uinv<-Uinv_fun(G=G,phi=phi,alpha=alpha,lam=lam)
      
      ###Updating the hyperparameters
      hyper<-hyper_fun(z,U,Uinv,x)
      V0<-update_V(a0=hyper$a0,a1=hyper$a1,a2=hyper$a2,a3=hyper$a3,a4=hyper$a4,x=x,z=z,U=U,Uinv=Uinv,d=d)
      
      ###Updating the parameters
      for (g in 1:G){
        obs<-which(z[,g]==1)
        if (length(obs)>1){
          Gamma[[g]]<-update_Gamma(hyper$a0,hyper$a3,g=g)
          ISig[[g]]<-update_ISig(hyper$a0,V0=V0,g=g,d=d)
          mu_beta<-update_MuBeta(hyper$a0,hyper$a1,hyper$a2,hyper$a3,hyper$a4,ISig,g=g)
          Mu[[g]]<-mu_beta$Mu
          Beta[[g]]<-matrix(mu_beta$Beta,ncol=d,nrow=1)
        } else {
          sam_com<-sample_common(up_p=up_p,V0=V0_co)
          Gamma[[g]]<-sam_com$Gamma
          ISig[[g]]<-sam_com$ISig
          Mu[[g]]<-sam_com$Mu
          Beta[[g]]<-sam_com$Beta
        }
      }
      
      ### Store parameters and the values for updating common hyperparameters
      gam_vec<-numeric(G)
      mu_isig_beta<-numeric(G)
      beta_isig<-matrix(nrow=G,ncol=d)
      mu_isig<-matrix(nrow=G,ncol=d)
      beta_isig_beta<-numeric(G)
      mu_isig_mu<-numeric(G)
      for (g in 1:G){
        # print(g)
        gam_vec[g]<-Gamma[[g]]
        mu_isig_beta[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
        beta_isig[g,]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]
        mu_isig[g,]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]
        beta_isig_beta[g]<-matrix(Beta[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Beta[[g]],ncol=d,nrow=1))
        mu_isig_mu[g]<-matrix(Mu[[g]],ncol=d,nrow=1)%*%ISig[[g]]%*%t(matrix(Mu[[g]],ncol=d,nrow=1))
      }
      
      ###Updating common hyperparameter using their posterior distributions
      b0_new=1/(1/b0-sum(log(pi_g))-sum(log(sapply(ISig,det))/2)-sum(gam_vec)+sum(mu_isig_beta))
      if(b0_new <= 0){b0_new=b0}
      b1_new=b1
      b2_new=b2
      c1_new=c1+(colSums(beta_isig))%*%b1
      c2_new=c2+(colSums(mu_isig))%*%b2
      b3_new=1/(1/b3+0.5*(sum(beta_isig_beta)+sum(gam_vec^2)))
      b4_new=1/(1/b4+0.5*(sum(mu_isig_mu)+G))
      
      a0<-rexp(1,rate=b0_new)
      a1<-as.matrix(rmvnorm(1,mean=c1_new,sigma=b1_new),nrow=1)
      a2<-as.matrix(rmvnorm(1,mean=c2_new,sigma=b2_new),nrow=1)
      a3<-rexp(1,rate=b3_new)
      a4<-rexp(1,rate=b4_new)
      
      nu0_new=nu0+G*a0
      L0_new=solve(solve(L0)+Reduce("+",ISig))
      
      up_p<-list(a0=a0,a1=a1,a2=a2,a3=a3,a4=a4)
      V0_co<-list(riwish(nu0_new,L0_new))
      
      if(verbs){
        cat("Chain:",multi_starts[chain],"start","\n")
        cat("Iteration", it, "\n")
        print(table(Class[,it]))
      }
      ###Calucluate the log likelihood
      likelihood <- matrix(0,nrow=n,ncol=G)
      for (g in 1:G){
        likelihood[,g] <- sapply(1:n,function(i){pi_g[g]*lik(mu=Mu[[g]],beta=Beta[[g]],gamma=Gamma[[g]],Sigma_inv=ISig[[g]],x=as.matrix(x[i,]))})
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
      ord_Beta <- list()
      ord_Gamma <- list()
      ord_ISig <- list()
      ord_class <- numeric(n)
      for (g in 1:G){
        order_ind <- which(orders == g)
        ord_Mu[[g]] <- Mu[[order_ind]]
        ord_Beta[[g]] <- Beta[[order_ind]]
        ord_Gamma[[g]] <- Gamma[[order_ind]]
        ord_ISig[[g]] <- ISig[[order_ind]]
        ord_class[which(Class[,it] == names(which(orders == g)))] = g
      }
      store_all_pars[[chain]][[it]] <- list(ord_Mu,ord_Beta,ord_Gamma,ord_ISig) 
      ord_Class_labels[[chain]][,(it+1-it_old)] <- ord_class
      
      it<-it+1
    }
    all_Class_matrix[[chain]] <- Class
    all_parameters[[chain]] <- list(phi,alpha,U,Uinv,up_p,V0_co,Gamma,Mu,Beta,ISig)
    all_classification[[chain]] <- list(z,labels,Ng,G,G_lab,d_alpha)
  }
  
  end_time <- Sys.time()
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
  
  # if(check.conv > 1.1){q()}
  
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
  find_diff <- function(x){
    a <- rep(1,3)
    if(any(is.na(x))){
      na.ind <- which(is.na(x))
      if(length(unique(x[-na.ind]))==1){
        a[na.ind]=0
      }else{a <- rep(0,3)}
    }else{
      if(x[1]!=x[2]&&x[2]==x[3]){a[1]=0}
      if(x[1]!=x[2]&&x[1]==x[3]){a[2]=0}
      if(x[1]==x[2]&&x[2]!=x[3]){a[3]=0}
      if(x[1]!=x[2]&&x[2]!=x[3]&&x[1]!=x[3]){a <- rep(0,3)}
    }
    return(a)
  }
  
  chain_conv_ind <- find_diff(all_final_G)
  convchains <- c(1,2,3)[which(chain_conv_ind==1)]
  # convchains <- c(1,2,3)[which(all_final_G==4)]
  if(length(convchains)==0){
    endlik <- c(all_loglik[[1]][it-1],all_loglik[[2]][it-1],all_loglik[[3]][it-1])
    convchains <- which.max(endlik)
  }
  
  ord_class_mat <- Reduce(cbind,ord_Class_labels[convchains])
  ord_final_class <- NULL
  for(i in 1:n){
    ord_final_class[i] <- names(which.max(table(ord_class_mat[i,])))
  }
  z<-as.matrix(unmap(ord_final_class))
  pi_g <- colMeans(z)
  ARI <- adjustedRandIndex(ord_final_class,true_lab)
  
  # summarize parameter estimation
  gamma_exp <-list()
  mu_exp <- list()
  beta_exp <- list()
  sig_exp <- list()
  class_labs <- list()
  loglik_approx_ls <- list()
  
  for (chain in convchains){
    len = length(all_loglik[[chain]])
    allparas <- store_all_pars[[chain]][(len-399):len]
    G = all_final_G[chain]
    loglik_approx_ls[[chain]] <- all_loglik[[chain]][(len-399):len]
    
    mu_exp[[chain]]<-matrix(nrow=400,ncol=d*G)
    beta_exp[[chain]]<-matrix(nrow=400,ncol=d*G)
    gamma_exp[[chain]]<-matrix(nrow=400,ncol=G)
    sig_exp[[chain]]<-matrix(nrow=400,ncol=d*d*G)
    
    for(i in 1:400){
      for(g in 1:G){
        mu.try <- try(allparas[[i]][[1]][[g]],silent=TRUE)
        if(is.numeric(mu.try)&!(class(mu.try) == "try-error")){
          mu_exp[[chain]][i,(1+d*(g-1)):(d*g)]<-mu.try
        }
        beta.try <- try(allparas[[i]][[2]][[g]],silent=TRUE)
        if(is.numeric(beta.try)&!(class(beta.try) == "try-error")){
          beta_exp[[chain]][i,(1+d*(g-1)):(d*g)]<-beta.try
        }
        gamma.try<-try(allparas[[i]][[3]][[g]],silent=TRUE)
        if(is.numeric(gamma.try)&!(class(gamma.try) == "try-error")){
          gamma_exp[[chain]][i,g]<-gamma.try
        }
        sig.try<-try(as.vector(solve(allparas[[i]][[4]][[g]])),silent=TRUE)
        if(is.numeric(sig.try)&!(class(sig.try) == "try-error")){
          sig_exp[[chain]][i,(1+d*d*(g-1)):(d*d*g)]<-sig.try
        }
      }
    }
  }
  
  mu <- apply(Reduce(rbind,mu_exp),2,function(x){mean(x,na.rm=TRUE)})
  beta <- apply(Reduce(rbind,beta_exp),2,function(x){mean(x,na.rm=TRUE)})
  gamma <- apply(Reduce(rbind,gamma_exp),2,function(x){mean(x,na.rm=TRUE)})
  sig <- apply(Reduce(rbind,sig_exp),2,function(x){mean(x,na.rm=TRUE)})
  G=unique(all_final_G[convchains])
  npar <- G+2*(d)*G+0.5*(d+1)*(d)*G+(G-1)
  BIC <- 2*loglik_approx - npar*log(n)
  
  ### get the empirical 0.025-qurtile and .975-quartile
  gam_l_var<-numeric(G)
  gam_u_var<-numeric(G)
  mu_l_var<-numeric(d*G)
  mu_u_var<-numeric(d*G)
  beta_l_var<-numeric(d*G)
  beta_u_var<-numeric(d*G)
  Sig_l_var<-numeric(d*d*G)
  Sig_u_var<-numeric(d*d*G)
  
  for(j in 1:G){
    num.ests <- sum(!is.na(Reduce(rbind,gamma_exp)[,j]))
    l.num <- floor(num.ests*0.025)
    u.num <- ceiling(num.ests*0.975)
    gam_l_var[j] <- Reduce(rbind,gamma_exp)[which(rank(Reduce(rbind,gamma_exp)[,j])==l.num),j]
    gam_u_var[j] <- Reduce(rbind,gamma_exp)[which(rank(Reduce(rbind,gamma_exp)[,j])==u.num),j]
  }
  
  for(j in 1:(d*G)){
    num.ests <- sum(!is.na(Reduce(rbind,mu_exp)[,j]))
    l.num <- floor(num.ests*0.025)
    u.num <- ceiling(num.ests*0.975)
    mu_l_var[j] <- Reduce(rbind,mu_exp)[which(rank(Reduce(rbind,mu_exp)[,j]) == l.num),j]
    mu_u_var[j] <- Reduce(rbind,mu_exp)[which(rank(Reduce(rbind,mu_exp)[,j]) == u.num),j]
  }
  
  for(j in 1:(d*G)){
    num.ests <- sum(!is.na(Reduce(rbind,beta_exp)[,j]))
    l.num <- floor(num.ests*0.025)
    u.num <- ceiling(num.ests*0.975)
    beta_l_var[j] <- Reduce(rbind,beta_exp)[which(rank(Reduce(rbind,beta_exp)[,j]) == l.num),j]
    beta_u_var[j] <- Reduce(rbind,beta_exp)[which(rank(Reduce(rbind,beta_exp)[,j]) == u.num),j]
  }
  
  for(j in 1:(d*d*G)){
    num.ests <- sum(!is.na(Reduce(rbind,sig_exp)[,j]))
    l.num <- floor(num.ests*0.025)
    u.num <- ceiling(num.ests*0.975)
    Sig_l_var[j] <- Reduce(rbind,sig_exp)[which(rank(Reduce(rbind,sig_exp)[,j]) == l.num),j]
    Sig_u_var[j] <- Reduce(rbind,sig_exp)[which(rank(Reduce(rbind,sig_exp)[,j]) == u.num),j]
  }
  
  Mu <- list()
  Beta <- list()
  Gamma <- list()
  Sigma <- list()
  
  for(g in 1:G){
    Mu[[g]] <- mu[(1+d*(g-1)):(d*g)]
    Beta[[g]] <- beta[(1+d*(g-1)):(d*g)]
    Gamma[[g]] <- gamma[g]
    Sigma[[g]] <- matrix(sig[(1+d*d*(g-1)):(d*d*g)],nrow=d,ncol=d)
  }
  ISig <- lapply(Sigma,solve)
  
  par.vec=list(mu=mu,beta=beta,gamma=gamma,sigma=sig,final.G=G,ARI=ARI)
  par.list=list(mu=Mu,beta=Beta,gamma=Gamma,sigma=Sigma,final.G=G,ARI=ARI)
  CIs=list(mu.ci=rbind(mu_l_var,mu_u_var),beta.ci=rbind(beta_l_var,beta_u_var),gamma.ci=rbind(gam_l_var,gam_u_var),sigma.ci=rbind(Sig_l_var,Sig_u_var))
 
  runtime <- end_time-start_time
###################################
##### Save all results ############
###################################
all_results <- c(gamma,mu,beta,sig,G,ARI,runtime)
gam_results <- c(gam_l_var,gam_u_var,gamma)
mu_results <- c(mu_l_var,mu_u_var,mu)
beta_results <- c(beta_l_var,beta_u_var,beta)
Sigma_results <- c(Sig_l_var,Sig_u_var,sig)

# save.image("/usr2/postdoc/yuanf/test/DPMNIG/fish_DPMNIG.RData")
save.image("C:/Users/echof/Dropbox/Work/DPMNIG_real/crabs_MING_new.RData")


# output.filename = paste("/usr2/postdoc/yuanf/test/DPMNIG/sim5/dalphaup/simulationOutput",".csv",sep='')
# write.table(t(all_results),file=output.filename,sep=",",row.names=paste("simulation",procid),col.names=F,append = TRUE)
# outputbic.filename = paste("/usr2/postdoc/yuanf/test/DPMNIG/sim5/dalphaup/GammaCI",".csv",sep='')
# write.table(t(gam_results),file=outputbic.filename,sep=",",row.names=paste("simulation",procid),col.names=F,append = TRUE)
# outputbic.filename = paste("/usr2/postdoc/yuanf/test/DPMNIG/sim5/dalphaup/MuCI",".csv",sep='')
# write.table(t(mu_results),file=outputbic.filename,sep=",",row.names=paste("simulation",procid),col.names=F,append = TRUE)
# outputbic.filename = paste("/usr2/postdoc/yuanf/test/DPMNIG/sim5/dalphaup/BetaCI",".csv",sep='')
# write.table(t(beta_results),file=outputbic.filename,sep=",",row.names=paste("simulation",procid),col.names=F,append = TRUE)
# outputbic.filename = paste("/usr2/postdoc/yuanf/test/DPMNIG/sim5/dalphaup/SigmaCI",".csv",sep='')
# write.table(t(Sigma_results),file=outputbic.filename,sep=",",row.names=paste("simulation",procid),col.names=F,append = TRUE)  