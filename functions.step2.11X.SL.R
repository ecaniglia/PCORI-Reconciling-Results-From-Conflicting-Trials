packages <- c("haven","dplyr","tidyr","expss","lubridate","geex","tidyverse", "SuperLearner", "ranger", "glmnet", "caret", "randomForest","xgboost",
              "nnet", "caret")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}
#11 covariates
#Weight estimation in 17p.step2.04212024.R
#Older weight generation procedure
generate_weights<-function(Smod,Amod, data){
  S1data<-subset(data, S==1)
  w_reg<-glm(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response") 
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response") 
  w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
  data$w<-w
  list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
  return(list)
}



#functions that do not account for estimation of the working models 
#most appropriate for bootstrap estimation of SEs
#Updated with additional covariates 3.26.24
Cov_set<-as.formula(paste("Y~X1+X22+X23+X24+X25+X3+X4_c+X5+X6+X7+X8+X9+X10+X12"))

OM_est<-function(data){
  S1data_A1<-subset(data, S==1 & A==1)
  OM1mod<-glm(formula=Cov_set, family="binomial", data=S1data_A1)
  p1<- predict(OM1mod,newdata=data, type="response") 
  data$p1<-p1
  S1data_A0<-subset(data, S==1 & A==0)
  OM0mod<-glm(formula=Cov_set, family="binomial", data=S1data_A0)
  p0<- predict(OM0mod, newdata=data, type="response") 
  data$p0<-p0
  S0sub<-subset(data, S==0)
  OM_1<-mean(S0sub$p1)
  OM_0<-mean(S0sub$p0)
  OM<-mean(S0sub$p1)-mean(S0sub$p0)
  list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
  return(list)
}

IOW2_est<-function(data){
  S0data<-subset(data, S==0)  
  S1data_A1<-subset(data, S==1 & A==1)
  IOW1mod<-glm(formula=Y~1, family="binomial", data=S1data_A1, weights=w)
  p1<- predict(IOW1mod,newdata=S0data, type="response") 
  S1data_A0<-subset(data, S==1 & A==0)
  IOW0mod<-glm(formula=Y~1, family="binomial", data=S1data_A0, weights=w)
  p0<- predict(IOW0mod,newdata=S0data, type="response") 
  IOW2_1<-mean(p1)
  IOW2_0<-mean(p0)
  IOW2<-mean(p1)-mean(p0)
  list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
  return(list)
}


DR2_est<-function(data){
  A<-data$A
  S<-data$S
  Y<-data$Y
  p1<-data$p1
  p0<-data$p0
  w<-data$w
  sum1_DR2<-sum(S*A*w*(Y-p1)) 
  sum0_DR2<-sum(S*(1-A)*w*(Y-p0)) 
  norm1<-(sum(S*A*w))^-1
  norm0<-(sum(S*(1-A)*w))^-1
  DR2_1<-norm1*sum1_DR2 + (sum(1-S)^-1)*sum((1-S)*p1)
  DR2_0<-norm0*sum0_DR2 + (sum(1-S)^-1)*sum((1-S)*p0)
  DR2<-DR2_1-DR2_0
  # se1_DR2<-sd(norm1*(S*A*w*(Y-p1)) + (sum(1-S)^-1)*((1-S)*p1))
  # se0_DR2<-sd(norm0*(S*(1-A)*w*(Y-p0)) + (sum(1-S)^-1)*((1-S)*p0))
  # se_DR2<-sd((norm1*(S*A*w*(Y-p1)) + (sum(1-S)^-1)*((1-S)*p1))-(norm0*(S*(1-A)*w*(Y-p0)) + (sum(1-S)^-1)*((1-S)*p0)))
  # se1_DR2<-sd(norm1*(S*A*w*(Y-p1)) + (sum(1-S)^-1)*((1-S)*p1))/sqrt(length((data)))
  # se0_DR2<-sd(norm0*(S*(1-A)*w*(Y-p0)) + (sum(1-S)^-1)*((1-S)*p0))/sqrt(length((data)))
  # se_DR2<-sd((norm1*(S*A*w*(Y-p1)) + (sum(1-S)^-1)*((1-S)*p1))-(norm0*(S*(1-A)*w*(Y-p0)) + (sum(1-S)^-1)*((1-S)*p0)))/sqrt(length((data)))
  # list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2, se_1=se1_DR2, se_0=se0_DR2, se=se_DR2)
  list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
  return(list)
}

#functions for M-estimation (geex)

OM_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X22, data$X23, data$X24, data$X25, data$X3, data$X4_c, data$X5, data$X6, data$X7, data$X8, data$X9, data$X10, data$X12) #update for covariate change
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  function(theta){ 
    #outcome model 
    beta<-theta[1:15] #update for covariate change
    alpha<-theta[16:30] #update for covariate change
    mu1<-theta[31] #update for covariate change
    mu0<-theta[32] #update for covariate change
    muate<-theta[33] #update for covariate change
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    #estimates
    mean1<-(1-S)*(m_A1-mu1) 
    mean0 <- (1-S)*(m_A0-mu0) 
    ate<-(1-S)*(m_A1-m_A0-muate) 
    c(ols_A1,ols_A0,mean1, mean0,ate)
  }
}

IOW2_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X22, data$X23, data$X24, data$X25, data$X3, data$X4_c, data$X5, data$X6, data$X7, data$X8, data$X9, data$X10, data$X12) #update for covariate change
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  function(theta){
    #participation model
    lp  <- X %*% theta[1:15] #update for covariate change
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[16:30] #update for covariate change
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model
    m_A1<-1 %*% theta[31] #update for covariate change
    m_A0<-1 %*% theta[32] #update for covariate change
    linear_eqns1<-crossprod(1, (S*A*w)*(Y -  m_A1) )
    linear_eqns0<-crossprod(1, (S*(1-A)*w)*(Y - 1 %*% theta[32]) ) #update for covariate change
    mu1<-theta[33] #update for covariate change
    mu0<-theta[34] #update for covariate change
    muate<-theta[35] #update for covariate change
    #estimates
    mean1 <- (1-S)*( m_A1 -mu1) 
    mean0 <- (1-S)*( m_A0 -mu0)
    ate <- (1-S)*(m_A1- m_A0 - muate)
    c(score_eqns,score_eqns2, linear_eqns1,linear_eqns0,  mean1,  mean0,ate)
  }
}

DR2_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X22, data$X23, data$X24, data$X25, data$X3, data$X4_c, data$X5, data$X6, data$X7, data$X8, data$X9, data$X10, data$X12) #update for covariate change
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  function(theta){
    #participation model
    lp  <- X %*% theta[1:15] #update for covariate change
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[16:30] #update for covariate change
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) ) 
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model 
    beta<-theta[31:45] #update for covariate change
    alpha<-theta[46:60] #update for covariate change
    mu1<-theta[64] #update for covariate change
    mu0<-theta[65] #update for covariate change
    mu<-theta[66] #update for covariate change
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    mu_S<-theta[61] #update for covariate change
    propS1<-S-mu_S
    one_over<-(1/(1-mu_S))
    mu_norm1<-theta[62] #update for covariate change
    norm1eq<-(A*S*w)-mu_norm1
    norm1<-1/mu_norm1
    mu_norm0<-theta[63] #update for covariate change
    norm0eq<-((1-A)*S*w)-mu_norm0
    norm0<-1/mu_norm0
    ey1<-norm1*((w*S*A*(Y-m_A1))) + one_over*((1-S)*m_A1)
    ey0<-norm0*((w*S*(1-A)*(Y-m_A0))) + one_over*((1-S)*m_A0)
    #estimates
    mean1<-ey1-mu1
    mean0<-ey0-mu0
    ate<-ey1-ey0-mu
    c(score_eqns,score_eqns2,ols_A1,ols_A0,propS1, norm1eq, norm0eq,mean1,mean0,ate)
  }
}
#Function to extract point estimate and SE from geex output
extractEST<-function(geex_output=OM_mest, est_name="m1",param_start=param_start_OM){
  param_num_EST<-match(est_name,names(param_start))
  EST<-geex_output@estimates[param_num_EST]
  sandwich_se <- diag(geex_output@vcov)^0.5 
  SE<-sandwich_se[param_num_EST]
  return(c(EST, SE=SE))
}
