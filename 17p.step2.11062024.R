#Author: Arti Virkud
#Code purpose: Replicate results of Prolong and Meis with 17P data
#Last updated:#03/03/2025 Subgroup analyses with smoking/drinking/substance use, PPROM in previous delivery,
#Secondary outcomes: Spontaneous delivery, PPROM current delivery
#01/01/2025 Add Rscript system arguments
#10/21/2024 Add bootstrapping to accommodate superlearner weights

packages <- c("haven","dplyr","tidyr","expss","lubridate","geex","tidyverse", "boot","caret", "SuperLearner", "ranger", "glmnet", "caret", "randomForest","xgboost",
              "nnet", "caret")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}

# Transport findings from PROLONG to Meis
# •	Use methods and R code described in Dahabreh et al. Statistics in Medicine 2020 to transport results using (1) standardization, (2) inverse odds weighting, and (3) doubly robust estimators.
# •	Measured effect modifiers: Those in Table 2, as well as race/ethnicity and country (US vs. non-US)
#Treating PROLONG as trial participants and Meis and non-participants
setwd("~/Library/CloudStorage/Box-Box/17P Trial Reanalysis/Data and Code/Code")
load("~/Library/CloudStorage/Box-Box/17P Trial Reanalysis/Data and Code/Code/mdf.R")
load("~/Library/CloudStorage/Box-Box/17P Trial Reanalysis/Data and Code/Code/pdf.R")

#Print column names
#labs <- sapply(df, var_lab)
#print(labs)

#Data setup
#Covariates X1: Maternal Age, X2: Maternal Race:
# 1=Black, 2=White, 3=Hispanic, 4=Asian, 5=Other
#X3: Gestational age at qualifying prior SPTB
#X4: >1 prior spontaneous PTB>1, X5: Marital status, X6: Prepregnancy BMI, X7: Years of education
#X8: Smoking during pregnancy, X9: Alcohol use during pregnancy, X10: Substance use during pregnancy
#X11: Region, X12: Missing data indicator
#Outcomes: Y=primary outcome, Y2=spontaneous delivery, Y3=PPROM delivery

#The first argument denotes whether we are transporting from Meis to PROLONG or PROLONG to Meis
#1: Meis to PROLONG
#0: PROLONG to Meis
#The second argument denotes which sensitivity analysis we are running
#0: Primary analysis
#1: US only
#2: >.05 P(Meis|Covariates) and >.1 P(PROLONG|Covariates)
#3: Individuals who have had at least 1 ptb
#4: Individuals who haven't drank/smoked/used substances during pregnancy
#5: Individuals who have drank/smoked/used substances during pregnancy
#6: Individuals who have had a previous PPROM
#7: Individuals who have not had a previous PPROM
#8: Secondary outcome: spontaneous delivery
#9: Secondary outcome: PPROM delivery

args <- commandArgs(trailingOnly = TRUE)
#args<-c(1,1)
setup1A <- mdf.all %>% rename(Y=PRE37, Y2=SPRE37, Y3=APROM37, X1=CXAGE, X3=BXGSTAGE, X4=PREVPT1, X4_c=PRETERM, X5=MARRIED, X6=BMI, X7=CXEDUC)  %>% 
                      mutate (S=as.numeric(args[1]), 
                              C_90=case_when(COMP>=90~1,.default=0), C_95=case_when(COMP>=95~1,.default=0),
                              A=case_when(TREAT=="17P"~1, TREAT=="PLACEBO"~0),
                              X2=case_when(CXRACE=="1"~1,CXRACE=="2"~2,CXRACE=="3"~3,
                                           CXRACE=="4"~4,CXRACE=="5"|CXRACE=="6"~5),
                              X8=case_when(CXSMOKE=="Y"~1,CXSMOKE=="N"~0),
                              X9=case_when(CXALCHOL=="Y"~1,CXALCHOL=="N"~0),
                              X10=case_when(CXDRUGS=="Y"~1,CXDRUGS=="N"~0),
                              X11=1, X12=case_when(is.na(X2)|is.na(X3)|is.na(X4)|is.na(X5)|
                                                     is.na(X6)|is.na(X7)|is.na(X8)|is.na(X9)|
                                                     is.na(X10)~1, .default = 0),
                              X13=case_when(DXDELCLS %in% c(3,4)~1,.default=0)) %>% 
          select(-PATID,-TREAT,-CXRACE,-CXSMOKE,-CXALCHOL,-CXDRUGS,-COMP,-FXDCLASS,-DXDELCLS,-CXINFEC) 
#setup1A %>% filter(!is.na(X6)) %>% summarise(median=median(X6))
#BMI mean 24.6

setup1 <- setup1A %>% mutate(X6=case_when(is.na(X6)~24.6, .default=X6))

setup2A <- pdf.all %>% rename(Y=outcome, Y2=outcome_spont, Y3=outcome_PPROM, X1=AGE, X3=GADEL, X4_c=PRETERM, X6=PREBMI,
                             X7=EDLEVEL, C_90=C10_90, C_95=C10_95)  %>% 
  mutate (S=(1-as.numeric(args[1])), A=case_when(ACTARM=="17P"~1, ACTARM=="Vehicle"~0),
          X2=case_when((RACE=="BLACK OR AFRICAN AMERICAN" & ETHNIC=="NOT HISPANIC OR LATINO")~1,
                       (RACE =="WHITE" & ETHNIC=="NOT HISPANIC OR LATINO")~2,
                       (ETHNIC=="HISPANIC OR LATINO")~3,
                       (RACE=="ASIAN" & ETHNIC=="NOT HISPANIC OR LATINO")~4,
                       (RACE=="AMERICAN INDIAN OR ALASKA NATIVE" & ETHNIC=="NOT HISPANIC OR LATINO")~5,
                       ((RACE=="OTHER" | RACE=="NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER"|RACE=="MIXED RACE") & ETHNIC=="NOT HISPANIC OR LATINO")~5),
          X4=case_when(preterm==">1 previous preterm delivery"~1, preterm=="1 previous preterm delivery"~0,is.na(preterm)~NA),
          X5=case_when(MARISTAT=="Married/living with partner"~1, 
                       MARISTAT=="Never married"|MARISTAT=="Divorced/widowed/separated"~0),
          X8=case_when(TOBACCO=="Y"~1,is.na(TOBACCO)~0),
          X9=case_when(ALCOHOL=="Y"~1,is.na(ALCOHOL)~0),
          X10=case_when(DRUG=="Y"~1,is.na(DRUG)~0),
          X11=case_when(REGION1=="US"~1,REGION1=="Non-US"~0),
          X12=case_when(is.na(X2)|is.na(X3)|is.na(X4)|is.na(X5)|
                          is.na(X6)|is.na(X7)|is.na(X8)|is.na(X9)|
                          is.na(X10)~1, .default = 0),
          X13=PPROM_prev) %>% 
  select(A,Y,Y2,Y3,S,C_90,C_95,X1,X2,X3,X4,X4_c,X5,X6,X7,X8,X9,X10,X11,X12,X13)

#setup2A %>% filter(!is.na(X3)) %>% summarise(median=median(X3))
#gestaional age median: 33
#setup2A %>% filter(!is.na(X6)) %>% summarise(median=median(X6))
#Prepregnancy BMI median 24.5
#setup2A %>% filter(!is.na(X7)) %>% summarise(median=median(X7))
#Years of education median: 13.0
setup2 <- setup2A %>% mutate(X3=case_when(is.na(X3)~33, .default=X3),
                             X4=case_when(is.na(X4)~0, .default=X4),
                             X4_c=case_when(is.na(X4_c)~1, .default=X4_c),
                             X6=case_when(is.na(X6)~23.0, .default=X6),
                             X7=case_when(is.na(X7)~13.0, .default=X7))
#Set up missing indicator
#Original analysis excludes missing outcome n=28: 24 from PROLONG and 4 from Meis
#Create missing indicator variable for remaining covariates
#Identify reference values

#Subgroup/Sensitivity Analyses:
if (as.numeric(args[2])==0) {
DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
  mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
         X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0))
covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")
}

if (as.numeric(args[2])==1) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X11==1)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")
}

  if (as.numeric(args[2])==2) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0))
  S<-glm(formula=S~X1+X22+X23+X24+X25+X3+X4_c+X5+X6+X7+X8+X9+X10+X12, family=binomial("logit"), data=DF)
  S.hist <- predict(S,newdata=DF, type="response")
  DF1 <- DF
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")
  if (as.numeric(args[1])==1){
  #When Meis S=1, and Prol S=0
  DF1$pMeis <- S.hist
  DF1$pProl <- 1-S.hist
  DF <- DF1 %>% filter(S==1 & pMeis>=.05|S==0) %>% select (-pMeis,-pProl)}
  if (as.numeric(args[1])==0){
  # When Meis S=0, and Prol S=1, not doing this anymore given the .05 threshold 11/05/24
  DF1$pMeis <- 1-S.hist
  DF1$pProl <- S.hist
  DF <- DF1 %>% filter(S==0 |S==1 & pProl>=.1) %>% select (-pMeis,-pProl)}}

if (as.numeric(args[2])==3) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X4==0)
  covlist <- c("X1","X22","X23","X24","X25","X3","X5","X6","X7","X8","X9","X10","X12")}

if (as.numeric(args[2])==4) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X8==0 & X9==0 & X10==0)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X12")}

if (as.numeric(args[2])==5) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X8==1 | X9==1 | X10==1)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")}

if (as.numeric(args[2])==6) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X13==1)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")}

if (as.numeric(args[2])==7) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% filter(X13==0)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")}

if (as.numeric(args[2])==8) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y2)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% select(-Y) %>% rename(Y=Y2)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")}

if (as.numeric(args[2])==9) {
  DF <-rbind(setup1,setup2) %>% filter(!is.na(Y3)) %>% 
    mutate(X22=ifelse(X2==2,1,0), X23=ifelse(X2==3,1,0),
           X24=ifelse(X2==4,1,0), X25=ifelse(X2==5,1,0)) %>% select(-Y) %>% rename(Y=Y3)
  covlist <- c("X1","X22","X23","X24","X25","X3","X4_c","X5","X6","X7","X8","X9","X10","X12")}

DF$Y[is.na(DF$Y)] <- 0
DF$A[is.na(DF$A)] <- 0
DF$epsilon<-NULL
DF$inter<-NULL
DF$d_pop<-NULL
DF$d_out<-NULL
DF$Y_1<-NULL
DF$Y_0<-NULL

#Updated everything with superlearner
#Set seed
set.seed(2024)

#Set up cross-fitting for datasets
  S1data<-subset(DF, S==1)
  S1d <- createDataPartition(S1data$Y, p=0.5, list=FALSE)
  S1d_cfit1 <- S1data[S1d,]
  S1d_cfit2 <- S1data[-S1d,]
  S0data<-subset(DF, S==0)
  S0d = createDataPartition(S0data$Y, p=0.5, list=FALSE)
  S0d_cfit1 <- S0data[S0d,]
  S0d_cfit2 <- S0data[-S0d,]
  Sd_cfit1 <- rbind(S1d_cfit1,S0d_cfit1)
  Sd_cfit2 <- rbind(S1d_cfit2,S0d_cfit2)

#Set up datsets to run superlearner
  data.X <- Sd_cfit1 %>% select(one_of(covlist))
  S1data.X <- S1d_cfit1 %>% select(one_of(covlist))
  data.X2 <- Sd_cfit2 %>% select(one_of(covlist))
  S1data.X2 <- S1d_cfit2 %>% select(one_of(covlist))
  
#Set superlearner library
  sl_lib = c("SL.xgboost", "SL.ranger", "SL.glmnet", "SL.nnet", "SL.glm","SL.gbm", "SL.gam")
#Train on cfit1, Predict on cfit2
  w_reg<-SuperLearner(Y = Sd_cfit1$S, X = data.X , family = binomial(), SL.library = sl_lib)
  ps1 <- predict(w_reg, newdata=data.X2)$pred
  w_reg2<-SuperLearner(Y = S1d_cfit1$A, X = S1data.X , family = binomial(), SL.library = sl_lib)
  pa1 <- predict(w_reg2, newdata=data.X2)$pred
  w1 <- (Sd_cfit2$A*Sd_cfit2$S*(1-ps1) )/(ps1*pa1) + ((1 -Sd_cfit2$A) *Sd_cfit2$S*(1-ps1) ) /(ps1*(1-pa1))
  Sd_cfit2$w<-w1
  
#Train on cfit2, Predict on cfit1
  w_reg3<-SuperLearner(Y = Sd_cfit2$S, X = data.X2 , family = binomial(), SL.library = sl_lib)
  ps2 <- predict(w_reg3, newdata=data.X)$pred
  w_reg4<-SuperLearner(Y = S1d_cfit2$A, X = S1data.X2 , family = binomial(), SL.library = sl_lib)
  pa2 <- predict(w_reg4, newdata=data.X)$pred
  w2 <- (Sd_cfit1$A*Sd_cfit1$S*(1-ps2) )/(ps2*pa2) + ((1 -Sd_cfit1$A) *Sd_cfit1$S*(1-ps2) ) /(ps2*(1-pa2))
  Sd_cfit1$w<-w2

#Outcome model
#OM fit on cfit1, OM predict on cfit2
S1data_A1fit1<-subset(S1d_cfit1, S==1 & A==1)
S1data_A1fit1.X2 <- S1data_A1fit1 %>% select(one_of(covlist))
OM1mod_c1<-SuperLearner(Y = S1data_A1fit1$Y, X = S1data_A1fit1.X2 , family = binomial(), SL.library = sl_lib)
p1_c1 <- predict(OM1mod_c1, newdata=data.X2)$pred
Sd_cfit2$p1<-p1_c1

S1data_A0fit1<-subset(S1d_cfit1, S==1 & A==0)
S1data_A0fit1.X2 <- S1data_A0fit1 %>% select(one_of(covlist))
OM0mod_c1<-SuperLearner(Y = S1data_A0fit1$Y, X = S1data_A0fit1.X2 , family = binomial(), SL.library = sl_lib)
p0_c1 <- predict(OM0mod_c1, newdata=data.X2)$pred
Sd_cfit2$p0<-p0_c1


#OM fit on cfit2, OM predict on cfit1
S1data_A1fit2<-subset(S1d_cfit2, S==1 & A==1)
S1data_A1fit2.X2 <- S1data_A1fit2 %>% select(one_of(covlist))
OM1mod_c2<-SuperLearner(Y = S1data_A1fit2$Y, X = S1data_A1fit2.X2 , family = binomial(), SL.library = sl_lib)
p1_c2 <- predict(OM1mod_c2, newdata=data.X)$pred
Sd_cfit1$p1<-p1_c2

S1data_A0fit2<-subset(S1d_cfit2, S==1 & A==0)
S1data_A0fit2.X2 <- S1data_A0fit2 %>% select(one_of(covlist))
OM0mod_c2<-SuperLearner(Y = S1data_A0fit2$Y, X = S1data_A0fit2.X2 , family = binomial(), SL.library = sl_lib)
p0_c2 <- predict(OM0mod_c2, newdata=data.X)$pred
Sd_cfit1$p0<-p0_c2

#S0sub_cfit2<-subset(Sd_cfit2, S==0)

data <- rbind(Sd_cfit1, Sd_cfit2)
S0sub_data<-subset(data, S==0)
OM_1_cfit2<-mean(S0sub_data$p1)
OM_0_cfit2<-mean(S0sub_data$p0)
OM_cfit2<-mean(S0sub_data$p1)-mean(S0sub_data$p0)
OM<-list(OM_1=OM_1_cfit2, OM_0=OM_0_cfit2, OM=OM_cfit2, p1=data$p1, p0=data$p0, OM1mod_c2=OM1mod_c2, OM0mod_c1=OM0mod_c1)

weights2<-list(dat=data, Smod=w_reg3, Amod=w_reg4)
DF2<-weights2$dat

DF2$p1<-OM$p1
DF2$p0<-OM$p0

#Inverse odds weighting
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
IOW2<-IOW2_est(data=DF2)

#Doubly robust estimation
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
  list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
  return(list)
}
DR2<-DR2_est(data=data)

print <- c(OM$OM_1,OM$OM_0,OM$OM,IOW2$IOW2_1,IOW2$IOW2_0,IOW2$IOW2,DR2$DR2_1,DR2$DR2_0,DR2$DR2)

# setting up the results directory
results_dir <- "PrimaryR/"
# check if results sub directory exists and create if not
if (!dir.exists(results_dir)){dir.create(results_dir)}
if (as.numeric(args[1])==1){
    save(DF, file=paste(results_dir,"DF_M2P_S",as.numeric(args[2]),".Rdata",sep=""))
  # Save the  results
  save(print, file=paste(results_dir, "Results_M2P_S", as.numeric(args[2]), ".Rdata", sep = ""))
  }
if (as.numeric(args[1])==0){
      save(DF, file=paste(results_dir,"DF_P2M_S",as.numeric(args[2]),".Rdata",sep=""))
  # Save the  results
  save(print, file=paste(results_dir, "Results_P2M_S", as.numeric(args[2]), ".Rdata", sep = ""))}

