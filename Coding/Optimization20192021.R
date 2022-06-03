
rm(list = ls())
library(quadprog)
library(xts)
library(zoo)
library(quantmod)
library(lubridate)
library(PerformanceAnalytics)
library(rpart)
library(writexl)
library(openxlsx)
library(Matrix)


setwd("C:/Users/lyand/Desktop/FE800/Optimization/Data")


################################################################################
## 1.Data Handling

#1) Investment Universe

# 11 Sector ETFs in tickers

tickers <- c("XLK","XLY","XLV","XLF","XLI","XLP","XLU","XLB","XLE","XLC","XLRE")

# Download 11 ETFs and SPY

getSymbols(tickers,from="2018-08-21", to="2021-12-30")
getSymbols("SPY",from="2018-08-21",to="2021-12-30")

P <- merge(XLK$XLK.Close,XLY$XLY.Close,XLV$XLV.Close,XLF$XLF.Close,XLI$XLI.Close,
           XLP$XLP.Close,XLU$XLU.Close,XLB$XLB.Close,XLE$XLE.Close,XLC$XLC.Close,
           XLRE$XLRE.Close)

colnames(P) <- tickers

# Make the function for Return Matrix

Rfun <- function(p) {
  
  na.omit((p/lag(p))-1)

  }

R <- Rfun(P)

rm(XLK,XLY,XLV,XLF,XLI,XLP,XLU,XLB,XLE,XLC,XLRE)
SPY <- SPY$SPY.Close

RSPY <- Rfun(SPY)
colnames(RSPY) <- c("SPY")

#2) Download French Fama 3factors

# Fa is % data

Fa <- read.csv("Fama_from180822to211229.csv")
Fa <- Fa[,-1]

# Fama model is percent data, so it should be divided by 100

Fa <- Fa/100


Date <- index(R)
rownames(Fa) <- Date
colnames(Fa) <- c("RMRF","SMB","HML","RF")
Fa <- xts(Fa,order.by = Date)


#3) Download Regime

RG <- read.xlsx("Label.xlsx",1)

# original z3 is 2nd column, predicted z3 is 3rd column
# original vix is 4th column, predicted vix is 4th column
RG <- RG[,2]


RG <- as.matrix(RG)
colnames(RG) <- c("Regime")
RG <- xts(RG,order.by = Date)


################################################################################
## 2.Factor model for the Whole Period

# Make the R, Fa DataFrame

R <- as.data.frame(R)
Fa <- as.data.frame(Fa)
FD <- Fa # FD is copy version Fama factors model

FD1 <- cbind(rep(1,length(FD[,1])),FD[,1:3])
colnames(FD1) <- c("1","RMRF","SMB","HML")

###### Setting for regression

# RiRf = Return of ETFs - RiskFree Rate(Daily)
# Coef = Coefficient of Factor model
# Beta = CAPM beta using SPY

RiRf <- R-FD[,"RF"]
Coef <- as.data.frame(matrix(0,nrow=11,ncol=4))
Beta <- as.data.frame(matrix(0,nrow=11,ncol=1))

colnames(Beta) <- c("Beta")
rownames(Beta) <- tickers

colnames(Coef) <- c("Alpha","RMRF","SMB","HML")
rownames(Coef) <- tickers

Coef_M <- Coef
Coef_Q <- Coef


  

################################################################################
## 3.ProcessData function


processdata <- function(LB_Mu,LB_Q,startday){
  
  # 1. Set the data&value for getting Mu and Q and Beta
  
  # Startday is the point
  # LBst_M, LBed_M are the start and end day of lookback period for Mu
  # LBst_Q, LBed_Q are the start and end day of lookback period for Q
  
  point <- which(Date == startday)
  LBst_M <- point - LB_Mu
  LBed_M <- point - 1
  LBst_Q <- point - LB_Q
  LBed_Q <- point - 1
  
  # Fama model's lookback period for Mu and Q 
  
  FD_M <- FD[LBst_M:LBed_M,]
  FD1_M <- FD1[LBst_M:LBed_M,]
  FD_Q <- FD[LBst_Q:LBed_Q,]
  FD1_Q <- FD1[LBst_Q:LBed_Q,]
  
  # Return of ETFs for Mu and Q
  
  R_M <- R[LBst_M:LBed_M,]
  R_Q <- R[LBst_Q:LBed_Q,]
  
  # RiRf_* is 'Return of ETFs - Riskfree rate' for lookback period for Mu and Q
  RiRf_M <- RiRf[LBst_M:LBed_M,]
  RiRf_Q <- RiRf[LBst_Q:LBed_Q,]

  # RSPY_* is Return of SPY for lookback period for Mu, Q 
  RSPY_M <- RSPY[LBst_M:LBed_M,]
  RSPY_Q <- RSPY[LBst_Q:LBed_Q,]
  
  
  #2 Regression of Fama 3 factor model & CAPM Beta
  
  for (i in 1:11) {
    
    Coef_M[i,] <- coefficients(lm(RiRf_M[,i] ~ FD_M[,"RMRF"]+FD_M[,"SMB"]+FD_M[,"HML"]))
    Beta[i,] <- coefficients(lm(RiRf_M[,i] ~ RSPY_M))[2]
    
    Coef_Q[i,] <- coefficients(lm(RiRf_Q[,i] ~ FD_Q[,"RMRF"]+FD_Q[,"SMB"]+FD_Q[,"HML"]))
    
  }
  
  
  #3 Calculate the Mu
  
  Rho <- as.data.frame(as.matrix(FD1_M) %*% t(as.matrix(Coef_M)))
  Mu <- sapply(Rho,mean) * 252

  #4 Calculate the Covariance Q
  
  Pred_R_Q <- as.matrix(FD1_Q) %*% t(as.matrix(Coef_Q))
  Idio <- diag(cov(R_Q - Pred_R_Q))
  D <- matrix(diag(Idio),ncol=11)
  Q <- as.matrix(Coef[,2:4]) %*% cov(as.matrix(FD[,1:3])) %*% t(as.matrix(Coef[,2:4])) * 252 + D
  
  save(Mu,Q,Beta,file="inputs.rData")
  
}

processdata(60,90,"2020-01-02")

load("inputs.rData")
Mu
Q
Beta


################################################################################
## 4.Detect Regime function

# Last 5 days, The mean of regime is higher than 0 it is crisis, the other is normal

DetectRegime <- function(startday){
  
  point <- which(Date == startday)
  LBst_V <- point - 5
  LBed_V <- point - 1
  
  V_Mean <- mean(RG[LBst_V:LBed_V])
  
  if (V_Mean > 0.21) {
    
    Regime <- 1
    
  } else 
    
    Regime <- 0
  
  return(Regime)

  }

DetectRegime("2020-02-26")

################################################################################
## 5. Optimizer function

# When the Regime is 0, find the optimal weight
# When the Regime is not 0, the all weight is 0 

# Optimizer is the Target Beta is not band (Static AA)
# Optimizer1 has the different objective function (Dynamic AA)

Optimizer <- function(lambda,T_Beta,wp) {
  
  load("Inputs.rData")
  
  if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
  } else {
      
      Dmat <- as.matrix(2*lambda*Q)
  } 
    
    dvec <- as.matrix(Mu)
    Amat <- as.matrix(cbind(rep(1,11),Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta,rep(0,11),rep(-0.33,11))  # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
    rownames(w) <- tickers
    
  return(w)

    }


Optimizer1 <- function(lambda,T_Beta,wp,Regime) {
  
  load("Inputs.rData")
  
  Mu
  
  if (Regime < 0.1) {
    
    if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
    } else {
      
      Dmat <- as.matrix(2*lambda*Q)
    } 
    
    dvec <- as.matrix(Mu) + as.matrix(2* lambda * Q) %*% as.matrix(wp)
    Amat <- as.matrix(cbind(rep(1,11),Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta,rep(0,11),rep(-0.33,11))  # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
    rownames(w) <- tickers
    
  } else {
    
    w <- matrix(c(rep(0,11)),ncol=1)
    rownames(w) <- tickers
    
  }
  
  return(w)
  
}

########### Excute if the code work well ###############

wp <- matrix(c(0.1,0.2,-0.3,0.1,0.2,-0.1,0.2,0.4,-0.2,0.5,-0.1),ncol=1)
wp <- Optimizer(1,1,wp,1)
wp
sum(wp)

wp <- Optimizer(1,1,wp,0)
wp
sum(wp)


###############################################################################
## 6. Realized Portfolio Return

# 1) Static Asset Allocation (No considering Regime)
# Using Optimizer

R <- xts(R,order.by=Date)

Port_Return1 <- function(startday,endday,LB_Mu,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex + 1   # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex + 1
  
  Nrebal <- PlayTime/5
  
  RR <- R[startindex:endindex,] # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/11,11),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  
  for (i in 1:Nrebal) {
    
  
    processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    
    wp <- Optimizer(lambda,T_Beta,wp) # Get new portfolio wp from the Beta, Mu, Q
    W <- rbind(W,t(wp)) # Record the new portfolio
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    
    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
  
    }
  
  W <- xts(W,order.by=index(RR))
  
  ret <- apply(W*RR,1,sum)
  ret <- as.matrix(ret)
  ret <- xts(ret,order.by = index(W))

  save(W,RR,ret,file="WRr.rData")
  
  return(ret)
  
  }

ret1 <- Port_Return1("2019-01-02","2021-12-29",60,90,1,1)

Port_Return1("2019-01-02","2021-12-29",60,90,1,1)
load("WRr.rData")
W1 <- W

Return.cumulative(ret1)
chart.CumReturns(ret1)



###############################################################################
# 2) Stop Investing Strategy

R <- xts(R,order.by=Date)

Port_Return2 <- function(startday,endday,LB_Mu,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex + 1   # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex + 1
  
  Nrebal <- PlayTime/5
  
  RR <- R[startindex:endindex,] # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/11,11),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  
  for (i in 1:Nrebal) {
    
    
    processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    
    Regime <- DetectRegime(startday)
    
    wp <- Optimizer1(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
    W <- rbind(W,t(wp)) # Record the new portfolio
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    W <- rbind(W,t(wp))
    
    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
    
  }
  
  W <- xts(W,order.by=index(RR))
  
  ret <- apply(W*RR,1,sum)
  ret <- as.matrix(ret)
  ret <- xts(ret,order.by = index(W))
  
  save(W,RR,ret,file="WRr.rData")
  
  return(ret)
  
}

ret2 <- Port_Return2("2019-01-02","2021-12-29",60,90,1,1)

Port_Return2("2019-01-02","2021-12-29",60,90,1,1)
load("WRr.rData")
W2 <- W

W1_W2 <- W1-W2

plot(W1_W2[,1])

K_Data <- cbind(ret1,ret2)


Return.cumulative(K_Data)
chart.CumReturns(K_Data)


###############################################################################
# 2-1) Investing GLD&TLT when regime is crisis

Port_Return2("2019-01-02","2021-12-29",60,90,1,1)

load("WRr.rData")

W <- apply(W,1,sum)
W <- xts(W,order.by=index(ret))

colnames(W) <- c("Weight")

ret0 <- matrix(rep(0,length(W)),ncol=1)
ret0 <- xts(ret0,order.by=index(W))
ret0

getSymbols("GLD",from="2018-12-31",to="2021-12-30")
getSymbols("TLT",from="2018-12-31",to="2021-12-30")

P_GT <- merge(GLD$GLD.Close,TLT$TLT.Close)
R_GT0 <- Rfun(P_GT)
R_GT <- 0.5*(R_GT0[,1] + R_GT0[,2])
R_GT


for (i in 1:length(W)) {
  
  if (W[i] == 0) {
    
    ret0[i] <- R_GT[i] 
    } 
  }

ret21 <- ret0 + ret

K_Data <- cbind(ret1,ret2,ret21)


Return.cumulative(K_Data)
chart.CumReturns(K_Data)


###############################################################################
# 3) Changing Lookback Period Strategy -> Drop!!

Optimizer3 <- function(lambda,T_Beta,wp,Regime) {
  
  load("Inputs.rData")
  
  if (Regime < 0.1) {
    
    if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
    } else {
      
      Dmat <- as.matrix(2*lambda*Q)
    } 
    
    dvec <- as.matrix(Mu) + as.matrix(2* lambda * Q) %*% as.matrix(wp)
    Amat <- as.matrix(cbind(rep(1,11),Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta,rep(0,11),rep(-0.33,11))  
    
    # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
    rownames(w) <- tickers
    
  } else {
    
    if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
    } else {
      
      Dmat <- as.matrix(2*lambda*Q)
    } 
    
    dvec <- as.matrix(Mu)
    Amat <- as.matrix(cbind(rep(1,11),Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta,rep(0,11),rep(-0.33,11))  # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
    rownames(w) <- tickers
    
  }
  
  return(w)
  
}


Port_Return3 <- function(startday,endday,LB_Mu,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex + 1   # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex + 1
  
  Nrebal <- PlayTime/5
  
  RR <- R[startindex:endindex,] # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/11,11),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  
  for (i in 1:Nrebal) {
    
    Regime <- DetectRegime(startday)
    
    if (Regime > 0.1) {
      
      LB_Mu <- LB_Mu / 3
      LB_Q <- LB_Q / 3
      
      processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
      
      wp <- Optimizer3(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
      LB_Mu <- LB_Mu * 3
      LB_Q <- LB_Q * 3
      
    } else {
      
      processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
      
      wp <- Optimizer3(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
    }
    
    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
    
  }
  
  W <- xts(W,order.by=index(RR))
  
  ret <- apply(W*RR,1,sum)
  ret <- as.matrix(ret)
  ret <- xts(ret,order.by = index(W))
  
  save(W,RR,ret,file="WRr.rData")
  
  return(ret)
  
}

ret3 <- Port_Return3("2019-01-02","2021-12-29",60,90,1,1)

Port_Return3("2019-01-02","2021-12-29",60,90,1,1)
load("WRr.rData")
W3 <- W

W1_W3 <- W1 - W3
plot(W1_W3)

K_Data <- cbind(ret1,ret2,ret21,ret3)


Return.cumulative(K_Data)
chart.CumReturns(K_Data)



###############################################################################
# 4) Beta Neutral Strategy (Short is allowed)

# Optimizer4 is the optimization when short is allowed

Optimizer4 <- function(lambda,T_Beta,wp,Regime) {
  
  load("Inputs.rData")
  
  if (Regime < 0.1) {
    
    if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
    } else {
      
      Dmat <- as.matrix(2*lambda*Q)
    } 
    
    dvec <- as.matrix(Mu) + as.matrix(2* lambda * Q) %*% as.matrix(wp)
    Amat <- as.matrix(cbind(rep(1,11),Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta,rep(0,11),rep(-0.33,11))  # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
    rownames(w) <- tickers
    
  } else {
    
    if (det(as.matrix(2*lambda*Q)) < 0.1) {
      
      Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
      
    } else {
      
      Dmat <- as.matrix(2*lambda*Q)
    } 
    
    dvec <- as.matrix(Mu)
    Amat <- as.matrix(cbind(rep(1,11),Beta,-Beta,diag(11),-diag(11)))
    bvec <- c(1,T_Beta-0.25,-T_Beta-0.25,rep(-0.66,11),rep(-0.66,11))  # bvec = c(sum of weight,target beta, wi lower bound, wi upper bound)
    
    w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=1)$sol,ncol=1)
    rownames(w) <- tickers
    
  }
  
  return(w)
  
}


Port_Return4 <- function(startday,endday,LB_Mu,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex + 1   # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex + 1
  
  Nrebal <- PlayTime/5
  
  RR <- R[startindex:endindex,] # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/11,11),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  
  for (i in 1:Nrebal) {
    
    
    processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    
    Regime <- DetectRegime(startday)
    
    if (Regime < 0.1){
      
      wp <- Optimizer4(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
    } else {
      
      temp_T_Beta <- T_Beta
      T_Beta <- 0
      
      wp <- Optimizer4(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
      T_Beta <- temp_T_Beta
      
    }

    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
    
  }
  
  W <- xts(W,order.by=index(RR))
  
  ret <- apply(W*RR,1,sum)
  ret <- as.matrix(ret)
  ret <- xts(ret,order.by = index(W))
  
  save(W,RR,ret,file="WRr.rData")
  
  return(ret)
  
}

ret4 <- Port_Return4("2019-01-02","2021-12-29",60,90,1,1)

Port_Return4("2019-01-02","2021-12-29",60,90,1,1)
load("WRr.rData")
W4 <- W

W1_W4 <- W1 - W4
plot(W1_W4)


Return.cumulative(ret4)
chart.CumReturns(ret4)

K_Data <- cbind(ret1,ret2,ret21,ret3,ret4)


Return.cumulative(K_Data)
chart.CumReturns(K_Data)

###############################################################################
# 5) Change Lambda


Port_Return5 <- function(startday,endday,LB_Mu,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex + 1   # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex + 1
  
  Nrebal <- PlayTime/5
  
  RR <- R[startindex:endindex,] # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/11,11),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  
  for (i in 1:Nrebal) {
    
    
    processdata(LB_Mu,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    
    Regime <- DetectRegime(startday)
    
    if (Regime < 0.1){
      
      wp <- Optimizer3(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
    } else {
      
      temp_lambda <- lambda
      lambda <- lambda * 10
      
      wp <- Optimizer3(lambda,T_Beta,wp,Regime) # Get new portfolio wp from the Beta, Mu, Q
      W <- rbind(W,t(wp)) # Record the new portfolio
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      W <- rbind(W,t(wp))
      
      lambda <- temp_lambda
      
    }
    
    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
    
  }
  
  W <- xts(W,order.by=index(RR))
  
  ret <- apply(W*RR,1,sum)
  ret <- as.matrix(ret)
  ret <- xts(ret,order.by = index(W))
  
  save(W,RR,ret,file="WRr.rData")
  
  return(ret)
  
}

ret5 <- Port_Return5("2019-01-02","2021-12-29",60,90,1,1)

Port_Return5("2019-01-02","2021-12-29",60,90,1,1)
load("WRr.rData")
W5 <- W

W1_W5 <- W1 - W5
plot(W1_W5)



Return.cumulative(ret5)
chart.CumReturns(ret5)

K_Data <- cbind(ret1,ret2,ret21,ret3,ret4,ret5)


Return.cumulative(K_Data)
chart.CumReturns(K_Data)



###############################################################################
# 5) Summary

#Hr1 <- hist(ret1,plot=FALSE)
#Hr2 <- hist(ret2,plot=FALSE)
#Hr21 <- hist(ret21,plot=FALSE)
#Hr3 <- hist(ret3,plot=FALSE)
#Hr4 <- hist(ret4,plot=FALSE)

#plot(Hr1, col="lightblue",main="Static AA vs Stop Inv")
#plot(Hr2, col="pink",add = TRUE)

#plot(Hr1, col="lightblue",main="Static AA vs GLD+TLT")
#plot(Hr21, col="pink",add = TRUE)

#plot(Hr1, col="lightblue",main="Static AA vs Beta Neutral")
#plot(Hr4, col="pink",add = TRUE)

colnames(K_Data) <- c("Static AA","Stop Invest","GLD+TLT St","Change Lookback","Beta Neutral","Change lambda") 

CumRet <- sapply(K_Data,Return.cumulative)
Mean_K <- sapply(K_Data,mean)*252 # Annualized Mean Return Before Crisis
Vol_K <- sapply(K_Data,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_K <- sapply(K_Data,skewness)
Kur_K <- sapply(K_Data,kurtosis)
SR_K <- sapply(K_Data,SharpeRatio.annualized)
DD_K <- sapply(K_Data,maxDrawdown)
#DD_K <- min(DD_K)
#data.frame(DD_K)
#DD_K <- sapply(DD_K,min)



Return.annualized(K_Data)
#chart.Boxplot(K_Data)
#chart.Correlation(K_Data)
chart.Drawdown(K_Data,wealth.index = TRUE,legend.loc="bottomleft")
#chart.RiskReturnScatter(K_Data)
#charts.PerformanceSummary(K_Data)

Table_K <- rbind(CumRet,Mean_K,Vol_K,Ske_K,Kur_K,SR_K,DD_K)
rownames(Table_K) <- c("Cum Return","Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio","Max Draw Down")


chart.CumReturns(K_Data,wealth.index = TRUE,legend.loc = "topleft")

Table_K <- as.data.frame(Table_K)
Table_K

setwd("C:/Users/lyand/Desktop/FE800/Optimization/Code")

write.xlsx(K_Data,"z3_covid.xlsx")
write.xlsx(Table_K,"Table_z3_covid.xlsx")







