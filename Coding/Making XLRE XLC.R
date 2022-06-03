

rm(list = ls())
library(quadprog)
library(xts)
library(zoo)
library(quantmod)
library(lubridate)
library(PerformanceAnalytics)
library(rpart)
library(writexl)
library(xlsx)
library(Matrix)

setwd("C:/Users/lyand/Desktop/FE800/Optimization/Data")

################################################################################
## 1.Data Handling

#1) DownLoad Data

ticker1 <- c("AMT","PLD","CCI","EQIX","PSA","SPG","DLR","SBAC","WELL","AVB")
ticker2 <- c("GOOGL","GOOG","NFLX","CMCSA","T","DIS","VZ")

getSymbols(ticker1,from="2006-08-22", to="2010-12-29")
getSymbols(ticker2,from="2006-08-22", to="2010-12-29")


#length(AMT["2007-01-03/2010-12-28",1])
#length(AMT["2006-08-22/2006-12-31",1])

PXLRE <- merge(AMT$AMT.Close,PLD$PLD.Close,CCI$CCI.Close,EQIX$EQIX.Close,PSA$PSA.Close,
           SPG$SPG.Close,DLR$DLR.Close,SBAC$SBAC.Close,WELL$WELL.Close,AVB$AVB.Close)

PXLC <- merge(GOOGL$GOOGL.Close,GOOG$GOOG.Close,NFLX$NFLX.Close,CMCSA$CMCSA.Close,
              DIS$DIS.Close,T$T.Close,VZ$VZ.Close)

rm(AMT,PLD,CCI,EQIX,PSA,SPG,DLR,SBAC,WELL,AVB)
rm(GOOGL,GOOG,NFLX,CMCSA,T,DIS,VZ)

colnames(PXLRE) <- ticker1
colnames(PXLC) <- ticker2

# Make the function for Return Matrix

Rfun <- function(p) {
  
  na.omit((p/lag(p))-1)
  
}

RXLRE <- Rfun(PXLRE)
RXLC <- Rfun(PXLC)

RXLRE <- 0.2*RXLRE[,1]+0.15*RXLRE[,2]+0.15*RXLRE[,3]+0.12*RXLRE[,4]+0.08*RXLRE[,5]+0.07*RXLRE[,6]+
         0.07*RXLRE[,7]+0.055*RXLRE[,8]+0.055*RXLRE[,9]+0.05*RXLRE[,10]

RXLC <- 0.25*RXLC[,1]+0.25*RXLC[,2]+0.1*RXLC[,3]+0.1*RXLC[,4]+0.1*RXLC[,5]+0.1*RXLC[,6]+0.1*RXLC[,7]

colnames(RXLRE) <-c("XLRE")
colnames(RXLC) <-c("XLC")

write.xlsx(RXLRE,"XLRE.xlsx")
write.xlsx(RXLC,"XLC.xlsx")
