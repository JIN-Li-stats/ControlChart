# Adaptive CUSUM
# Article title: CUSUM Charts for Signalling Varying Location Shifts
# DOI: 10.1080/00224065.2000.11979987
# ARL0=400

rm(list = ls())

n <- 1
j <- 1
loop <- 1000

lambda <- 0.05
hz <- 0.9959    #lambda=0.05
#hz <- 0.9896    #lambda=0.10

qu_min <- 0.5
ql_max <- -0.5

qu_mean <- 1
ql_mean <- -1

rl <- vector(mode = "numeric", length = loop)
output <-matrix(1:81, nrow=9, ncol=9)
colnames(output)<-c("delta","ARL","SDRL","0.05","0.10","0.25","0.5","0.75","0.95")


for (shift in c(0.00,0.25,0.50,0.75,1.00,1.50,2.00,2.50,3.00)) {
  
  for (i in 1:loop) {
    CUSUM_U <- 0
    CUSUM_L <- 0
    
    qu <- qu_mean
    ql <- ql_mean
    
    count <- 0
    
    repeat{

      count <- count + 1
      
      z=rnorm(1,mean=0,sd=1) + shift
      
      #ARL=400
      hU=-5 + 37.69862/sqrt(1+15.405*qu-1.6150*qu^2) + 0.4978428*qu - 1.251752*qu^2 + 1.32498*qu^3 - 0.6796388*qu^4 + 0.169046*qu^5 - 0.01645198*qu^6
      hL=-5 + 37.69862/sqrt(1+15.405*(-ql)-1.6150*(-ql)^2) + 0.4978428*(-ql) - 1.251752*(-ql)^2 + 1.32498*(-ql)^3 - 0.6796388*(-ql)^4 + 0.169046*(-ql)^5 - 0.01645198*(-ql)^6
      
      CUSUM_L <- min(0, CUSUM_L+(z-ql/2)/hL)
      CUSUM_U <- max(0, CUSUM_U+(z-qu/2)/hU)
      
      qu <- max(lambda*z+(1-lambda)*qu, qu_min)
      ql <- min(lambda*z+(1-lambda)*ql, ql_max)
      
      if(CUSUM_U==0){
        qu <- qu_mean
      }
      if(CUSUM_L==0){
        ql <- ql_mean
      }
      
      if(CUSUM_L < -hz | CUSUM_U > hz | count==10000){
        rl[i] <- count
        break
      }
    }
  }
  arl <- mean(rl)
  sdrl <- sd(rl)
  quan <- quantile(rl,c(0.05,0.10,0.25,0.5,0.75,0.95))
  quan <- as.numeric(quan)
  
  output[j,] <- c(shift,arl,sdrl,quan)
  j <- j+1
}

print(output)
write.table(output, file = "D:\\OneDrive\\SPC\\Code\\Adaptive CUSUM Program\\output.csv",sep = ",", row.names = FALSE, col.names = TRUE)

