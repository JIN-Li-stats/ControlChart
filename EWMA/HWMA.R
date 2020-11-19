# program HWMA
# Homogeneously weighted moving average control chart with an application in substrate manufacturing process
# DOI: 10.1016/j.cie.2018.05.009
# R version 4.0.3
# 2020/10/03

rm(list=ls())
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(magrittr)

n <- 1; loop <- 10000

l <- 2.272; r <- 0.03

mu0 <- 0; sigma0 <- 1

j <- 1
output <-matrix(1:80, nrow=10, ncol=8)
rl <- vector(mode = "numeric", length = loop)

t1=proc.time()
# for (shift in c(seq(0.00, 0.75, 0.25))){
for (shift in c(seq(0.00, 1.00, 0.25), seq(1.50, 3.00, 0.50), 5.00)){
  for (i in 1:loop){
    k <- 0
    hwma <- 0
    count <- 0
    
    repeat{   
     count <- count + 1
      x <- rnorm(n)+shift
      xb <- mean(x)
      if(count==1){
        hwma <- r*xb
        k <- k+xb
        UCL <- mu0 + l*sqrt(r^2*sigma0^2/n)
        LCL <- mu0 - l*sqrt(r^2*sigma0^2/n)
      }else{
        hwma <-r*xb+((1.0-r)*k)/(count-1.0)
        k <- k+xb
        UCL <- mu0 + l*sqrt((r^2)*(sigma0^2)/n + ((1.0-r)^2)*(sigma0^2)/(n*(count-1.0)))
        LCL <- mu0 - l*sqrt((r^2)*(sigma0^2)/n + ((1.0-r)^2)*(sigma0^2)/(n*(count-1.0)))
      }
      
      if(hwma > UCL | hwma < LCL | count == 8000){
        rl[i] <- count
        break
      }
    }
  }
  
  arl <- mean(rl)
  sdrl <- sd(rl)
  quan <- quantile(rl,c(0.05,0.25,0.5,0.75,0.95)) %>% as.numeric()
  output[j,] <- c(shift,arl,sdrl,quan)
  j <- j+1
}
print(output)

t2=proc.time()
print(t2-t1)
