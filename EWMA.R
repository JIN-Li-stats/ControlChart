# EWMA program (Steady state)
# ARL0 = 370

rm(list = ls())
n <- 1
j <- 1
l <- 2.302
r <- 0.03
loop <- 10000
h <- l * sqrt(r / (2 - r))


rl <- vector(mode = "numeric", length = loop)
output <-matrix(1:81, nrow=9, ncol=9)
colnames(output)<-c("shift","ARL","SDRL","0.05","0.10","0.25","0.5","0.75","0.95")

for (shift in seq(0, 2, 0.25)) {
   if (shift == 0.0) {
      w <- 0
   }else {
      w <- 0
   }

   for (i in 1:loop) {
      ewma <- 0
      count <- 0
      
      while (ewma < h & ewma > -h & count < 8000) {
         x <- rnorm(n)
         count <- count + 1
         if (count > w) {
            x <- x + shift
         }
         y <- sum(x) / n
         ewma  <- (1 - r) * ewma + r * y
      }
      if (count > w) {
         rl[i] <- count - w
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
#write.table(output, file = "D:\\OneDrive\\SPC\\Code\\output.csv",sep = ",", row.names = FALSE, col.names = TRUE)
