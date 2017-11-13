funkyduck <- function(i. = i){
  library(RWeka)
  library(caret)
  library(edgeR)
  
  
  T.cpm <- read.csv(paste(i.,"Training.csv",sep=""))
  T.cpm <- T.cpm[,-1]
  
  
  smo <- SMO(Status ~ ., data = as.data.frame(T.cpm),control = Weka_control(C = 1.0, L = 0.001, P = 1.0E-12, N = 1, V = -1, W = 1, K = list("weka.classifiers.functions.supportVector.PolyKernel -E 1.0 -C 250007"), calibrator = list("weka.classifiers.functions.Logistic -R 1.0E-8 -M -1 -num-decimal-places 4")))
  
  
  T.cpm <- read.csv(paste(i.,"Validation.csv",sep=""))
  T.cpm <- T.cpm[,-1]
  
  
  pred <- predict(smo,newdata = as.data.frame(T.cpm))
  p <- table(pred, T.cpm[,ncol(T.cpm)])
  #save(pred ,file = paste(i., " pred.txt", sep=""))
  write.csv(as.matrix(p) ,file = paste(i., " CM.csv", sep=""))
  rm(T.cpm,T,keep,lrt,fit,DEnames,smo,e)
  gc()
  
}


library(foreach)
library(doParallel)
library(parallel)

registerDoParallel(cores = 5)
noOfSplits <- 10
foreach(i = 1:noOfSplits) %dopar% {
  system.time(funkyduck(i. = i))
}






for(i in 1:10){
  t <- read.csv(paste(i," CM.csv",sep=""))
  a <- t[1,2]/(t[1,2] +t[2,2])
  b <- (t[1,2] + t[2,3])/(sum(t[,2:3]))
  k[i,] <- cbind(b,a)
}
write.csv(k,file = "SMO.csv")
