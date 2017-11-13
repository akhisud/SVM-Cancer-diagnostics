funkyduck <- function(i.,T2,Part, Factors, minval = 5, mincount = 1, FDRVal = 0.00005, logCPMval = 3, cpmVal = 1, cpmCount = 5){
  library(edgeR)
  library(RWeka)
  T <- T2[,Part]
  # Filtering of genes
  keep <- rowSums(as.matrix(T)>=minval) >= mincount
  keep2 <- rowSums(cpm(T)>cpmVal) >= cpmCount
  keep <- keep&keep2
  rm(keep2)
  #Removes the Y chromosome genes
  #YList <- read.csv("/home/manesh/Manesh/R/Untitled Folder/GeneListY")    ###
  #YList <- YList[,1]
  for(j in YList){
    keep[rownames(T)==j] <- FALSE
  }
  rm(j)
  T <- DGEList(T, group = factor(Factors[Part]))  ###
  T <- T[keep, , keep.lib.sizes=FALSE]
  T <- calcNormFactors(T)
  rm(keep)
  # Create model
  design <- model.matrix(~factor(Factors[Part]))
  rownames(design) <- colnames(T)
  # Estimate dispersions, fit model and perform GLM likelihood ratio test
  T <- estimateGLMCommonDisp(T,design = design)
  T <- estimateGLMTrendedDisp(T,design = design)
  T <- estimateGLMTagwiseDisp(T,design = design)
  fit <- glmFit(T, design = design)
  lrt <- glmLRT(fit)
  # Second round of selection for DEG
  FDR <- p.adjust(lrt$table$PValue, method="BH")
  keep <- FDR < FDRVal
  keep2 <- lrt$table$logCPM>=logCPMval
  keep <- keep&keep2
  rm(keep2)
  keep3 <- lrt$table$logFC>log(1.2)     ###
  keep4 <- lrt$table$logFC<log(0.8)     ###
  keep3 <- (keep4|keep3)
  keep <- keep&keep3
  rm(keep3)
  rm(keep4)
  #dt <- decideTestsDGE(lrt[keep,], lfc=1)      # Can use logFC as a selection measure
  # Exporting the csv files for use
  DEnames <- rownames(T)[keep]
  T.cpm <- cpm(T, normalized.lib.size=TRUE)
  T.cpm <- t(T.cpm[DEnames,])
  T.cpm <- cbind(T.cpm,Factors[Part])
  colnames(T.cpm)[ncol(T.cpm)]<-"Status"
  write.csv(T.cpm,file = paste(i., "Training.csv", sep=""))
  T <- T2[,-Part]
  T.cpm <- cpm(T, normalized.lib.size=TRUE)
  T.cpm <- t(T.cpm[DEnames,])
  T.cpm <- cbind(T.cpm,Factors[-Part])
  colnames(T.cpm)[ncol(T.cpm)]<-"Status"
  write.csv(T.cpm,file = paste(i., "Validation.csv", sep=""))
  rm(T.cpm,T,keep,lrt,fit,DEnames,smo,e)
  #### Insert code over here for the relevant function
}
library(foreach)
library(doParallel)
library(parallel)
noOfSplits <- 10
registerDoParallel(cores = 5)
#T <- read.table("/home/manesh/Manesh/GSE68086_TEP_data_matrix.txt")  ###
foreach(i = 1:noOfSplits) %dopar% {
  set.seed(i)
  library(caret)
  Part <- createDataPartition(y = Factors.bi,times = 1,p = 0.7,list = FALSE)
  system.time(funkyduck(i.=i,T2 = T.,Part = Part, Factors = Factors.bi))
}
stopImplicitCluster()