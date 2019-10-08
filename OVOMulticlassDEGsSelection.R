# Loading data (previously obtained from GeneExpressionIntegration R script)
load("FULL_DATASET_NORMALIZED_FINAL_27_SERIES_28_BATCHES.RData")
colnames(data_normalized) <- sub("\\.CEL$","",colnames(data_normalized))

# ***************************************************************************************************************************
# Sample separation by classes (each skin pathological state (SPS))
# ***************************************************************************************************************************
NSK.samples <- data_normalized[,which(labels.dataset.series %in% 1)]
NEV.samples <- data_normalized[,which(labels.dataset.series %in% 2)]
BCC.samples <- data_normalized[,which(labels.dataset.series %in% 3)]
ISCC.samples <- data_normalized[,which(labels.dataset.series %in% 4)]
PMCC.samples <- data_normalized[,which(labels.dataset.series %in% 5)]
MMCC.samples <- data_normalized[,which(labels.dataset.series %in% 6)]
PRIMEL.samples <- data_normalized[,which(labels.dataset.series %in% 7)]
METMEL.samples <- data_normalized[,which(labels.dataset.series %in% 8)]
AK.samples <- data_normalized[,which(labels.dataset.series %in% 9)]
PS.samples <- data_normalized[,which(labels.dataset.series %in% 10)]

# ***************************************************************************************************************************
# Random shuffle for KFOLDs (in our case, 10-CV), ensuring similar representativeness of each SPS within each KFOLD
# ***************************************************************************************************************************
require(Biobase)
KFOLDS = 10
size.each.class <- c(dim(NSK.samples)[2],dim(NEV.samples)[2],dim(BCC.samples)[2],dim(ISCC.samples)[2],dim(PMCC.samples)[2],
                     dim(MMCC.samples)[2],dim(PRIMEL.samples)[2],dim(METMEL.samples)[2],dim(AK.samples)[2],dim(PS.samples)[2])

random.classes.folds <- list()
matrix.size.each.fold <- as.data.frame(matrix(0, nrow = length(size.each.class), ncol = KFOLDS))
for(i in 1:length(size.each.class)){
  distribution.samples.folds <- split(sample(1:size.each.class[i]),rep(1:KFOLDS,ceiling(size.each.class[i]/KFOLDS)))
  size.each.fold <- listLen(distribution.samples.folds) 
  distribution.samples.folds <- list(distribution.samples.folds)
  random.classes.folds <- c(random.classes.folds,distribution.samples.folds)
  matrix.size.each.fold[i,] <- size.each.fold
}

save(random.classes.folds,matrix.size.each.fold, file = "PARTITIONS_FOLDS_28_BATCHES.RData")

# ***************************************************************************************************************************
# Differential Expression Analysis within 10-CV procedure, tuning LFC and NMAX
# ***************************************************************************************************************************
require(limma)
for(NUMFOLD in seq(1,10,1)){     ### NUMFOLD represents each experiment (in our case, M = 10 experiments)
  training_data <- cbind(NSK.samples[,-as.numeric(unlist(random.classes.folds[[1]][NUMFOLD]))],
                         NEV.samples[,-as.numeric(unlist(random.classes.folds[[2]][NUMFOLD]))],
                         BCC.samples[,-as.numeric(unlist(random.classes.folds[[3]][NUMFOLD]))],
                         ISCC.samples[,-as.numeric(unlist(random.classes.folds[[4]][NUMFOLD]))],
                         PMCC.samples[,-as.numeric(unlist(random.classes.folds[[5]][NUMFOLD]))],
                         MMCC.samples[,-as.numeric(unlist(random.classes.folds[[6]][NUMFOLD]))],
                         PRIMEL.samples[,-as.numeric(unlist(random.classes.folds[[7]][NUMFOLD]))],
                         METMEL.samples[,-as.numeric(unlist(random.classes.folds[[8]][NUMFOLD]))],
                         AK.samples[,-as.numeric(unlist(random.classes.folds[[9]][NUMFOLD]))],
                         PS.samples[,-as.numeric(unlist(random.classes.folds[[10]][NUMFOLD]))])
  
  training_labels <- c(rep(1,size.each.class[1]-matrix.size.each.fold[1,NUMFOLD]),
                       rep(2,size.each.class[2]-matrix.size.each.fold[2,NUMFOLD]),
                       rep(3,size.each.class[3]-matrix.size.each.fold[3,NUMFOLD]),
                       rep(4,size.each.class[4]-matrix.size.each.fold[4,NUMFOLD]),
                       rep(5,size.each.class[5]-matrix.size.each.fold[5,NUMFOLD]),
                       rep(6,size.each.class[6]-matrix.size.each.fold[6,NUMFOLD]),
                       rep(7,size.each.class[7]-matrix.size.each.fold[7,NUMFOLD]),
                       rep(8,size.each.class[8]-matrix.size.each.fold[8,NUMFOLD]),
                       rep(9,size.each.class[9]-matrix.size.each.fold[9,NUMFOLD]),
                       rep(10,size.each.class[10]-matrix.size.each.fold[10,NUMFOLD]))
  
  training_labels <- as.character(training_labels)
  
  test_data <- cbind(NSK.samples[,as.numeric(unlist(random.classes.folds[[1]][NUMFOLD]))],
                     NEV.samples[,as.numeric(unlist(random.classes.folds[[2]][NUMFOLD]))],
                     BCC.samples[,as.numeric(unlist(random.classes.folds[[3]][NUMFOLD]))],
                     ISCC.samples[,as.numeric(unlist(random.classes.folds[[4]][NUMFOLD]))],
                     PMCC.samples[,as.numeric(unlist(random.classes.folds[[5]][NUMFOLD]))],
                     MMCC.samples[,as.numeric(unlist(random.classes.folds[[6]][NUMFOLD]))],
                     PRIMEL.samples[,as.numeric(unlist(random.classes.folds[[7]][NUMFOLD]))],
                     METMEL.samples[,as.numeric(unlist(random.classes.folds[[8]][NUMFOLD]))],
                     AK.samples[,as.numeric(unlist(random.classes.folds[[9]][NUMFOLD]))],
                     PS.samples[,as.numeric(unlist(random.classes.folds[[10]][NUMFOLD]))])
  
  test_labels <- c(rep(1,matrix.size.each.fold[1,NUMFOLD]),
                   rep(2,matrix.size.each.fold[2,NUMFOLD]),
                   rep(3,matrix.size.each.fold[3,NUMFOLD]),
                   rep(4,matrix.size.each.fold[4,NUMFOLD]),
                   rep(5,matrix.size.each.fold[5,NUMFOLD]),
                   rep(6,matrix.size.each.fold[6,NUMFOLD]),
                   rep(7,matrix.size.each.fold[7,NUMFOLD]),
                   rep(8,matrix.size.each.fold[8,NUMFOLD]),
                   rep(9,matrix.size.each.fold[9,NUMFOLD]),
                   rep(10,matrix.size.each.fold[10,NUMFOLD]))
  
  test_labels <- as.character(test_labels)
  
  # Differential Expression Analysis by FOLD
  DEA.by.fold <- rbind(training_data,training_labels)
  
  designMulti.by.fold <- model.matrix(~0+DEA.by.fold[dim(DEA.by.fold)[1],])
  colnames(designMulti.by.fold) = c("NSK","PS","NEV","BCC","ISCC","PMCC","MMCC","PRIMEL","METMEL","AK") # OJO CON EL DESIGNMULTI!
  fit.by.fold <- lmFit(training_data, designMulti.by.fold)
  
  cont.matrix.by.fold = makeContrasts(
    NSKvsNEV = NSK-NEV,
    NSKvsBCC = NSK-BCC,
    NSKvsISCC = NSK-ISCC,
    NSKvsPMCC = NSK-PMCC,
    NSKvsMMCC = NSK-MMCC,
    NSKvsPRIMEL = NSK-PRIMEL,
    NSKvsMETMEL = NSK-METMEL,
    NSKvsAK = NSK-AK,
    NSKvsPS = NSK-PS,
    NEVvsBCC = NEV-BCC,
    NEVvsISCC = NEV-ISCC,
    NEVvsPMCC = NEV-PMCC,
    NEVvsMMCC = NEV-MMCC,
    NEVvsPRIMEL = NEV-PRIMEL,
    NEVvsMETMEL = NEV-METMEL,
    NEVvsAK = NEV-AK,
    NEVvsPS = NEV-PS,
    BCCvsISCC = BCC-ISCC,
    BCCvsPMCC = BCC-PMCC,
    BCCvsMMCC = BCC-MMCC,
    BCCvsPRIMEL = BCC-PRIMEL,
    BCCvsMETMEL = BCC-METMEL,
    BCCvsAK = BCC-AK,
    BCCvsPS = BCC-PS,
    ISCCvsPMCC = ISCC-PMCC,
    ISCCvsMMCC = ISCC-MMCC,
    ISCCvsPRIMEL = ISCC-PRIMEL,
    ISCCvsMETMEL = ISCC-METMEL,
    ISCCvsAK = ISCC-AK,
    ISCCvsPS = ISCC-PS,
    PMCCvsMMCC = PMCC-MMCC,
    PMCCvsPRIMEL = PMCC-PRIMEL,
    PMCCvsMETMEL = PMCC-METMEL,
    PMCCvsAK = PMCC-AK,
    PMCCvsPS = PMCC-PS,
    MMCCvsPRIMEL = MMCC-PRIMEL,
    MMCCvsMETMEL = MMCC-METMEL,
    MMCCvsAK = MMCC-AK,
    MMCCvsPS = MMCC-PS,
    PRIMELvsMETMEL = PRIMEL-METMEL,
    PRIMELvsAK = PRIMEL-AK,
    PRIMELvsPS = PRIMEL-PS,
    METMELvsAK = METMEL-AK,
    METMELvsPS = METMEL-PS,
    AKvsPS = AK-PS,
    levels = designMulti.by.fold)
  
  fitmulti.by.fold = contrasts.fit(fit.by.fold,cont.matrix.by.fold)
  fitmulti.by.fold <- eBayes(fitmulti.by.fold)
  
  for(LFCMIN in seq(1,5,0.5)){     ### Change this loop to LFCMIN = 1 as it was selected for our final analysis
    res.val <- decideTests(fitmulti.by.fold,p.value=0.001,lfc = LFCMIN)
    ind.val <- which(apply(res.val,1,function(x) {length(which(x != 0))>0}) == T)
    lfcIndmatrix.sig <- fitmulti.by.fold$coefficients[ind.val,]
    lfcIndmatrix.abs <- abs(lfcIndmatrix.sig)
    
    lfcs.sig <- as.data.frame(lfcIndmatrix.sig)
    lfcs.sig <- cbind(rownames(lfcs.sig),lfcs.sig)
    
    for(numMaxLFC in seq(1,10,1)){ ### Change this loop to numMaxLFC = 1 as it was selected for our final analysis
      
      if(dim(lfcIndmatrix.abs)[1] >= numMaxLFC){
        genesSeveralMaxLFC <- rownames(lfcIndmatrix.abs)[apply(lfcIndmatrix.abs,2, order, decreasing = T)]
        genesSeveralMaxLFC <- data.frame(matrix(unlist(genesSeveralMaxLFC), nrow = dim(lfcIndmatrix.abs)[1], byrow = F), stringsAsFactors = F)
        genesSeveralMaxLFC <- genesSeveralMaxLFC[1:numMaxLFC,]
        
        genesFilteredByLFC <- data.frame(matrix(character(),dim(genesSeveralMaxLFC)[1],dim(genesSeveralMaxLFC)[2]),stringsAsFactors = F)
        for (i in 1:dim(genesSeveralMaxLFC)[1]){
          for (j in 1:dim(genesSeveralMaxLFC)[2]){
            if(lfcIndmatrix.abs[genesSeveralMaxLFC[i,j],j] >= LFCMIN){
              genesFilteredByLFC[i,j] <- as.character(genesSeveralMaxLFC[i,j])
            }
          }
        }
        genesFilteredByLFC <- as.data.frame(genesFilteredByLFC)
        colnames(genesFilteredByLFC) <- colnames(lfcIndmatrix.abs)
        genesSeveralMaxLFC <- unique(unlist(genesFilteredByLFC))
        genesSeveralMaxLFC <- genesSeveralMaxLFC[!is.na(genesSeveralMaxLFC)]
        TOTALDEGS <- length(genesSeveralMaxLFC)
        filename <- paste("GENES_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,"_MDEGS_",TOTALDEGS,".RData",sep = "")
        save(genesSeveralMaxLFC, file = filename)
        
        # TRAIN
        DEGsMultiClass_train <- as.data.frame(training_data[genesSeveralMaxLFC,])
        if(length(genesSeveralMaxLFC)==1){
          colnames(DEGsMultiClass_train) <- genesSeveralMaxLFC
        }else{
          DEGsMultiClass_train <- t(DEGsMultiClass_train)
        }
        DEGsMultiClass_train <- cbind(DEGsMultiClass_train,training_labels)
        finalDataSelected <- as.matrix(DEGsMultiClass_train)
        matrix.train <- matrix(as.numeric(finalDataSelected), nrow = dim(finalDataSelected)[1], ncol = dim(finalDataSelected)[2])
        
        filename <- paste("DATASET_TRAIN_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
        save(matrix.train, file = filename)
        
        # TEST
        DEGSMulticlass_test <- as.data.frame(test_data[genesSeveralMaxLFC,])
        if(length(genesSeveralMaxLFC)==1){
          colnames(DEGSMulticlass_test) <- genesSeveralMaxLFC
        }else{
          DEGSMulticlass_test <- t(DEGSMulticlass_test)
        }
        DEGSMulticlass_test <- cbind(DEGSMulticlass_test,test_labels)
        finalDataSelected <- as.matrix(DEGSMulticlass_test)
        matrix.test <- matrix(as.numeric(finalDataSelected), nrow = dim(finalDataSelected)[1], ncol = dim(finalDataSelected)[2])
        
        filename <- paste("DATASET_TEST_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
        save(matrix.test, file = filename)
        
        # LOG
        print(paste(NUMFOLD,";",LFCMIN,";",numMaxLFC,";",TOTALDEGS,sep = ""))
        
      } else { 
        next   
      }
    }
  }
}

# ***************************************************************************************************************************
# Identification of common genes matching all the M experiments for each parameter combination (LFC and NMAX)
# ***************************************************************************************************************************
for(numMaxLFC in seq(1,10,1)){ # Change this loop to numMaxLFC = 1 as it was selected for our final analysis
  for(LFCMIN in seq(1,5,0.5)){ # Change this loop to LFCMIN = 1 as it was selected for our final analysis
    pattern1 <- paste("GENES_KFOLD_1_LFC_",LFCMIN,"_NMAX_",numMaxLFC,"_",sep = "")
    filename <- list.files(path = ".", pattern = pattern1)
    load(filename)
    commonGenes <- genesSeveralMaxLFC
    for(NUMFOLD in seq(2,10,1)){
      pattern <- paste("GENES_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,"_",sep = "")
      filename <- list.files(path = ".", pattern = pattern)
      load(filename)
      commonGenes <- intersect(commonGenes,genesSeveralMaxLFC)
    }
    final_filename <- paste("COMMON_GENES_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
    save(commonGenes, file = final_filename)
    rm(final_filename,filename)
  }
}

# ***************************************************************************************************************************
# Selection of common genes matching all the M experiments for each parameter combination (LFC and NMAX)
# ***************************************************************************************************************************
for(NUMFOLD in seq(1,10,1)){ 
  training_data <- cbind(NSK.samples[,-as.numeric(unlist(random.classes.folds[[1]][NUMFOLD]))],
                         NEV.samples[,-as.numeric(unlist(random.classes.folds[[2]][NUMFOLD]))],
                         BCC.samples[,-as.numeric(unlist(random.classes.folds[[3]][NUMFOLD]))],
                         ISCC.samples[,-as.numeric(unlist(random.classes.folds[[4]][NUMFOLD]))],
                         PMCC.samples[,-as.numeric(unlist(random.classes.folds[[5]][NUMFOLD]))],
                         MMCC.samples[,-as.numeric(unlist(random.classes.folds[[6]][NUMFOLD]))],
                         PRIMEL.samples[,-as.numeric(unlist(random.classes.folds[[7]][NUMFOLD]))],
                         METMEL.samples[,-as.numeric(unlist(random.classes.folds[[8]][NUMFOLD]))],
                         AK.samples[,-as.numeric(unlist(random.classes.folds[[9]][NUMFOLD]))],
                         PS.samples[,-as.numeric(unlist(random.classes.folds[[10]][NUMFOLD]))])
  
  training_labels <- c(rep(1,size.each.class[1]-matrix.size.each.fold[1,NUMFOLD]),
                       rep(2,size.each.class[2]-matrix.size.each.fold[2,NUMFOLD]),
                       rep(3,size.each.class[3]-matrix.size.each.fold[3,NUMFOLD]),
                       rep(4,size.each.class[4]-matrix.size.each.fold[4,NUMFOLD]),
                       rep(5,size.each.class[5]-matrix.size.each.fold[5,NUMFOLD]),
                       rep(6,size.each.class[6]-matrix.size.each.fold[6,NUMFOLD]),
                       rep(7,size.each.class[7]-matrix.size.each.fold[7,NUMFOLD]),
                       rep(8,size.each.class[8]-matrix.size.each.fold[8,NUMFOLD]),
                       rep(9,size.each.class[9]-matrix.size.each.fold[9,NUMFOLD]),
                       rep(10,size.each.class[10]-matrix.size.each.fold[10,NUMFOLD]))
  
  training_labels <- as.character(training_labels)
  
  test_data <- cbind(NSK.samples[,as.numeric(unlist(random.classes.folds[[1]][NUMFOLD]))],
                     NEV.samples[,as.numeric(unlist(random.classes.folds[[2]][NUMFOLD]))],
                     BCC.samples[,as.numeric(unlist(random.classes.folds[[3]][NUMFOLD]))],
                     ISCC.samples[,as.numeric(unlist(random.classes.folds[[4]][NUMFOLD]))],
                     PMCC.samples[,as.numeric(unlist(random.classes.folds[[5]][NUMFOLD]))],
                     MMCC.samples[,as.numeric(unlist(random.classes.folds[[6]][NUMFOLD]))],
                     PRIMEL.samples[,as.numeric(unlist(random.classes.folds[[7]][NUMFOLD]))],
                     METMEL.samples[,as.numeric(unlist(random.classes.folds[[8]][NUMFOLD]))],
                     AK.samples[,as.numeric(unlist(random.classes.folds[[9]][NUMFOLD]))],
                     PS.samples[,as.numeric(unlist(random.classes.folds[[10]][NUMFOLD]))])
  
  test_labels <- c(rep(1,matrix.size.each.fold[1,NUMFOLD]),
                   rep(2,matrix.size.each.fold[2,NUMFOLD]),
                   rep(3,matrix.size.each.fold[3,NUMFOLD]),
                   rep(4,matrix.size.each.fold[4,NUMFOLD]),
                   rep(5,matrix.size.each.fold[5,NUMFOLD]),
                   rep(6,matrix.size.each.fold[6,NUMFOLD]),
                   rep(7,matrix.size.each.fold[7,NUMFOLD]),
                   rep(8,matrix.size.each.fold[8,NUMFOLD]),
                   rep(9,matrix.size.each.fold[9,NUMFOLD]),
                   rep(10,matrix.size.each.fold[10,NUMFOLD]))
  
  test_labels <- as.character(test_labels)
  
  for(LFCMIN in seq(1,5,0.5)){     # Change this loop to LFCMIN = 1 as it was selected for our final analysis
    for(numMaxLFC in seq(1,10,1)){ # Change this loop to numMaxLFC = 1 as it was selected for our final analysis#1,10,1
      filename <- paste("COMMON_GENES_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
      load(filename)
      # TRAIN
      DEGsMultiClass_train <- as.data.frame(training_data[commonGenes,])
      if(length(commonGenes)==1){
        colnames(DEGsMultiClass_train) <- commonGenes
      }else{
        DEGsMultiClass_train <- t(DEGsMultiClass_train)
      }
      DEGsMultiClass_train <- as.matrix(DEGsMultiClass_train)
      filename <- paste("DATASET_TRAIN_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
      save(DEGsMultiClass_train, training_labels, file = filename)
      
      # TEST
      DEGSMulticlass_test <- as.data.frame(test_data[commonGenes,])
      if(length(commonGenes)==1){
        colnames(DEGSMulticlass_test) <- commonGenes
      }else{
        DEGSMulticlass_test <- t(DEGSMulticlass_test)
      }
      DEGSMulticlass_test <- as.matrix(DEGSMulticlass_test)
      filename <- paste("DATASET_TEST_KFOLD_",NUMFOLD,"_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
      save(DEGSMulticlass_test, test_labels, file = filename)

      # LOG
      print(paste(NUMFOLD,";",LFCMIN,";",numMaxLFC,";",length(commonGenes),sep = "")) 
    }
  }
  
}

# ***************************************************************************************************************************
# mRMR ranking of each parameter combination has to be calculated from the whole integrated dataset
# ***************************************************************************************************************************
for(LFCMIN in seq(1,5,0.5)){     # Change this loop to LFCMIN = 1 as it was selected for our final analysis
  for(numMaxLFC in seq(1,10,1)){ # Change this loop to numMaxLFC = 1 as it was selected for our final analysis
    filename <- paste("COMMON_GENES_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
    load(filename)
    
    DEGsMultiClass_whole <- as.data.frame(data_normalized[commonGenes,])
    if(length(commonGenes)==1){
      colnames(DEGsMultiClass_whole) <- commonGenes
    }else{
      DEGsMultiClass_whole <- t(DEGsMultiClass_whole)
    }
    DEGsMultiClass_whole <- as.matrix(DEGsMultiClass_whole)
    filename <- paste("DATASET_WHOLE_LFC_",LFCMIN,"_NMAX_",numMaxLFC,".RData",sep = "")
    save(DEGsMultiClass_whole, labels.dataset.series, file = filename)

    # LOG
    print(paste(LFCMIN,";",numMaxLFC,";",length(commonGenes),sep = "")) 
  }
}