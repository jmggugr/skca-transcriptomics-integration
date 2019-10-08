# Loading data (previously obtained from RawDataAcquisition_and_Preprocessing R script)
load("EDATA_RAW_27_SERIES_28_BATCHES_.RData")
load("LABELS_27_SERIES_28_BATCHES.RData")

#######################################################################################################################################
# Single value summarization for each gene transcribing the same identifier
# ¡! This loop is not automated, each iteration has to be separately executed
files <- load("EDATA_RAW_27_SERIES_28_BATCHES.RData")
for(i in 1:length(files)){
  edata <- get(files[i])
  colnames.edata <- colnames(edata)
  edata <- cbind(rownames(edata),edata)
  edata <- as.data.frame(edata)
  edata <- edata[order(edata$V1),]
  unique.genes.edata <- unique(edata$V1)                      
  length(unique.genes.edata)
  summarized.edata <- data.frame()
  for(j in 1:length(unique.genes.edata)){
    print(j)
    subset.edata <- edata[which(edata$V1 == unique.genes.edata[j]),]
    subset.edata <- as.matrix(subset.edata[,-1])
    dims.subset <- dim(subset.edata)
    subset.edata <- as.double(subset.edata)
    dim(subset.edata) <- dims.subset
    subset.edata <- colMeans(subset.edata)                    # Summary by means (it can be replaced by colMedians)
    summarized.edata <- rbind(summarized.edata,subset.edata)
  }
  colnames(summarized.edata) <- colnames.edata
  rownames(summarized.edata) <- unique.genes.edata
  
  # Here, each series has to be separately saved:
  # NAME.EDATA <- summarized.edata (Example: edata02503.trans <- summarized.edata)
  # edata2503.rna.trans <- summarized.edata
}

save(edata02503.trans,edata03189.trans,edata06710.trans,edata07553.trans,edata13355.trans,edata14905.trans,edata15605.trans,
     edata30999.trans,edata32407.trans,edata32628.trans,edata32924.trans,edata36150.trans,edata39612.trans,edata42109.trans,
     edata42677a2.trans,edata42677plus2.trans,edata45216.trans,edata46517.trans,edata50451.trans,edata52471.trans,
     edata53223.trans,edata53462.trans,edata82105.trans,edata5678.ae.trans,edata54456.rna.trans,edata67785.rna.trans,
     edata84293.rna.trans,edata98394.rna.trans, file = "EDATA_TRANS_27_SERIES_28_BATCHES.RData")

# ***************************************************************************************************************************
# Logarithmic transformation (if it is needed for specific series)
# ***************************************************************************************************************************
plot.new()
boxplot(log2(edata02503.trans), at=1:12,  xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=F,las=3, col="blue", yaxt="n")  # 18.44 -> 20 (HGU133A)
boxplot(log2(edata03189.trans), at=13:34, xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="red", yaxt="n")   # 16.75 -> 20 (HGU133A)
boxplot(edata06710.trans, at=35:47, xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="orange", yaxt="n")      # 15.46 -> 20 (HGU133A)
boxplot(edata07553.trans, at=48:91, xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="pink", yaxt="n")        # 14.56 -> 16 (HGU133PLUS2)
boxplot(edata13355.trans, at=92:208,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="brown", yaxt="n")       # 14.74 -> 16 (HGU133PLUS2)
boxplot(edata14905.trans, at=209:255,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="purple", yaxt="n")     # 14.40 -> 16 (HGU133PLUS2)
boxplot(edata15605.trans, at=256:300,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="yellow", yaxt="n")     # 14.23 -> 16 (HGU133PLUS2)
boxplot(edata30999.trans, at=301:373,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="grey", yaxt="n")       # 14.61 -> 16 (HGU133PLUS2)
boxplot(edata32407.trans, at=374:383,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="darkgreen", yaxt="n")  # 14.42 -> 16 (HGU133A2)
boxplot(edata32628.trans, at=384:410,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="white", yaxt="n")      # 15.76 -> 16 (HUMANALL)
boxplot(edata32924.trans, at=411:417,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="darkblue", yaxt="n")   # 14.85 -> 16 (HGU133PLUS2)
boxplot(edata36150.trans, at=418:427,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="darkred", yaxt="n")    # 11.63 -> 12 (HUEX1.0ST)
boxplot(edata39612.trans, at=428:448,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="darkgrey", yaxt="n")   # 14.28 -> 16 (HGU133PLUS2)
boxplot(edata42109.trans, at=449:458,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="khaki", yaxt="n")      # 14.50 -> 16 (HGU133A2)
boxplot(edata42677a2.trans, at=459:468,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="bisque", yaxt="n")   # 14.82 -> 16 (HGU133A2)
boxplot(edata42677plus2.trans, at=469:477,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="sienna", yaxt="n")# 14.38 -> 16 (HGU133PLUS2)
boxplot(edata45216.trans, at=478:512,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="orchid", yaxt="n")     # 14.20 -> 16 (HGU133PLUS2)
boxplot(edata46517.trans, at=513:580,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="royalblue2", yaxt="n") # 14.10 -> 20 (HGU133A)
boxplot(edata50451.trans, at=581:602,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="coral", yaxt="n")      # 14.48 -> 16 (HGU133PLUS2)
boxplot(edata52471.trans, at=603:626,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="chocolate", yaxt="n")  # 14.47 -> 16 (HGU133A2)
boxplot(edata53223.trans, at=627:639,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="chartreuse", yaxt="n") # 15.02 -> 16 (HGU133PLUS2)
boxplot(edata53462.trans, at=640:660,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="grey", yaxt="n")       # 14.19 -> 16 (HUMANALL)
boxplot(edata82105.trans, at=661:666,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="hotpink", yaxt="n")    # 14.27 -> 16 (HGU133PLUS2)
boxplot(edata5678.ae.trans, at=667:688,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="turquoise", yaxt="n")   # 22.73 -> 24 (HISEQ2000)
boxplot(edata54456.rna.trans, at=689:857,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="springgreen", yaxt="n") # 21.75 -> 22 (GA)
boxplot(edata67785.rna.trans, at=858:871,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="gold", yaxt="n")        # 19.56 -> 20 (GAIIX)
boxplot(edata84293.rna.trans, at=872:890,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="cyan", yaxt="n")        # 21.57 -> 24 (HISEQ2000)
boxplot(edata98394.rna.trans, at=891:968,xlim=c(0.5,968.5),ylim=c(1,24),log="y",add=T,las=3, col="firebrick", yaxt="n")   # 20.44 -> 22 (HISEQ2500)

axis(side = 2,at = c(1,2,4,6,8,10,12,14,16,18,20,22,24))

# ***************************************************************************************************************************
# 16-bit depth homogeneization
# ***************************************************************************************************************************
plot.new()
boxplot(16/20*log2(edata02503.trans), at=1:12,  xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=F,las=3, col="blue", yaxt="n")  # 18.44 -> 20 (HGU133A)
boxplot(16/20*log2(edata03189.trans), at=13:34, xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="red", yaxt="n")   # 16.75 -> 20 (HGU133A)
boxplot(16/20*edata06710.trans, at=35:47, xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="orange", yaxt="n")      # 15.46 -> 20 (HGU133A)
boxplot(edata07553.trans, at=48:91, xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="pink", yaxt="n")        # 14.56 -> 16 (HGU133PLUS2)
boxplot(edata13355.trans, at=92:208,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="brown", yaxt="n")       # 14.74 -> 16 (HGU133PLUS2)
boxplot(edata14905.trans, at=209:255,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="purple", yaxt="n")     # 14.40 -> 16 (HGU133PLUS2)
boxplot(edata15605.trans, at=256:300,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="yellow", yaxt="n")     # 14.23 -> 16 (HGU133PLUS2)
boxplot(edata30999.trans, at=301:373,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="grey", yaxt="n")       # 14.61 -> 16 (HGU133PLUS2)
boxplot(edata32407.trans, at=374:383,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="darkgreen", yaxt="n")  # 14.42 -> 16 (HGU133A2)
boxplot(edata32628.trans, at=384:410,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="white", yaxt="n")      # 15.76 -> 16 (HUMANALL)
boxplot(edata32924.trans, at=411:417,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="darkblue", yaxt="n")   # 14.85 -> 16 (HGU133PLUS2)
boxplot(16/12*edata36150.trans, at=418:427,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="darkred", yaxt="n")    # 11.63 -> 12 (HUEX1.0ST)
boxplot(edata39612.trans, at=428:448,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="darkgrey", yaxt="n")   # 14.28 -> 16 (HGU133PLUS2)
boxplot(edata42109.trans, at=449:458,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="khaki", yaxt="n")      # 14.50 -> 16 (HGU133A2)
boxplot(edata42677a2.trans, at=459:468,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="bisque", yaxt="n")   # 14.82 -> 16 (HGU133A2)
boxplot(edata42677plus2.trans, at=469:477,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="sienna", yaxt="n")# 14.38 -> 16 (HGU133PLUS2)
boxplot(edata45216.trans, at=478:512,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="orchid", yaxt="n")     # 14.20 -> 16 (HGU133PLUS2)
boxplot(16/20*edata46517.trans, at=513:580,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="royalblue2", yaxt="n") # 14.10 -> 20 (HGU133A)
boxplot(edata50451.trans, at=581:602,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="coral", yaxt="n")      # 14.48 -> 16 (HGU133PLUS2)
boxplot(edata52471.trans, at=603:626,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="chocolate", yaxt="n")  # 14.47 -> 16 (HGU133A2)
boxplot(edata53223.trans, at=627:639,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="chartreuse", yaxt="n") # 15.02 -> 16 (HGU133PLUS2)
boxplot(edata53462.trans, at=640:660,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="grey", yaxt="n")       # 14.19 -> 16 (HUMANALL)
boxplot(edata82105.trans, at=661:666,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="hotpink", yaxt="n")    # 14.27 -> 16 (HGU133PLUS2)
boxplot(16/24*edata5678.ae.trans, at=667:688,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="turquoise", yaxt="n")   # 22.73 -> 24 (HISEQ2000)
boxplot(16/22*edata54456.rna.trans, at=689:857,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="springgreen", yaxt="n") # 21.75 -> 22 (GA)
boxplot(16/20*edata67785.rna.trans, at=858:871,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="gold", yaxt="n")        # 19.56 -> 20 (GAIIX)
boxplot(16/24*edata84293.rna.trans, at=872:890,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="cyan", yaxt="n")        # 21.57 -> 24 (HISEQ2000)
boxplot(16/22*edata98394.rna.trans, at=891:968,xlim=c(0.5,968.5),ylim=c(1,16),log="y",add=T,las=3, col="firebrick", yaxt="n")   # 20.44 -> 22 (HISEQ2500)

axis(side = 2,at = c(1,2,4,6,8,10,12,14,16))

edata02503.bdc <- 16/20*log2(edata02503.trans)
edata03189.bdc <- 16/20*log2(edata03189.trans)
edata06710.bdc <- 16/20*edata06710.trans
edata07553.bdc <- edata07553.trans
edata13355.bdc <- edata13355.trans
edata14905.bdc <- edata14905.trans
edata15605.bdc <- edata15605.trans
edata30999.bdc <- edata30999.trans
edata32407.bdc <- edata32407.trans
edata32628.bdc <- edata32628.trans
edata32924.bdc <- edata32924.trans
edata36150.bdc <- 16/12*edata36150.trans
edata39612.bdc <- edata39612.trans
edata42109.bdc <- edata42109.trans
edata42677a2.bdc <- edata42677a2.trans
edata42677plus2.bdc <- edata42677plus2.trans
edata45216.bdc <- edata45216.trans
edata46517.bdc <- 16/20*edata46517.trans
edata50451.bdc <- edata50451.trans
edata52471.bdc <- edata52471.trans
edata53223.bdc <- edata53223.trans
edata53462.bdc <- edata53462.trans
edata82105.bdc <- edata82105.trans
edata5678.ae.bdc <- 16/24*edata5678.ae.trans
edata54456.bdc <- 16/22*edata54456.rna.trans
edata67785.bdc <- 16/20*edata67785.rna.trans
edata84293.bdc <- 16/24*edata84293.rna.trans
edata98394.bdc <- 16/22*edata98394.rna.trans

save(edata02503.bdc,edata03189.bdc,edata06710.bdc,edata07553.bdc,edata13355.bdc,edata14905.bdc,edata15605.bdc,
     edata30999.bdc,edata32407.bdc,edata32628.bdc,edata32924.bdc,edata36150.bdc,edata39612.bdc,edata42109.bdc,
     edata42677a2.bdc,edata42677plus2.bdc,edata45216.bdc,edata46517.bdc,edata50451.bdc,edata52471.bdc,
     edata53223.bdc,edata53462.bdc,edata82105.bdc,edata5678.ae.bdc,edata54456.bdc,
     edata67785.bdc,edata84293.bdc,edata98394.bdc, file = "EDATA_BDC_27_SERIES_28_BATCHES.RData")

# ***************************************************************************************************************************
# Complete cases selection
# ***************************************************************************************************************************
files <- load("EDATA_BDC_27_SERIES_28_BATCHES.RData")

unique.genes.all.series <- unique(c(rownames(edata02503.bdc),rownames(edata03189.bdc),rownames(edata06710.bdc),rownames(edata07553.bdc),
                                    rownames(edata13355.bdc),rownames(edata14905.bdc),rownames(edata15605.bdc),rownames(edata30999.bdc),
                                    rownames(edata32407.bdc),rownames(edata32628.bdc),rownames(edata32924.bdc),rownames(edata36150.bdc),
                                    rownames(edata39612.bdc),rownames(edata42109.bdc),rownames(edata42677a2.bdc),rownames(edata42677plus2.bdc),
                                    rownames(edata45216.bdc),rownames(edata46517.bdc),rownames(edata50451.bdc),rownames(edata52471.bdc),
                                    rownames(edata53223.bdc),rownames(edata53462.bdc),rownames(edata82105.bdc),rownames(edata5678.ae.bdc),
                                    rownames(edata54456.bdc),rownames(edata67785.bdc),
                                    rownames(edata84293.bdc),rownames(edata98394.bdc)))

samples.full.dataset <- c(colnames(edata02503.bdc),colnames(edata03189.bdc),colnames(edata06710.bdc),colnames(edata07553.bdc),
                          colnames(edata13355.bdc),colnames(edata14905.bdc),colnames(edata15605.bdc),colnames(edata30999.bdc),
                          colnames(edata32407.bdc),colnames(edata32628.bdc),colnames(edata32924.bdc),colnames(edata36150.bdc),
                          colnames(edata39612.bdc),colnames(edata42109.bdc),colnames(edata42677a2.bdc),colnames(edata42677plus2.bdc),
                          colnames(edata45216.bdc),colnames(edata46517.bdc),colnames(edata50451.bdc),colnames(edata52471.bdc),
                          colnames(edata53223.bdc),colnames(edata53462.bdc),colnames(edata82105.bdc),colnames(edata5678.ae.bdc),
                          colnames(edata54456.bdc),colnames(edata67785.bdc),
                          colnames(edata84293.bdc),colnames(edata98394.bdc))

samples.full.dataset <- sub("\\.CEL$","",samples.full.dataset)

full.dataset.series.bdc.final <- 0        
for(i in 1:length(files)){
  edata <- get(files[i])
  dim.edata <- dim(edata)
  edata.genes <- rownames(edata)
  different.genes <- setdiff(unique.genes.all.series,edata.genes)
  aggregate.NAs <- matrix(NA, nrow = length(different.genes), ncol = dim.edata[2])
  aggregate.NAs <- as.data.frame(aggregate.NAs)
  colnames(aggregate.NAs) <- colnames(edata)
  complete.edata <- rbind(edata,aggregate.NAs)
  lim.inf <- dim(edata)[1]+1
  lim.sup <- dim(complete.edata)[1]
  rownames(complete.edata)[lim.inf:lim.sup] <- different.genes
  complete.edata <- cbind(rownames(complete.edata),complete.edata)
  colnames(complete.edata)[1] <- "Gene"
  complete.edata <- complete.edata[order(complete.edata$Gene),]
  complete.edata <- complete.edata[,-1]
  full.dataset.series.bdc.final <- cbind(full.dataset.series.bdc.final,complete.edata)
}
full.dataset.series.bdc.final  <- full.dataset.series.bdc.final[,-1]

save(full.dataset.series.bdc.final, file = "FULL_DATASET_BDC_WITH_NAS_27_SERIES_28_BATCHES.RData")

data_selected <- full.dataset.series.bdc.final[complete.cases(full.dataset.series.bdc.final),]

# ***************************************************************************************************************************
# Batch effect correction 
# ***************************************************************************************************************************
require(sva)
targets <- c(GSE02503Labels,GSE03189Labels,GSE06710Labels,GSE07553Labels,GSE13355Labels,GSE14905Labels,GSE15605Labels,
             GSE30999Labels,GSE32407Labels,GSE32628Labels,GSE32924Labels,GSE36150Labels,GSE39612Labels,GSE42109Labels,
             GSE42677a2Labels,GSE42677plus2Labels,GSE45216Labels,GSE46517Labels,GSE50451Labels,GSE52471Labels,GSE53223Labels,
             GSE53462Labels,GSE82105Labels,AE5678Labels,GSE54456Labels,GSE67785Labels,GSE84293Labels,GSE98394Labels)

datasets <- c(rep("GSE02503",12), rep("GSE03189",22), rep("GSE06710",13), rep("GSE07553",44), rep("GSE13355",117), rep("GSE14905",47),
              rep("GSE15605",45), rep("GSE30999",73), rep("GSE32407",10), rep("GSE32628",27), rep("GSE32924",7),   rep("GSE36150",10),
              rep("GSE39612",21), rep("GSE42109",10), rep("GSE42677a2",10), rep("GSE42677plus2",9), rep("GSE45216",35), rep("GSE46517",68),
              rep("GSE50541",22), rep("GSE52471",24), rep("GSE53223",13), rep("GSE53462",21), rep("GSE82105",6),   rep("AE5678",22),
              rep("GSE54456",169), rep("GSE67785",14), rep("GSE84293",19),  rep("GSE98394",78))

batch <- as.vector(datasets)
datasets_No <- max(as.numeric(factor(datasets)))
f.model <- factor(targets, levels = levels(factor(targets)))
mod <- model.matrix(~f.model)
data <- data_selected

data_combat <- ComBat(dat = as.matrix(data), batch = batch, mod = mod)
boxplot(data_combat)
save(data_combat, file = "FULL_DATASET_COMBAT_27_SERIES_28_BATCHES.RData")

# ***************************************************************************************************************************
# Inter-array normalization
# ***************************************************************************************************************************
require(limma)
data_normalized <- normalizeBetweenArrays(data_combat)
boxplot(data_normalized)
max(data_normalized)
min(data_normalized)

#############################################################################################################################
#############################################################################################################################
# Saving integrated dataset (matrix and labels)
#############################################################################################################################
#############################################################################################################################
labels.dataset.series <- gsub("NSK",1,targets)
labels.dataset.series <- gsub("NEV",2,labels.dataset.series)
labels.dataset.series <- gsub("BCC",3,labels.dataset.series)
labels.dataset.series <- gsub("ISCC",4,labels.dataset.series)
labels.dataset.series <- gsub("SCC",4,labels.dataset.series)
labels.dataset.series <- gsub("PMCC",5,labels.dataset.series)
labels.dataset.series <- gsub("MMCC",6,labels.dataset.series)
labels.dataset.series <- gsub("PRIMEL",7,labels.dataset.series)
labels.dataset.series <- gsub("METMEL",8,labels.dataset.series)
labels.dataset.series <- gsub("AK",9,labels.dataset.series)
labels.dataset.series <- gsub("PS",10,labels.dataset.series)

save(data_normalized,labels.dataset.series, file = "FULL_DATASET_NORMALIZED_FINAL_27_SERIES_28_BATCHES.RData")