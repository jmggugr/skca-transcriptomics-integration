pwd <- getwd()

source("http://bioconductor.org/biocLite.R")

suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(lumi))
suppressPackageStartupMessages(library(arrayQualityMetrics))
suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(cqn))
suppressPackageStartupMessages(library(NOISeq))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))

#############################################################################################################################
#############################################################################################################################
# Raw data acquisition and preprocessing of microarray series
#############################################################################################################################
#############################################################################################################################

# ***************************************************************************************************************************
# GSE2503 AFFY Study
# ***************************************************************************************************************************
GSE02503 <- getGEO("GSE2503", destdir = "ReferenceFiles/Samples/Microarray/GSE2503/")[[1]]     # No hay .CEL's
outlier <- c("GSM47616")                                                                       # La muestra de UAK
GSE02503 <- GSE02503[,setdiff(sampleNames(GSE02503),outlier)]
arrayQualityMetrics(GSE02503, outdir = "GeneratedFiles/QualityAnalysis/GSE02503", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM47614","GSM47621")                                                       
GSE02503 <- GSE02503[,setdiff(sampleNames(GSE02503),outlier)]
arrayQualityMetrics(GSE02503, outdir = "GeneratedFiles/QualityAnalysis/GSE02503_2", force = TRUE, do.logtransform = TRUE)

GSE02503Norm <- GSE02503
phenoGSE02503 <- pData(GSE02503Norm)
edataGSE02503 <- exprs(GSE02503Norm)
GSE02503Labels <- c(rep("AK",3),rep("NSK",4),rep("ISCC",5))

require("hgu133a.db", character.only = TRUE)
genesIDGSE02503 <- rownames(edataGSE02503)
SymbolsGSE02503 <- getSYMBOL(genesIDGSE02503, "hgu133a.db")
EntrezGSE02503 <- getEG(genesIDGSE02503, "hgu133a.db")

NAs <- which(is.na(SymbolsGSE02503))
SymbolsGSE02503 <- SymbolsGSE02503[-NAs]
edataGSE02503 <- edataGSE02503[-NAs,]

rownames(edataGSE02503) <- SymbolsGSE02503

# SUMMARY: 4 NSK, 5 ISCC, 3 AK

# ***************************************************************************************************************************
# GSE3189 AFFY Study
# ***************************************************************************************************************************
GSE03189 <- getGEO("GSE3189", destdir = "ReferenceFiles/Samples/Microarray/GSE3189/")[[1]]     # No hay .CEL's
selection <- match(c(paste0("GSM716",71:95)),sampleNames(GSE03189)) 
GSE03189 <- GSE03189[,selection]
arrayQualityMetrics(GSE03189, outdir = "GeneratedFiles/QualityAnalysis/GSE03189", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM71675","GSM71694","GSM71695")
GSE03189 <- GSE03189[,setdiff(sampleNames(GSE03189),outlier)]
arrayQualityMetrics(GSE03189, outdir = "GeneratedFiles/QualityAnalysis/GSE03189_2", force = TRUE, do.logtransform = TRUE)

GSE03189Norm <- GSE03189
phenoGSE03189 <- pData(GSE03189Norm)
edataGSE03189 <- exprs(GSE03189Norm)
GSE03189Labels <- c(rep("NSK",6),rep("NEV",16))

require("hgu133a.db", character.only = TRUE)
genesIDGSE03189 <- rownames(edataGSE03189)
SymbolsGSE03189 <- getSYMBOL(genesIDGSE03189, "hgu133a.db")
EntrezGSE03189 <- getEG(genesIDGSE03189, "hgu133a.db")

NAs <- which(is.na(SymbolsGSE03189))
SymbolsGSE03189 <- SymbolsGSE03189[-NAs]
edataGSE03189 <- edataGSE03189[-NAs,]

rownames(edataGSE03189) <- SymbolsGSE03189

# SUMMARY: 6 NSK, 16 NEV

# ***************************************************************************************************************************
# GSE6710 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE6710")
GSE06710 <- read.celfiles(list.celfiles()) 
setwd(pwd)
arrayQualityMetrics(GSE06710, outdir = "GeneratedFiles/QualityAnalysis/GSE06710", force = TRUE, do.logtransform = TRUE)

GSE06710Norm <- rma(GSE06710)
phenoGSE06710 <- pData(GSE06710Norm)
edataGSE06710 <- exprs(GSE06710Norm)
GSE06710Labels <- rep("PS",13)

require("hgu133a.db", character.only = TRUE)
genesIDGSE06710 <- rownames(edataGSE06710)
SymbolsGSE06710 <- getSYMBOL(genesIDGSE06710, "hgu133a.db")
EntrezGSE06710 <- getEG(genesIDGSE06710, "hgu133a.db")

NAs <- which(is.na(SymbolsGSE06710))
SymbolsGSE06710 <- SymbolsGSE06710[-NAs]
edataGSE06710 <- edataGSE06710[-NAs,]

rownames(edataGSE06710) <- SymbolsGSE06710

# SUMMARY: 13 PS

# ***************************************************************************************************************************
# GSE7553 AFFY Study
# ***************************************************************************************************************************
GSE07553 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE7553/") 
arrayQualityMetrics(GSE07553, outdir = "GeneratedFiles/QualityAnalysis/GSE07553", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM183265.CEL","GSM183304.CEL")
GSE07553 <- GSE07553[,setdiff(sampleNames(GSE07553),outlier)]
arrayQualityMetrics(GSE07553, outdir = "GeneratedFiles/QualityAnalysis/GSE07553_2", force = TRUE, do.logtransform = TRUE)

GSE07553Norm <- affy::rma(GSE07553)
phenoGSE07553 <- pData(GSE07553Norm)
edataGSE07553 <- exprs(GSE07553Norm)
GSE07553Labels <- c(rep("BCC",4),rep("PRIMEL",4),"NSK","PRIMEL",rep("SCC",5),rep("BCC",11),rep("PRIMEL",8),rep("SCC",6),rep("NSK",3),"PRIMEL")

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE07553 <- rownames(edataGSE07553)
SymbolsGSE07553 <- getSYMBOL(genesIDGSE07553, "hgu133plus2.db")
EntrezGSE07553 <- getEG(genesIDGSE07553, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE07553))
SymbolsGSE07553 <- SymbolsGSE07553[-NAs]
edataGSE07553 <- edataGSE07553[-NAs,]

rownames(edataGSE07553) <- SymbolsGSE07553

# SUMMARY: 4 NSK, 15 BCC, 11 SCC, 14 PRIMEL

# ***************************************************************************************************************************
# GSE13355 AFFY Study
# ***************************************************************************************************************************
GSE13355 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE13355") 
sampleNames(GSE13355) <- sub("\\.CEL$", "", sampleNames(GSE13355));
outlier <- c(paste0("GSM3372",61:99),paste0("GSM33730",0:9),paste0("GSM33731",0:8)) 
GSE13355 <- GSE13355[,setdiff(sampleNames(GSE13355),outlier)]
arrayQualityMetrics(GSE13355, outdir = "GeneratedFiles/QualityAnalysis/GSE13355", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM337240","GSM337243","GSM337250","GSM337332","GSM337360")
GSE13355 <- GSE13355[,setdiff(sampleNames(GSE13355),outlier)]
arrayQualityMetrics(GSE13355, outdir = "GeneratedFiles/QualityAnalysis/GSE13355_2", force = TRUE, do.logtransform = TRUE)

GSE13355Norm <- affy::rma(GSE13355)
phenoGSE13355 <- pData(GSE13355Norm)
edataGSE13355 <- exprs(GSE13355Norm)

GSE13355Labels <- c(rep("NSK",61),rep("PS",56))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE13355 <- rownames(edataGSE13355)
SymbolsGSE13355 <- getSYMBOL(genesIDGSE13355, "hgu133plus2.db")
EntrezGSE13355 <- getEG(genesIDGSE13355, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE13355))
SymbolsGSE13355 <- SymbolsGSE13355[-NAs]
edataGSE13355 <- edataGSE13355[-NAs,]

rownames(edataGSE13355) <- SymbolsGSE13355

# SUMMARY: 61 NSK, 56 PS

# ***************************************************************************************************************************
# GSE14905 AFFY Study
# ***************************************************************************************************************************
GSE14905 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE14905") 
sampleNames(GSE14905) <- sub("\\.CEL$", "", sampleNames(GSE14905));
arrayQualityMetrics(GSE14905, outdir = "GeneratedFiles/QualityAnalysis/GSE14905", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM372287","GSM372294","GSM372295","GSM372350")
GSE14905 <- GSE14905[,setdiff(sampleNames(GSE14905),outlier)]
arrayQualityMetrics(GSE14905, outdir = "GeneratedFiles/QualityAnalysis/GSE14905_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM372289","GSM372348")
GSE14905 <- GSE14905[,setdiff(sampleNames(GSE14905),outlier)]
arrayQualityMetrics(GSE14905, outdir = "GeneratedFiles/QualityAnalysis/GSE14905_3", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM372286")
GSE14905 <- GSE14905[,setdiff(sampleNames(GSE14905),outlier)]
arrayQualityMetrics(GSE14905, outdir = "GeneratedFiles/QualityAnalysis/GSE14905_4", force = TRUE, do.logtransform = TRUE)

GSE14905Norm <- affy::rma(GSE14905)
phenoGSE14905 <- pData(GSE14905Norm)
edataGSE14905 <- exprs(GSE14905Norm)

GSE14905Labels <- c(rep("NSK",16),rep("PS",31))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE14905 <- rownames(edataGSE14905)
SymbolsGSE14905 <- getSYMBOL(genesIDGSE14905, "hgu133plus2.db")
EntrezGSE14905 <- getEG(genesIDGSE14905, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE14905))
SymbolsGSE14905 <- SymbolsGSE14905[-NAs]
edataGSE14905 <- edataGSE14905[-NAs,]

rownames(edataGSE14905) <- SymbolsGSE14905

# SUMMARY: 16 NSK, 31 PS

# ***************************************************************************************************************************
# GSE15605 AFFY Study
# ***************************************************************************************************************************
GSE15605 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE15605") 
sampleNames(GSE15605) <- sub("\\.CEL$", "", sampleNames(GSE15605));
arrayQualityMetrics(GSE15605, outdir = "GeneratedFiles/QualityAnalysis/GSE15605", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM390210","GSM390222","GSM390223","GSM390224","GSM390225","GSM390229","GSM390232","GSM390235",
             "GSM390238","GSM390240","GSM390246","GSM390252","GSM390254","GSM390257","GSM390259","GSM390264",
             "GSM390265","GSM390267")

GSE15605 <- GSE15605[,setdiff(sampleNames(GSE15605),outlier)]
arrayQualityMetrics(GSE15605, outdir = "GeneratedFiles/QualityAnalysis/GSE15605_2", force = TRUE, do.logtransform = TRUE)

GSE15605Norm <- affy::rma(GSE15605)
phenoGSE15605 <- pData(GSE15605Norm)
edataGSE15605 <- exprs(GSE15605Norm)

GSE15605Labels <- c(rep("NSK",13),rep("PRIMEL",30),rep("METMEL",2))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE15605 <- rownames(edataGSE15605)
SymbolsGSE15605 <- getSYMBOL(genesIDGSE15605, "hgu133plus2.db")
EntrezGSE15605 <- getEG(genesIDGSE15605, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE15605))
SymbolsGSE15605 <- SymbolsGSE15605[-NAs]
edataGSE15605 <- edataGSE15605[-NAs,]

rownames(edataGSE15605) <- SymbolsGSE15605

# SUMMARY: 13 NSK, 30 PRIMEL, 2 METMEL

# ***************************************************************************************************************************
# GSE30999 AFFY Study
# ***************************************************************************************************************************
GSE30999 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE30999") 
sampleNames(GSE30999) <- sub("\\.CEL$", "", sampleNames(GSE30999));
arrayQualityMetrics(GSE30999, outdir = "GeneratedFiles/QualityAnalysis/GSE30999", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM768038","GSM768063","GSM768089","GSM768103","GSM768107","GSM768115","GSM768125")
GSE30999 <- GSE30999[,setdiff(sampleNames(GSE30999),outlier)]
arrayQualityMetrics(GSE30999, outdir = "GeneratedFiles/QualityAnalysis/GSE30999_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM768026","GSM768073","GSM768119")
GSE30999 <- GSE30999[,setdiff(sampleNames(GSE30999),outlier)]
arrayQualityMetrics(GSE30999, outdir = "GeneratedFiles/QualityAnalysis/GSE30999_3", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM768028","GSM768131")
GSE30999 <- GSE30999[,setdiff(sampleNames(GSE30999),outlier)]
arrayQualityMetrics(GSE30999, outdir = "GeneratedFiles/QualityAnalysis/GSE30999_4", force = TRUE, do.logtransform = TRUE)

GSE30999Norm <- affy::rma(GSE30999)
phenoGSE30999 <- pData(GSE30999Norm)
edataGSE30999 <- exprs(GSE30999Norm)

GSE30999Labels <- rep("PS",73)

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE30999 <- rownames(edataGSE30999)
SymbolsGSE30999 <- getSYMBOL(genesIDGSE30999, "hgu133plus2.db")
EntrezGSE30999 <- getEG(genesIDGSE30999, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE30999))
SymbolsGSE30999 <- SymbolsGSE30999[-NAs]
edataGSE30999 <- edataGSE30999[-NAs,]

rownames(edataGSE30999) <- SymbolsGSE30999

# SUMMARY: 73 PS

# ***************************************************************************************************************************
# GSE32407 AFFY Study
# ***************************************************************************************************************************
GSE32407 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE32407")
sampleNames(GSE32407) <- sub("\\.CEL$", "", sampleNames(GSE32407));
arrayQualityMetrics(GSE32407, outdir = "GeneratedFiles/QualityAnalysis/GSE32407", force = TRUE, do.logtransform = TRUE)

GSE32407Norm <- affy::rma(GSE32407)
phenoGSE32407 <- pData(GSE32407Norm)
edataGSE32407 <- exprs(GSE32407Norm)

GSE32407Labels <- rep("NSK",10)

require("hgu133a2.db", character.only = TRUE)
genesIDGSE32407 <- rownames(edataGSE32407)
SymbolsGSE32407 <- getSYMBOL(genesIDGSE32407, "hgu133a2.db")
EntrezGSE32407 <- getEG(genesIDGSE32407, "hgu133a2.db")

NAs <- which(is.na(SymbolsGSE32407))
SymbolsGSE32407 <- SymbolsGSE32407[-NAs]
edataGSE32407 <- edataGSE32407[-NAs,]

rownames(edataGSE32407) <- SymbolsGSE32407

# SUMMARY: 10 NSK

# ***************************************************************************************************************************
# GSE32628 ILLU Study
# ***************************************************************************************************************************
GSE32628 <- lumiR(file = "ReferenceFiles/Samples/Microarray/GSE32628/GSE32628_non_normalized.txt", sep = "\t",
                  lib.mapping = "lumiHumanIDMapping")

# First, identify uninvolved actinic keratosis sample and exclude it
GSE32628 <- GSE32628[,setdiff(sampleNames(GSE32628),outlier)]
arrayQualityMetrics(GSE32628, outdir = "GeneratedFiles/QualityAnalysis/GSE32628", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM808779","GSM808808")
GSE32628 <- GSE32628[,setdiff(sampleNames(GSE32628),outlier)]
arrayQualityMetrics(GSE32628, outdir = "GeneratedFiles/QualityAnalysis/GSE32628_2", force = TRUE, do.logtransform = TRUE)

GSE32628Norm <- lumiExpresso(GSE32628)
phenoGSE32628 <- pData(GSE32628Norm)
edataGSE32628 <- exprs(GSE32628Norm)

# Later, identify and label the final selected samples for this series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32628
GSE32628Labels <- as.vector("SCC","SCC","AK","SCC","AK","SCC","AK","SCC","AK","SCC","AK","SCC","AK","SCC",
                            "AK","SCC","AK","SCC","AK","AK","SCC","AK","SCC","AK","SCC","AK","SCC")

require("lumiHumanAll.db", character.only = TRUE)
genesIDGSE32628 <- rownames(edataGSE32628)
SymbolsGSE32628 <- getSYMBOL(genesIDGSE32628, "lumiHumanAll.db")
EntrezGSE32628 <- getEG(genesIDGSE32628, "lumiHumanAll.db")

NAs <- which(is.na(SymbolsGSE32628))
SymbolsGSE32628 <- SymbolsGSE32628[-NAs]
edataGSE32628 <- edataGSE32628[-NAs,]

rownames(edataGSE32628) <- SymbolsGSE32628

# SUMMARY: 14 ISCC, 13 AK 

# ***************************************************************************************************************************
# GSE32924 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE32924")
file.rename(list.files(pattern = "GSM815.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM815.*.CEL")))
setwd(pwd)
GSE32924 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE32924")
sampleNames(GSE32924) <- sub("\\.CEL$", "", sampleNames(GSE32924));
arrayQualityMetrics(GSE32924, outdir = "GeneratedFiles/QualityAnalysis/GSE32924", force = TRUE, do.logtransform = TRUE)

outlier <- "GSM815458"
GSE32924 <- GSE32924[,setdiff(sampleNames(GSE32924),outlier)]
arrayQualityMetrics(GSE32924, outdir = "GeneratedFiles/QualityAnalysis/GSE32924_2", force = TRUE, do.logtransform = TRUE)

GSE32924Norm <- affy::rma(GSE32924)
phenoGSE32924 <- pData(GSE32924Norm)
edataGSE32924 <- exprs(GSE32924Norm)

GSE32924Labels <- rep("NSK",7)

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE32924 <- rownames(edataGSE32924)
SymbolsGSE32924 <- getSYMBOL(genesIDGSE32924, "hgu133plus2.db")
EntrezGSE32924 <- getEG(genesIDGSE32924, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE32924))
SymbolsGSE32924 <- SymbolsGSE32924[-NAs]
edataGSE32924 <- edataGSE32924[-NAs,]

rownames(edataGSE32924) <- SymbolsGSE32924

# SUMMARY: 7 NSK

# ***************************************************************************************************************************
# GSE36150 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE36150")
file.rename(list.files(pattern = "GSM8818.*.CEL"), paste0("GSM8818", 44:58,".CEL"))
GSE36150 <- read.celfiles(list.files())
sampleNames(GSE36150) <- sub("\\.CEL$", "", sampleNames(GSE36150));
setwd(pwd)
arrayQualityMetrics(GSE36150, outdir = "GeneratedFiles/QualityAnalysis/GSE36150", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM881845","GSM881847","GSM881848","GSM881856","GSM881858")
GSE36150 <- GSE36150[,setdiff(sampleNames(GSE36150),outlier)]

GSE36150Norm <- rma(GSE36150)
phenoGSE36150 <- pData(GSE36150Norm)
edataGSE36150 <- exprs(GSE36150Norm)

GSE36150Labels <- c(rep("PMCC",5),rep("MMCC",5))

require("huex10sttranscriptcluster.db", character.only = TRUE)
genesIDGSE36150 <- rownames(edataGSE36150)
SymbolsGSE36150 <- getSYMBOL(genesIDGSE36150, "huex10sttranscriptcluster.db")
EntrezGSE36150 <- getEG(genesIDGSE36150, "huex10sttranscriptcluster.db")

NAs <- which(is.na(SymbolsGSE36150))
SymbolsGSE36150 <- SymbolsGSE36150[-NAs]
edataGSE36150 <- edataGSE36150[-NAs,]

rownames(edataGSE36150) <- SymbolsGSE36150

# SUMMARY: 5 PMCC, 5 MMCC

# ***************************************************************************************************************************
# GSE39612 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE39612")
file.rename(list.files(pattern = "GSM97.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM97.*.CEL")))
setwd(pwd)
GSE39612 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE39612")
sampleNames(GSE39612) <- sub("\\.CEL$", "", sampleNames(GSE39612));
arrayQualityMetrics(GSE39612, outdir = "GeneratedFiles/QualityAnalysis/GSE39612", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM972996","GSM972997","GSM973010","GSM973016","GSM973019","GSM973026","GSM973029","GSM973030","GSM973036","GSM973039","GSM973042")
GSE39612 <- GSE39612[,setdiff(sampleNames(GSE39612),outlier)]
arrayQualityMetrics(GSE39612, outdir = "GeneratedFiles/QualityAnalysis/GSE39612_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM973003","GSM973008","GSM973022","GSM973040")
GSE39612 <- GSE39612[,setdiff(sampleNames(GSE39612),outlier)]
arrayQualityMetrics(GSE39612, outdir = "GeneratedFiles/QualityAnalysis/GSE39612_3", force = TRUE, do.logtransform = TRUE)

GSE39612Norm <- affy::rma(GSE39612)
phenoGSE39612 <- pData(GSE39612Norm)
edataGSE39612 <- exprs(GSE39612Norm)

GSE39612Labels <- c(rep("ISCC",2),rep("BCC",2),rep("MMCC",2),rep("PMCC",10),"MMCC",rep("PMCC",2),rep("MMCC",2))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE39612 <- rownames(edataGSE39612)
SymbolsGSE39612 <- getSYMBOL(genesIDGSE39612, "hgu133plus2.db")
EntrezGSE39612 <- getEG(genesIDGSE39612, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE39612))
SymbolsGSE39612 <- SymbolsGSE39612[-NAs]
edataGSE39612 <- edataGSE39612[-NAs,]

rownames(edataGSE39612) <- SymbolsGSE39612

# SUMMARY: 2 BCC, 2 ISCC, 12 PMCC, 5 MMCC

# ***************************************************************************************************************************
# GSE42109 AFFY Study
# ***************************************************************************************************************************
GSE42109 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE42109") 
setwd("./ReferenceFiles/Samples/Microarray/GSE42109")
file.rename(list.files(pattern = "GSM103.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM103.*.CEL")))
setwd(pwd)
sampleNames(GSE42109) <- sub("\\.CEL$", "", sampleNames(GSE42109));
arrayQualityMetrics(GSE42109, outdir = "GeneratedFiles/QualityAnalysis/GSE42109", force = TRUE, do.logtransform = TRUE)

outlier <- "GSM1032618"
GSE42109 <- GSE42109[,setdiff(sampleNames(GSE42109),outlier)]
arrayQualityMetrics(GSE42109, outdir = "GeneratedFiles/QualityAnalysis/GSE42109_2", force = TRUE, do.logtransform = TRUE)

GSE42109Norm <- affy::rma(GSE42109)
phenoGSE42109 <- pData(GSE42109Norm)
edataGSE42109 <- exprs(GSE42109Norm)

GSE42109Labels <- rep("BCC",10)

require("hgu133a2.db", character.only = TRUE)
genesIDGSE42109 <- rownames(edataGSE42109)
SymbolsGSE42109 <- getSYMBOL(genesIDGSE42109, "hgu133a2.db")
EntrezGSE42109 <- getEG(genesIDGSE42109, "hgu133a2.db")

NAs <- which(is.na(SymbolsGSE42109))
SymbolsGSE42109 <- SymbolsGSE42109[-NAs]
edataGSE42109 <- edataGSE42109[-NAs,]

rownames(edataGSE42109) <- SymbolsGSE42109

# SUMMARY: 10 BCC

# ***************************************************************************************************************************
# GSE42677 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE42677")
file.rename(list.files(pattern = "GSM104.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM104.*.CEL")))

# HGU133A2
setwd("./hgu133a2/")
GSE42677a2 <- read.celfiles(list.files())
sampleNames(GSE42677a2) <- sub("\\.CEL$", "", sampleNames(GSE42677a2));
setwd(pwd)
arrayQualityMetrics(GSE42677a2, outdir = "GeneratedFiles/QualityAnalysis/GSE42677a2", force = TRUE, do.logtransform = TRUE)

GSE42677a2Norm <- rma(GSE42677a2)
phenoGSE42677a2 <- pData(GSE42677a2Norm)
edataGSE42677a2 <- exprs(GSE42677a2Norm)

GSE42677a2Labels <- c(rep("AK",5), rep("ISCC",5))

require("hgu133a2.db", character.only = TRUE)
genesIDGSE42677a2 <- rownames(edataGSE42677a2)
SymbolsGSE42677a2 <- getSYMBOL(genesIDGSE42677a2, "hgu133a2.db")
EntrezGSE42677a2 <- getEG(genesIDGSE42677a2, "hgu133a2.db")

NAs <- which(is.na(SymbolsGSE42677a2))
SymbolsGSE42677a2 <- SymbolsGSE42677a2[-NAs]
edataGSE42677a2 <- edataGSE42677a2[-NAs,]

rownames(edataGSE42677a2) <- SymbolsGSE42677a2

# HGU133A2 SUMMARY: 5 AK, 5 ISCC 

# HGU133plus2
setwd("./ReferenceFiles/Samples/Microarray/GSE42677/hgu133plus2/")
GSE42677plus2 <- read.celfiles(list.files())
sampleNames(GSE42677plus2) <- sub("\\.CEL$", "", sampleNames(GSE42677plus2));
setwd(pwd)
arrayQualityMetrics(GSE42677plus2, outdir = "GeneratedFiles/QualityAnalysis/GSE42677plus2", force = TRUE, do.logtransform = TRUE)

outlier <- "GSM1047855"
GSE42677plus2 <- GSE42677plus2[,setdiff(sampleNames(GSE42677plus2),outlier)]
arrayQualityMetrics(GSE42677plus2, outdir = "GeneratedFiles/QualityAnalysis/GSE42677plu2_2", force = TRUE, do.logtransform = TRUE)

GSE42677plus2Norm <- rma(GSE42677plus2)
phenoGSE42677plus2 <- pData(GSE42677plus2Norm)
edataGSE42677plus2 <- exprs(GSE42677plus2Norm)

GSE42677plus2Labels <- rep("NSK",9)

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE42677plus2 <- rownames(edataGSE42677plus2)
SymbolsGSE42677plus2 <- getSYMBOL(genesIDGSE42677plus2, "hgu133plus2.db")
EntrezGSE42677plus2 <- getEG(genesIDGSE42677plus2, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE42677plus2))
SymbolsGSE42677plus2 <- SymbolsGSE42677plus2[-NAs]
edataGSE42677plus2 <- edataGSE42677plus2[-NAs,]

rownames(edataGSE42677plus2) <- SymbolsGSE42677plus2

# HGU133PLUS2 SUMMARY: 9 NSK

# SUMMARY: 5 AK, 5 ISCC, 9 NSK

# ***************************************************************************************************************************
# GSE45216 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE45216")
file.rename(list.files(pattern = "GSM109.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM109.*.CEL")))
setwd(pwd)
GSE45216 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE45216")
sampleNames(GSE45216) <- sub("\\.CEL$", "", sampleNames(GSE45216));
arrayQualityMetrics(GSE45216, outdir = "GeneratedFiles/QualityAnalysis/GSE45216", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1099229","GSM1099240","GSM1099263")
GSE45216 <- GSE45216[,setdiff(sampleNames(GSE45216),outlier)]
arrayQualityMetrics(GSE45216, outdir = "GeneratedFiles/QualityAnalysis/GSE45216_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1099241","GSM1099264")
GSE45216 <- GSE45216[,setdiff(sampleNames(GSE45216),outlier)]
arrayQualityMetrics(GSE45216, outdir = "GeneratedFiles/QualityAnalysis/GSE45216_3", force = TRUE, do.logtransform = TRUE)

GSE45216Norm <- affy::rma(GSE45216)
phenoGSE45216 <- pData(GSE45216Norm)
edataGSE45216 <- exprs(GSE45216Norm)

GSE45216Labels <- c(rep("ISCC",27),rep("AK",8))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE45216 <- rownames(edataGSE45216)
SymbolsGSE45216 <- getSYMBOL(genesIDGSE45216, "hgu133plus2.db")
EntrezGSE45216 <- getEG(genesIDGSE45216, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE45216))
SymbolsGSE45216 <- SymbolsGSE45216[-NAs]
edataGSE45216 <- edataGSE45216[-NAs,]

rownames(edataGSE45216) <- SymbolsGSE45216

# SUMMARY: 27 ISCC, 8 AK

# ***************************************************************************************************************************
# GSE46517 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE46517")
file.rename(list.files(pattern = "GSM113.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM113.*.CEL")))
GSE46517 <- read.celfiles(list.celfiles()) 
sampleNames(GSE46517) <- sub("\\.CEL$", "", sampleNames(GSE46517));
setwd(pwd)
arrayQualityMetrics(GSE46517, outdir = "GeneratedFiles/QualityAnalysis/GSE46517", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1131592","GSM1131593","GSM1131596","GSM1131634","GSM1131637","GSM1131638","GSM1131642","GSM1131657","GSM1131673","GSM1131676")
GSE46517 <- GSE46517[,setdiff(sampleNames(GSE46517),outlier)]
arrayQualityMetrics(GSE46517, outdir = "GeneratedFiles/QualityAnalysis/GSE46517_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1131588","GSM1131597","GSM1131609","GSM1131644","GSM1131651","GSM1131660","GSM1131590","GSM1131653","GSM1131670","GSM1131679")
GSE46517 <- GSE46517[,setdiff(sampleNames(GSE46517),outlier)]
arrayQualityMetrics(GSE46517, outdir = "GeneratedFiles/QualityAnalysis/GSE46517_3", force = TRUE, do.logtransform = TRUE)

GSE46517Norm <- rma(GSE46517)
phenoGSE46517 <- pData(GSE46517Norm)
edataGSE46517 <- exprs(GSE46517Norm)
GSE46517Labels <- c(rep("METMEL",31),rep("PRIMEL",25),rep("NEV",6),rep("NSK",6))

require("hgu133a.db", character.only = TRUE)
genesIDGSE46517 <- rownames(edataGSE46517)
SymbolsGSE46517 <- getSYMBOL(genesIDGSE46517, "hgu133a.db")
EntrezGSE46517 <- getEG(genesIDGSE46517, "hgu133a.db")

NAs <- which(is.na(SymbolsGSE46517))
SymbolsGSE46517 <- SymbolsGSE46517[-NAs]
edataGSE46517 <- edataGSE46517[-NAs,]

rownames(edataGSE46517) <- SymbolsGSE46517

# SUMMARY: 6 NSK, 6 NEV, 25 PRIMEL, 31 METMEL

# ***************************************************************************************************************************
# GSE50451 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE50451")
file.rename(list.files(pattern = "GSM12194.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM12194.*.CEL")))
setwd(pwd)
GSE50451 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE50451")
sampleNames(GSE50451) <- sub("\\.CEL$", "", sampleNames(GSE50451));
arrayQualityMetrics(GSE50451, outdir = "GeneratedFiles/QualityAnalysis/GSE50451", force = TRUE, do.logtransform = TRUE)

outlier <- "GSM1219455"
GSE50451 <- GSE50451[,setdiff(sampleNames(GSE50451),outlier)]

GSE50451Norm <- affy::rma(GSE50451)
phenoGSE50451 <- pData(GSE50451Norm)
edataGSE50451 <- exprs(GSE50451Norm)

GSE50451Labels <- c(rep("MMCC",3),"PMCC",rep("MMCC",2),rep("PMCC",2),rep("MMCC",3),"PMCC",rep("MMCC",3),"PMCC","MMCC",rep("PMCC",4),"MMCC")

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE50451 <- rownames(edataGSE50451)
SymbolsGSE50451 <- getSYMBOL(genesIDGSE50451, "hgu133plus2.db")
EntrezGSE50451 <- getEG(genesIDGSE50451, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE50451))
SymbolsGSE50451 <- SymbolsGSE50451[-NAs]
edataGSE50451 <- edataGSE50451[-NAs,]

rownames(edataGSE50451) <- SymbolsGSE50451

# SUMMARY: 9 PMCC, 13 MMCC

# ***************************************************************************************************************************
# GSE52471 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE52471")
file.rename(list.files(pattern = "GSM12674.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM12674.*.CEL")))
setwd(pwd)
GSE52471 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE52471")
sampleNames(GSE52471) <- sub("\\.CEL$", "", sampleNames(GSE52471));
arrayQualityMetrics(GSE52471, outdir = "GeneratedFiles/QualityAnalysis/GSE52471", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1267418","GSM1267444","GSM1267445","GSM1267446")
GSE52471 <- GSE52471[,setdiff(sampleNames(GSE52471),outlier)]
arrayQualityMetrics(GSE52471, outdir = "GeneratedFiles/QualityAnalysis/GSE52471_2", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1267447","GSM1267448","GSM1267449")
GSE52471 <- GSE52471[,setdiff(sampleNames(GSE52471),outlier)]
arrayQualityMetrics(GSE52471, outdir = "GeneratedFiles/QualityAnalysis/GSE52471_3", force = TRUE, do.logtransform = TRUE)

GSE52471Norm <- affy::rma(GSE52471)
phenoGSE52471 <- pData(GSE52471Norm)
edataGSE52471 <- exprs(GSE52471Norm)

GSE52471Labels <- c(rep("PS",14),rep("NSK",10))

require("hgu133a2.db", character.only = TRUE)
genesIDGSE52471 <- rownames(edataGSE52471)
SymbolsGSE52471 <- getSYMBOL(genesIDGSE52471, "hgu133a2.db")
EntrezGSE52471 <- getEG(genesIDGSE52471, "hgu133a2.db")

NAs <- which(is.na(SymbolsGSE52471))
SymbolsGSE52471 <- SymbolsGSE52471[-NAs]
edataGSE52471 <- edataGSE52471[-NAs,]

rownames(edataGSE52471) <- SymbolsGSE52471

# SUMMARY: 10 NSK, 14 PS

# ***************************************************************************************************************************
# GSE53223 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE53223")
file.rename(list.files(pattern = "GSM12879.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM12879.*.CEL")))
setwd(pwd)
GSE53223 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE53223") 
sampleNames(GSE53223) <- sub("\\.CEL$", "", sampleNames(GSE53223));
arrayQualityMetrics(GSE53223, outdir = "GeneratedFiles/QualityAnalysis/GSE53223", force = TRUE, do.logtransform = TRUE)

outlier <- c("GSM1287903","GSM1287904","GSM1287910","GSM1287911","GSM1287913")
GSE53223 <- GSE53223[,setdiff(sampleNames(GSE53223),outlier)]
arrayQualityMetrics(GSE53223, outdir = "GeneratedFiles/QualityAnalysis/GSE53223_2", force = TRUE, do.logtransform = TRUE)

GSE53223Norm <- affy::rma(GSE53223)
phenoGSE53223 <- pData(GSE53223Norm)
edataGSE53223 <- exprs(GSE53223Norm)

GSE53223Labels <- c(rep("NEV",8),rep("NSK",5))

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE53223 <- rownames(edataGSE53223)
SymbolsGSE53223 <- getSYMBOL(genesIDGSE53223, "hgu133plus2.db")
EntrezGSE53223 <- getEG(genesIDGSE53223, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE53223))
SymbolsGSE53223 <- SymbolsGSE53223[-NAs]
edataGSE53223 <- edataGSE53223[-NAs,]

rownames(edataGSE53223) <- SymbolsGSE53223

# SUMMARY: 8 NEV, 5 NSK

# ***************************************************************************************************************************
# GSE53462 LUMI Study
# ***************************************************************************************************************************
GSE53462 <- lumiR(file = "ReferenceFiles/Samples/Microarray/GSE53462/GSE53462_non_normalized.txt", 
                  lib.mapping="lumiHumanIDMapping")

outlier <- c(paste0("GSM12940",33:37))
GSE53462 <- GSE53462[,setdiff(sampleNames(GSE53462),outlier)]
arrayQualityMetrics(GSE53462, outdir = "GeneratedFiles/QualityAnalysis/GSE53462_1", force = TRUE, do.logtransform = TRUE)

GSE53462Norm <- lumiExpresso(GSE53462)
phenoGSE53462 <- pData(GSE53462Norm)
edataGSE53462 <- exprs(GSE53462Norm)


GSE53462Labels <- as.vector("BCC","BCC","BCC","BCC","BCC","SCC","SCC","BCC","BCC","BCC","BCC",
                            "BCC","BCC","BCC","SCC","BCC","BCC","BCC","SCC","BCC","SCC")

require("lumiHumanAll.db", character.only = TRUE)
genesIDGSE53462 <- rownames(edataGSE53462)
SymbolsGSE53462 <- getSYMBOL(genesIDGSE53462, "lumiHumanAll.db")
EntrezGSE53462 <- getEG(genesIDGSE53462, "lumiHumanAll.db")

NAs <- which(is.na(SymbolsGSE53462))
SymbolsGSE53462 <- SymbolsGSE53462[-NAs]
edataGSE53462 <- edataGSE53462[-NAs,]

rownames(edataGSE53462) <- SymbolsGSE53462

# SUMMARY: 16 BCC, 5 ISCC

# ***************************************************************************************************************************
# GSE82105 AFFY Study
# ***************************************************************************************************************************
setwd("./ReferenceFiles/Samples/Microarray/GSE82105")
file.rename(list.files(pattern = "GSM2183.*.CEL"), sub("\\_.*", ".CEL",list.files(pattern = "GSM2183.*.CEL")))
setwd(pwd)
GSE82105 <- ReadAffy(celfile.path="ReferenceFiles/Samples/Microarray/GSE82105")
sampleNames(GSE82105) <- sub("\\.CEL$", "", sampleNames(GSE82105));
arrayQualityMetrics(GSE82105, outdir = "GeneratedFiles/QualityAnalysis/GSE82105", force = TRUE, do.logtransform = TRUE)

GSE82105Norm <- affy::rma(GSE82105)
phenoGSE82105 <- pData(GSE82105Norm)
edataGSE82105 <- exprs(GSE82105Norm)

GSE82105Labels <- rep("METMEL",6)

require("hgu133plus2.db", character.only = TRUE)
genesIDGSE82105 <- rownames(edataGSE82105)
SymbolsGSE82105 <- getSYMBOL(genesIDGSE82105, "hgu133plus2.db")
EntrezGSE82105 <- getEG(genesIDGSE82105, "hgu133plus2.db")

NAs <- which(is.na(SymbolsGSE82105))
SymbolsGSE82105 <- SymbolsGSE82105[-NAs]
edataGSE82105 <- edataGSE82105[-NAs,]

rownames(edataGSE82105) <- SymbolsGSE82105

# SUMMARY: 6 METMEL

#############################################################################################################################
#############################################################################################################################
# Raw data acquisition and preprocessing of RNA-Seq series
#############################################################################################################################
#############################################################################################################################
setwd(pwd)

# ***************************************************************************************************************************
# E-MTAB-5678 RNA-seq Study
# ***************************************************************************************************************************
i = 1
countf <- vector(mode="character", length=0)
countPath <- "ReferenceFiles/Samples/RNAseq/CountFiles/E-MTAB-5678"
countFiles <- list.files(countPath)

counts = readDGE(countFiles, path = countPath)$counts
noint  =  rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms  =  cpm(counts)
keep  =  rowSums(cpms  > 1) >= 2 & !noint
counts  =  counts[keep,]

write.table(counts, file="ReferenceFiles/counts_AE5678.txt", row.names=TRUE, col.names=TRUE)
AE5678 <- read.table("ReferenceFiles/counts_AE5678.txt")
factorsRnaSeq <- read.delim("ReferenceFiles/design_AE5678.txt", sep = "\t")

# Gene annotation from Biomart --------------------------------------------
biomartHomoSapiens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
atributos = listAttributes(biomartHomoSapiens)
atributos[grep("gc", atributos$description, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id","external_gene_name", "percentage_gene_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(AE5678), mart=biomartHomoSapiens)
head(myannot)

# Length info
geneLength = read.csv("ReferenceFiles/Genes_length_Homo_Sapiens.csv", header = TRUE, sep = ",")
geneLength = geneLength[1:63677,]
geneLength = geneLength[,-2]
geneLength = geneLength[,-2]

AE5678Data = readData(data = AE5678,
                      factors = factorsRnaSeq,
                      gc = myannot[,c(1,2)],
                      length = geneLength[,1:2])

# LOW COUNTS
myCounts = dat(AE5678Data, type = "countsbio")
pdf("AE5678_Plot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()

# Normalization -----------------------------------------------------------
# Correcting GC content bias with cqn library
myGCannot = myannot$percentage_gene_gc_content
names(myGCannot) = myannot$ensembl_gene_id
myGCannot = myGCannot[rownames(AE5678)]
summary(myGCannot)  # NA's!!

mygenes = intersect(rownames(AE5678), myannot$ensembl_gene_id)
mylength = setNames(geneLength[match(mygenes, geneLength$Gene_stable_ID), 2],
                    nm = mygenes)
geneNames <- rownames(AE5678)
AE5678 <- sapply(AE5678,as.numeric,2)
rownames(AE5678) <- geneNames
AE5678 <- AE5678[mygenes,]
mycqn <- cqn(as.matrix(AE5678), lengths = mylength,
             x = myGCannot[mygenes], sizeFactors = apply(AE5678, 2, sum),
             verbose = TRUE)

# Computing the normalized values
rnaseqCQN = mycqn$y + mycqn$offset
boxplot(rnaseqCQN)
min(rnaseqCQN)
AE5678Matrix = rnaseqCQN - min(rnaseqCQN) + 1
rownames(AE5678Matrix) = myannot$external_gene_name[which(myannot$ensembl_gene_id == rownames(AE5678Matrix))]

AE5678Labels <- c(rep("AK",13),rep("NSK",4),rep("ISCC",5))

# ***************************************************************************************************************************
# GSE54456 RNA-seq Study
# ***************************************************************************************************************************
# Duplicated and removed samples: GSM1315754 (PS), GSM1315768 (NSK), GSM1315618 (PS), GSM1315621 (NSK), GSM1315637 (PS)
i = 1
countf <- vector(mode="character", length=0)
countPath <- "ReferenceFiles/Samples/RNAseq/CountFiles/GSE54456"
countFiles <- list.files(countPath)

counts = readDGE(countFiles, path = countPath)$counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms  > 1) >= 2 & !noint
counts = counts[keep,]

write.table(counts, file="ReferenceFiles/counts_GSE54456.txt", row.names=TRUE, col.names=TRUE)
GSE54456 <- read.table("ReferenceFiles/counts_GSE54456.txt")
factorsRnaSeq <- read.delim("ReferenceFiles/design_GSE54456.txt",sep = "\t")

# Gene annotation from Biomart --------------------------------------------
biomartHomoSapiens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
atributos = listAttributes(biomartHomoSapiens)
atributos[grep("gc", atributos$description, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id","external_gene_name", "percentage_gene_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(GSE54456), mart=biomartHomoSapiens)
head(myannot)

# Length info
geneLength = read.csv("ReferenceFiles/Genes_length_Homo_Sapiens.csv", header = TRUE, sep = ",")
geneLength = geneLength[1:63677,]
geneLength = geneLength[,-2]
geneLength = geneLength[,-2]

GSE54456Data = readData(data = GSE54456,
                        factors = factorsRnaSeq,
                        gc = myannot[,c(1,2)],
                        length = geneLength[,1:2])


# LOW COUNTS
myCounts = dat(GSE54456Data, type = "countsbio")
pdf("GSE54456_Plot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()

# Normalization -----------------------------------------------------------
### Correcting GC content bias with cqn library
myGCannot = myannot$percentage_gene_gc_content
names(myGCannot) = myannot$ensembl_gene_id
myGCannot = myGCannot[rownames(GSE54456)]
summary(myGCannot)  # NA's!!

mygenes = intersect(rownames(GSE54456), myannot$ensembl_gene_id)
mylength = setNames(geneLength[match(mygenes, geneLength$Gene_stable_ID), 2],
                    nm = mygenes)
mylengthFiltered = mylength
mygenesFiltered = mygenes
myGCannotFiltered = myGCannot[mygenesFiltered]
mycqn <- cqn(GSE54456[mygenesFiltered,], lengths = mylengthFiltered,
             x = myGCannotFiltered, sizeFactors = apply(GSE54456, 2, sum),
             verbose = TRUE)

# Computing the normalized values
rnaseqCQN = mycqn$y + mycqn$offset
boxplot(rnaseqCQN)
min(rnaseqCQN)
GSE54456Matrix = rnaseqCQN - min(rnaseqCQN) + 1
rownames(GSE54456Matrix) = myannot$external_gene_name[which(myannot$ensembl_gene_id == rownames(GSE54456Matrix))]

GSE54456Labels <- as.vector(factorsRnaSeq$Condition)

# ***************************************************************************************************************************
# GSE67785 RNA-seq Study
# ***************************************************************************************************************************
i = 1
countf <- vector(mode="character", length=0)
countPath <- "ReferenceFiles/Samples/RNAseq/CountFiles/GSE67785"
countFiles <- list.files(countPath)

counts = readDGE(countFiles, path = countPath)$counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms  > 1) >= 2 & !noint
counts = counts[keep,]

write.table(counts, file="ReferenceFiles/counts_GSE67785.txt", row.names=TRUE, col.names=TRUE)
GSE67785 <- read.table("ReferenceFiles/counts_GSE67785.txt")
factorsRnaSeq <- read.delim("ReferenceFiles/design_GSE67785.txt",sep = "\t")

# Gene annotation from Biomart --------------------------------------------
biomartHomoSapiens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
atributos = listAttributes(biomartHomoSapiens)
atributos[grep("gc", atributos$description, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id","external_gene_name", "percentage_gene_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(GSE67785), mart=biomartHomoSapiens)
head(myannot)

# Length info
geneLength = read.csv("ReferenceFiles/Genes_length_Homo_Sapiens.csv", header = TRUE, sep = ",")
geneLength = geneLength[1:63677,]
geneLength = geneLength[,-2]
geneLength = geneLength[,-2]

GSE67785Data = readData(data = GSE67785,
                        factors = factorsRnaSeq,
                        gc = myannot[,c(1,2)],
                        length = geneLength[,1:2])

# LOW COUNTS
myCounts = dat(GSE67785Data, type = "countsbio")
pdf("GSE67785_Plot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()

# Normalization -----------------------------------------------------------
### Correcting GC content bias with cqn library
myGCannot = myannot$percentage_gene_gc_content
names(myGCannot) = myannot$ensembl_gene_id
myGCannot = myGCannot[rownames(GSE67785)]
summary(myGCannot)  # NA's!!

mygenes = intersect(rownames(GSE67785), myannot$ensembl_gene_id)
mylength = setNames(geneLength[match(mygenes, geneLength$Gene_stable_ID), 2],
                    nm = mygenes)
mylengthFiltered = mylength
mygenesFiltered = mygenes
myGCannotFiltered = myGCannot[mygenesFiltered]
mycqn <- cqn(GSE67785[mygenesFiltered,], lengths = mylengthFiltered,
             x = myGCannotFiltered, sizeFactors = apply(GSE67785, 2, sum),
             verbose = TRUE)

# Computing the normalized values
rnaseqCQN = mycqn$y + mycqn$offset
boxplot(rnaseqCQN)
min(rnaseqCQN)
GSE67785Matrix = rnaseqCQN - min(rnaseqCQN) + 1
rownames(GSE67785Matrix) = myannot$external_gene_name[which(myannot$ensembl_gene_id == rownames(GSE67785Matrix))]

GSE67785Labels <- as.vector(factorsRnaSeq$Condition)

# ***************************************************************************************************************************
# GSE84293 RNA-seq Study
# ***************************************************************************************************************************
i = 1
countf <- vector(mode="character", length=0)
countPath <- "ReferenceFiles/Samples/RNAseq/CountFiles/GSE84293"
countFiles <- list.files(countPath)

counts = readDGE(countFiles, path = countPath)$counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms  > 1) >= 2 & !noint
counts = counts[keep,]

write.table(counts, file="ReferenceFiles/counts_GSE84293.txt", row.names=TRUE, col.names=TRUE)
GSE84293 <- read.table("ReferenceFiles/counts_GSE84293.txt")
factorsRnaSeq <- read.delim("ReferenceFiles/design_GSE84293.txt",sep = "\t")

# Gene annotation from Biomart --------------------------------------------
biomartHomoSapiens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
atributos = listAttributes(biomartHomoSapiens)
atributos[grep("gc", atributos$description, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id","external_gene_name", "percentage_gene_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(GSE84293), mart=biomartHomoSapiens)
head(myannot)

# Length info
geneLength = read.csv("ReferenceFiles/Genes_length_Homo_Sapiens.csv", header = TRUE, sep = ",")
geneLength = geneLength[1:63677,]
geneLength = geneLength[,-2]
geneLength = geneLength[,-2]

GSE84293Data = readData(data = GSE84293,
                        factors = factorsRnaSeq,
                        gc = myannot[,c(1,2)],
                        length = geneLength[,1:2])


# LOW COUNTS
myCounts = dat(GSE84293Data, type = "countsbio")
pdf("GSE84293_Plot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()

# Normalization -----------------------------------------------------------
### Correcting GC content bias with cqn library
myGCannot = myannot$percentage_gene_gc_content
names(myGCannot) = myannot$ensembl_gene_id
myGCannot = myGCannot[rownames(GSE84293)]
summary(myGCannot)  # NA's!!

mygenes = intersect(rownames(GSE84293), myannot$ensembl_gene_id)
mylength = setNames(geneLength[match(mygenes, geneLength$Gene_stable_ID), 2],
                    nm = mygenes)
mylengthFiltered = mylength
mygenesFiltered = mygenes
myGCannotFiltered = myGCannot[mygenesFiltered]
mycqn <- cqn(GSE84293[mygenesFiltered,], lengths = mylengthFiltered,
             x = myGCannotFiltered, sizeFactors = apply(GSE84293, 2, sum),
             verbose = TRUE)

# Computing the normalized values
rnaseqCQN = mycqn$y + mycqn$offset
boxplot(rnaseqCQN)
min(rnaseqCQN)
GSE84293Matrix = rnaseqCQN - min(rnaseqCQN) + 1
rownames(GSE84293Matrix) = myannot$external_gene_name[which(myannot$ensembl_gene_id == rownames(GSE84293Matrix))]

GSE84293Labels <- as.vector(factorsRnaSeq$Condition)

# ***************************************************************************************************************************
# GSE98394 RNA-seq Study
# ***************************************************************************************************************************
i = 1
countf <- vector(mode="character", length=0)
countPath <- "ReferenceFiles/Samples/RNAseq/CountFiles/GSE98394"
countFiles <- list.files(countPath)

counts = readDGE(countFiles, path = countPath)$counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms  > 1) >= 2 & !noint
counts = counts[keep,]

write.table(counts, file="ReferenceFiles/counts_GSE98394.txt", row.names=TRUE, col.names=TRUE)
GSE98394 <- read.table("ReferenceFiles/counts_GSE98394.txt")
factorsRnaSeq <- read.delim("ReferenceFiles/design_GSE98394.txt",sep = "\t")

# Gene annotation from Biomart --------------------------------------------
biomartHomoSapiens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
atributos = listAttributes(biomartHomoSapiens)
atributos[grep("gc", atributos$description, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id","external_gene_name", "percentage_gene_gc_content", "gene_biotype"),
                filters = "ensembl_gene_id", values=rownames(GSE98394), mart=biomartHomoSapiens)
head(myannot)

# Length info
geneLength = read.csv("ReferenceFiles/Genes_length_Homo_Sapiens.csv", header = TRUE, sep = ",")
geneLength = geneLength[1:63677,]
geneLength = geneLength[,-2]
geneLength = geneLength[,-2]

GSE98394Data = readData(data = GSE98394,
                        factors = factorsRnaSeq,
                        gc = myannot[,c(1,2)],
                        length = geneLength[,1:2])

# LOW COUNTS
myCounts = dat(GSE98394Data, type = "countsbio")
pdf("GSE98394_Plot_distribution_lowcounts.pdf", width = 7, height = 7)
explo.plot(myCounts, plottype = "barplot", samples = NULL)
dev.off()

# Normalization -----------------------------------------------------------
### Correcting GC content bias with cqn library
myGCannot = myannot$percentage_gene_gc_content
names(myGCannot) = myannot$ensembl_gene_id
myGCannot = myGCannot[rownames(GSE98394)]
summary(myGCannot)  # NA's!!

mygenes = intersect(rownames(GSE98394), myannot$ensembl_gene_id)
mylength = setNames(geneLength[match(mygenes, geneLength$Gene_stable_ID), 2],
                    nm = mygenes)
mylengthFiltered = mylength
mygenesFiltered = mygenes
myGCannotFiltered = myGCannot[mygenesFiltered]
mycqn <- cqn(GSE98394[mygenesFiltered,], lengths = mylengthFiltered,
             x = myGCannotFiltered, sizeFactors = apply(GSE98394, 2, sum),
             verbose = TRUE)

# Computing the normalized values
rnaseqCQN = mycqn$y + mycqn$offset
boxplot(rnaseqCQN)
min(rnaseqCQN)
GSE98394Matrix = rnaseqCQN - min(rnaseqCQN) + 1
rownames(GSE98394Matrix) = myannot$external_gene_name[which(myannot$ensembl_gene_id == rownames(GSE98394Matrix))]

GSE98394Labels <- as.vector(factorsRnaSeq$Condition)

#############################################################################################################################
#############################################################################################################################
# Saving preprocessed datasets (matrices and labels)
#############################################################################################################################
#############################################################################################################################
save(edataGSE02503,edataGSE03189,edataGSE06710,edataGSE07553,edataGSE13355,edataGSE14905,edataGSE15605,
     edataGSE30999,edataGSE32407,edataGSE32628,edataGSE32924,edataGSE36150,edataGSE39612,edataGSE42109,
     edataGSE42677a2,edataGSE42677plus2,edataGSE45216,edataGSE46517,edataGSE50451,edataGSE52471,
     edataGSE53223,edataGSE53462,edataGSE82105,AE5678Matrix,GSE54456Matrix,GSE67785Matrix,GSE84293Matrix,
     GSE98394Matrix, file = "EDATA_RAW_27_SERIES_28_BATCHES.RData")

save(GSE02503Labels,GSE03189Labels,GSE06710Labels,GSE07553Labels,GSE13355Labels,GSE14905Labels,GSE15605Labels,
     GSE30999Labels,GSE32407Labels,GSE32628Labels,GSE32924Labels,GSE36150Labels,GSE39612Labels,GSE42109Labels,
     GSE42677a2Labels,GSE42677plus2Labels,GSE45216Labels,GSE46517Labels,GSE50451Labels,GSE52471Labels,
     GSE53223Labels,GSE53462Labels,GSE82105Labels,AE5678Labels,GSE54456Labels,GSE67785Labels,GSE84293Labels,
     GSE98394Labels, file = "LABELS_27_SERIES_28_BATCHES.RData")
