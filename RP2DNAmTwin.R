#Load dependencies

library(missMethyl)
library(limma)
library(minfi)
library(shinyMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)


# Set directory
baseDir <- "~/Desktop/RP2/twin_methylation"
list.files(baseDir)

# Load in target file and DNAm data
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)

#Quality control
#summary <- shinySummarize(RGset)
#runShinyMethyl(summary)


# Subset quantile within array normalisation (SWAN)
mSet <- preprocessRaw(RGset)
mSetSWAN <- SWAN(mSet, verbose = TRUE)


# Filter out poor quality probes
detP <- detectionP(RGset)
keep <- rowSums(detP < 0.01) == ncol(RGset)
mSetSWAN <- mSetSWAN[keep,]


#Extract Beta & M values
meth <- getMeth(mSetSWAN)
unmeth <- getUnmeth(mSetSWAN)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSWAN)


# Testing for differential nethylation with LIMMA
#group <- factor(targets$Sample_Group)
#gender <- factor(targets$gender)
SibShip <- factor(targets$Family)
Status <- factor(targets$Sample_Group, level = c("healthy_cotwin", "RA"))

design <- model.matrix(~SibShip+Status)
#design
#dupcor <- duplicateCorrelation(Mval, design, block=targets$Family)
#dupcor$consensus.correlation  # 0.364 

fit.reduced <- lmFit(Mval, design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))

# Example of one SibShip and the Disease status overall
#           SibShip9446 StatusRA
# Down          6688        2
# NotSig      445396   481500
# Up           29425        7


# Look at top values differentially methylated between Healthy/Unhealthy
topRA <-topTable(fit.reduced, coef = "StatusRA")
topRA

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

ann450kSub <- ann450k[match(rownames(Mval),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs_Status <- topTable(fit.reduced, coef = "StatusRA", num=Inf,  genelist=ann450kSub)
head(DMPs)


write.csv(DMPs_Status,"~/Desktop/RP2/twin_methylation/RP2DNAmTwin_DMPsDiseaseStatus.csv", row.names = FALSE)






