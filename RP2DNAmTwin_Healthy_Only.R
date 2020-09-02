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

Gender <- factor(targets$gender)
design <- model.matrix(~Gender)
fit.reduced <- lmFit(Mval, design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))

# Gender results of diff DNAm in RA 
#         (Intercept) Gender1
# Down        213789    8345
# NotSig        6049  471818
# Up          263615    3290



# Look at top values differentially methylated between Healthy/Unhealthy
topRA <-topTable(fit.reduced, coef = 'Gender1')
topRA

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#head(ann450k)

ann450kSub <- ann450k[match(rownames(Mval),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit.reduced, coef = 'Gender1', num=Inf,  genelist=ann450kSub)
#head(DMPs)

write.csv(DMPs,"~/Desktop/RP2/RP2DNAmTwin_DMPsGender_Healthy_Only.csv", row.names = FALSE)
