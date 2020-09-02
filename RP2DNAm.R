#Load dependencies

library(missMethyl)
library(limma)
library(minfi)
library(shinyMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)

# Set directory

baseDir <- "~/Desktop/RP2/DNAmInd/ind_methylation"
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
#par(mfrow=c(1,2), cex=1.25)
#densityByProbeType(mSet[,1], main = "Raw")
#densityByProbeType(mSetSWAN[,1], main = "SWAN")


# Filter out poor quality probes
detP <- detectionP(RGset)
keep <- rowSums(detP < 0.01) == ncol(RGset)
mSetSWAN <- mSetSWAN[keep,]


#Extract Beta & M values

meth <- getMeth(mSetSWAN)
unmeth <- getUnmeth(mSetSWAN)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSWAN)

write.csv(Mval, "RP2DNAmInd_Female_Mval.csv")


# Multidimensional scaling plots based on "Gender", "Sample_Group"
#par(mfrow=c(1,2), cex=1.25)
#plotMDS(Mval, labels=targets$sm, col=as.integer(factor(targets$Gender)))
#legend("topright",legend=c("Female","Male"),pch=16,cex=1.2,col=1:2)  #
#plotMDS(Mval, labels=targets$sm, col=as.integer(factor(targets$Sample_Group)))
#legend("topright",legend=c("Poor","Good"),pch=16,cex=1.2,col=1:2) #Is poor/good the right way around?


# Testing for differential nethylation with LIMMA
group <- factor(targets$Gender)
id <- factor(targets$sm)
design <- model.matrix(~group)
design

fit.reduced <- lmFit(Mval, design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))

#(Intercept) groupmale
#Down        203804      6106
#NotSig        5229    465515
#Up          265248      2660

# Look at top values differentially methylated between gender
top<-topTable(fit.reduced)
#top
# logFC    AveExpr         t      P.Value    adj.P.Val        B
# cg06489418 -3.825631 -0.6311413 -53.34186 1.959701e-62 4.354686e-57 131.7673
# cg01461762 -3.764818 -0.7652439 -53.07538 2.849861e-62 4.354686e-57 131.4059
# cg09790289 -4.440151 -1.0570355 -52.97308 3.292058e-62 4.354686e-57 131.2667
# cg07867687 -4.614025 -1.2897248 -52.89561 3.672663e-62 4.354686e-57 131.1611
# cg13311758 -4.014251 -1.3061818 -52.72338 4.686480e-62 4.445417e-57 130.9257
# cg09720515 -4.342337 -1.0725416 -52.13442 1.084851e-61 8.575404e-57 130.1148
# cg21860846 -3.024446 -0.5289122 -51.79832 1.758472e-61 1.191442e-56 129.6477
# cg10088372 -4.695783 -1.0325861 -51.08285 4.965846e-61 2.944008e-56 128.6430
# cg10912974 -4.032036 -0.8639436 -50.36602 1.424720e-60 7.004014e-56 127.6218
# cg17831869 -3.811476 -1.1054879 -50.34178 1.476765e-60 7.004014e-56 127.5870


cpgs <- rownames(top)
par(mfrow=c(2,2))
for(i in 1:4){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("Female","Male"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgs[i],cex.main=1.5)
}

#---------------------------------------------------------------------------------------#
#Remove unwanted variation (RUV)

# get M-values for ALL probes
meth_RUV <- getMeth(mSet)
unmeth_RUV <- getUnmeth(mSet)
Mval_RUV <- log2((meth + 100)/(unmeth + 100))
# setup the factor of interest
group_RUV <- factor(targets$Gender)
# extract Illumina negative control data
INCs <- getINCs(RGset)

#head(INCs)

# add negative control data to M-values
Mc <- rbind(Mval_RUV,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)

table(ctl1)
ctl1
#FALSE   TRUE 
#474281    613 

# Rank CpG's to their association with factor of interest (Gender)
rfit1 <- RUVfit(Y = Mc, X = group_RUV, ctl = ctl1) # Stage 1 analysis
rfit2 <- RUVadj(Y = Mc, fit = rfit1)

# Designate CpGs least associated with factor of interest as empirical control probes (ECPs)
top1 <- topRUV(rfit2, num=Inf, p.BH = 1)

head(top1)
# F.p F.p.BH p_X1.male p.BH_X1.male b_X1.male     sigma2 var.b_X1.male fit.ctl       mean
# cg06489418 1e-24  1e-24     1e-24        1e-24 -3.874967 0.04516863   0.005024522   FALSE -0.6311413
# cg05086798 1e-24  1e-24     1e-24        1e-24 -3.350616 0.03954369   0.004398809   FALSE -1.0036306
# cg08560117 1e-24  1e-24     1e-24        1e-24 -4.047844 0.05991521   0.006664920   FALSE -0.5651509
# cg01461762 1e-24  1e-24     1e-24        1e-24 -3.739259 0.05206580   0.005791759   FALSE -0.7652439
# cg14350469 1e-24  1e-24     1e-24        1e-24 -4.519226 0.07729009   0.008597687   FALSE -0.7663017
# cg05130312 1e-24  1e-24     1e-24        1e-24  3.238823 0.04033337   0.004486652   FALSE  0.5886610

ctl2 <- rownames(Mval_RUV) %in% rownames(top1[top1$p.BH_X1.male > 0.5,])
table(ctl2)
# ctl2
# FALSE   TRUE 
# 22836 451445 

# Use ECPs to perfrom diff-meth with RUV-inverse, which is adjusted for unnwanted variation
rfit3 <- RUVfit(Y = Mval_RUV, X = group_RUV, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = Mval_RUV, fit = rfit3)

# Look at table of top results
topRUV(rfit4)
# F.p F.p.BH p_X1.male p.BH_X1.male b_X1.male     sigma2 var.b_X1.male fit.ctl        mean
# cg03028216 1e-24  1e-24     1e-24        1e-24 -3.449531 0.02458445   0.002356247   FALSE -0.55984257
# cg24863802 1e-24  1e-24     1e-24        1e-24 -3.425848 0.02551638   0.002445566   FALSE -0.44713634
# cg05204037 1e-24  1e-24     1e-24        1e-24 -3.422250 0.02578959   0.002471752   FALSE -0.27955154
# cg06144999 1e-24  1e-24     1e-24        1e-24 -3.858157 0.03283041   0.003146564   FALSE -0.50916412
# cg08656747 1e-24  1e-24     1e-24        1e-24 -4.457238 0.04398627   0.004215775   FALSE -0.89416966
# cg13516940 1e-24  1e-24     1e-24        1e-24 -3.954188 0.03467446   0.003323304   FALSE -0.50638016
# cg19586382 1e-24  1e-24     1e-24        1e-24 -3.820285 0.03248260   0.003113230   FALSE -0.72507905
# cg08841290 1e-24  1e-24     1e-24        1e-24 -4.026385 0.03746576   0.003590830   FALSE -0.50618013
# cg18091964 1e-24  1e-24     1e-24        1e-24 -3.196472 0.02471522   0.002368781   FALSE -0.01296737
# cg01079126 1e-24  1e-24     1e-24        1e-24 -3.575420 0.03100792   0.002971892   FALSE -1.11944145

# Visualising effect of RUVm adjustment
Madj <- getAdj(Mval_RUV, rfit3) # get adjusted values

par(mfrow=c(1,2))
plotMDS(Mval_RUV, labels=targets$sm, col=as.integer(factor(targets$Gender)),
        main="Unadjusted", gene.selection = "common")
legend("topleft",legend=c("Male","Female"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=targets$sm, col=as.integer(factor(targets$Gender)),
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topleft",legend=c("Male","Female"),pch=16,cex=1,col=1:2)

# Section of adjusting K val of RUV akipped for time being


# Differential variability of Gender
fitvarSWAN <- varFit(Mval, design = design)
summary(decideTests(fitvarSWAN))

#           (Intercept) groupmale
# Down             0        10  #HypoVariable sites??
# NotSig           2    473754
# Up          474279       517  # Hypervariable sites??

topDV <- topVar(fitvarSWAN)
# topDV
# SampleVar LogVarRatio DiffLevene         t      P.Value  Adj.P.Value
# cg21341170 1.9646530    3.583466  0.8612092 11.511076 1.525993e-20 7.237496e-15
# cg09566331 0.4765442    2.851231  0.8969998 10.138681 2.067914e-17 4.903861e-12
# cg23936476 2.8401731    3.130588  0.6110022  9.300022 1.698817e-15 2.644056e-10
# cg03823092 1.6665711    2.876075  0.7109942  9.241599 2.307277e-15 2.644056e-10
# cg12239365 1.0815379    3.527840  0.5528459  9.205502 2.787436e-15 2.644056e-10
# cg21782309 0.5693286    2.475173  0.9185402  9.006469 7.893803e-15 6.239802e-10
# cg01191832 0.4277516    2.518402  0.6222734  8.726434 3.397918e-14 2.302240e-09
# cg15035437 1.8883738    3.070685  0.5079440  8.580855 7.237406e-14 4.290705e-09
# cg09975171 1.2833770    2.463439  0.4149373  8.458011 1.367564e-13 7.206776e-09
# cg23685102 2.7042062    2.579782  0.5097771  8.259170 3.816931e-13 1.810298e-08


fitvarRUV <- varFit(Madj, design = design)
#Warning message:
#  791 very small variances detected, have been offset away from zero 
summary(decideTests(fitvarRUV))
# (       Intercept) groupmale
# Down             0    392020
# NotSig        2018     45340
# Up          472263     36921
# Doesn't semm right... RUV method not applicable to diffVAR?


# Beta values for top four differentially variable CpGs
cpgsDV <- rownames(topDV)
par(mfrow=c(2,2))
for(i in 1:4){
  stripchart(beta[rownames(beta)==cpgsDV[i],]~design[,2],method="jitter",
             group.names=c("Females","Males"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgsDV[i],cex.main=1.5)
}

#-----------------------------------------------------------------------------------#

# CpG Annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

ann450kSub <- ann450k[match(rownames(Mval),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit.reduced, coef='groupmale', num=Inf,  genelist=ann450kSub)
head(DMPs)

colnames(DMPs)[colnames(DMPs) == "adj.P.Val"] <- "AdjP"

autosomal <- DMPs[!(DMPs$chr=="chrX" | DMPs$chr=="chrY"),] # 463577 observatiosn
autosomal <- autosomal[!(autosomal$AdjP > 0.05),]  # 458969 sig observations

sex <- DMPs[(DMPs$chr=="chrX" | DMPs$chr=="chrY"),] # 10704 observations
sex <- sex[!(sex$AdjP > 0.05),] # 10083 sig observations
# This produces a different list than printing topTable(fit.reduced) ????

# Top named gene locations:
# ZMYM6 - Plays a role in the regulation of cell morphology and cytoskeletal organization.
# DDX4 - This gene encodes a DEAD box protein, which is a homolog of VASA proteins in Drosophila
#       and several other species. The gene is specifically expressed in the germ cell lineage in both sexes and functions in germ cell development
# CD37 - CD37 expression is restricted to cells of the immune system, with highest abundance on mature B cells
# TM9SF2 - The protein may play a role in small molecule transport or act as an ion channel. 
#          A pseudogene associated with this gene is located on the X chromosome.


#None of these genes have specific RA literature
#   A couple mentions in microarray tables etc.


# Removing coef = 1 in DMP assignment
DMPs_nocoef <- topTable(fit.reduced, num=Inf,  genelist=ann450kSub)
head(DMPs_nocoef)

# Now the same as toptable(fit.reduced)
# All top genes now on X chromosome too!

# Look into autosomal and sex loci separately
colnames(DMPs_nocoef)[colnames(DMPs_nocoef) == "adj.P.Val"] <- "AdjP"

autosomal <- DMPs_nocoef[!(DMPs_nocoef$chr=="chrX" | DMPs_nocoef$chr=="chrY"),] # 463577 observatiosn
autosomal <- autosomal[!(autosomal$AdjP > 0.05),]  # 928 sig observations

sex <- DMPs_nocoef[(DMPs_nocoef$chr=="chrX" | DMPs_nocoef$chr=="chrY"),] # 10704 observations
sex <- sex[!(sex$AdjP >0.05),] # 7838sig observations

Xchr <- DMPs_nocoef[(DMPs_nocoef$chr=="chrX"),] # 10664 observations
Xchr <- Xchr[!(Xchr$AdjP >0.05),] # 7818 sig observations

Ychr <- DMPs_nocoef[(DMPs_nocoef$chr=="chrY"),] # 40 observations
Ychr <- Ychr[!(Ychr$AdjP >0.05),] # 20 sig observations


# Looking into breakdown of methylation profiles

AutoOverMales <- autosomal[(autosomal$logFC > 0),] # 398
AutoUnderMales <- autosomal[(autosomal$logFC < 0),] # 530

# Is there chromosomal correlations to differential methylation?

AutoOverTab <- prop.table(table(AutoOverMales$chr)) * 100 
AutoOverTab
# chr1      chr10      chr11      chr12      chr13      chr14      chr15      chr16      chr17      chr18      chr19       chr2      chr20 
# 12.5628141  5.5276382  5.0251256  7.7889447  3.5175879  3.0150754  4.5226131  4.7738693  4.5226131  1.0050251  9.5477387  5.0251256  2.7638191 
# chr21      chr22       chr3       chr4       chr5       chr6       chr7       chr8       chr9 
# 1.2562814  0.5025126  5.2763819  4.2713568  4.2713568  5.7788945  2.5125628  3.2663317  3.2663317 

AutoUnderTab <- prop.table(table(AutoUnderMales$chr)) * 100
AutoUnderTab

# chr1      chr10      chr11      chr12      chr13      chr14      chr15      chr16      chr17      chr18      chr19       chr2      chr20 
# 9.2452830  4.3396226  5.0943396  4.7169811  1.5094340  3.2075472  3.0188679  4.3396226  5.8490566  0.9433962  6.4150943  8.3018868  3.9622642 
# chr21      chr22       chr3       chr4       chr5       chr6       chr7       chr8       chr9 
# 0.5660377  1.1320755  5.0943396  1.6981132  5.4716981 10.1886792  7.1698113  4.3396226  3.3962264 

#Seems some chromosomes are targeted for methylation more than others
# May just be size based - chr1 is 8% of total DNA


# Percentage point difference between each chromosome for Over and Under expression
AutoChrDiff <- AutoOverTab - AutoUnderTab 
AutoChrDiff
# chr1       chr10       chr11       chr12       chr13       chr14       chr15       chr16       chr17       chr18       chr19        chr2 
# 3.31753105  1.18801555 -0.06921399  3.07196359  2.00815398 -0.19247179  1.50374514  0.43424671 -1.32644354  0.06162890  3.13264435 -3.27676116 
# chr20       chr21       chr22        chr3        chr4        chr5        chr6        chr7        chr8        chr9 
# -1.19844506  0.69024367 -0.62956291  0.18204229  2.57324358 -1.20034133 -4.40978477 -4.65724851 -1.07329098 -0.12989476 

#-----------------------------------------------------------------------------------#
SexOverMales <- sex[(sex$logFC >0),] # 2262
SexUnderMales <- sex[(sex$logFC < 0),] # 5576

# What percentage of male Over/Under methylated profiles are on y/x
prop.table(table(SexOverMales$chr))      
# chrX        chrY 
# 0.992926614 0.007073386  
prop.table(table(SexUnderMales$chr))
# chrX         chrY 
# 0.9992826399 0.0007173601 
#-------------------------------------------------------------------------------------#
# Looking into genes associated with methylation 

Autogenes <- prop.table(table(autosomal$UCSC_RefGene_Name))
head(Autogenes)
# 34% don't have reference gene names


Sexgenes <- prop.table(table(sex$UCSC_RefGene_Name))
head(Sexgenes)
# 15% don't have reference gene names
# Multiple rows for the same gene/CpG location 


#------------------------------------------------------------------------------------#
# Write to table
fwrite(DMPs, "RP2DNAmInd_Gender.csv")
