# Load in dependencies
library(limma)
setwd("~/Desktop/RP2/DiffExp")

# First attempt which fails at background correction step
#--------------------------------------------------------------------------------#

# Load in data, including detection values

# x <- read.ilmn(files = "SampleProbeProfile.txt", ctrlfiles="ControlProbeProfile.txt",
#                other.columns = "Detection")
# 
# # Don't have "targets file" mentioned in vignettes
# # Presume this would have gender/age details etc?
# 
# # Description of probes
# table(x$genes$Status)
# 
# # Quick look at expression values for each probe
# options(digits=3) 
# head(x$E)
# boxplot(log2(x$E),range=0,ylab="log2 intensity")
# 
# # Properly expressed probes
# pe <- propexpr(x)
# pe # RAMS_BL and RAMS_fup: Are these before/after treatment??
# dim(pe) <- c(2,90)
# dimnames(pe) <- list(Treatment=c("Before MTX", "After MTX"),Patient=c(1:90))
# pe
# 
# min(pe) #0   Faulty probes???
# mean(pe) # 0.45
# max(pe) #0.535
# 
# 
# dim(x) # 47423, 180
# 
# # Background correction/Normalisation
# y <- neqc(x)
# # Returns following error:
# # Error in normexp.signal(normexp.par[i, ], x$E[, i]) : 
# #  sigma must be positive
# 
# # Searching this error only brings up the source code!!!

#----------------------------------------------------------------------------#

# Checking to see if this error is caused by probes with zero values
# RAMS23003_fup, RAMS04013_BL, RAMS28006_fup have 0 properly expressed probes
# Removed these to make new files (in python)
# Also removed RAMS19013, contained in samples, but not in targets
# Now have 86 patients in all files

# 
#  x <- read.ilmn(files = "SampleEdit.txt", ctrlfiles="ControlEdit.txt",
#                 other.columns = "Detection")
#  
#  targets = readTargets('TargetsEdit.txt')

# Read in Sorted datasets 
  x <- read.ilmn(files = "SampleSorted.txt", ctrlfiles="ControlSorted.txt",
                 other.columns = "Detection")
  
  targets = readTargets('TargetsEdit.txt')

# Don't have "targets file" mentioned in vignettes
# Presume this would have gender/age details etc?

# Description of probes
#table(x$genes$Status)

# Quick look at expression values for each probe
#options(digits=3)
#head(x$E)
#boxplot(log2(x$E),range=0,ylab="log2 intensity")

# # Properly expressed probes
# pe <- propexpr(x)
# pe # 
# dim(pe) <- c(2,86)
# dimnames(pe) <- list(Treatment=c("Before MTX", "After MTX"),Patient=c(1:86))
# pe
# 
# min(pe) #0.232 #Seems low given mean and max values
# mean(pe) # 0.458
# max(pe) #0.535
# 
# 
# dim(x) # 47423, 172

# Background correction/Normalisation
 y <- neqc(x)
# This now works!!
 
# Dimensions
dim(y) #46548, 172

# Filter unexpressed probes
expressed <- rowSums(y$other$Detection < 0.05) >= 3
y <- y[expressed,] 
dim(y)  #35786, 172

# Need targets data to carry on with analysis

#plotMDS(y, label = targets$time_point)
#plotMDS(y, label = targets$gender) # Gender 2 = Female, due to number of cases


# Correlation between donors (before/after treatment)
ct <- factor(targets$time_point)
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
dupcor <- duplicateCorrelation(y, design, block=targets$sample_name)
dupcor$consensus.correlation  # 0.219 

# Differential expression between treatment pre/post
fit <- lmFit(y,design,block=targets$sample_name, 
             correlation=dupcor$consensus.correlation)

contrasts <- makeContrasts("Baseline-FourWeeks", levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, method="global"))
results <- decideTests(fit2, method = 'global')
# Baseline-FourWeeks
# Down                    9
# NotSig              35751
# Up                     26

DiffExpTreatment <- topTable(fit2, number = 35, coef=1)

TreatmentDown <-DiffExpTreatment[(DiffExpTreatment$logFC < 0),]
TreatmentDown[,"SYMBOL"]
TreatmentUp <-DiffExpTreatment[(DiffExpTreatment$logFC > 0),]
TreatmentUp[,"SYMBOL"]

