library(limma)

xmal <- read.ilmn(files = "MaleSample.txt", ctrlfiles="MaleControl.txt",
               other.columns = "Detection")

targetsmal = readTargets('TargetsMale.txt')

xfem <- read.ilmn(files = "FemaleSample.txt", ctrlfiles="FemaleControl.txt",
                  other.columns = "Detection")

targetsfem = readTargets('TargetsFemale.txt')


# Background correction/Normalisation
ymal <- neqc(xmal)
yfem <- neqc(xfem)


# Filter unexpressed probes
expressedmal <- rowSums(ymal$other$Detection < 0.05) >= 3
ymal <- ymal[expressedmal,] 

expressedfem <- rowSums(yfem$other$Detection < 0.05) >= 3
yfem <- yfem[expressedfem,] 


# Correlation between donors (before/after treatment)

# Males
ctmal <- factor(targetsmal$time_point)
designmal <- model.matrix(~0+ctmal)
colnames(designmal) <- levels(ctmal)
dupcormal <- duplicateCorrelation(ymal, designmal, block=targetsmal$sample_name)
dupcormal$consensus.correlation  # 0.319

# Females
ctfem <- factor(targetsfem$time_point)
designfem <- model.matrix(~0+ctfem)
colnames(designfem) <- levels(ctfem)
dupcorfem <- duplicateCorrelation(yfem, designfem, block=targetsfem$sample_name)
dupcorfem$consensus.correlation  # 0.226


# Differential expression between males and females

# Treatment differences in males
fitmal <- lmFit(ymal,designmal,block=targetsmal$sample_name, 
             correlation=dupcormal$consensus.correlation)

contrastsmal <- makeContrasts("Baseline-FourWeeks", levels=designmal)
fit2mal <- contrasts.fit(fitmal, contrastsmal)
fit2mal <- eBayes(fit2mal, trend=TRUE)
summary(decideTests(fit2mal, method="global"))
resultsmal <- decideTests(fit2mal, method = 'global')

# Baseline-FourWeeks
# Down                    0
# NotSig              26277
# Up                      0

DiffExpTreatmentMales <- topTable(fit2mal, number = 50000,  coef=1)
DiffExpTreatmentMales[(DiffExpTreatmentMales$SYMBOL == 'SSU72'),]

# NS due to Adj P val/Small sample size (20 males)... Lowest raw P value: 0.000005


# Treatment difference in females
fitfem <- lmFit(yfem,designfem,block=targetsfem$sample_name, 
                correlation=dupcorfem$consensus.correlation)

contrastsfem <- makeContrasts("Baseline-FourWeeks", levels=designfem)
fit2fem <- contrasts.fit(fitfem, contrastsfem)
fit2fem <- eBayes(fit2fem, trend=TRUE)
summary(decideTests(fit2fem, method="global"))
resultsfem <- decideTests(fit2fem, method = 'global')

# Baseline-FourWeeks
# Down                    3
# NotSig              34127
# Up                     10

DiffExpTreatmentFemales <- topTable(fit2fem, number = 50000, coef=1)
DiffExpTreatmentFemales[(DiffExpTreatmentFemales$SYMBOL == 'RHOA'),]

TreatmentDownFem <-DiffExpTreatmentFemales[(DiffExpTreatmentFemales$logFC < 0),]
TreatmentDownFem[,"SYMBOL"]
TreatmentUpFem <-DiffExpTreatmentFemales[(DiffExpTreatmentFemales$logFC > 0),]
TreatmentUpFem[,"SYMBOL"]

#The three occurences of doewnregulation are: RNU1-5, RNU1G2, RNU1-3. 










