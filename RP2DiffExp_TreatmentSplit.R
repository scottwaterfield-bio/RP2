library(limma)
library(data.table)

xpre <- read.ilmn(files = "SamplePre.txt", ctrlfiles="ControlPre.txt",
                  other.columns = "Detection")

targetspre = readTargets('TargetsPre.txt')

xpost <- read.ilmn(files = "SamplePost.txt", ctrlfiles="ControlPost.txt",
                  other.columns = "Detection")

targetspost = readTargets('TargetsPost.txt')


# Background correction/Norpreisation
ypre <- neqc(xpre)
ypost <- neqc(xpost)


# Filter unexpressed probes
expressedpre <- rowSums(ypre$other$Detection < 0.05) >= 3
ypre <- ypre[expressedpre,] 

expressedpost <- rowSums(ypost$other$Detection < 0.05) >= 3
ypost <- ypost[expressedpost,] 


# Model fitting

# Gene difference between gender in pretreatment
ctpre <- factor(targetspre$gender)
designpre <- model.matrix(~0+ctpre)
fitpre <- lmFit(ypre, designpre)

contrastspre <- makeContrasts("ctpreMale-ctpreFemale", levels=designpre)
fit2pre <- contrasts.fit(fitpre, contrastspre)
fit2pre <- eBayes(fit2pre, trend=TRUE)
summary(decideTests(fit2pre, method="global"))
resultspre <- decideTests(fit2pre, method = 'global')

# ctpreMale-ctpreFemale              
# Down                       9
# NotSig                 30907
# Up                        16

DiffExpTreatmentPre <- topTable(fit2pre, number = 50000, coef=1)
#DiffExpTreatmentPre[,"SYMBOL"]
DiffExpTreatmentPre[(DiffExpTreatmentPre$SYMBOL == 'TNXB'),]
#fwrite(DiffExpTreatmentPre, 'RP2PreMTXDiffExp.csv')

# [1] "RPS4Y1"       "LOC100133662" "EIF1AY"       "RPS4Y2"       "EIF1AY"       "JARID1D"      "CYorf15A"     "PRKY"         "XIST"         "UTY"         
# [11] "ZFY"          "CYorf15B"     "TMSB4Y"       "USP9Y"        "CYorf15B"     "RPS4X"        "RPS4X"        "LOC391777"    "UTX"          "PRKX"        
# [21] "MED8"         "ZFX"          "SEPT6"        "PALLD"        "EIF1AX" 

PreTreatmentDown <-DiffExpTreatmentPre[(DiffExpTreatmentPre$logFC < 0),]
#PreTreatmentDown[,"SYMBOL"]

PreTreatmentUp <-DiffExpTreatmentPre[(DiffExpTreatmentPre$logFC > 0),]
#PreTreatmentUp[,"SYMBOL"]

# Gene difference between gender in posttreatment
ctpost <- factor(targetspost$gender)
designpost <- model.matrix(~0+ctpost)
fitpost <- lmFit(ypost, designpost)

contrastspost <- makeContrasts("ctpostMale-ctpostFemale", levels=designpost)
fit2post <- contrasts.fit(fitpost, contrastspost)
fit2post <- eBayes(fit2post, trend=TRUE)
summary(decideTests(fit2post, method="global"))
resultspost <- decideTests(fit2post, method = 'global')

# ctpostMale-ctpostFemale
# Down                         7
# NotSig                   31253
# Up                          15


DiffExpTreatmentpost <- topTable(fit2post, number = 50000, coef=1)
DiffExpTreatmentpost[(DiffExpTreatmentpost$SYMBOL == 'TNXB'),]
#fwrite(DiffExpTreatmentpost, 'RP2PostMTXDiffExp.csv')

# [1] "LOC100133662" "RPS4Y1"       "EIF1AY"       "CYorf15A"     "RPS4Y2"       "EIF1AY"       "JARID1D"      "PRKY"         "UTY"          "XIST"        
# [11] "ZFY"          "CYorf15B"     "USP9Y"        "TMSB4Y"       "CYorf15B"     "UTX"          "RPS4X"        "PRKX"         "ZFX"          "DDIT4"       
# [21] "LOC391777"    "RPS4X" 


PostTreatmentDown <-DiffExpTreatmentpost[(DiffExpTreatmentpost$logFC < 0),]
#PostTreatmentDown[,"SYMBOL"]

PostTreatmentUp <-DiffExpTreatmentpost[(DiffExpTreatmentpost$logFC > 0),]
#PostTreatmentUp[,"SYMBOL"]
