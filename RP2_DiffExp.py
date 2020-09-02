
import pandas as pd

'''
# Load in data
ControlExp = pd.read_csv('ControlProbeProfile.txt',sep = '\t', header = 0)
SampleExp = pd.read_csv('SampleProbeProfile.txt',sep = '\t', header = 0)
targetsExp = pd.read_csv('Targets_Exp.txt', sep = '\t', header = 0)

# Break up Data into 
dropfiles = ['RAMS23003', 'RAMS04013', 'RAMS28006']  # zero functional probes
for i in dropfiles:
    Sampledrops = [col for col in SampleExp.columns if i in col]
    SampleExp = SampleExp.drop(Sampledrops, axis=1)

    Controldrops = [col for col in ControlExp.columns if i in col]
    ControlExp = ControlExp.drop(Controldrops, axis=1)


for i in dropfiles: 
    targetsExp = targetsExp[~targetsExp.sample_name.str.contains(i)]

# Sample and Control now have 87 patients, targets has 86
# Find the remaining patient to drop 

targetIDs = list(targetsExp["sample_name"])

# drop any columns which dont contain these values from control/sample
for i in targetIDs:
    Sampledrops = [col for col in SampleExp.columns if i in col]
    SampleExp = SampleExp.drop(Sampledrops, axis=1)

    Controldrops = [col for col in ControlExp.columns if i in col]
    ControlExp = ControlExp.drop(Controldrops, axis=1)
    
# This dropped all the matching columns (Obviously...) Oops. 
# Does identify RAMS19013 as remaining mismatch though!
'''

'''
# Break up Data into 
dropfiles = ['RAMS23003', 'RAMS04013', 'RAMS28006', 'RAMS19013']  # 3 zero functional probes, and 1 mismatch
for i in dropfiles:
    Sampledrops = [col for col in SampleExp.columns if i in col]
    SampleExp = SampleExp.drop(Sampledrops, axis=1)

    Controldrops = [col for col in ControlExp.columns if i in col]
    ControlExp = ControlExp.drop(Controldrops, axis=1)


for i in dropfiles: 
    targetsExp = targetsExp[~targetsExp.sample_name.str.contains(i)]
    


ControlExp.to_csv('ControlEdit.txt', sep = '\t', index=False)
SampleExp.to_csv('SampleEdit.txt', sep = '\t', index=False)
targetsExp.to_csv('TargetsEdit.txt', sep = '\t', index=False)


#----------------------------------------------------------------------------#
# Make new control/sample files with alphanumeric order columns

Control = pd.read_csv('ControlEdit.txt',sep = '\t', header = 0)
Sample = pd.read_csv('SampleEdit.txt',sep = '\t', header = 0)
targetsExp = pd.read_csv('TargetsEdit.txt', sep = '\t', header = 0)

Control = Control.reindex(sorted(Control.columns), axis=1)
Sample = Sample.reindex(sorted(Sample.columns), axis=1)

Control.to_csv('ControlSorted.txt', sep = '\t', index=False)
Sample.to_csv('SampleSorted.txt', sep = '\t', index=False)
targetsExp.to_csv('TargetsEdit.txt', sep = '\t', index=False)

#----------------------------------------------------------------------------#
# Remove _fup and _BL from targets sample_name, as this makes in donor comparision 
# difficult to carry out. 
# There is still a column displaying Treatment: before/after

Targets = pd.read_csv('TargetsEdit.txt', sep = '\t', header = 0)

Targets['sample_name'] = Targets['sample_name'].str.rstrip("_BL")
Targets['sample_name'] = Targets['sample_name'].str.rstrip("_fup")

Targets.to_csv('TargetsEdit.txt', sep = '\t', index=False)

#----------------------------------------------------------------------------#
# Change tartgets 4weeks to fourweeks (integer causing issues)

d = {'4weeks':'FourWeeks'}
Targets = Targets.replace(d)
Targets.to_csv('TargetsEdit.txt', sep = '\t', index=False)
'''


#----------------------------------------------------------------------------#
# Alter datasets to before drug/after drug, and male/female

Control = pd.read_csv('ControlSorted.txt', sep = '\t', header=0)
Sample = pd.read_csv('SampleSorted.txt', sep = '\t', header =0)
Targets = pd.read_csv('TargetsEdit.txt', sep = '\t', header = 0)


# Male targets file and female targets file
TargetsMale = Targets[Targets.gender == 1 ]
TargetsFemale = Targets[Targets.gender == 2 ]

# Now to Match the control and sample files to these four groups...

Males = list(TargetsMale['sample_name'])
Females = list(TargetsFemale['sample_name'])

FemaleSample = Sample
FemaleControl = Control

for i in Males:
    Sampledrops = [col for col in FemaleSample.columns if i in col]
    FemaleSample = FemaleSample.drop(Sampledrops, axis=1)

    Controldrops = [col for col in FemaleControl.columns if i in col]
    FemaleControl = FemaleControl.drop(Controldrops, axis=1)
    
    
MaleSample = Sample
MaleControl = Control  

for i in Females:
    Sampledrops = [col for col in MaleSample.columns if i in col]
    MaleSample = MaleSample.drop(Sampledrops, axis=1)

    Controldrops = [col for col in MaleControl.columns if i in col]
    MaleControl = MaleControl.drop(Controldrops, axis=1)



# Output this to files
TargetsMale.to_csv('TargetsMale.txt', sep = '\t', index=False)
TargetsFemale.to_csv('TargetsFemale.txt', sep = '\t', index=False)

FemaleSample.to_csv('FemaleSample.txt', sep = '\t', index=False)
FemaleControl.to_csv('FemaleControl.txt', sep = '\t', index=False)
      
MaleSample.to_csv('MaleSample.txt', sep = '\t', index=False)
MaleControl.to_csv('MaleControl.txt', sep = '\t', index=False) 



# Make the Before/After treatment treatment dataframes

# Before/after treatment targets files 
TargetsPre = Targets[Targets.time_point == 'Baseline' ]
TargetsPost = Targets[Targets.time_point == 'FourWeeks' ]


#Drop _fup columns from pre-MTX files
Sampledrops = [col for col in Sample.columns if '_fup' in col]
SamplePre = Sample.drop(Sampledrops, axis=1)

Controldrops = [col for col in Control.columns if '_fup' in col]
ControlPre = Control.drop(Controldrops, axis=1)


# Drop _BL columns from post-MTX files
Sampledrops = [col for col in Sample.columns if '_BL' in col]
SamplePost = Sample.drop(Sampledrops, axis=1)

Controldrops = [col for col in Control.columns if '_BL' in col]
ControlPost = Control.drop(Controldrops, axis=1)


# Change gender to male/female (currently 1/2)

d = {1:'Male', 2:'Female'}
TargetsPre = TargetsPre.replace(d)
TargetsPost = TargetsPost.replace(d)
# Export to csv

SamplePre.to_csv('SamplePre.txt', sep = '\t', index=False)
ControlPre.to_csv('ControlPre.txt', sep = '\t', index=False)
      
SamplePost.to_csv('SamplePost.txt', sep = '\t', index=False)
ControlPost.to_csv('ControlPost.txt', sep = '\t', index=False) 


TargetsPre.to_csv('TargetsPre.txt', sep = '\t', index= False)
TargetsPost.to_csv('TargetsPost.txt', sep = '\t', index= False)
















