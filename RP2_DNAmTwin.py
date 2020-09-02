#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 16:21:22 2020

@author: sw1906
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
cwd = os.getcwd()

'''
targets  = pd.read_csv("SampleSheet.csv", header = 7)


# Remove controls and nontwins
targets = targets[~targets.Sample_Group.str.contains("control")]
targets = targets[~targets.Sample_Group.str.contains("nontwin")]

#remove repeats 
targets = targets[~targets.Sample_Group.str.contains("repeat")]


# Set All RA and healthy patients to have same value names

targets = targets.replace(['PBL_RA'],'RA')
targets = targets.replace(['PBL_healthy_cotwin',],'healthy_cotwin')

targets.to_csv('SampleSheetEdit.csv', sep = '\t', index= False)
'''

'''
# Edit targets files into relevant groups
targets  = pd.read_csv("SampleSheetEdit.csv", header = 7)


femtargets = targets[targets['gender'] == 0] 
maltargets = targets[targets['gender'] == 1] 

Healthy = targets[targets['Sample_Group'] == 'healthy_cotwin']
RA = targets[targets['Sample_Group'] == 'RA']


femtargets.to_csv('SampleSheetFem.csv', sep = '\t', index= False)
maltargets.to_csv('SampleSheetMal.csv', sep = '\t', index= False)

Healthy.to_csv('SampleSheetHealthy.csv', sep = '\t', index= False)
RA.to_csv('SampleSheetRA.csv', sep = '\t', index= False)
'''

'''
#----------------------------------------------------------------------------#
# Load in Differential DNAm results

cols= ['chr', 'pos', 'Name', 'UCSC_RefGene_Name', 'logFC',
                           'P.Value', 'adj.P.Val' ]

TwinAllData = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatus.csv',header = 0)
TwinAllData = TwinAllData[cols]

TwinFemData = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusFemales.csv',header = 0)
TwinFemData = TwinFemData[cols]

TwinMalData = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusMales.csv', header = 0)
TwinMalData = TwinMalData[cols]


SigAllData = TwinAllData[TwinAllData['adj.P.Val'] < 0.05] 
SigFemData = TwinFemData[TwinFemData['adj.P.Val'] < 0.05] 

#----------------------------------------------------------------------------#
# Looking at alldata CPGs in females and males

AllDataNames = list(SigAllData["Name"])

AllDataCpg_FemVals = TwinFemData[TwinFemData['Name'].isin(AllDataNames)]
AllDataCpg_MalVals = TwinMalData[TwinMalData['Name'].isin(AllDataNames)]

# Sort df's by Cg name to allow row by row comparison 
SigAllData = SigAllData.sort_values(by = 'Name')
AllDataCpg_FemVals = AllDataCpg_FemVals.sort_values(by = 'Name')
AllDataCpg_MalVals = AllDataCpg_MalVals.sort_values(by = 'Name')


# Compare logFC values in Sig All data CPGs to males and females
AllDatavsFemlogFC = list(SigAllData['logFC'] - list(AllDataCpg_FemVals['logFC']))
AllDatavsMallogFC = list(SigAllData['logFC'] - list(AllDataCpg_MalVals['logFC']))

print('Alldata Vs. Female data logFC differences: ', AllDatavsFemlogFC)
print('Alldata Vs. Male data logFC differences: ', AllDatavsMallogFC)

# cg09838169 (INPP5A) only cpg with logFC diff greater than 0.1 
# AllData logFC: 0.209, Male logFC: 0.412
# Potential that this is stronger in males?
# The males are actually pushing this into adjusted signifcance, 
# not sig in just females (logFC = 0.17, adj P = 0.18)

#----------------------------------------------------------------------------#
# Looking at significant female CPGs in males (and all data)

FemDataNames = list(SigFemData["Name"])

FemCpg_AllDataVals = TwinAllData[TwinAllData['Name'].isin(FemDataNames)]
FemCpg_MalVals = TwinMalData[TwinMalData['Name'].isin(FemDataNames)]


# Sort df's by Cg name to allow row by row comparison 
SigFemData = SigFemData.sort_values(by = 'Name')
FemCpg_AllVals = FemCpg_AllDataVals.sort_values(by = 'Name')
FemCpg_MalVals = FemCpg_MalVals.sort_values(by = 'Name')


# Compare logFC values in Sig All data CPGs to males and females
FemDatavsAlllogFC = list(SigFemData['logFC'] - list(FemCpg_AllVals['logFC']))
FemDatavsMallogFC = list(SigFemData['logFC'] - list(FemCpg_MalVals['logFC']))

print('Femdata Vs.  Alldata logFC differences: ', FemDatavsAlllogFC)
print('Femdata Vs. Male data logFC differences: ', FemDatavsMallogFC)


# Cpgs with greater than 0.1 logFC change in males vs females

# cg07090025 (SSU72) isn't differentially methylated in males
# Males logFC: -0.0045, Females logFC: 0.20

# cg20045320 (No gene) isn't differentially methylated in males either
# Males logFC: -0.049, Females logFC: -0.156

#----------------------------------------------------------------------------#

# Looking into the X chromosome

FemX = TwinFemData[TwinFemData['chr'] == 'chrX'] # Lowest Adj P Values ~ 0.2
MalX = TwinMalData[TwinMalData['chr'] == 'chrX'] # Lowest Adj P values ~ 0.8


# Looking into the X chromosome psuedoautosomal gene regions (PARs)
PAR1Start = 10001 	
PAR1Stop = 2781479
PAR2Start = 155701383 	
PAR2Stop = 156030895

FemPAR1 =  FemX[FemX['pos'].between(PAR1Start, PAR1Stop)]
FemPAR2 = FemX[FemX['pos'].between(PAR2Start, PAR2Stop )]

# Only CPGs identified here are PAR1: GYG2 and XG 
# 7 GYG2 sites, lowest adj P = 0.3
# No differences in 5 XG sites either  
'''

#----------------------------------------------------------------------------#

# Looking into the Gender DNAm patterns in healthy people and RA patients


cols= ['chr', 'pos','strand', 'Name', 'UCSC_RefGene_Name', 'logFC', 'Relation_to_Island',
       'Enhancer', 'Regulatory_Feature_Group', 'P.Value', 'adj.P.Val' ]


Twin = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatus.csv',header = 0)
Twin = Twin[cols]
TwinFem = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusFemales.csv',header = 0)
TwinFem = TwinFem[cols]
TwinMal = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusMales.csv',header = 0)
TwinMal = TwinMal[cols]

TwinHealthy = pd.read_csv('RP2DNAmTwin_DMPsGender_Healthy_Only.csv',header = 0)
TwinHealthy = TwinHealthy[cols]

TwinRA = pd.read_csv('RP2DNAmTwin_DMPsGender_RA_Only.csv',header = 0)
TwinRA = TwinRA[cols]

#IndRA = pd.read_csv('RP2DNAmInd_Gender.csv')
#IndRA = IndRA[cols]

#----------------------------------------------------------------------------#
'''
# Below contains all IndDNA data analysis
# Comment out/Save in individual script afer use
SigRA = TwinRA[TwinRA['adj.P.Val'] < 0.05] # 15340
SigIndRA = IndRA[IndRA['adj.P.Val'] < 0.05] # 8782
SigBoth = list(set(SigRA['Name']) & set(SigIndRA['Name'])) #7995

SigBothTwin = TwinRA.loc[TwinRA['Name'].isin(SigBoth)]
SigBothInd = IndRA.loc[IndRA['Name'].isin(SigBoth)]

Sexchr = ['chrY', 'chrX']
SigBothTwinAuto =  SigBothTwin[~SigBothTwin['chr'].isin(Sexchr)] # 379/7995
SigBothIndAuto =  SigBothInd[~SigBothInd['chr'].isin(Sexchr)]


# Identify the ~800 insignificant CpGs in TwinRA
SiginTwin = list(set(SigRA['Name']).difference(SigIndRA['Name'])) #7345
SiginInd = list(set(SigIndRA['Name']).difference(SigRA['Name'])) #787
# Look to see if these are merely close to non significant in Twin/vice versa
TwindataSigInd = TwinRA.loc[TwinRA['Name'].isin(SiginInd)] #773, 14 not in data
AlmostSig = TwindataSigInd[TwindataSigInd['adj.P.Val'] <0.15] #285 results

InddataSigRA = IndRA.loc[IndRA['Name'].isin(SiginTwin)] #6817, ~500 not in data
AlmostSigInd = InddataSigRA[InddataSigRA['adj.P.Val'] <0.15] # Only 285
NotSigInd = InddataSigRA[InddataSigRA['adj.P.Val'] > 0.5] # 5702
'''
#----------------------------------------------------------------------------#





SigHealthy = TwinHealthy[TwinHealthy['adj.P.Val'] < 0.05] #11636
SigRA = TwinRA[TwinRA['adj.P.Val'] < 0.05] # 15340

# LogFC values are for men (compared to women)
OverMethMenHealthy = SigHealthy[SigHealthy['logFC'] > 0 ] # 3290
UnderMethMenHealthy = SigHealthy[SigHealthy['logFC'] < 0 ] # 8345

OverMethMenRA = SigRA[SigRA['logFC'] > 0 ] # 4594
UnderMethMenRA = SigRA[SigRA['logFC'] < 0 ] # 10746


# Looking for Gender similarities and differences in healthy and RA patients

# CpGs that are over/under methylated in both healthy and unhealthy patients
OverBoth = list(set(OverMethMenHealthy['Name']) & set(OverMethMenRA['Name'])) #2581
UnderBoth = list(set(UnderMethMenHealthy['Name']) & set(UnderMethMenRA['Name'])) # 6548


# Cpgs that are only significant in either healthy or unhealthy
OverinHealthy = list(set(OverMethMenHealthy['Name']).difference(SigRA['Name'])) #709
UnderinHealthy = list(set(UnderMethMenHealthy['Name']).difference(SigRA['Name'])) #1797
OverinRA = list(set(OverMethMenRA['Name']).difference(SigHealthy['Name'])) #2013
UnderinRA = list(set(UnderMethMenRA['Name']).difference(SigHealthy['Name'])) #4198


# No CpGs undermethylated in healthy/RA and the reverse in RA/Healthy
OverRAUnderHealthy = list(set(UnderMethMenHealthy['Name']) & set(OverMethMenRA['Name']))
UnderRAOverHealthy = list(set(OverMethMenHealthy['Name']) & set(UnderMethMenRA['Name']))

#----------------------------------------------------------------------------#
# Autosomal/Sex CpGs in these lists
Sexchr = ['chrX', 'chrY']

# CpGs common to both Healthy and RA
OverBothHealthy = SigHealthy.loc[SigHealthy['Name'].isin(OverBoth)]
AutoOverBothHealthy = OverBothHealthy[~OverBothHealthy['chr'].isin(Sexchr)] #240
SexOverBothHealthy = OverBothHealthy[OverBothHealthy['chr'].isin(Sexchr)] #2341


OverBothRA = SigRA.loc[SigRA['Name'].isin(OverBoth)]
AutoOverBothRA = OverBothRA[~OverBothRA['chr'].isin(Sexchr)] #240
SexOverBothRA = OverBothRA[OverBothRA['chr'].isin(Sexchr)] #2341


UnderBothHealthy = SigHealthy.loc[SigHealthy['Name'].isin(UnderBoth)]
AutoUnderBothHealthy = UnderBothHealthy[~UnderBothHealthy['chr'].isin(Sexchr)] #616
SexUnderBothHealthy = UnderBothHealthy[UnderBothHealthy['chr'].isin(Sexchr)] #5932


UnderBothRA = SigRA.loc[SigRA['Name'].isin(UnderBoth)]
AutoUnderBothRA = UnderBothRA[~UnderBothRA['chr'].isin(Sexchr)] #616
SexUnderBothRA = UnderBothRA[UnderBothRA['chr'].isin(Sexchr)] #5932


#----------------------------------------------------------------------------#

#CpGs seen to be significant in only one of the two states
OverHealthy = SigHealthy.loc[SigHealthy['Name'].isin(OverinHealthy)]
AutoOverHealthy = OverHealthy[~OverHealthy['chr'].isin(Sexchr)] #392
SexOverHealthy = OverHealthy[OverHealthy['chr'].isin(Sexchr)] # 317

OverRA = SigRA.loc[SigRA['Name'].isin(OverinRA)]
AutoOverRA = OverRA[~OverRA['chr'].isin(Sexchr)] #1750
SexOverRA = OverRA[OverRA['chr'].isin(Sexchr)] # 263

UnderHealthy = SigHealthy.loc[SigHealthy['Name'].isin(UnderinHealthy)]
AutoUnderHealthy = UnderHealthy[~UnderHealthy['chr'].isin(Sexchr)] #1693
SexUnderHealthy = UnderHealthy[UnderHealthy['chr'].isin(Sexchr)] #104

UnderRA = SigRA.loc[SigRA['Name'].isin(UnderinRA)]
AutoUnderRA = UnderRA[~UnderRA['chr'].isin(Sexchr)] #3834
SexUnderRA = UnderRA[UnderRA['chr'].isin(Sexchr)] #364


# Evaluating if these are truly only significant in one State (healthy/RA)
# Choose a sensible Adj. P to show 'Unsignifcant' results
# A (adjusted) coin flip?
UnsigHealthy = TwinHealthy[TwinHealthy['adj.P.Val'] > 0.5] #230957
UnsigRA = TwinRA[TwinRA['adj.P.Val'] > 0.5] # 217,693

# looking for CpGs which have true RA/disease gender differences 

SigAutoOverHealthy = pd.merge(AutoOverHealthy, UnsigRA, how='inner', on=['Name']) #172, was 392
SigSexOverHealthy = pd.merge(SexOverHealthy, UnsigRA, how='inner', on=['Name']) # 43, was 317

SigAutoUnderHealthy = pd.merge(AutoUnderHealthy, UnsigRA, how='inner', on=['Name']) # 387, was 1693
SigSexUnderHealthy = pd.merge(SexUnderHealthy, UnsigRA, how='inner', on=['Name']) # 3, was 104

SigAutoOverRA = pd.merge(AutoOverRA, UnsigHealthy, how='inner', on=['Name']) # 1241, was 1750
SigSexOverRA = pd.merge(SexOverRA, UnsigHealthy, how='inner', on=['Name']) # 57, was 263

SigAutoUnderRA = pd.merge(AutoUnderRA, UnsigHealthy, how='inner', on=['Name']) # 2038, was 3834
SigSexUnderRA = pd.merge(SexUnderRA, UnsigHealthy, how='inner', on=['Name']) # 93, was 364

# Analysing TrueSig Autosomal RA datasets
SigAutoOverRADupes = SigAutoOverRA[SigAutoOverRA.duplicated(subset=['UCSC_RefGene_Name_x'], keep=False)]
SigAutoOverRADupes = SigAutoOverRADupes[SigAutoOverRADupes['UCSC_RefGene_Name_x'].notna()]
len(SigAutoOverRADupes['UCSC_RefGene_Name_x'].unique().tolist()) #50 genes with at least two CpGs



SigAutoUnderRADupes = SigAutoUnderRA[SigAutoUnderRA.duplicated(subset=['UCSC_RefGene_Name_x'], keep=False)]
SigAutoUnderRADupes = SigAutoUnderRADupes[SigAutoUnderRADupes['UCSC_RefGene_Name_x'].notna()]
len(SigAutoUnderRADupes['UCSC_RefGene_Name_x'].unique().tolist()) #128 genes
SigAutoUnderRADupes = SigAutoUnderRADupes[SigAutoUnderRADupes['Regulatory_Feature_Group_x'] == 'Promoter_Associated']



# Looking into identified genes of interest

# add disease status identifier to datasets
TwinHealthy['Disease'] = 'Healthy'
TwinRA['Disease'] = 'RA'


# HOXC4 analysis
HoxC4 = ['ARID5B' ]
HOXHealthy = TwinHealthy[TwinHealthy['UCSC_RefGene_Name'].isin(HoxC4)]
HOXRA = TwinRA[TwinRA['UCSC_RefGene_Name'].isin(HoxC4)]
HoxBoth = pd.concat([HOXRA, HOXHealthy], ignore_index=True)

sns.scatterplot(x = 'pos', y= 'adj.P.Val' , data = HoxBoth, hue = 'Disease')
ax = plt.gca()
ax.set_title("HOXC4 CpG site Differential Methylation between Males and Females (Adj. P value)")


sns.scatterplot(x = 'pos', y= 'logFC' , data = HoxBoth, hue = 'Disease')
ax = plt.gca()
ax.set_title("HOXC4 CpG site Differential Methylation between Males and Females (logFC) ")

# sns.lmplot(x = 'pos', y = 'adj.P.Val', data=HoxBoth, hue = 'Disease') # Nice graph


# Look into RA-related sex chromosomes differences first
HOXC4 = ['PRDM16;PRDM16']
HOXHealthy = TwinHealthy[TwinHealthy['UCSC_RefGene_Name'].isin(HOXC4)]
HOXRA = TwinRA[TwinRA['UCSC_RefGene_Name'].isin(HOXC4)]
HoxBoth = pd.concat([ HOXHealthy,HOXRA], ignore_index=True)

sns.scatterplot(x = 'pos', y= 'adj.P.Val' , data = HoxBoth, hue = 'Disease')
ax = plt.gca()
ax.set_title("HOXC4 CpG site Differential Methylation between Males and Females (Adj. P value)")
plt.clf()

sns.scatterplot(x = 'pos', y= 'logFC' , data = HoxBoth, hue = 'Disease')
ax = plt.gca()
ax.set_title("HOXC4 CpG site Differential Methylation between Males and Females (Adj. P value)")

# Plot function
def plotter(Healthydf, RAdf, genename): 
    
    
    HealthySubset = Healthydf[Healthydf['UCSC_RefGene_Name'] == genename]
    RASubset = RAdf[RAdf['UCSC_RefGene_Name'] == genename]
    Combineddf = pd.concat([HealthySubset, RASubset])
    
    sns.scatterplot(x = 'pos', y= 'adj.P.Val', data = Combineddf, hue = 'Disease')
    ax = plt.gca()
    ax.set_title(genename + " CpG site Differential Methylation between Males and Females (Adj. P value)")
    plt.savefig(genename + '_PVal_Plot')
    plt.clf()
    
    sns.scatterplot(x = 'pos', y= 'logFC', data = Combineddf, hue = 'Disease')
    ax = plt.gca()
    ax.set_title(genename + " CpG site Differential Methylation between Males and Females (logFC)")
    plt.savefig(genename + '_logFC_Plot')
    plt.clf()
    
    
# Get to plotting 
listofgenes = ['CDH8']
for i in listofgenes:
    plotter(TwinHealthy, TwinRA, i)

SanityCheck = ['PRDM16;PRDM16', 'TNXB', 'ARID5B', 'COL11A2;COL11A2;COL11A2']
SanityHealthy = TwinHealthy[TwinHealthy['UCSC_RefGene_Name'].isin(SanityCheck)]
SanityRA = TwinRA[TwinRA['UCSC_RefGene_Name'].isin(SanityCheck)]
#----------------------------------------------------------------------------#
# Searching for 'AID' in the dataset (HOXC4 pathway)
'''
AIDAlias = [ 'AICDA', 'AID', 'ARP2', 'CDA2', 'HEL-S-284', 'HIGM2']

AIDHealthy = TwinHealthy[TwinHealthy['UCSC_RefGene_Name'].isin(AIDAlias)]
AIDRA = TwinRA[TwinRA['UCSC_RefGene_Name'].isin(AIDAlias)]

'''

#----------------------------------------------------------------------------#
# Loooking at specific genes in the data 
GOIHealthy = TwinHealthy[TwinHealthy['UCSC_RefGene_Name'] == 'RTKN']
GOIRA = TwinRA[TwinRA['UCSC_RefGene_Name'] == 'RTKN']

#----------------------------------------------------------------------------#







