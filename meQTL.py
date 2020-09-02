#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 11:27:19 2020

@author: sw1906
"""

import pandas as pd

TwinRA = pd.read_csv('RP2DNAmTwin_DMPsGender_RA_Only.csv',header = 0)
Mval_Fem = pd.read_csv('RP2DNAmInd_Female_Mval.csv', header = 0)
Mval_Mal = pd.read_csv('RP2DNAmInd_Male_Mval.csv', header = 0)
SNPs = pd.read_csv('meQTL/meQTL_SNPs.txt', sep = ' ', header = 0)
Males = pd.read_csv('SampleSheetMales.csv', header=0, skiprows = 7)
Females = pd.read_csv('SampleSheetFemales.csv', header=0, skiprows = 7)

# list DNAm patients
MalesDNAm = list(Males['sm']) 
FemalesDNAm = list(Females['sm'])

# Only inlcude DNAm patients in SNP data
SNPs_Mal = SNPs.loc[SNPs['FID'].isin(MalesDNAm)]
SNPs_Fem = SNPs.loc[SNPs['FID'].isin(FemalesDNAm)]


# list of cpg ranges of interest 
lol = [['chr10', 63800880, 63809170], 
       ['chr6', 33173278, 33173501],
       ['chr12', 54441450, 54446576],
       ['chr6', 32729358, 32729808],
       ['chr2', 74669349, 74669516],
       ['chr6', 18122965, 18123232],
       ['chr16', 62067937, 62071193]]
       

CpGlistname = []
df1 = pd.DataFrame() # Blank df for CpG dataframe
for i in lol:
    df = TwinRA
    df = df[df['chr'] == i[0]]
    df = df[df['pos'] >= i[1]]
    df = df[df['pos'] <= i[2]]
    df1 = df1.append(df) # Add Data to CpG dataframe
    CpGlistname += list(df['Name'])
        
# filter out irrelevant CpG values

Mval_Fem = Mval_Fem.loc[Mval_Fem['Unnamed: 0'].isin(CpGlistname)]
Mval_Mal = Mval_Mal.loc[Mval_Mal['Unnamed: 0'].isin(CpGlistname)]

# Set CpGs as index (as was initally the case)
Mval_Fem = Mval_Fem.set_index(['Unnamed: 0'])
Mval_Mal = Mval_Mal.set_index(['Unnamed: 0'])

#Alter columns to pe patient name
Mval_Mal.columns = ['BIOP006001', 'BIOP006021', 'BIOP006035', 'BIOP007012', 'BIOP007070', 'BIOP009002', 'BIOP009011', 'BIOP009020', 'BIOP009036', 'BIOP010013', 'BIOP011001', 'BIOP012014', 'BIOP020040', 'BIOP025006', 'BIOP025007']
Mval_Fem.columns = ['BIOP002015', 'BIOP003004', 'BIOP003009', 'BIOP003010', 'BIOP003014', 'BIOP003015', 'BIOP003039', 'BIOP004008', 'BIOP004011', 'BIOP004015', 'BIOP004017', 'BIOP005012', 'BIOP005017', 'BIOP005027', 'BIOP006005', 'BIOP006015', 'BIOP006016', 'BIOP006028', 'BIOP007014', 'BIOP007046', 'BIOP007069', 'BIOP007072', 'BIOP007082', 'BIOP008006', 'BIOP009008', 'BIOP009009', 'BIOP009010', 'BIOP009018', 'BIOP009019', 'BIOP009032', 'BIOP010003', 'BIOP010033', 'BIOP010036', 'BIOP010048', 'BIOP010050', 'BIOP010051', 'BIOP010053', 'BIOP011003', 'BIOP011004', 'BIOP017003', 'BIOP017007', 'BIOP017009', 'BIOP017010', 'BIOP018002', 'BIOP018013', 'BIOP018015', 'BIOP020008', 'BIOP022006', 'BIOP024018', 'BIOP026007', 'BIOP028005', 'BIOP033018', 'BIOP033039', 'BIOP033042', 'BIOP037007', 'BIOP039006', 'BIOP042001']

# List SNP patients 
MalesSNP = list(SNPs_Mal['FID'])
FemalesSNP = list(SNPs_Fem['FID'])

# Only inclide SNP patients in DNAm data
Mval_Mal = Mval_Mal[MalesSNP]
Mval_Fem = Mval_Fem[FemalesSNP]


# Reorganise SNP data to fit needs of Matrix e-QTL

SNPs_Mal = SNPs_Mal.drop(['IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1)
SNPs_Mal = SNPs_Mal.set_index(['FID'])
SNPs_Mal = SNPs_Mal.T

SNPs_Fem = SNPs_Fem.drop(['IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1)
SNPs_Fem = SNPs_Fem.set_index(['FID'])
SNPs_Fem = SNPs_Fem.T

# Make a CpG location df
CpGsLoc = df1[['Name', 'chr', 'pos']]
CpGsLoc['s2'] = CpGsLoc['pos']
CpGsLoc.rename(columns = {'pos':'s1'}, inplace = True) 

# Change indexes to columns
SNPs_Fem = SNPs_Fem.reset_index()
SNPs_Fem.rename(columns = {'index':'snp'}, inplace = True)

SNPs_Mal = SNPs_Mal.reset_index()
SNPs_Mal.rename(columns = {'index':'snp'}, inplace= True)

Mval_Fem = Mval_Fem.reset_index()
Mval_Fem.rename(columns = {'Unnamed: 0':'CpG'}, inplace = True)

Mval_Mal = Mval_Mal.reset_index()
Mval_Mal.rename(columns = {'Unnamed: 0':'CpG'}, inplace = True)


# Output these files for Matrix e-QTL work

#CpGsLoc.to_csv('meQTL_CPGsloc.csv', sep = '\t', index= False)

#SNPs_Fem.to_csv('meQTL_FemSNPs.csv', sep = '\t', index= False)
#SNPs_Mal.to_csv('meQTL_MalSNPs.csv', sep = '\t', index= False)

#Mval_Fem.to_csv('meQTL_FemCpGsMval.csv', sep = '\t', index = False)
#Mval_Mal.to_csv('meQTL_MalCpGsMval.csv', sep = '\t', index = False)

# Edit the .txt SNPlocs to .csv
#SNPlocsfile = pd.read_csv('meQTL_SNPLocs.txt', sep = ' ', header = None)
#SNPlocsfile.rename(columns = {0:'snp', 1:'chr', 2:'pos'}, inplace = True)
#SNPlocsfile.to_csv('meQTL_SNPsLoc.csv', sep = '\t', index = False)

# Sort all data by SNP/CpG name
SNPs_Mal = pd.read_csv('meQTL/meQTL_MalSNPs.csv', header= 0, sep = '\t')
SNPs_Fem = pd.read_csv('meQTL/meQTL_FemSNPs.csv', header= 0, sep = '\t')
SNPs_locs = pd.read_csv('meQTL/meQTL_SNPsLoc.csv', header= 0, sep = '\t')

Mval_Mal = pd.read_csv('meQTL/meQTL_MalCpGsMval.csv', header= 0, sep = '\t')
Mval_Fem = pd.read_csv('meQTL/meQTL_FemCpGsMval.csv', header= 0, sep = '\t')
CpGsloc = pd.read_csv('meQTL/meQTL_CPGsloc.csv', header= 0, sep = '\t')


SNPs_locs = SNPs_locs.drop_duplicates('snp', keep='first')

# Sort files to have same order
SNPs_Mal = SNPs_Mal.sort_values(by=['snp'])
SNPs_Fem = SNPs_Fem.sort_values(by=['snp'])
SNPs_locs = SNPs_locs.sort_values(by=['snp'])

SNPs_locs.insert(0, 'snps', range(1, 1 + len(SNPs_locs)))
SNPs_locs = SNPs_locs[['snps', 'chr', 'pos']]

SNPs_Mal.insert(0, 'snps', range(1, 1 + len(SNPs_Mal)))
SNPs_Mal = SNPs_Mal.drop('snp', axis = 1)

SNPs_Fem.insert(0, 'snps', range(1, 1 + len(SNPs_Fem)))
SNPs_Fem = SNPs_Fem.drop('snp', axis = 1)


Mval_Fem = Mval_Fem.sort_values(by=['CpG'])
Mval_Mal = Mval_Mal.sort_values(by=['CpG'])
CpGsloc = CpGsloc.sort_values(by=['Name'])



# Output these files for Matrix e-QTL work
'''
CpGsloc.to_csv('meQTL_CPGsloc.csv', sep = '\t', index= False)

SNPs_Fem.to_csv('meQTL_FemSNPs.csv', sep = '\t', index= False)
SNPs_Mal.to_csv('meQTL_MalSNPs.csv', sep = '\t', index= False)

Mval_Fem.to_csv('meQTL_FemCpGsMval.csv', sep = '\t', index = False)
Mval_Mal.to_csv('meQTL_MalCpGsMval.csv', sep = '\t', index = False)


SNPs_locs.to_csv('meQTL_SNPsLoc.csv', sep = '\t', index = False)

'''

#----------------------------------------------------------------------------#
# Analysis of meQTL results

import pandas as pd

# Reference for SNP number to SNP name
SNPlocsfile = pd.read_csv('meQTL/meQTL_SNPLocs.txt', sep = ' ', header = None)
SNPlocsfile.rename(columns = {0:'snp', 1:'chr', 2:'pos'}, inplace = True)
SNPlocsfile = SNPlocsfile.sort_values(by=['snp'])
SNPlocsfile = SNPlocsfile.drop_duplicates('snp', keep='first')
SNPlocsfile.insert(0, 'snpNumber', range(1, 1 + len(SNPlocsfile)))


CpGs_loc = pd.read_csv('meQTL/meQTL_CPGsloc.csv', header = 0, sep = '\t')


meQTL_Mal = pd.read_csv('meQTL/meQTL_CisMal_All.csv', sep = ',', header=0,
                        index_col =0)


meQTL_Fem = pd.read_csv('meQTL/meQTL_CisFem.csv', sep = ',', header=0,
                        index_col =0)


cols= ['chr', 'pos','strand', 'Name', 'UCSC_RefGene_Name', 'logFC', 'Relation_to_Island',
       'Enhancer', 'Regulatory_Feature_Group', 'P.Value', 'adj.P.Val' ]
TwinRA = pd.read_csv('RP2DNAmTwin_DMPsGender_RA_Only.csv',header = 0)
TwinRA = TwinRA[cols]

lol = [['chr10', 63800880, 63809170], 
       ['chr6', 33173278, 33173501],
       ['chr12', 54441450, 54446576],
       ['chr6', 32729358, 32729808],
       ['chr2', 74669349, 74669516],
       ['chr6', 18122965, 18123232],
       ['chr16', 62067937, 62071193]]

CpGs = pd.DataFrame() # Blank df for CpG dataframe
for i in lol:
    df = TwinRA
    df = df[df['chr'] == i[0]]
    df = df[df['pos'] >= i[1]]
    df = df[df['pos'] <= i[2]]
    CpGs = CpGs.append(df) # Add Data to CpG dataframe

# Use CpGs df and CpGs_loc df to identify interest SNPs/CpGs

# 10 sig, all linked to one CpG (HOXC4 island)
Fem_Sig = meQTL_Fem[meQTL_Fem['FDR'] < 0.05 ]
Fem_Sig_SNPs = list(Fem_Sig['snps'])
Fem_Sig_SNPs_df = SNPlocsfile.loc[SNPlocsfile['snpNumber'].isin(Fem_Sig_SNPs)]


# Looking at these locations in males
Mal_SNPs = meQTL_Mal.loc[meQTL_Mal['snps'].isin(Fem_Sig_SNPs)]
Mal_SNPs = Mal_SNPs.loc[Mal_SNPs['gene'] == 	'cg08465346']














