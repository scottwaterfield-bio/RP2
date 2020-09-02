#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:51:20 2020

@author: sw1906
"""
import pandas as pd

SNPs = pd.read_csv('meQTL/HOXC4_SNPs.txt', sep = ' ', header = 0)
#SNPlocs = pd.read_csv('meQTL/HOXC4_SNPlocs.txt', sep = ' ', header = None)
#SNPlocs.rename(columns = {0:'snp', 1:'chr', 2:'pos'}, inplace = True)
#SNPlocs.to_csv('meQTL/eQTL_HOXC4SNPsLoc.csv', sep = '\t', index = False)


# Sort SNPs by gender
MalSNPs = SNPs.loc[SNPs['SEX'] == 1]
FemSNPs = SNPs.loc[SNPs['SEX'] == 2]


# make sure both groups are same patients
Exp_Pre = pd.read_csv('meQTL/PreExpression.csv', header = 0)
Exp_Pre = Exp_Pre.loc[Exp_Pre['probes.SYMBOL'] == 'HOXC4']
Exp_Pre =  Exp_Pre.set_index(['Unnamed: 0'])
Exp_Pre = Exp_Pre.iloc[:, :-3]
Exp_Pre.columns = Exp_Pre.columns.str.lstrip('exprs.')
Exp_Pre.columns = Exp_Pre.columns.str.rstrip('_BL')

MalesSNP = list(MalSNPs['FID'])
FemalesSNP = list(FemSNPs['FID'])
Pre_Males = Exp_Pre[Exp_Pre.columns.intersection(MalesSNP)]
Pre_Females = Exp_Pre[Exp_Pre.columns.intersection(FemalesSNP)]


FemaleExp = list(Pre_Females.columns)
MaleExp = list(Pre_Males.columns)
SNPs_Fem = FemSNPs.loc[FemSNPs['FID'].isin(FemaleExp)]
SNPs_Mal = MalSNPs.loc[MalSNPs['FID'].isin(MaleExp)]


# sort all data files
SNPs_locs = pd.read_csv('meQTL/eQTL_HOXC4SNPsLoc.csv', sep = '\t', header = 0)
SNPs_locs = SNPs_locs.drop_duplicates('snp', keep='first')

# Sort files to have same order
SNPs_Mal = SNPs_Mal.drop(['IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1)
SNPs_Mal = SNPs_Mal.set_index(['FID'])
SNPs_Mal = SNPs_Mal.T

SNPs_Fem = SNPs_Fem.drop(['IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1)
SNPs_Fem = SNPs_Fem.set_index(['FID'])
SNPs_Fem = SNPs_Fem.T


SNPs_Fem = SNPs_Fem.reset_index()
SNPs_Fem.rename(columns = {'index':'snp'}, inplace = True)

SNPs_Mal = SNPs_Mal.reset_index()
SNPs_Mal.rename(columns = {'index':'snp'}, inplace= True)

SNPs_Mal = SNPs_Mal.sort_values(by=['snp'])
SNPs_Fem = SNPs_Fem.sort_values(by=['snp'])
SNPs_locs = SNPs_locs.sort_values(by=['snp'])

SNPs_locs.insert(0, 'snps', range(1, 1 + len(SNPs_locs)))
SNPs_locs = SNPs_locs[['snps', 'chr', 'pos']]

SNPs_Mal.insert(0, 'snps', range(1, 1 + len(SNPs_Mal)))
SNPs_Mal = SNPs_Mal.drop('snp', axis = 1)

SNPs_Fem.insert(0, 'snps', range(1, 1 + len(SNPs_Fem)))
SNPs_Fem = SNPs_Fem.drop('snp', axis = 1)


# Organise SNP data to match Expression data patient order
SNPs_Fem = SNPs_Fem.reindex(sorted(SNPs_Fem.columns), axis=1)
cols = list(SNPs_Fem.columns)
cols = [cols[-1]] + cols[:-1]
SNPs_Fem = SNPs_Fem[cols]

SNPs_Mal = SNPs_Mal.reindex(sorted(SNPs_Mal.columns), axis = 1)
cols = list(SNPs_Mal.columns)
cols = [cols[-1]] + cols[:-1]
SNPs_Mal = SNPs_Mal[cols]

# Save expression data
Pre_Males.reset_index(level=0, inplace=True)
Pre_Males.rename(columns = {'Unnamed: 0':'Gene'}, inplace = True)

Pre_Females.reset_index(level=0, inplace=True)
Pre_Females.rename(columns = {'Unnamed: 0':'Gene'}, inplace = True)


# Output files 
#SNPs_Fem.to_csv('meQTL/eQTL_FemSNPs.csv', sep = '\t', index= False)
#SNPs_Mal.to_csv('meQTL/eQTL_MalSNPs.csv', sep = '\t', index= False)
#SNPs_locs.to_csv('meQTL/eQTL_SNPsLoc.csv', sep = '\t', index = False)
#Pre_Females.to_csv('meQTL/eQTL_FemExp.csv', sep = '\t', index = False)
#Pre_Males.to_csv('meQTL/eQTL_MalExp.csv', sep = '\t', index = False)




# Look into the results of cis eQTL's

CisFem = pd.read_csv('meQTL/eQTL_CisFem.csv', header = 0, index_col=0)
CisMal = pd.read_csv('meQTL/eQTL_CisMal.csv', header = 0, index_col=0)

# File for comparing snp number to snp name
SNPs_locs = pd.read_csv('eQTL_HOXC4SNPsLoc.csv', sep = '\t', header = 0)
SNPs_locs = SNPs_locs.sort_values(by=['snp'])
SNPs_locs.insert(0, 'snps', range(1, 1 + len(SNPs_locs)))























