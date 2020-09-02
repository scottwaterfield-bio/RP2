#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 19:57:13 2020

@author: sw1906
"""

import pandas as pd

# Load in data


cols= ['chr', 'pos','strand', 'Name', 'UCSC_RefGene_Name', 'logFC', 'Relation_to_Island',
       'Enhancer', 'Regulatory_Feature_Group', 'P.Value', 'adj.P.Val' ]


TwinRA = pd.read_csv('RP2DNAmTwin_DMPsGender_RA_Only.csv',header = 0)
TwinRA = TwinRA[cols]
IndRA = pd.read_csv('RP2DNAmInd_Gender.csv')
IndRA = IndRA[cols]

PreMTX = pd.read_csv('RP2PreMTXDiffExp_Gender.csv', header=0)
PostMTX = pd.read_csv('RP2PostMTXDiffExp_Gender.csv', header=0)

#Drop NaN Refgene rows for this purpose
TwinRA = TwinRA.dropna(axis=0, subset=['UCSC_RefGene_Name'])
IndRA = IndRA.dropna(axis=0, subset=['UCSC_RefGene_Name'])
#----------------------------------------------------------------------------#

# Looking at DNAm results of significant autosomal DiffExp results
 

#PALLD and MED8 significant upregulated in males before MTX treatment
PALLDTwin = TwinRA[TwinRA['UCSC_RefGene_Name'].str.contains('PALLD')] #79
PALLDInd = IndRA[IndRA['UCSC_RefGene_Name'].str.contains('PALLD')] #76


MED8Twin = TwinRA[TwinRA['UCSC_RefGene_Name'].str.contains('MED8')]
MED8Ind = IndRA[IndRA['UCSC_RefGene_Name'].str.contains('MED8')]
MED8Twin = MED8Twin[~MED8Twin['UCSC_RefGene_Name'].str.contains('TMED')] #24
MED8Ind = MED8Ind[~MED8Ind['UCSC_RefGene_Name'].str.contains('TMED')]#23


#DDIT4 expression sigdiff between gender after MTX treatment
DDIT4Twin = TwinRA[TwinRA['UCSC_RefGene_Name'].str.contains('DDIT4')]
DDIT4Ind = IndRA[IndRA['UCSC_RefGene_Name'].str.contains('DDIT4')]
DDIT4Twin = DDIT4Twin[~DDIT4Twin['UCSC_RefGene_Name'].str.contains('DDIT4L')] #19
DDIT4Ind = DDIT4Ind[~DDIT4Ind['UCSC_RefGene_Name'].str.contains('DDIT4L')] #19


# Looking at 'Potentially significant' (Adj. P <0.25) CpGs in these genes
DDIT4 = DDIT4Twin[DDIT4Twin['adj.P.Val'] < 0.25] #8/19 
PALLD = PALLDTwin[PALLDTwin['adj.P.Val'] < 0.25] #18/79 All but one negative logFC
MED8 = MED8Twin[MED8Twin['adj.P.Val'] < 0.25] #4/24

#----------------------------------------------------------------------------#
#Autosomal Genes significant in both DNAm Twin and Ind
# Haven't used control TwinHealthy to rule out normal gender differences yet

AutoGenes = ['TFDP1', 'TLE1', 'TIMM8B', 'UBE2QP1', 'SSC5D', 'PPT2']
PreMTX[PreMTX['SYMBOL'].isin(AutoGenes)]
PostMTX[PostMTX['SYMBOL'].isin(AutoGenes)]




