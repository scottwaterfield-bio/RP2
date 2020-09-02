#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 15:11:12 2020

@author: sw1906
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# All SNPs
SNPs = pd.read_csv('meQTL/HOXC4GWAS.raw', sep = ' ', header = 0)

# Patients with known gender (removes 2)
SNPs = SNPs[SNPs.SEX != 0]

# check gender split
SNPs[SNPs.SEX == 1] #459
SNPs[SNPs.SEX == 2] #1511

# Remove irrelavant columns
SNPs = SNPs.drop(['FID','IID', 'PAT', 'MAT', 'PHENOTYPE'], axis=1)
maldf = SNPs.loc[SNPs['SEX'] == 1]
femdf = SNPs.loc[SNPs['SEX'] == 2]

SNPName = SNPs.filter(regex='rs538217839')


# Gender proportion of SNPs of interest

# rs12814447
(maldf['rs12814447_G'].value_counts()/maldf['rs12814447_G'].count())*100
#0    91.067538
#1     8.932462

(femdf['rs12814447_G'].value_counts()/femdf['rs12814447_G'].count())*100
#0    92.058240
#1     7.610854
#2     0.330907
#-----------------------------------------------------------------------------#
#rs182062166
(maldf['rs182062166_T'].value_counts()/maldf['rs182062166_T'].count())*100
#0    100

(femdf['rs182062166_T'].value_counts()/femdf['rs182062166_T'].count())*100
#0    99.735275
#1     0.264725
#-----------------------------------------------------------------------------#
#rs189988893
(maldf['rs189988893_T'].value_counts()/maldf['rs189988893_T'].count())*100
#0    100

(femdf['rs189988893_T'].value_counts()/femdf['rs189988893_T'].count())*100
#0    99.801456
#1     0.198544
#-----------------------------------------------------------------------------#
#rs34448539
(maldf['rs34448539_A'].value_counts()/maldf['rs34448539_A'].count())*100
#0    91.067538
#1     8.932462

(femdf['rs34448539_A'].value_counts()/femdf['rs34448539_A'].count())*100
#0    92.190602
#1     7.478491
#2     0.330907
#----------------------------------------------------------------------------#
#rs538217839
(maldf['rs538217839_A'].value_counts()/maldf['rs538217839_A'].count())*100
#0    100

(femdf['rs538217839_A'].value_counts()/femdf['rs538217839_A'].count())*100
#0    99.735275
#1     0.264725






























