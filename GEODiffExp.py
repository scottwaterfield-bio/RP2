#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 16:37:25 2020

@author: sw1906
"""

import pandas as pd

DiffExpFem = pd.read_csv('GEODiffExpFemales.csv',sep = '\t', header = 0)
DiffExpFem = DiffExpFem[DiffExpFem['adj.P.Val'] < 0.05] #274 sig results 
OverRA = DiffExpFem[DiffExpFem['logFC'] > 0 ] # 92 over in RA, 78 named genes
UnderRA = DiffExpFem[DiffExpFem['logFC'] < 0 ] # 182 under in RA, 161 named 


OverRAgenes = list(OverRA['Gene.symbol'])
OverRAgenes
UnderRAgenes = list(UnderRA['Gene.symbol'])

[print(i) for i in OverRAgenes]
[print(i) for i in UnderRAgenes]
