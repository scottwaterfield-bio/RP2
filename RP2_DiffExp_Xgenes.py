#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 11:02:07 2020

@author: sw1906
"""

import pandas as pd

XCI = pd.read_excel('NIHMS905235-supplement-supp_table1.xlsx', index_col=None)


# SEPT6 was converted to a date!!
Xgenes = ["XIST","RPS4X", "KDM6A" ,"PRKX","ZFX", "SEPT6", "EIF1AX"  ]

XCI_interest = XCI[XCI['Gene name'].isin(Xgenes)]
                   
XCI_interest.to_csv('DiffExpXgenes.csv', sep = '\t', index= False)
