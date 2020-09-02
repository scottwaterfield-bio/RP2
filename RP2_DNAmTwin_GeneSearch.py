#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 18:48:07 2020

@author: sw1906
"""

import pandas as pd

cols= ['chr', 'pos', 'Name', 'UCSC_RefGene_Name', 'logFC',
                           'P.Value', 'adj.P.Val' ]


TwinFemData = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusFemales.csv',header = 0)
TwinFemData = TwinFemData[cols]

TwinMalData = pd.read_csv('RP2DNAmTwin_DMPsDiseaseStatusMales.csv', header = 0)
TwinMalData = TwinMalData[cols]

SSU72Fem = TwinFemData[TwinFemData['UCSC_RefGene_Name'] == 'SOCS3']
SSU72Mal = TwinMalData[TwinMalData['UCSC_RefGene_Name'] == 'SOCS3']
