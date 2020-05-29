# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 18:01:45 2020

@author: mra013
"""

import os
os.chdir('..\\fsf')

import fsf_engine

def Call_Case(a):
    if a==1:
        path='..\\inputs\Static_Test_3.txt'
           
    elif a==2:
        path='..\\inputs\Steady_Test.txt'
    
    elif a==3:
        path='..\\inputs\Transient_Test_A_2.txt'
            
    R1=fsf_engine.Calc_FSF(path)
    R2=R1[0,:,0]
    return R2


# =============================================================================
# # Call the cases from the command below
# =============================================================================
Result=Call_Case(3)

