#!/usr/bin/env python3

import numpy as np 

from local_Support_mM import *
import hash_tools as ht

output_file = "hash.dat"

ELEMENTS = np.loadtxt('ELEMENTS.dat')
IND_mask = np.loadtxt('IND_mask.dat').astype(int)
IND_mask_tot = np.loadtxt('IND_mask_tot.dat').astype(int)
U1 = np.loadtxt('U1.dat')
U2 = np.loadtxt('U2.dat')
U3 = np.loadtxt('U3.dat')

local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, 'standard')


with open(output_file, 'w') as file:
    file.write(ht.hash(local_support)+'\n')
    file.write(ht.hash(BF_support)+'\n')
    file.write(ht.hash(IND_mask_active)+'\n')
