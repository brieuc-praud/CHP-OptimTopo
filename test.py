
import numpy as np 
from local_Support_mM import *

ELEMENTS = np.loadtxt('ELEMENTS.dat')
IND_mask = np.loadtxt('IND_mask.dat',dtype=int)
IND_mask_tot = np.loadtxt('IND_mask_tot.dat',dtype=int)
U1 = np.loadtxt('U1.dat')
U2 = np.loadtxt('U2.dat')
U3 = np.loadtxt('U3.dat')

local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, 'standard')
np.savetxt('Bf_support.dat',BF_support)

