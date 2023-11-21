
import numpy as np 
from local_Support_mM import *

size=30
filename_prefix="./NPZ/"
filename_extension=".npz"

size_str=str(size)
filename= 2*(size_str + "x") + str(size)

data = np.load(filename_prefix+filename+filename_extension)
ELEMENTS = data['ELEMENTS']
IND_mask = data['IND_mask']
IND_mask_tot = data['IND_mask_tot']
U1 = data['U1']
U2 = data['U2']
U3 = data['U3']

local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, 'standard')
#np.savetxt('Bf_support.dat',BF_support)

