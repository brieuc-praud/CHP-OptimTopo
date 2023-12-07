#!/usr/bin/env python3

import numpy as np 
import requests
import os
import warnings
import re

import hash_tools as ht

# Set test size
size=60


settings_file = "Problem_Setting.py"

size_str=str(size)
filename=2*(size_str + "x")+ size_str + ".npz"
hashfile="hash/hash"+size_str+".dat"

if (size == 30):
    url = "https://drive.google.com/uc?export=download&id=1x2chwlaU1CAVvV_2n44eAsBAoPXRtMTH&confirm=t"
    n1, n2, n3 = 27, 27, 27
elif (size == 60):
    url = "https://drive.google.com/uc?export=download&id=1oZcU0eCD9naiIyWGoX2pUNK5OaR4QyWy&confirm=t"
    n1, n2, n3 = 39, 39, 39
elif (size == 100):
    url = "https://drive.google.com/uc?export=download&id=1i0pPRIxCvnr76J8Z57A322bWX0bxOoZ6&confirm=t"
    n1, n2, n3 = 69, 69, 69
else:
    raise ValueError("Unavailable size " + size_str)

# Modify the setting file according to the test case
with open(settings_file, 'r') as file:
    filedata = file.read()
filedata = re.sub("n1=\d+; n2=\d+; n3=\d+;",
        "n1={n1:d}; n2={n2:d}; n3={n3:d};".format(n1=n1,n2=n2,n3=n3), filedata)
with open(settings_file, 'w') as file:
    file.write(filedata)

# Get test data
if not os.path.exists(filename):

    print("Downloading", filename)
    print("...")
    request = requests.get(url, allow_redirects="True")
    with open(filename, 'wb') as file:
        file.write(request.content)
    print("Done.")


# Load data
data = np.load(filename)
ELEMENTS = data['ELEMENTS']
IND_mask = data['IND_mask']
IND_mask_tot = data['IND_mask_tot']
U1 = data['U1']
U2 = data['U2']
U3 = data['U3']


from local_Support_mM import *
# Do the computation
local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, 'standard')


# Validation
print("=== VALIDATION:")
local_support_hash = ht.hash(local_support)
BF_support_hash = ht.hash(BF_support)
IND_mask_active_hash = ht.hash(IND_mask_active)

if os.path.exists(hashfile):
    with open(hashfile, 'r') as file:
        for array_hash, varname in zip([local_support_hash, BF_support_hash, IND_mask_active_hash], ["local_support", "BF_support", "IND_mask_active"]):
            print( array_hash )
            # Test if one hash fits
            validated = False
            for h in file.readline().strip().split():
                validated |= (array_hash == h)
            if (validated):
                print(varname, "->", "\033[92mPASS\033[0m")
            else:
                print(varname, "->", "\033[91mFAIL\033[0m")
else:
    warnings.warn("No test available for size "+size_str)
    for array_hash, varname in zip([local_support_hash, BF_support_hash, IND_mask_active_hash], ["local_support", "BF_support", "IND_mask_active"]):
        print(varname)
        print(array_hash)
print("===")
