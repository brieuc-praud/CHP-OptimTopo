import numpy as np 

ELEMENTS = np.loadtxt('ELEMENTS.dat')
IND_mask = np.loadtxt('IND_mask.dat').astype(int)
IND_mask_tot = np.loadtxt('IND_mask_tot.dat').astype(int)
U1 = np.loadtxt('U1.dat')
U2 = np.loadtxt('U2.dat')
U3 = np.loadtxt('U3.dat')

ELEMENTS.tofile('ELEMENTS.bin')
IND_mask.tofile('IND_mask.bin')
IND_mask_tot.tofile('IND_mask_tot.bin')
U1.tofile('U1.bin')
U2.tofile('U2.bin')
U3.tofile('U3.bin')

np.savez_compressed("30x30x30",
        ELEMENTS=ELEMENTS,
        IND_mask=IND_mask,
        IND_mask_tot=IND_mask_tot,
        U1=U1,
        U2=U2,
        U3=U3)
