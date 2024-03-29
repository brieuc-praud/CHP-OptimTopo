#!/usr/bin/env python3

import os
import sys

#bibliothèque chronomètrage
from time import perf_counter

#---------------------------------------------

# Obtenir le chemin du dossier courant
current_dir = os.getcwd()
main_dir = os.path.join(current_dir, "../gradients")

sys.path.append(main_dir)
os.chdir(main_dir)

#---------------------------------------------

from compliance_mM import *
from volume_mM import *
from local_Support_mM import *
from local_Support_mM_CSR import * 
from volume_mM_CSR import *
from compliance_mM_CSR import *

#---------------------------------------------

ELEMENTS = np.loadtxt(os.path.join(main_dir, 'ELEMENTS.dat'))
rho_e = np.loadtxt(os.path.join(main_dir, 'rho_e_test.dat')).reshape(-1, 1)
IND_mask = np.loadtxt(os.path.join(main_dir, 'IND_mask_test.dat'), dtype=int)
IND_mask_tot = np.loadtxt(os.path.join(main_dir, 'IND_mask_tot_test.dat'), dtype=int)
U1 = np.loadtxt(os.path.join(main_dir, 'U1.dat'))
U2 = np.loadtxt(os.path.join(main_dir, 'U2.dat'))
U3 = np.loadtxt(os.path.join(main_dir, 'U3.dat'))

P_rho = np.loadtxt(os.path.join(main_dir, 'P_rho_test.dat')).reshape(28, 28, 28)
W = np.loadtxt(os.path.join(main_dir, 'W_test.dat')).reshape(28, 28, 28)

#---------------------------------------------
t_start = perf_counter()
local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, scale)

grad_c_T = compliance_grad_fun(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, scale)
grad_v = volume_grad_fun(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, scale)
t_finish = perf_counter()

#---------------------------------------------

# Calcul des résultats avec les fonctions CSR
t_start_csr = perf_counter()
local_support_csr, BF_support_csr, IND_mask_active_csr = local_support_fun_csr(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, scale)

grad_c_T_csr = compliance_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support_csr, BF_support_csr, IND_mask_tot, IND_mask_active_csr, scale)#mettre '_csr' après les trois arg quand compliance csr faite
grad_v_csr = volume_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support_csr, BF_support_csr, IND_mask_tot, IND_mask_active_csr, scale)
t_finish_csr = perf_counter()


#---------------------------------------------

# Comparaison résultats
if np.allclose(grad_c_T, grad_c_T_csr):
    print("TEST VALIDE : Les résultats de grad_c_T et grad_c_T_csr sont identiques.")
else:
    print("TEST NON VALIDE : Les résultats de grad_c_T et grad_c_T_csr ne sont pas identiques.")

if np.allclose(grad_v, grad_v_csr):
    print("TEST VALIDE : Les résultats de grad_v et grad_v_csr sont identiques.")
else:
    print("TEST NON VALIDE : Les résultats de grad_v et grad_v_csr ne sont pas identiques.")


print("TEMPS DE CALCUL ( en secondes ):")
print("Code d'origine ( non-optimisé ) :", t_finish - t_start)
print("Code optimisé ( lil/csr ) :", t_finish_csr - t_start_csr)
print("Différence de temps entre les deux versions : temps_original - temps_csr )", t_finish - t_start - ( t_finish_csr - t_start_csr ))
#---------------------------------------------

os.chdir(current_dir)

sys.path.remove(main_dir)
