import os
import sys

# Obtenir le chemin du dossier courant
current_dir = os.getcwd()

main_dir = os.path.join(current_dir, "..")

sys.path.append(main_dir)

from compliance_mM import *
from volume_mM import *
from local_Support_mM import *
from local_Support_mM_CSR import * 
from volume_mM_CSR import *
from compliance_mM_CSR import *

test_dir = os.path.join(main_dir, "test")
os.chdir(test_dir)

ELEMENTS = np.loadtxt(os.path.join(main_dir, 'ELEMENTS.dat'))
rho_e = np.loadtxt(os.path.join(main_dir, 'rho_e_test.dat')).reshape(-1, 1)
IND_mask = np.loadtxt(os.path.join(main_dir, 'IND_mask_test.dat'), dtype=int)
IND_mask_tot = np.loadtxt(os.path.join(main_dir, 'IND_mask_tot_test.dat'), dtype=int)
U1 = np.loadtxt(os.path.join(main_dir, 'U1.dat'))
U2 = np.loadtxt(os.path.join(main_dir, 'U2.dat'))
U3 = np.loadtxt(os.path.join(main_dir, 'U3.dat'))

P_rho = np.loadtxt(os.path.join(main_dir, 'P_rho_test.dat')).reshape(28, 28, 28)
W = np.loadtxt(os.path.join(main_dir, 'W_test.dat')).reshape(28, 28, 28)

local_support, BF_support, IND_mask_active = local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, scale)

grad_c_T = compliance_grad_fun(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, scale)
grad_v = volume_grad_fun(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, scale)

# Calculez les résultats avec les fonctions CSR
local_support_csr, BF_support_csr, IND_mask_active_csr = local_support_fun_csr(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, scale)
grad_c_T_csr = compliance_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, scale)#mettre '_csr' après les trois arg quand compliance csr faite
grad_v_csr = volume_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support_csr, BF_support_csr, IND_mask_tot, IND_mask_active_csr, scale)



# Comparaison résultats
if np.allclose(grad_c_T, grad_c_T_csr):
    print("TEST VALIDE : Les résultats de grad_c_T et grad_c_T_csr sont identiques.")
else:
    print("TEST NON VALIDE : Les résultats de grad_c_T et grad_c_T_csr ne sont pas identiques.")

if np.allclose(grad_v, grad_v_csr):
    print("TEST VALIDE : Les résultats de grad_v et grad_v_csr sont identiques.")
else:
    print("TEST NON VALIDE : Les résultats de grad_v et grad_v_csr ne sont pas identiques.")

#if np.allclose(local_support, local_support_csr) and np.allclose(BF_support, BF_support_csr.toarray()) and np.allclose(IND_mask_active, IND_mask_active_csr) :
    #print("TEST VALIDE : Les résultats de local_support_fun sont identiques.")
#else:
    #print("TEST NON VALIDE : Les résultats de local_support_fun ne sont pas identiques.")


os.chdir(current_dir)

sys.path.remove(main_dir)
