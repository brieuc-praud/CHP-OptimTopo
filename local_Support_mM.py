'''--------------------------------------------------------------------------
# - Routine name: local_support.m
# - Routine description: Script for the calculation of the local support 
# - of the active CP, of the Basis Function matrix 
# --------------------------------------------------------------------------
#- INPUT
#- ELEMENTS: elements features matrix [shape (n_elem,15) for 3D and shape (n_elem,9) for 2D]
#- IND_mask: matrix of the indexes of the active CP stocked into the CP matrix [shape (n_cp_active,3)]
#- IND_mask_tot: matrix of the indexes of the active and inactive CP [shape (n_cp_tot,3)]
#- U1: KV uniformly distributed between [0,1] [shape (n1+p1+2,)]
#- U2: KV uniformly distributed between [0,1] [shape (n2+p2+2,)]
#- U3: KV uniformly distributed between [0,1] [shape (n3+p3+2,)] (optional)
#- flag_scale: parameter used to define the type of scale of the problem at hand, defined in Problem_Setting
#--------------------------------------------------------------------------
#- OUTPUT
#- local_Support: list of lists of indexes of elements contained in each CP's 
#-                local support [shape (n_CP_tot,)] 
#- BF_Support: matrix of Basis Function evaluated in the each parameter u1, u2, u3 
#-            (related to the element center of gravity) [shape (n_elem,n_CP_tot)]
#- IND_mask_active: list of indexes of active CP [shape (n_CP_active,)] 
# --------------------------------------------------------------------------'''
#General modules import:
from numba import njit
from Problem_Setting import *
import numpy as np
from SurfacePoint_numba import *
from HyperSurfacePoint_numba import *

def local_support_fun(ELEMENTS, IND_mask, IND_mask_tot, U1, U2, U3, flag_scale):
    #- Local support definition
    if flag_scale == 'micro_macro':
        p1_temp = p1M;
        p2_temp = p2M;
        p3_temp = p3M;
        n1_temp = n1M;
        n2_temp = n2M;
        n3_temp = n3M;
    else:
        p1_temp = p1;
        p2_temp = p2;
        p3_temp = p3;
        n1_temp = n1;
        n2_temp = n2;
        n3_temp = n3;

    if DIM == 2:
        u1 = ELEMENTS[:, 11]
        u2 = ELEMENTS[:, 12]
        local_Support, BF_Support, IND_mask_active = ls_2d_numba(IND_mask_tot, IND_mask, u1, u2, U1, U2, p1_temp, p2_temp, n1_temp, n2_temp)
            
    elif DIM == 3:
        u1 = ELEMENTS[:, 14]
        u2 = ELEMENTS[:, 15]
        u3 = ELEMENTS[:, 16]
        local_Support, BF_Support, IND_mask_active = ls_3d_numba(IND_mask_tot, IND_mask, u1, u2, u3, U1, U2, U3, p1_temp, p2_temp, p3_temp, n1_temp, n2_temp, n3_temp)

    return local_Support, BF_Support, IND_mask_active

@njit
def ls_2d_numba(IND_mask_tot, IND_mask, u1, u2, U1, U2, p1_temp, p2_temp, n1_temp, n2_temp):
    IND_mask_temp = [[IND_mask[i,j] for j in range(len(IND_mask[i]))] for i in range(len(IND_mask))]
    BF_Support = np.zeros((len(u1), len(IND_mask_tot)))
    IND_mask_active = []
    local_Support = []

    # - For each CP of the model (active and inactive):
    # - stock in local_support list the indexes of the elements which are into the LS
    # - stock in IND_mask_active the index related to the active CP, considering the totality of CP
    for k in range(len(IND_mask_tot)):
        IND_mask_tot_temp = [IND_mask_tot[k,i] for i in range(len(IND_mask_tot[k]))]
        b1 = np.logical_or(np.logical_and(u1 < U1[IND_mask_tot[k, 0] + p1_temp + 1], u1 >= U1[IND_mask_tot[k, 0]]), u1 == U1[IND_mask_tot[k, 0] + p1_temp + 1])
        b2 = np.logical_or(np.logical_and(u2 < U2[IND_mask_tot[k, 1] + p2_temp + 1], u2 >= U2[IND_mask_tot[k, 1]]), u2 == U2[IND_mask_tot[k, 1] + p2_temp + 1])
        contr = np.logical_and(b1, b2)
        ind = np.where(contr == np.array(True, np.bool_))[0]

        if IND_mask_tot_temp in IND_mask_temp:
            IND_mask_active.append(k)
            local_Support.append(list(ind))

        # Evaluate Basis Function product of elements belonging to the local support (try on all the elements belonging to LS)
        P_rho_aux = np.zeros((n1_temp + 1, n2_temp + 1))
        P_rho_aux[IND_mask_tot[k, 0], IND_mask_tot[k, 1]] = 1
        u = u1[ind].reshape((len(ind), 1))
        v = u2[ind].reshape((len(ind), 1))
        w = np.ones((n1_temp + 1, n2_temp + 1))
        BF_Support[ind, k] = SurfacePoint_fun_numba(u, n1_temp, p1_temp, U1, v, n2_temp, p2_temp, U2, P_rho_aux, w, 0)
    return local_Support, BF_Support, IND_mask_active

@njit
def ls_3d_numba(IND_mask_tot, IND_mask, u1, u2, u3, U1, U2, U3, p1_temp, p2_temp, p3_temp, n1_temp, n2_temp, n3_temp):
    IND_mask_temp = [[IND_mask[i,j] for j in range(len(IND_mask[i]))] for i in range(len(IND_mask))]
    BF_Support = np.zeros((len(u1), len(IND_mask_tot)))
    IND_mask_active = []
    local_Support = []

    # - For each CP of the model (active and inactive):
    # - stock in local_support list the indexes of the elements which are into the LS
    # - stock in IND_mask_active the index related to the active CP, considering the totality of CP
    for k in range(len(IND_mask_tot)):
        IND_mask_tot_temp = [IND_mask_tot[k,i] for i in range(len(IND_mask_tot[k]))]
        b1 = np.logical_or(np.logical_and(u1 < U1[IND_mask_tot[k, 0] + p1_temp + 1],u1 >= U1[IND_mask_tot[k, 0]]),u1 == U1[IND_mask_tot[k, 0] + p1_temp + 1])
        b2 = np.logical_or(np.logical_and(u2 < U2[IND_mask_tot[k, 1] + p2_temp + 1],u2 >= U2[IND_mask_tot[k, 1]]),u2 == U2[IND_mask_tot[k, 1] + p2_temp + 1])
        b3 = np.logical_or(np.logical_and(u3 < U3[IND_mask_tot[k, 2] + p3_temp + 1],u3 >= U3[IND_mask_tot[k, 2]]),u3 == U3[IND_mask_tot[k, 2] + p3_temp + 1])
        #TODO
        contr = np.logical_and(np.logical_and(b1,b2),b3)
        #contr = np.logical_and(b1, b2, b3)
        ind = np.where(contr == np.array(True, np.bool_))[0]

        if IND_mask_tot_temp in IND_mask_temp:
            IND_mask_active.append(k)
            local_Support.append(list(ind))

        # Evaluate Basis Function product of elements belonging to the local support (try on all the elements belonging to LS)
        P_rho_aux = np.zeros((n1_temp + 1, n2_temp + 1, n3_temp + 1))
        P_rho_aux[IND_mask_tot[k, 0], IND_mask_tot[k, 1], IND_mask_tot[k, 2]] = 1
        u = u1[ind].reshape((len(ind), 1, 1))
        v = u2[ind].reshape((len(ind), 1, 1))
        w = u3[ind].reshape((len(ind), 1, 1))
        w_t = np.ones((n1_temp + 1, n2_temp + 1, n3_temp + 1))
        BF_Support[ind, k] = HyperSurfacePoint_fun_numba(n1_temp, p1_temp, U1, n2_temp, p2_temp, U2, n3_temp, p3_temp, U3, P_rho_aux, w_t, u, v, w, 0)
    return local_Support, BF_Support, IND_mask_active