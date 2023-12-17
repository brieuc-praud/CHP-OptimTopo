"""
Created on Wed Oct 21 16:15:59 2020

@author: gbertolino
"""
from Problem_Setting import *

def der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e):
    
    if DIM == 2:
        
        local_support_flat = [item for sublist in local_support for item in sublist]
        BF_mask = list(set(local_support_flat))
        BF_support_temp = BF_support[BF_mask,:]
        BF_support_temp = BF_support_temp[:,IND_mask_active]
        W_temp = W[IND_mask[:,0],IND_mask[:,1]].reshape((len(IND_mask),1))
        W_S_temp = W[IND_mask_tot[:,0],IND_mask_tot[:,1]].reshape((len(IND_mask_tot),1))
        P_temp = P_rho[IND_mask[:,0],IND_mask[:,1]].reshape((len(IND_mask),1))
        
        Nij_w = BF_support_temp * W_temp.T
        S_w = (BF_support.dot(W_S_temp))[BF_mask]
        der_CP = Nij_w/S_w
        
        Nij_P = (BF_support_temp * P_temp.T)
        Nij_nurbs = (BF_support_temp * rho_e[BF_mask])
        der_W = Nij_P/S_w - Nij_nurbs/S_w
    
    elif DIM == 3:
        local_support_flat = [item for sublist in local_support for item in sublist]
        BF_mask = list(set(local_support_flat))
        BF_support_temp = BF_support[BF_mask,:]
        BF_support_temp = BF_support_temp[:,IND_mask_active]
        W_temp = W[IND_mask[:,0],IND_mask[:,1],IND_mask[:,2]].reshape((len(IND_mask),1))
        W_S_temp = W[IND_mask_tot[:,0],IND_mask_tot[:,1],IND_mask_tot[:,2]].reshape((len(IND_mask_tot),1))
        P_temp = P_rho[IND_mask[:,0],IND_mask[:,1],IND_mask[:,2]].reshape((len(IND_mask),1))
        
        Nij_w = BF_support_temp * W_temp.T
        S_w = (BF_support.dot(W_S_temp))[BF_mask]
        der_CP = Nij_w/S_w
        
        Nij_P = (BF_support_temp * P_temp.T)
        Nij_nurbs = (BF_support_temp * rho_e[BF_mask])
        der_W = Nij_P/S_w - Nij_nurbs/S_w
        
    return der_CP, der_W, BF_mask

def der_BSPLINE(IND_mask_active,BF_support):
    
    BF_support_temp = BF_support[:,IND_mask_active]
    
    return BF_support_temp