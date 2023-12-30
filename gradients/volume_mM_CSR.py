'''all TESTED
#--------------------------------------------------------------------------
Routine name: volume.py

Routine description: volume objective/constraint function
Evaluation of the volume of the structure. The respective derivative is 
computed as well (the derivative is computed with respect to the control 
points and/or to the weights, depending on the problem at hand).

Symmetries are handle by varialble "sym_coef" calculated in Problem_Setting.py
#--------------------------------------------------------------------------
volume_fun
STANDARD INPUT ARGUMENT 
    rho_e, ELEMENTS
OPTIONAL INPUT ARGUMENT   

OUTPUT
    v
#--------------------------------------------------------------------------
volume_grad_fun
STANDARD INPUT ARGUMENT
    rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active
OPTIONAL INPUT ARGUMENT   
    args = flag_scale -> when scale is 0, 'micro', 'micro_macro'
OUTPUT
    grad_v_T
#--------------------------------------------------------------------------
'''
#Numpy import:
import numpy as np

#Global variables setting up the study:
from Problem_Setting import *

#Surface and Hypersurface computation import:
from derTopo import der_NURBS, der_BSPLINE
#--------------------------------------------------------------------------
def volume_fun(rho_e, ELEMENTS, *args):
    #--------------------------------------------------------------------------
    flag_scale = args[0]
    if flag_scale == 'micro':
        v_NDR = 0
    else:
        from Problem_Setting import v_NDR

    if DIM==2: 
        vol_temp = (ELEMENTS[:,8]*ELEMENTS[:,10]).reshape((len(ELEMENTS),1))                   
        v=np.sum(vol_temp*rho_e) 
        
    elif DIM==3:
        vol_temp = (ELEMENTS[:,12]).reshape((len(ELEMENTS),1))
        v=np.sum(vol_temp*rho_e) 
    #--------------------------------------------------------------------------
    return v+v_NDR

def volume_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, *args):#BF_support is a lil/csr matrix
    #--------------------------------------------------------------------------
    if DERIVATIVES ==1:

        flag_scale = args[0]
        if flag_scale == 'micro_macro':   #MACROSCALE 
            sym_coef_temp = sym_coefM
            
        elif flag_scale != 'micro_macro':   #MICROSCALE / Standard scale
            sym_coef_temp = sym_coef
        
        #--------------------------------------------------------------------------
        if DIM==2:            
            #--------------------------------------------------------------------------
            if NURBS==1: #NURBS Surfaces
                # Derivatives of compliance respecting to micro scale design variables 
                der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)#all outputs are lil/csr

                #changing them to numpy because the lil/sar format does not wotrh it for these :
                der_CP = der_CP.toarray()
                der_W = der_W.toarray()
                BF_mask = BF_mask.toarray()
                               
                vol_temp  = (ELEMENTS[BF_mask,8]*ELEMENTS[BF_mask,10]).reshape((len(BF_mask),1))
                grad_v_cp = np.sum(sym_coef_temp*der_CP*vol_temp,axis=0).reshape((len(IND_mask),1)) 
                grad_v_w = np.sum(sym_coef_temp*der_W*vol_temp,axis=0).reshape((len(IND_mask),1))
                grad_v= np.concatenate((grad_v_cp,grad_v_w)) 

            #--------------------------------------------------------------------------        
            else: #BSPLINE Surface
                # Derivatives of compliance respecting to micro scale design variables 
                BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)#we get there a lil/csr matrix
                if flag_shepard == 1:
                    BF_support_temp = shepard_support(BF_support_temp,shepard_dist)
                
                vol_temp  = (ELEMENTS[:,8]*ELEMENTS[:,10]).reshape((len(ELEMENTS),1))
                grad_v = np.sum(sym_coef_temp*BF_support_temp*vol_temp,axis=0).reshape((len(IND_mask),1))#np.sum outputs an np matrix even though lil/csr in argument, we want grad_v to not be lil/csr anymore because small enough  
        
        elif DIM==3:
            #--------------------------------------------------------------------------
            if NURBS==1: #NURBS Surfaces                  
                # Derivatives of compliance respecting to micro scale design variables 
                der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)#all outputs are lil/csr

                #changing them to numpy because the lil/sar format does not wotrh it for these :
                der_CP = der_CP.toarray()
                der_W = der_W.toarray()
                BF_mask = BF_mask.toarray()
                
                vol_temp = (ELEMENTS[BF_mask,12]).reshape((len(BF_mask),1))
                grad_v_cp = np.sum(sym_coef_temp*der_CP*vol_temp,axis=0).reshape((len(IND_mask),1))
                grad_v_w = np.sum(sym_coef_temp*der_W*vol_temp,axis=0).reshape((len(IND_mask),1))
                grad_v= np.concatenate((grad_v_cp,grad_v_w)) 

            #--------------------------------------------------------------------------
            else: #BSPLINE Surface                
                # Derivatives of compliance respecting to micro scale design variables 
                BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)
                if flag_shepard == 1:
                    BF_support_temp = shepard_support(BF_support_temp,shepard_dist)
                
                vol_temp = (ELEMENTS[:,12]).reshape((len(ELEMENTS),1))
                grad_v= np.sum(sym_coef_temp*BF_support_temp*vol_temp,axis=0).reshape((len(IND_mask),1)) 
    #--------------------------------------------------------------------------
    else:
        grad_v=False
    #--------------------------------------------------------------------------
    return grad_v
