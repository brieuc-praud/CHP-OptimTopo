'''all TESTED
#--------------------------------------------------------------------------
Routine name: compliance_mM.py
Routine description: compliance objective/constraint function
Evaluation of the compliance of the structure. The respective derivative  
is computed as well (the derivative is computed with respect to the 
control points and/or to the weights, depending on the problem at hand). 
Symmetries are handle by varialble "sym_coef" calculated in Problem_Setting.py
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
compliance_fun
STANDARD INPUT ARGUMENT 
    rho_e: vector of the pseudo density evaluated at the center of gravity of each element [shape (len(ELEMENTS),1)]
    ELEMENTS: elements features matrix [shape (n_elem,15) for 3D and shape (n_elem,9) for 2D]
#--------------------------------------------------------------------------
OPTIONAL INPUT ARGUMENT   
    flag_scale: parameter used to define the type of scale of the problem at hand, defined in Problem_Setting  (always)
    rho_e_M: vector of the pseudo density evaluated at the center of gravity of each element of the macroscopic scale [shape (len(ELEMENTS_macro),1)] (only if the problem is multiscale - double scale)
    ELEMENTS_macro: elements features matrix [shape (n_elem,15) for 3D and shape (n_elem,9) for 2D] (only if the problem is multiscale - double scale)
#--------------------------------------------------------------------------
OUTPUT
    c: objective function value (scalar)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
compliance_grad_fun
STANDARD INPUT ARGUMENT
    rho_e: vector of the pseudo density evaluated at the center of gravity of each element [shape (len(ELEMENTS),1)]
    P_rho: matrix of CP coordinates along the 3rd or 4th direction (defining the density concentrated at the CP) [shape (n1+1,n2+1,n3+1)]
    W: matrix of Weights defining the weight of the pseudo density concentrated at the CP [shape (n1+1,n2+1,n3+1)]
    ELEMENTS: elements features matrix [shape (n_elem,15) for 3D and shape (n_elem,9) for 2D]
    IND_mask: matrix of the indexes of the active CP stocked into the CP matrix [shape (n_CP_active,3)]
    local_support: list of lists of indexes of elements contained in each CP's local support [shape (n_CP_tot,)] 
    BF_support: matrix of Basis Function evaluated in the each parameter u1, u2, u3 (related to the element center of gravity) [shape (n_elem,n_CP_tot)] 
    IND_mask_tot: list of indexes of all the CP [shape (n_CP_tot,)] 
    IND_mask_active: list of indexes of active CP [shape (n_CP_active,)] 
#--------------------------------------------------------------------------
OPTIONAL INPUT ARGUMENT   
    flag_scale: parameter used to define the type of scale of the problem at hand, defined in Problem_Setting (always)
    ELEMENTS_macro: elements features matrix of the macroscopic scale [shape (n_elem_macro,15) for 3D and shape (n_elem_macro,9) for 2D] (only if the problem is multiscale - lower/double scale)
    rho_e_M: vector of the pseudo density evaluated at the center of gravity of each element of the macroscopic scale [shape (len(ELEMENTS_macro),1)] (only if the problem is multiscale - double scale)
    P_rho_M: matrix of CP coordinates along the 3rd or 4th direction (defining the density concentrated at the CP) [shape (n1M+1,n2M+1,n3M+1)] (only if the problem is multiscale - double scale)
    W_M: matrix of Weights defining the weight of the pseudo density concentrated at the CP [shape (n1M+1,n2M+1,n3M+1)] (only if the problem is multiscale - double scale)
    IND_mask_M: matrix of the indexes of the active CP stocked into the CP matrix [shape (n_CP_M_active,3)] (only if the problem is multiscale - double scale)
    local_support_M: list of lists of indexes of elements contained in each CP's local support [shape (n_CP_M_tot,)] (only if the problem is multiscale - double scale)
    BF_support_M: matrix of Basis Function evaluated in the each parameter u1, u2, u3 (related to the element center of gravity) [shape (n_elem_macr,n_CP_M_tot)] (only if the problem is multiscale - double scale)
    IND_mask_M_tot: list of indexes of all the CP [shape (n_CP_M_tot,)] (only if the problem is multiscale - double scale)
    IND_mask_active_M: list of indexes of active CP [shape (n_CP_M_active,)] (only if the problem is multiscale - double scale)
#--------------------------------------------------------------------------
OUTPUT
    grad_c_T: vector of gradients of the objective function respecting to the design variable vector [shape (len(x),1)] 
#--------------------------------------------------------------------------
'''
#Numpy, OS import:
import numpy as np
from numpy import linalg as LA
import os
from time import process_time

#Scipy.sparse import:
import scipy.sparse as sc


#Global variables setting up the study:
from Problem_Setting import *

#Surface and Hypersurface computation import:
from derTopo import der_NURBS, der_BSPLINE

#--------------------------------------------------------------------------

def compliance_fun (rho_e, ELEMENTS, *args):
        global c_vec
            
        flag_scale = args[0]
        if flag_scale == 'micro_macro':
            rho_e_M = args[1]
            ELEMENTS_macro = args[2]

        os.chdir(''.join(op_dir))
        # Writing the rho_e file for ANSYS computation for the microscale/standard analysis
        mask = np.concatenate((ELEMENTS[:,0].reshape((len(ELEMENTS),1)),rho_e),axis = 1)            
        fid = open('rho_e.txt','w')
        for i in range (len(ELEMENTS[:,0])):
                fid.write("%23.15e %23.15e\n"  % (mask[i,0], mask[i,1]))
        fid.close()
        del (fid)

        if flag_scale == 'micro_macro':
            # Writing the rho_e file for ANSYS computation for the macroscale analysis
            mask = np.concatenate((ELEMENTS_macro[:,0].reshape((len(ELEMENTS_macro),1)),rho_e_M),axis = 1)            
            fid = open('rho_e_M.txt','w')
            for i in range (len(ELEMENTS_macro[:,0])):
                fid.write("%23.15e %23.15e \n" % (mask[i, 0], mask[i, 1]))
            fid.close()
            del (fid)
        
        ANSYS_stat = 100
        
        while np.logical_and(np.logical_and(ANSYS_stat != 0, ANSYS_stat != 1),ANSYS_stat!=8):
            ANSYS_stat=os.system(ANSYS_comline_opt)
            os.system(del_string)

        if 'compliance' in con_fun:
            c_index = args[1]
            if c_index == 0:
                fid = open('FEM_SIMP_RES.txt')
            elif c_index == 1:
                fid=open('FEM_SIMP_RES_mom.txt')
        else:
            fid = open('FEM_SIMP_RES.txt')
            #TODO
            if flag_BCs == 'mixed-w' or flag_BCs == 'mixed-c':
                c = np.loadtxt(fid)
            else:                
                c_vec = np.loadtxt(fid)*2
                c = np.sum(c_vec)
                
        fid.close()
        del (fid)

        
        os.chdir(WD)
        #--------------------------------------------------------------------------      
        return c


def compliance_grad_fun_csr(rho_e, P_rho, W, ELEMENTS, IND_mask, local_support, BF_support, IND_mask_tot, IND_mask_active, *args):
        #--------------------------------------------------------------------------
        if DERIVATIVES == 1:

            flag_scale = args[0]            

            if flag_scale == 'micro_macro':
                rho_e_M = args[1]
                P_rho_M = args[2]
                W_M = args[3]
                ELEMENTS_macro = args[4]
                IND_mask_M = args[5]
                local_support_M = args[6]
                BF_support_M = args[7]
                IND_mask_M_tot = args[8]
                IND_mask_active_M = args[9]
            elif flag_scale == 'micro':                
                ELEMENTS_macro = args[1]                               

            # Reading files
            os.chdir(''.join(op_dir))
            if flag_scale == 'micro_macro' or flag_scale == 'micro':
                # Reading micro strain energy
                fid=open('FEM_SIMP_RES_micro.txt')
                str_en_vec = np.loadtxt(fid)
                c_vec_micro = 2*str_en_vec
                fid.close()
                del (fid)
                # Reading macro strains
                fid=open('STRAIN_ELEMENTS.txt')
                strain_M = np.loadtxt(fid)
                strain_M = strain_M[np.int_(ELEMENTS_macro[:,0])-1]
                fid.close()
                del (fid)
            else:
                if flag_BCs == 'mixed-w':
                    fid = open('C_diff_elements.txt')
                    c_vec = np.loadtxt(fid)
                #TODO
                elif flag_BCs == 'mixed-c':
                    fid = open('C_diff_elements.txt')
                    c_vec = np.loadtxt(fid)*2
                else:
                    fid = open('FEM_SIMP_RES.txt')
                    c_vec = np.loadtxt(fid)*2
                if 'compliance' in con_fun:
                    c_index = args[1]
                    os.chdir(''.join(op_dir))
                    if c_index == 0:
                        fid = open('FEM_SIMP_RES.txt')
                    elif c_index == 1:
                        fid = open('FEM_SIMP_RES_mom.txt')
                    c_vec = np.loadtxt(fid) * 2
                fid.close()
                del (fid)
            os.chdir(WD)
            #--------------------------------------------------------------------------
            if DIM == 2:
                #--------------------------------------------------------------------------
                if NURBS == 1:                    
                    if flag_scale == 'micro_macro' or flag_scale == 'micro':                        
                        # 1 - Derivatives of topology descriptor
                        der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)
                        
                        ### Ajout - Passage des variables CSR suivantes en format array
                        der_CP = der_CP.toarray()
                        der_W = der_W.toarray()
                        BF_mask = BF_mask.toarray()

                        # 2 - Derivatives of Stiffness Matrix coefficients (anisotropic, orthotropic)
                        grad_C_coef_cp = np.zeros((size_grad_C,len(IND_mask)),'f')  
                        grad_C_coef_w = np.zeros((size_grad_C,len(IND_mask)),'f')
                        #for a in range(3):
                            #c_vec_micro_temp = c_vec_micro[BF_mask,a].reshape((len(BF_mask),1))
                            #grad_C_coef_cp[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*der_CP/(a1*a2*rho_e[BF_mask]),axis=0)
                            #grad_C_coef_w[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*der_W/(a1*a2*rho_e[BF_mask]),axis=0)                        
                        #cont = 1
                        #for i in range(i_range_end):
                            #for j in range(i+1,j_range_end):
                                #comp_temp = (c_vec_micro[BF_mask,a+cont]-c_vec_micro[BF_mask,i]-c_vec_micro[BF_mask,j]).reshape((len(BF_mask),1))
                                #grad_C_coef_cp[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*der_CP/(2*a1*a2*rho_e[BF_mask]),axis=0)   
                                #grad_C_coef_w[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*der_W/(2*a1*a2*rho_e[BF_mask]),axis=0)
                                #cont += 1

                         # Reshape c_vec_micro for easier indexing
                        c_vec_micro_reshaped = c_vec_micro[BF_mask].reshape(len(BF_mask), 3, 1)

                        for a in range(3):
                            c_vec_micro_temp = c_vec_micro_reshaped[:, a, :]
                            common_term = sym_coef * p_c * c_vec_micro_temp / (a1 * a2 * rho_e[BF_mask])

                            grad_C_coef_cp[a, :] = np.sum(common_term * der_CP, axis=0)
                            grad_C_coef_w[a, :] = np.sum(common_term * der_W, axis=0)

                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i + 1, j_range_end):
                                comp_temp = (c_vec_micro_reshaped[:, a + cont] - c_vec_micro_reshaped[:, i] - c_vec_micro_reshaped[:, j])
                                grad_C_coef_cp[a + cont, :] = np.sum(sym_coef * p_c * comp_temp / (2 * a1 * a2 * rho_e[BF_mask]), axis=0)
                                grad_C_coef_w[a + cont, :] = np.sum(sym_coef * p_c * comp_temp / (2 * a1 * a2 * rho_e[BF_mask]), axis=0)
                                cont += 1

                        
                        # 3 - Derivatives of the Macro scale Compliance respecting to micro scale design variables                                                             
                        if flag_scale == 'micro':
                            # Topology descriptor is defined just at micro scale  
                            rho_e_M_temp = 1
                        else:
                            # Topology descriptor is defined also at the macro scale 
                            rho_e_M_temp = rho_e_M.reshape((len(rho_e_M),))
                            strain_M = strain_M[np.int_(ELEMENTS_macro[:, 0]) - 1]
                        
                        # Dependencies on the macroscopic strains and on the derivatives of the Stiffness Matrix coefficients                           
                        grad_c_cp_matrix = np.zeros((len(strain_M),len(IND_mask)))
                        grad_c_w_matrix = np.zeros((len(strain_M),len(IND_mask)))
                        
                        # Evaluating the contributes of the gradient of the diagonal Stiffness Matrix coefficients 
                        for a in range(3): 
                            # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                            # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                            grad_C_coef_cp_matrix = grad_C_coef_cp[a,:]*np.ones((len(strain_M),1))
                            grad_C_coef_w_matrix = grad_C_coef_w[a,:]*np.ones((len(strain_M),1))
                            strain_e = (strain_M[:,a]**2).reshape((len(strain_M),1))
                            
                            grad_c_cp_matrix = grad_c_cp_matrix + grad_C_coef_cp_matrix*strain_e
                            grad_c_w_matrix = grad_c_w_matrix + grad_C_coef_w_matrix*strain_e
                        
                        # Evaluating the contributes of the gradient of the out of diagonal Stiffness Matrix coefficients
                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                # Step to obtain a matrix (n_elem * n_ind_mask) of the same row (1,n_ind_mask) repeted n_elem times, to multiply each element contribute, of strain by penalised volume, by the same quantity of derivate
                                grad_C_coef_cp_matrix = grad_C_coef_cp[a+cont,:]*np.ones((len(strain_M),1))
                                grad_C_coef_w_matrix = grad_C_coef_w[a+cont,:]*np.ones((len(strain_M),1))
                                strain_e = (2*strain_M[:,i]*strain_M[:,j]).reshape((len(strain_M),1))
                                
                                grad_c_cp_matrix = grad_c_cp_matrix + strain_e*grad_C_coef_cp_matrix
                                grad_c_w_matrix = grad_c_w_matrix + strain_e*grad_C_coef_w_matrix
                                cont += 1 
                        
                        # Total gradient
                        grad_c_cp = np.sum(-grad_c_cp_matrix*(strain_M[:,-1]*rho_e_M_temp**p_c).reshape((len(strain_M),1)),axis=0).reshape(len(IND_mask),1)
                        grad_c_w = np.sum(-grad_c_w_matrix*(strain_M[:,-1]*rho_e_M_temp**p_c).reshape((len(strain_M),1)),axis=0).reshape(len(IND_mask),1)
                                                                      
                        if flag_scale == 'micro':
                            grad_c_T = np.concatenate((grad_c_cp, grad_c_w))
                        else: 
                            # 4 - Derivatives of the Macro scale Compliance respecting to macro scale design variables
                            der_CP, der_W, BF_mask = der_NURBS(local_support_M,BF_support_M,IND_mask_active_M,IND_mask_M,IND_mask_M_tot,P_rho_M,W_M,rho_e_M) 
                            
                            c_vec_macro_temp = c_vec_macro[BF_mask].reshape((len(BF_mask),1))
                            grad_c_M_cp = np.sum(-(sym_coefM * p_c * c_vec_macro_temp * der_CP)/rho_e_M[BF_mask],axis=0).reshape(len(IND_mask_M),1)
                            grad_c_M_w = np.sum(-(sym_coefM * p_c * c_vec_macro_temp * der_W)/rho_e_M[BF_mask],axis=0).reshape(len(IND_mask_M),1)

                            #print('micro cp','{:.2E}'.format(LA.norm(grad_c_cp)),'micro w','{:.2E}'.format(LA.norm(grad_c_w)),'macro cp','{:.2E}'.format(LA.norm(grad_c_M_cp)),'macro w','{:.2E}'.format(LA.norm(grad_c_M_w)))
                            grad_c_T = np.concatenate((grad_c_cp, grad_c_w, grad_c_M_cp, grad_c_M_w))
                    else:
                        # Derivatives of compliance respecting to standard scale design variables                        
                        der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)                       
                        

                        #c_vec_temp = c_vec[BF_mask].reshape((len(BF_mask),1))
                        c_vec_temp = c_vec[np.int_(ELEMENTS[:,0])-1].reshape((len(BF_mask),1))
                        #print('c_vec_temp dimension',c_vec_temp.shape)
                        #stoppa
                        #TODO
                        if flag_penalisation == 'SIMP':
                            if flag_BCs == 'mixed-w':                                                    
                                grad_c_cp = np.sum((sym_coef * p_c * c_vec_temp * der_CP)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum((sym_coef * p_c * c_vec_temp * der_W)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                            else:                         
                                grad_c_cp = np.sum(-(sym_coef * p_c * c_vec_temp * der_CP)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum(-(sym_coef * p_c * c_vec_temp * der_W)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                        else:
                            if flag_BCs == 'mixed-w':                                                    
                                grad_c_cp = np.sum((sym_coef * c_vec_temp * der_CP),axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum((sym_coef * c_vec_temp * der_W),axis=0).reshape(len(IND_mask),1)
                            else:                         
                                grad_c_cp = np.sum(-(sym_coef * c_vec_temp * der_CP),axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum(-(sym_coef * c_vec_temp * der_W),axis=0).reshape(len(IND_mask),1)
                            
                        grad_c_T = np.concatenate((grad_c_cp, grad_c_w))   
                                        
                else: # B-spline surfaces        
                    if flag_scale == 'micro_macro' or flag_scale == 'micro':

                        # 1 - Derivatives of topology descriptor  
                        BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)
                        if flag_shepard == 1:
                            BF_support_temp = shepard_support(BF_support_temp,shepard_dist)

                        # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                        # it contributes to the c_vec total but the punctual contribute has to be neglected
                        c_vec_micro = c_vec_micro[np.int_(ELEMENTS[:,0])-1]

                        # 2 - Derivatives of Stiffness Matrix coefficients (anisotropic, orthotropic)
                        grad_C_coef = np.zeros((size_grad_C,len(IND_mask)),'f')

                        ### Ajout - Construction intermédiaire du produit terme à terme
                        for a in range(3):
                            c_vec_micro_temp = c_vec_micro[:,a].reshape((len(c_vec_micro),1))

                            lignes = np.shape(BF_support_temp)[0]
                            colonnes = np.shape(BF_support_temp)[1]
                            intermediar = sc.lil_matrix((lignes, colonnes))
                            for i in range(lignes):
                                intermediar[i,:] = c_vec_micro_temp[i,0] * BF_support_temp[i,:]

                            grad_C_coef[a,:] = np.sum(sym_coef*p_c*intermediar/(a1*a2*rho_e),axis=0)
                            #grad_C_coef[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*BF_support_temp/(a1*a2*rho_e),axis=0)


                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                comp_temp = (c_vec_micro[:,a+cont]-c_vec_micro[:,i]-c_vec_micro[:,j]).reshape((len(c_vec_micro),1))

                                ### Ajout - Construction intermédiaire du produit terme à terme
                                lignes = np.shape(BF_support_temp)[0]
                                colonnes = np.shape(BF_support_temp)[1]
                                intermediar = sc.lil_matrix((lignes, colonnes))

                                for k in range(lignes):
                                    intermediar[k,:] = comp_temp[k,0] * BF_support_temp[k,:]

                                grad_C_coef[a+cont,:] = np.sum(sym_coef*p_c*intermediar/(2*a1*a2*rho_e),axis=0)   
                                #grad_C_coef[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*BF_support_temp/(2*a1*a2*rho_e),axis=0)
                                cont += 1

                        # 3 - Derivatives of the Macro scale Compliance               
                        if flag_scale == 'micro':
                            # Topology descriptor is defined just at micro scale
                            rho_e_M_temp = 1
                        else:
                            # Topology descriptor is defined also at macro scale
                            rho_e_M_temp = rho_e_M.reshape((len(rho_e_M),))
                            strain_M = strain_M[np.int_(ELEMENTS_macro[:, 0]) - 1]

                        # Dependencies on the macroscopic strains and on the derivatives of the Stiffness Matrix coefficients
                        grad_c_matrix = np.zeros((len(strain_M),len(IND_mask)))

                        # Evaluating the contributes of the gradient of the diagonal Stiffness Matrix coefficients
                        for a in range(3): 
                            # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                            # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                            grad_C_coef_matrix = grad_C_coef[a,:]*np.ones((len(strain_M),1))
                            strain_e = (strain_M[:,a]**2).reshape((len(strain_M),1))
                            grad_c_matrix = grad_c_matrix + grad_C_coef_matrix*strain_e

                        # Evaluating the contributes of the gradient of the out of diagonal Stiffness Matrix coefficients
                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                                # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                                grad_C_coef_matrix = grad_C_coef[a+cont,:]*np.ones((len(strain_M),1))
                                strain_e = (2*strain_M[:,i]*strain_M[:,j]).reshape((len(strain_M),1))
                                grad_c_matrix = grad_c_matrix + strain_e*grad_C_coef_matrix
                                cont += 1

                        # Total gradient                
                        grad_c = np.sum(-grad_c_matrix*(strain_M[:,-1]*rho_e_M_temp**p_c).reshape((len(strain_M),1)),axis=0).reshape(len(IND_mask),1)

                        if flag_scale == 'micro':
                            grad_c_T = grad_c
                        else:                                           
                            # 4 - Derivatives of macroscale compliance respecting to macro scale design variables
                            BF_support_temp = der_BSPLINE(IND_mask_active_M,BF_support_M)
                            if flag_shepard == 1:
                                BF_support_temp = shepard_support(BF_support_temp,shepard_dist)
                            
                            # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                            # it contributes to the c_vec total but the punctual contribute has to be neglected
                            c_vec_macro = c_vec_macro[np.int_(ELEMENTS_macro[:,0])-1]                            
                            c_vec_macro_temp = c_vec_macro.reshape((len(c_vec_macro),1))

                            ### Ajout - Construction intermédiaire du produit terme à terme
                            lignes = np.shape(BF_support_temp)[0]
                            colonnes = np.shape(BF_support_temp)[1]
                            intermediar = sc.lil_matrix((lignes, colonnes))

                            for i in range(lignes):
                                    intermediar[i,:] = c_vec_macro_temp[i,0] * BF_support_temp[i,:]

                            grad_c_M = np.sum(-sym_coefM*p_c*intermediar/rho_e_M,axis=0).reshape(len(IND_mask_M),1)    
                            #grad_c_M = np.sum(-sym_coefM*p_c*c_vec_macro_temp*BF_support_temp/rho_e_M,axis=0).reshape(len(IND_mask_M),1)  

                            #print('micro','{:.2E}'.format(LA.norm(grad_c)),'macro','{:.2E}'.format(LA.norm(grad_c_M)))                              
                            grad_c_T = np.concatenate((grad_c,grad_c_M))
                    else:
                        # Derivatives of compliance respecting to standard scale design variables
                        BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)
                        if flag_shepard == 1:
                            BF_support_temp = shepard_support(BF_support_temp,shepard_dist)
                                                           
                        # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                        # it contributes to the c_vec total but the punctual contribute has to be neglected
                        c_vec = c_vec[np.int_(ELEMENTS[:,0])-1]                        
                        c_vec_temp = c_vec.reshape((len(c_vec),1))
                        #print('c_vec_temp dimension',c_vec_temp.shape)
                        #stoppa
                        #TODO

                        ### Ajout - Construction intermédiaire du produit terme à terme
                        lignes = np.shape(BF_support_temp)[0]
                        colonnes = np.shape(BF_support_temp)[1]
                        intermediar = sc.lil_matrix((lignes, colonnes))

                        for i in range(lignes):
                                intermediar[i,:] = c_vec_temp[i,0] * BF_support_temp[i,:]


                        if flag_penalisation == 'SIMP':
                            if flag_BCs == 'mixed-w':
                                grad_c = np.sum(sym_coef*p_c*intermediar/rho_e,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(sym_coef*p_c*c_vec_temp*BF_support_temp/rho_e,axis=0).reshape(len(IND_mask),1)
                            else: 
                                grad_c = np.sum(-sym_coef*p_c*intermediar/rho_e,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(-sym_coef*p_c*c_vec_temp*BF_support_temp/rho_e,axis=0).reshape(len(IND_mask),1)
                        else:
                            if flag_BCs == 'mixed-w':
                                grad_c = np.sum(sym_coef*intermediar,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(sym_coef*c_vec_temp*BF_support_temp,axis=0).reshape(len(IND_mask),1)
                            else: 
                                grad_c = np.sum(-sym_coef*intermediar,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(-sym_coef*c_vec_temp*BF_support_temp,axis=0).reshape(len(IND_mask),1)
                        grad_c_T = grad_c
            #--------------------------------------------------------------------------                                                
            elif DIM == 3:
                #--------------------------------------------------------------------------
                if NURBS == 1:                 
                    if flag_scale == 'micro_macro' or flag_scale == 'micro':
                        # 1 - Derivatives of topology descriptor  
                        der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)
                        
                        ### Ajout - Passage des variables CSR suivantes en format array
                        der_CP = der_CP.toarray()
                        der_W = der_W.toarray()
                        BF_mask = BF_mask.toarray()
                        
                        # 2 - Derivatives of Stiffness Matrix coefficients (anisotropic, orthotropic)                        
                        grad_C_coef_cp = np.zeros((size_grad_C,len(IND_mask)),'f')  
                        grad_C_coef_w = np.zeros((size_grad_C,len(IND_mask)),'f')                                                                       
                        for a in range(6):
                            c_vec_micro_temp = c_vec_micro[BF_mask,a].reshape((len(BF_mask),1))
                            grad_C_coef_cp[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*der_CP/(a1*a2*a3*rho_e[BF_mask]),axis=0)
                            grad_C_coef_w[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*der_W/(a1*a2*a3*rho_e[BF_mask]),axis=0)
                        
                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                comp_temp = (c_vec_micro[BF_mask,a+cont]-c_vec_micro[BF_mask,i]-c_vec_micro[BF_mask,j]).reshape((len(BF_mask),1))
                                grad_C_coef_cp[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*der_CP/(2*a1*a2*a3*rho_e[BF_mask]),axis=0)   
                                grad_C_coef_w[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*der_W/(2*a1*a2*a3*rho_e[BF_mask]),axis=0)
                                cont += 1      

                        # 3 - Derivatives of the Macro scale Compliance                        
                        if flag_scale == 'micro':
                            # Topology descriptor is defined just at micro scale
                            rho_e_M_temp = 1
                        else:
                            # Topology descriptor is defined also at macro scale
                            rho_e_M_temp = rho_e_M.reshape((len(rho_e_M),))
                            strain_M = strain_M[np.int_(ELEMENTS_macro[:, 0]) - 1]

                        # Dependencies on the macroscopic strains and on the derivatives of the Stiffness Matrix coefficients
                        grad_c_cp_matrix = np.zeros((len(strain_M),len(IND_mask)))
                        grad_c_w_matrix = np.zeros((len(strain_M),len(IND_mask)))
                        
                        # Evaluating the contributes of the gradient of the diagonal Stiffness Matrix coefficients
                        for a in range(6): 
                            # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                            # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                            grad_C_coef_cp_matrix = grad_C_coef_cp[a,:]*np.ones((len(strain_M),1))
                            grad_C_coef_w_matrix = grad_C_coef_w[a,:]*np.ones((len(strain_M),1))
                            strain_e = (strain_M[:,a]**2).reshape((len(strain_M),1))
                            grad_c_cp_matrix = grad_c_cp_matrix + grad_C_coef_cp_matrix*strain_e
                            grad_c_w_matrix = grad_c_w_matrix + grad_C_coef_w_matrix*strain_e
                        

                        
                        # Evaluating the contributes of the gradient of the out of diagonal Stiffness Matrix coefficients
                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                                # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                                grad_C_coef_cp_matrix = grad_C_coef_cp[a+cont,:]*np.ones((len(strain_M),1))
                                grad_C_coef_w_matrix = grad_C_coef_w[a+cont,:]*np.ones((len(strain_M),1))
                                strain_e = (2*strain_M[:,i]*strain_M[:,j]).reshape((len(strain_M),1))
                                
                                grad_c_cp_matrix = grad_c_cp_matrix + strain_e*grad_C_coef_cp_matrix
                                grad_c_w_matrix = grad_c_w_matrix + strain_e*grad_C_coef_w_matrix
                                cont += 1
                                
                        # Total gradient                 
                        vol_pen = (strain_M[:,-1]*rho_e_M_temp**p_c).reshape((len(strain_M),1))
                        grad_c_cp = np.sum(-grad_c_cp_matrix*vol_pen,axis=0).reshape(len(IND_mask),1)
                        grad_c_w = np.sum(-grad_c_w_matrix*vol_pen,axis=0).reshape(len(IND_mask),1)                                                       
                            
                        if flag_scale == 'micro':
                            grad_c_T = np.concatenate((grad_c_cp, grad_c_w))
                        else:                             
                            # Derivatives of compliance respecting to macro scale design variables
                            der_CP, der_W, BF_mask = der_NURBS(local_support_M,BF_support_M,IND_mask_active_M,IND_mask_M,IND_mask_M_tot,P_rho_M,W_M,rho_e_M)
                                                       
                            c_vec_macro_temp = c_vec_macro[BF_mask].reshape((len(BF_mask),1))
                            grad_c_M_cp = np.sum(-(sym_coefM * p_c * c_vec_macro_temp * der_CP)/rho_e_M[BF_mask],axis=0).reshape(len(IND_mask_M),1)
                            grad_c_M_w = np.sum(-(sym_coefM * p_c * c_vec_macro_temp * der_W)/rho_e_M[BF_mask],axis=0).reshape(len(IND_mask_M),1)
                            
                            #print('{:.2E}      {:.2E}      {:.2E}      {:.2E}'.format(LA.norm(grad_c_cp),LA.norm(grad_c_w),LA.norm(grad_c_M_cp),LA.norm(grad_c_M_w)))          
                            grad_c_T = np.concatenate((grad_c_cp, grad_c_w, grad_c_M_cp, grad_c_M_w))     

                    else:
                        # Derivatives of compliance respecting to standard scale design variables
                        der_CP, der_W, BF_mask = der_NURBS(local_support,BF_support,IND_mask_active,IND_mask,IND_mask_tot,P_rho,W,rho_e)
                                           
                        c_vec_temp = c_vec[BF_mask].reshape((len(BF_mask),1)) 
                        #TODO
                        if flag_penalisation == 'SIMP':
                            if flag_BCs == 'mixed-w':                                                       
                                grad_c_cp = np.sum((sym_coef * p_c * c_vec_temp * der_CP)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum((sym_coef * p_c * c_vec_temp * der_W)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                            else:                         
                                grad_c_cp = np.sum(-(sym_coef * p_c * c_vec_temp * der_CP)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum(-(sym_coef * p_c * c_vec_temp * der_W)/rho_e[BF_mask],axis=0).reshape(len(IND_mask),1)
                        else:
                            if flag_BCs == 'mixed-w':                                                       
                                grad_c_cp = np.sum((sym_coef * c_vec_temp * der_CP),axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum((sym_coef * c_vec_temp * der_W),axis=0).reshape(len(IND_mask),1)
                            else:                         
                                grad_c_cp = np.sum(-(sym_coef * c_vec_temp * der_CP),axis=0).reshape(len(IND_mask),1)
                                grad_c_w = np.sum(-(sym_coef * c_vec_temp * der_W),axis=0).reshape(len(IND_mask),1)                                  
                        grad_c_T = np.concatenate((grad_c_cp, grad_c_w))                                                                                 
                         
                else: #B-spline hypersurfaces  
                    if flag_scale == 'micro_macro' or flag_scale == 'micro':
                        # 1 - Derivatives of topology descriptor 
                        BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)
                        

                        # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                        # it contributes to the c_vec total but the punctual contribute has to be neglected
                        c_vec_micro = c_vec_micro[np.int_(ELEMENTS[:,0])-1]

                        # 2 - Derivatives of Stiffness Matrix coefficients (anisotropic, orthotropic)
                        grad_C_coef = np.zeros((size_grad_C,len(IND_mask)),'f')  

                        ### Ajout - Construction intermédiaire du produit terme à terme
                        for a in range(6):
                            c_vec_micro_temp = c_vec_micro[:,a].reshape((len(c_vec_micro),1))

                            lignes = np.shape(BF_support_temp)[0]
                            colonnes = np.shape(BF_support_temp)[1]
                            intermediar = sc.lil_matrix((lignes, colonnes))

                            for i in range(lignes):
                                intermediar[i,:] = c_vec_micro_temp[i,0] * BF_support_temp[i,:]


                            grad_C_coef[a,:] = np.sum(sym_coef*p_c*intermediar/(a1*a2*a3*rho_e),axis=0)
                            #grad_C_coef[a,:] = np.sum(sym_coef*p_c*c_vec_micro_temp*BF_support_temp/(a1*a2*a3*rho_e),axis=0)


                        cont = 1 
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                comp_temp = (c_vec_micro[:,a+cont]-c_vec_micro[:,i]-c_vec_micro[:,j]).reshape((len(c_vec_micro),1))

                                ### Ajout - Construction intermédiaire du produit terme à terme
                                lignes = np.shape(BF_support_temp)[0]
                                colonnes = np.shape(BF_support_temp)[1]
                                intermediar = sc.lil_matrix((lignes, colonnes))

                                for k in range(lignes):
                                    intermediar[k,:] = comp_temp[k,0] * BF_support_temp[k,:]

                                grad_C_coef[a+cont,:] = np.sum(sym_coef*p_c*intermediar/(2*a1*a2*a3*rho_e),axis=0)
                                #grad_C_coef[a+cont,:] = np.sum(sym_coef*p_c*comp_temp*BF_support_temp/(2*a1*a2*a3*rho_e),axis=0)     
                                cont += 1



                        # 3 - Derivatives of the Macro scale Compliance
                        grad_c = np.zeros((len(IND_mask),1),'f')                    
                        if flag_scale == 'micro':
                            # Topology descriptor is defined just at micro scale
                            rho_e_M_temp = 1
                        else:
                            # Topology descriptor is defined also at macro scale
                            rho_e_M_temp = rho_e_M.reshape((len(rho_e_M),))
                            strain_M = strain_M[np.int_(ELEMENTS_macro[:, 0]) - 1]
                        
                        # Dependencies on the macroscopic strains and on the derivatives of the Stiffness Matrix coefficients
                        grad_c_matrix = np.zeros((len(strain_M),len(IND_mask)))
                        
                        # Evaluating the contributes of the gradient of the diagonal Stiffness Matrix coefficients
                        for a in range(6): 
                            # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                            # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                            grad_C_coef_matrix = grad_C_coef[a,:]*np.ones((len(strain_M),1))
                            strain_e = (strain_M[:,a]**2).reshape((len(strain_M),1))
                            grad_c_matrix = grad_c_matrix + grad_C_coef_matrix*strain_e
                         
                        # Evaluating the contributes of the gradient of the out of diagonal Stiffness Matrix coefficients
                        cont = 1
                        for i in range(i_range_end):
                            for j in range(i+1,j_range_end):
                                # Step to obtain a matrix (n_elem,n_ind_mask) of the row vector gradient of the Stiffness Matrix coefficients (1,n_ind_mask)  
                                # repeted n_elem times, to multiply each column vector strain (n_elem,1) contribute, by the same quantity of derivate
                                grad_C_coef_matrix = grad_C_coef[a+cont,:]*np.ones((len(strain_M),1))
                                strain_e = (2*strain_M[:,i]*strain_M[:,j]).reshape((len(strain_M),1))
                                grad_c_matrix = grad_c_matrix + strain_e*grad_C_coef_matrix
                                cont += 1
                                
                        # Total gradient
                        vol_pen = (strain_M[:,-1]*rho_e_M_temp**p_c).reshape((len(strain_M),1))
                        grad_c = np.sum(-grad_c_matrix*vol_pen,axis=0).reshape(len(IND_mask),1)
                            
                        if flag_scale == 'micro':
                            grad_c_T = grad_c
                            
                        else: 
                            # Derivatives of compliance respecting to macro scale design variables
                            BF_support_temp = der_BSPLINE(IND_mask_active_M,BF_support_M)
                            
                            # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                            # it contributes to the c_vec total but the punctual contribute has to be neglected
                            c_vec_macro = c_vec_macro[np.int_(ELEMENTS_macro[:,0])-1] 
                            c_vec_macro_temp = c_vec_macro.reshape((len(c_vec_macro),1))

                            ### Ajout - Construction intermédiaire du produit terme à terme
                            lignes = np.shape(BF_support_temp)[0]
                            colonnes = np.shape(BF_support_temp)[1]
                            intermediar = sc.lil_matrix((lignes, colonnes))

                            for i in range(lignes):
                                    intermediar[i,:] = c_vec_macro_temp[i,0] * BF_support_temp[i,:]

                            grad_c_M = np.sum(-(sym_coefM * p_c * intermediar)/rho_e_M,axis=0).reshape(len(IND_mask_M),1)
                            #grad_c_M = np.sum(-(sym_coefM * p_c * c_vec_macro_temp*BF_support_temp)/rho_e_M,axis=0).reshape(len(IND_mask_M),1)
                            
                            #print('{:.2E}      {:.2E}'.format(LA.norm(grad_c),LA.norm(grad_c_M)))           
                            grad_c_T = np.concatenate((grad_c,grad_c_M))
                            
                    else:
                        # Derivatives of compliance respecting to standard scale design variables
                        BF_support_temp = der_BSPLINE(IND_mask_active,BF_support)
                        
                            
                        # If we have non design region the compliance of the elements belonging to that region are stocked in c_vec, 
                        # it contributes to the c_vec total but the punctual contribute has to be neglected
                        c_vec = c_vec[np.int_(ELEMENTS[:,0])-1] 
                        c_vec_temp = c_vec.reshape((len(c_vec),1))
                        #TODO

                        ### Ajout - Construction intermédiaire du produit terme à terme
                        lignes = np.shape(BF_support_temp)[0]
                        colonnes = np.shape(BF_support_temp)[1]
                        intermediar = sc.lil_matrix((lignes, colonnes))

                        for i in range(lignes):
                                intermediar[i,:] = c_vec_temp[i,0] * BF_support_temp[i,:]


                        if flag_penalisation == 'SIMP':
                            if flag_BCs == 'mixed-w':
                                print('******************')
                                print(c_vec_temp*BF_support_temp)
                                grad_c = np.sum(sym_coef*p_c*intermediar/rho_e,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(sym_coef*p_c*c_vec_temp*BF_support_temp/rho_e,axis=0).reshape(len(IND_mask),1)
                            else: 
                                grad_c = np.sum(-sym_coef*p_c*intermediar/rho_e,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(-sym_coef*p_c*c_vec_temp*BF_support_temp/rho_e,axis=0).reshape(len(IND_mask),1)  
                        else:
                            if flag_BCs == 'mixed-w':
                                grad_c = np.sum(sym_coef*intermediar,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(sym_coef*c_vec_temp*BF_support_temp,axis=0).reshape(len(IND_mask),1)
                            else: 
                                grad_c = np.sum(-sym_coef*intermediar,axis=0).reshape(len(IND_mask),1)
                                #grad_c = np.sum(-sym_coef*c_vec_temp*BF_support_temp,axis=0).reshape(len(IND_mask),1)                                              
                        grad_c_T = grad_c               
        #--------------------------------------------------------------------------
        else:
            grad_c_T = False
        #--------------------------------------------------------------------------
        
        
        return grad_c_T
