# --------------------------------------------------------------------------
# - Routine name: Problem_Setting.py
# -
# - Routine description: This script constitutes the only file to be
# - modified according to the problem at hand. It should be carefully filled
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Numpy, OS import:
from math import pi
import numpy as np
import os

# Computation of symmetry coefficient:
from symmetry_coef import symmetry_coef_fun
from anisotropy import anisotropy_fun
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# - Definition of the optimisation problem:
# --------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!! USER-DEPENDENT GLOBAL PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!TO BE MODIFY BEFORE RUNNING OPTIMIZATION!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# User variables defined below global variables in the other routines:
global problem_name, res_dir, DIM, NURBS, NDR, obj_fun, n_constr, con_fun
global lb_nlc, ub_nlc, INIT, x0, lb, ub
global DERIVATIVES, optimoptions, method_optim
global WD, OP, ANSYS_stat, ANSYS_comline_start, ANSYS_comline_opt, ANSYS_comline_post
global obj_count, constr_count, obj_count_grad, constr_count_grad
global scale
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - problem_name: string with problem name (the same of the folder)
problem_name = '3D_cube_DN_C'
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - res_dir: specify the directory name for the folder collecting rqesults
res_dir = '3D_cube_DN_C'
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - DIM: 3 or 2 for 3D or 2D problems
DIM = 3
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - NURBS: 0for BSpline entities, 1 for NURBS entities
NURBS = 0
flag_nurbs = NURBS + 1
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - Boundary conditions (BCs) type for the loaf case (LC): 
# - 'forces' (only forces are applied);
# - 'displacements' (only displacements are applied);
# - 'mixed-w' (both displacements and forces are applied in the same LC, the work of external forces is taken as objective function)
# - 'mixed-c' (both displacements and forces are applied in the same LC, the compliance is taken as objective function)
flag_BCs = 'mixed-c'
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - scale: 'micro', 'micro_macro' or 0 if microscale analysis, macro and microscale optimisation or standard optimisation
scale = 'standard'
anisotropy = 'none'
i_range_end, j_range_end, size_grad_C = anisotropy_fun(DIM, anisotropy)
# - import variables relatives to the micro and macro scale
if scale == 'micro_macro':
    from Problem_Setting_mM import *
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - INITIALISATION for the gradient-based algorithm
# - Initialisation is automatically provided by the algorithm for some
# - standard cases. In particular, an uniform initialisation, consistent
# - with the volume/mass constraint or with the compliance constraint is
# - available. If the user wants to provide customised initialisation, it
# - will be possible in this section:
# - INIT = 0 uniform initialisation
# - INIT = 3 starting point definition from file (set index value)
# - INIT = 4 random starting point definition
#- INIT = 5 starting perforated topology (set n_holes value = 1/2/4 holes 2D - 1/4/8 holes 3D)
INIT = 0
index = 0
n_holes = 0
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - obj_fun: char variable
# - 'compliance': minimum compliance problem
# - 'volume': minimum volume problem
# - 'mass': minimum mass problem
# - 'buckling': maximum buckling problem
obj_fun = 'compliance'
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - n_constr number of topology optimisation constrained
# - quantities
n_constr = 1
# - con_fun: cell variable
# - 'compliance' or 'compliance_macro': compliance constraint at the macroscale or microscale
# - 'volume' or 'volume_macro': volume constraint at the microscale or macroscale
# - 'mass' or 'mass_macro': mass constraint at the microscale or macroscale
con_fun = ['volume']
# - constraint of the optimisation problem
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - lb_nlc: lower bounds for non linear constraints; cell variable with
# - either numerical or 'none' values.
if 'moduli_micro' in con_fun:
    C = fromMod_toC_fun(8888.89, 8888.89, 8888.89, 4583.34, 4583.34, 4583.34, 0.28, 0.28, 0.28, 3)
    C = C[0:3]
elif 'GUEST' in con_fun:
    # - Dmax values
    r_maxi = 12.5
    psi = 0.05
    eta = 1.05
    qq = 14
    mesh_dim = 6.25
# - lb_nlc.append('...') should be used to add every
# - lower bounds of the optimisation problem
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - ub_nlc: upper bounds for non linear constraints; cell variable with
# - either numerical or 'none' values.
length=200	  	  # Size along x [mm]
height=200	  # Size along y [mm]
thick=200		      # Plate thickness [mm]

V_B3=thick*length*height
gamma=0.4

ub_nlc = []
ub_nlc.append((gamma)*V_B3)

lb_nlc = []
lb_nlc.append('none')
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - NDR: 0 or 1 if non-design regions are not or are active
NDR = 0
m_NDR = 0
v_NDR = 0
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Reference setting option
REF = 0
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
#TODO
#Penalisation scheme
#flag_penalisation = 'SIMP' in the case of the SIMP method
#flag_penalisation = 'User_defined' for an arbitrary penalisation scheme
flag_penalisation = 'SIMP'
# - p_c: penalisation factor, used only for the SIMP method
p_c = 3
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
#- a1, a2, a3: domain weight (x direction), height (y direction), thickness (z direction)
#- if normal analysis, a1, a2 and a3 are the dimensions of the envelope volume.
#- if multiscale analysis, a1, a2 and a3 dimensions of the VER and a1_glob, a2_glob and a3_glob the dimensions of the structure at the macroscale
a1 = length
a2 = height
a3 = thick

a1_glob = 0 #is equal to 0 if non multi-scale analysis
a2_glob = 0 #is equal to 0 if non multi-scale analysis
a3_glob = 0 #is equal to 0 if non multi-scale analysis

#- check: if DIM=2, then a3 is forced to be 0 and a warning appears.
#- check: if DIM=3, then a3 must be set through a warning message.
if DIM==2 and a3 != 0:
    print('Two dimensions problem: a3 should be 0. Please check data.')
    a3=0
elif DIM==3 and a3==0:
    a3=input('Three dimensions problem and a3=0, please set a3=');
    a3=int(a3)######Not sure!!!!! Be careful!!!!

#TODO
#Options for plots during post-processing phase
#--------------------------------------------------------------------------
#- x_div, y_div, z_div: discretization in x, y, z direction for NURBS
#- representation. Normally, these quantities are
#- unrelated to the underlying mesh. HOWEVER they must
#- coincide with the mapped mesh when the 'member
#- size' constraint is active.
x_div = 200
y_div = 200
z_div = 200

if DIM==2:
    if scale == 'micro_macro':
        label_x = '$x_1^M [mm]$'
        label_y = '$x_2^M [mm]$'
        label_z = '$rho^M $'
    else:
        label_x = '$x_1 [mm]$'
        label_y = '$x_2 [mm]$'
        label_z = '$rho$'
else:
    if scale == 'micro_macro':
        label_x = r'$x_1^M [mm]$'
        label_y = r'$x_2^M [mm]$'
        label_z = r'$x_3^M [mm]$'
        label_rho = r'$rho^M $'
    else:
        label_x = r'$y_1^m [mm]$'
        label_y = r'$y_2^m [mm]$'
        label_z = r'$y_3^m [mm]$'
        label_rho = r'$rho^m  $'


#- check: if DIM=2, then z_div is forced to be 0 and a warning appears.
#- check: if DIM=3, then z_div must be set through a warning message.
if DIM==2 and z_div != 0:
    print('Two dimensions problem: z_div should be 0. Please check data.')
    z_div=0

elif DIM==3 and z_div==0:
    z_div=input('Three dimensions problem and z_div=0, please set z_div=')
    z_div=int(z_div)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#- p1, p2, p3: NURBS degrees
p1=2; p2=2; p3=2;
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# n1, n2, n3: NURBS control points maximum indeces
n1=27; n2=27; n3=27;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#- SYMMETRY_xy: activation of the symmetry condition with respect to
#- z=a3/2 (only available if DIM==3)
SYMMETRY_xy=0
#- SYMMETRY_xz: activation of the symmetry condition with respect to
#- y=a2/2
SYMMETRY_xz=0
#- SYMMETRY_yz: activation of the symmetry condition with respect to
#- x=a1/2
SYMMETRY_yz=0
if DIM==2 and SYMMETRY_xy==1:
    print('2D problem with in plane symmetry.')
    print('Please choose a symmetry with respect to other planes.')
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------    
#- NURBS_sys: the NURBS local coordinate system corresponds to
#- the global coordinte system if NURBS_sys=0.
NURBS_sys=0
OO_star=[]
P_e_1=np.array([1,2,3])
P_e_2=np.array([3,2,1])
if NURBS_sys==1 and (OO_star==[] or P_e_1==[] or P_e_2==[]):
    print('Local and global coordinates systems do not coincide but either')
    print('the local system centre or unit vectors are not provided.')
    print('Please check data')
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#- DERIVATIVES: 0 or 1 if the gradient for obj and constr are not
#- provided or are provided, respectively
DERIVATIVES = 1
method_optim = 'GCMMA'

#- Optimization Settings
if DERIVATIVES==1: #- Derivatives are provided for the computation
    if method_optim == 'GCMMA':
        optimoptions = {'GradObj':'on','GradConstr':'on','save_x':'on','tol_p':10**-6,'tol_f':10**-6,'max_it':250}
    elif method_optim == 'MMA':
        optimoptions = {'GradObj':'on','GradConstr':'on','save_x':'on','tol_p':10**-16,'tol_f':10**-16,'max_it':300}	        
    elif method_optim == 'SLSQP':
        optimoptions = {'ftol': 1e-6, 'eps' : 1e-6, 'maxiter' :10000, 'disp': True}
        
else: #- Derivatives are not provided for the computation
    if method_optim == 'GCMMA':
        optimoptions = {'GradObj':'off','GradConstr':'off','save_x':'on','tol_p':10**-8,'tol_f':10**-8,'max_it':6000}
    elif method_optim == 'MMA':
        optimoptions = {'GradObj':'off','GradConstr':'off','save_x':'on','tol_p':10**-8,'tol_f':10**-8,'max_it':6000}          
    elif method_optim == 'SLSQP':
        optimoptions = {'ftol': 1e-4, 'eps' : 1e-4, 'maxiter' :1000, 'disp': True}
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!! NO USER-DEPENDENT GLOBAL PARAMETERS !!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!DO NOT MODIFY SCRIPT BELOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# WD = Working Directory
WD = os.getcwd()
OP = 'OP'
RESULTS = 'RESULTS'
ITER = 'ITER'
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# - ANSYS_stat: check variable for ANSYS
# - it is set on the licence error code by default
ANSYS_stat = 100
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# serveur de calcul utilisé: HOME, I2M ou COFFA

op_dir = os.path.join(WD, OP, problem_name)

# --------------------------------------------------------------------------
# Computation of Symmetry coefficient:
sym_coef = symmetry_coef_fun (DIM, SYMMETRY_xy , SYMMETRY_xz, SYMMETRY_yz)
#--------------------------------------------------------------------------
