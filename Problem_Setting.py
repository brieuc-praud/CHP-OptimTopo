# --------------------------------------------------------------------------
# - Routine name: Problem_Setting.py
# -
# - Routine description: This script constitutes the only file to be
# - modified according to the problem at hand. It should be carefully filled
# -------------------------------------------------------------------------
# --------------------------------------------------------------------------
# User variables defined below global variables in the other routines:
global problem_name, res_dir, DIM, NURBS, NDR

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
length=200	  	  # Size along x [mm]
height=200	  # Size along y [mm]
thick=200		      # Plate thickness [mm]
# --------------------------------------------------------------------------

#- a1, a2, a3: domain width (x direction), height (y direction), thickness (z direction)
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
#--------------------------------------------------------------------------
#- x_div, y_div, z_div: discretization in x, y, z direction for NURBS
#- representation. Normally, these quantities are
#- unrelated to the underlying mesh. HOWEVER they must
#- coincide with the mapped mesh when the 'member
#- size' constraint is active.
x_div = 200
y_div = 200
z_div = 200
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

