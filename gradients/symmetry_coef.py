'''--------------------------------------------------------------------------
#- Routine name: Symmetry_coef.py
#- 
#- Routine description: This script computes symmetry coefficient according to
#- symmetry set up inside problem setting
#- Value taken could be : 1, 2, 4 in 2D and 1, 2, 4, 8 in 3D
#--------------------------------------------------------------------------'''

def symmetry_coef_fun (DIM, SYMMETRY_xy, SYMMETRY_xz, SYMMETRY_yz):
    if DIM ==2:
        if SYMMETRY_xz == 0 and SYMMETRY_yz == 0:
            #- Non-Symmetric NURBS surfaces for 2D structures
            sym_coef = 1

        elif (SYMMETRY_xz == 1 and SYMMETRY_yz==0) or (SYMMETRY_xz==0 and SYMMETRY_yz==1):
            #- 2D Symmetric NURBS surfaces with respect to the xz plane
            #- or with respect to the yz plane
            sym_coef = 2
            
        elif SYMMETRY_xz==1 and SYMMETRY_yz==1:
            #- 2D Symmetric NURBS surfaces with respect to both
            #- or xz plane and to the yz plane
            sym_coef = 4

    if DIM ==3:
        if  SYMMETRY_xy==0 and SYMMETRY_xz == 0 and SYMMETRY_yz == 0:
            #- Non-Symmetric NURBS hypersurfaces for 3D structures
            sym_coef = 1
            
        elif ((SYMMETRY_xy==1 and SYMMETRY_xz==0 and SYMMETRY_yz==0)or
              (SYMMETRY_xy==0 and SYMMETRY_xz==1 and SYMMETRY_yz==0) or
              (SYMMETRY_xy==0 and SYMMETRY_xz==0 and SYMMETRY_yz==1)):
            #-3D Symmetric NURBS hypersurfaces with respect to one
            #-among the planes x=a1/2, y=a2/2, z=a3/2
            sym_coef = 2
            
        elif ((SYMMETRY_xy==1 and SYMMETRY_xz==1 and SYMMETRY_yz==0)or
              (SYMMETRY_xy==1 and SYMMETRY_xz==0 and SYMMETRY_yz==1) or
              (SYMMETRY_xy==0 and SYMMETRY_xz==1 and SYMMETRY_yz==1)):
            #-3D Symmetric NURBS hypersurfaces with respect to two
            #-among the planes x=a1/2, y=a2/2, z=a3/2
            sym_coef = 4
            
        elif (SYMMETRY_xy==1 and SYMMETRY_xz==1 and SYMMETRY_yz==1):
            #-3D Symmetric NURBS hypersurfaces with respect to two
            #-among the planes x=a1/2, y=a2/2, z=a3/2
            sym_coef = 8
            
    return sym_coef
