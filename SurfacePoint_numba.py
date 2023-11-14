# - Routine name: SurfacePoint_fun_numba.m
# -
# - Routine description: modification of the routine found in NURBS Book, Piegl and Tiller
# - This algorithm takes as input:
# - u        the complete vector of the parameters in the 1st par axis (vector (n_elem,1) if the routine is called from pre processing, grid (x_div,y_div) if the routine is called from post processing)
# - n        max index for control points for the B-Spline along the parametric axis u
# - p        the degree of the B-Spline polynomial along the u parametric axis
# - U        the knot vector for the u parametric axis
# - v        the complete vector of the parameters in the 2nd par axis (vector (n_elem,1) if the routine is called from pre processing, grid (x_div,y_div) if the routine is called from post processing)
# - m        max index for control points for the B-Spline along the parametric axis v
# - q        the degree of the B-Spline polynomial along the v parametric axis
# - V        the knot vector for the v parametric axis
# - P        matrix of components of control points along the two assigned directions (dimensions [n+1,m+1])
# - w        matrix of components of weights along the two assigned directions (dimensions [n+1,m+1])
# - The output is S, the vector or grid of B-Spline/NURBS surface points according to "The NURBS book".
import numpy as np
from numba import njit

#from Problem_Setting import *

@njit
def SurfacePoint_fun_numba(u, n, p, U, v, m, q, V, P, w, NURBS):
    
    def find_span_basis(pt, p, n, U):
        #N = np.zeros(p + 1)
        #left = np.zeros(p + 1)
        #right = np.zeros(p + 1)
        N = np.zeros((p + 1), dtype=np.float32) #TODO
        left = np.zeros((p + 1), dtype=np.float32) #TODO
        right = np.zeros((p + 1), dtype=np.float32) #TODO
        N[0] = 1.0
        span = 0
        if pt == U[n + 1]:
            span = n

        elif pt <= U[p + 1]:
            span = p

        else:
            low = p
            high = n + 1
            mid = int((high + low) / 2)

            while (pt < U[mid] or pt >= U[mid + 1]):
                if pt < U[mid]:
                    high = mid
                else:
                    low = mid
                mid = int((low + high) / 2)
            span = mid

        for j in range(1, p + 1):
            left[j] = pt - U[span + 1 - j]
            right[j] = U[span + j] - pt
            saved = 0.0

            for r in range(j):
                temp = N[r] / (right[r + 1] + left[j - r])
                N[r] = saved + right[r + 1] * temp
                saved = left[j - r] * temp
            N[j] = saved
        return span, N

    def SurfacePoint_fun2(u, v, n, m, p, q, U, V, P, uspan, vspan, Nu, Nv):
        uind = int(uspan) - p
        surf = 0.0
        for l in range(q + 1):
            temp = 0.0
            vind = int(vspan) - q + l
            for k in range(p + 1):
                temp = temp + Nu[k] * P[uind + k, vind]
            surf = surf + Nv[l] * temp
        return surf

    size_u_r, size_u_c = np.shape(u)
    if size_u_c == 1:  # when u is passed as a vector
        surf = np.zeros(size_u_r)
        
        if NURBS==1:#np.all(w != 1.):  # NURBS
            p_w = np.zeros((n + 1, m + 1))
            p_w = P * w
            #print('NURBS')
            for i in range(size_u_r):
                uspan, nu = find_span_basis(u[i, 0], p, n, U)
                vspan, nv = find_span_basis(v[i, 0], q, m, V)
                surf_cp = SurfacePoint_fun2(u[i, 0], v[i, 0], n, m, p, q, U, V, p_w, uspan, vspan, nu, nv)
                surf_w = SurfacePoint_fun2(u[i, 0], v[i, 0], n, m, p, q, U, V, w, uspan, vspan, nu, nv)
                surf[i] = surf_cp / surf_w
            #print(surf)
        else:  # BSPLINE
            for i in range(size_u_r):
                uspan, nu = find_span_basis(u[i, 0], p, n, U)
                vspan, nv = find_span_basis(v[i, 0], q, m, V)
                surf[i] = SurfacePoint_fun2(u[i, 0], v[i, 0], n, m, p, q, U, V, P, uspan, vspan, nu, nv)
    
    else:  # when u is passed as a grid
        surf_temp = np.zeros((size_u_r, size_u_c))

        if NURBS==1:#np.all(w != 1.):  # NURBS
            p_w = P * w
            for i in range(size_u_r):
                for j in range(size_u_c):
                    uspan, nu = find_span_basis(u[i, j], p, n, U)
                    vspan, nv = find_span_basis(v[i, j], q, m, V)
                    c_w = SurfacePoint_fun2(u[i, j], v[i, j], n, m, p, q, U, V, w, uspan, vspan, nu, nv)
                    surf_temp[i, j] = SurfacePoint_fun2(u[i, j], v[i, j], n, m, p, q, U, V, p_w, uspan, vspan, nu, nv) / c_w
            surf = np.reshape(surf_temp, (size_u_r * size_u_c,))

        else:  # BSPLINE
            for i in range(size_u_r):
                for j in range(size_u_c):
                    uspan, nu = find_span_basis(u[i, j], p, n, U)
                    vspan, nv = find_span_basis(v[i, j], q, m, V)
                    surf_temp[i, j] = SurfacePoint_fun2(u[i, j], v[i, j], n, m, p, q, U, V, P, uspan, vspan, nu, nv)
            surf = np.reshape(surf_temp, (size_u_r * size_u_c,))

    return surf