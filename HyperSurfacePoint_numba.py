# Routine name: HyperSurfacePoint_fun_numba.m
# - Routine description:
# - This algorithm takes as input: 
# - u1, u2, u3         the value of dimensionless parameters (vector (n_elem,1) if the routine is called from pre processing, grid (x_div,y_div,z_div) if the routine is called from post processing)
# - p1, p2, p3         the degrees of the NURBS polynomial
# - U1, U2, U3         the knot vectors  
# - n1, n2, n3         max index for control points for the NURBS 
# - P        matrix of components of control points along the three assigned
# -          directions (dimensions [n1+1,n2+1,n3+1])
# - W        matrix of components of weights along the three assigned
# -          directions (dimensions [n1+1,n2+1,n3+1])
# - The output is H the vector or grid of B-Spline/NURBS hypersurface points
import numpy as np
from numba import njit

#from Problem_Setting import *

@njit
def HyperSurfacePoint_fun_numba(n1, p, U1, n2, q, U2, n3, r, U3, P, w, u1, u2, u3, NURBS):
    def find_span_basis(pt, p, n, U):
        #bf = np.zeros(p + 1)
        #left = np.zeros(p + 1)
        #right = np.zeros(p + 1)
        bf = np.zeros((p + 1), dtype=np.float32) #TODO
        left = np.zeros((p + 1), dtype=np.float32) #TODO
        right = np.zeros((p + 1), dtype=np.float32) #TODO
        bf[0] = 1.0
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
                temp = bf[r] / (right[r + 1] + left[j - r])
                bf[r] = saved + right[r + 1] * temp
                saved = left[j - r] * temp
            bf[j] = saved
        return span, bf

    def HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, P, u1, u2, u3, u1span, u2span, u3span, Nu1, Nu2, Nu3):
        hypersurf = 0
        u1ind = u1span - p

        for i3 in range(r + 1):
            S = 0
            u3ind = u3span - r + i3
            for i2 in range(q + 1):
                temp = 0
                u2ind = u2span - q + i2
                for i1 in range(p + 1):
                    temp = temp + Nu1[i1] * P[u1ind + i1, u2ind, u3ind]
                S = S + Nu2[i2] * temp
            hypersurf = hypersurf + Nu3[i3] * S
        return hypersurf

    size_u_r, size_u_c, size_u_t = np.shape(u1)
    if size_u_c == 1 and size_u_t == 1:  # when u is passed as a vector
        hypersurf = np.zeros(size_u_r)

        if NURBS==1:#np.all(w != 1.):  # NURBS
            P_w = P * w
            for i in range(size_u_r):
                u1span, nu1 = find_span_basis(u1[i, 0, 0], p, n1, U1)
                u2span, nu2 = find_span_basis(u2[i, 0, 0], q, n2, U2)
                u3span, nu3 = find_span_basis(u3[i, 0, 0], r, n3, U3)
                S_cp = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, P_w, u1[i, 0, 0], u2[i, 0, 0], u3[i, 0, 0],
                                          u1span, u2span, u3span, nu1, nu2, nu3)
                hs_w = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, w, u1[i, 0, 0], u2[i, 0, 0], u3[i, 0, 0],
                                          u1span, u2span, u3span, nu1, nu2, nu3)
                hypersurf[i] = S_cp / hs_w

        else:  # BSPLINE
            for i in range(size_u_r):
                u1span, nu1 = find_span_basis(u1[i, 0, 0], p, n1, U1)
                u2span, nu2 = find_span_basis(u2[i, 0, 0], q, n2, U2)
                u3span, nu3 = find_span_basis(u3[i, 0, 0], r, n3, U3)
                hypersurf[i] = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, P, u1[i, 0, 0], u2[i, 0, 0],
                                                  u3[i, 0, 0], u1span, u2span, u3span, nu1, nu2, nu3)

    else:  # when u is passed as a grid
        hypersurf_temp = np.zeros((size_u_r, size_u_c, size_u_t))

        if NURBS==1:#np.all(w != 1.):  # NURBS
            P_w = P * w
            for i in range(size_u_r):
                for j in range(size_u_c):
                    for k in range(size_u_t):
                        u1span, nu1 = find_span_basis(u1[i, j, k], p, n1, U1)
                        u2span, nu2 = find_span_basis(u2[i, j, k], q, n2, U2)
                        u3span, nu3 = find_span_basis(u3[i, j, k], r, n3, U3)
                        hs_w = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, w, u1[i, j, k], u2[i, j, k],
                                                  u3[i, j, k], u1span, u2span, u3span, nu1, nu2, nu3)
                        hypersurf_temp[i, j, k] = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, P_w, u1[i, j, k],
                                                                     u2[i, j, k], u3[i, j, k], u1span, u2span, u3span,
                                                                     nu1, nu2, nu3) / hs_w
            hypersurf = np.reshape(hypersurf_temp, (size_u_r * size_u_c * size_u_t,))

        else:  # BSPLINE
            for i in range(size_u_r):
                for j in range(size_u_c):
                    for k in range(size_u_t):
                        u1span, nu1 = find_span_basis(u1[i, j, k], p, n1, U1)
                        u2span, nu2 = find_span_basis(u2[i, j, k], q, n2, U2)
                        u3span, nu3 = find_span_basis(u3[i, j, k], r, n3, U3)
                        hypersurf_temp[i, j, k] = HSurfacePoint_fun2(n1, p, U1, n2, q, U2, n3, r, U3, P, u1[i, j, k],
                                                                     u2[i, j, k], u3[i, j, k], u1span, u2span, u3span,
                                                                     nu1, nu2, nu3)
            hypersurf = np.reshape(hypersurf_temp, (size_u_r * size_u_c * size_u_t,))

    return hypersurf
