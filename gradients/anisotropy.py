# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:15:41 2020

@author: gbertolino
"""

def anisotropy_fun(DIM,anisotropy):
    if DIM == 2:
        if anisotropy == 'ortho':
            i_range_end = 1
            j_range_end = 2
            size_grad_C = 4
        else:
            i_range_end = 2
            j_range_end = 3
            size_grad_C = 6
            
    elif DIM == 3:
        if anisotropy == 'ortho':
            i_range_end =  2
            j_range_end = 3
            size_grad_C = 9
        else:
            i_range_end = 5
            j_range_end = 6
            size_grad_C = 21
            
    return i_range_end, j_range_end, size_grad_C