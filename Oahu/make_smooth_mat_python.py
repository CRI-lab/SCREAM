# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:25:54 2023
converted from "Make_smooth_mat.m" by Tiffany Anderson.
Original Description: 
% MAKE_SMOOTH_MAT creates the sparse matrix S used to smooth with 
% a length 5 moving average defined as (1/13)*[1 3 5 3 1].  
% Smooths at ends following Ayesha's implementation in smoothst.m file.  
%
% S = SMOOTH_MAT_QUARTER(L) returns the sparse matrix S used 
%   to smooth. Size of S is L x L. Smooths at ends. 
%
% L: Length of column vector to smooth 

% Tiffany Anderson
% Created: 3/27/13

@author: rum
"""

import numpy as np
from scipy.sparse import csr_matrix

def make_smooth_mat(L):
    if L == 1:
        return csr_matrix([1])
    elif L == 2:
        return csr_matrix([[5/8, 3/8], [3/8, 5/8]])
    elif L == 3:
        return csr_matrix([[5/9, 3/9, 1/9], [3/11, 5/11, 3/11], [1/9, 3/9, 5/9]])
    elif L == 4:
        return csr_matrix([[5/9, 3/9, 1/9, 0], [3/12, 5/12, 3/12, 1/12],
                           [1/12, 3/12, 5/12, 3/12], [0, 1/9, 3/9, 5/9]])
    else:
        row1 = np.array([5/9, 3/9, 1/9])
        row2 = np.array([3/12, 5/12, 3/12, 1/12])
        rowL2 = np.flip(row2)
        rowL = np.flip(row1)

        smooth_vec = np.array([1/13, 3/13, 5/13, 3/13, 1/13])

        i = np.concatenate([np.ones(3), 
                            2 * np.ones(4), 
                            np.repeat(np.arange(3, L-1), 5),
                            (L-1)*np.ones(4), 
                            L * np.ones(3)])

        j = np.concatenate([np.arange(1, 4),
                            np.arange(1, 5),
                            np.reshape(np.arange(1, L-3).repeat(5) + np.tile(np.arange(5), L-4), (L-4)*5, order='F'),
                            np.arange(L-3, L+1),
                            np.arange(L-2, L+1)
        ])

        s = np.concatenate([row1, row2, np.tile(smooth_vec, L - 4), rowL2, rowL])

        return csr_matrix((s, (i.astype(int) - 1, j.astype(int) - 1)), shape=(L, L))
