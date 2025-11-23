# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:42:27 2022

@author: Hannes Grimm-Strele
"""

import numpy as np

A_random = np.random.rand(3,3)
w,R = np.linalg.eig(A_random)
Rt  = np.transpose(R)
w   = w.clip(min=0.01)
w  /= np.linalg.norm(w, 1)
wDiagonal = np.diag(w)

A_positiveDefinite = np.real(np.dot(np.dot(R, wDiagonal), Rt))
print(A_positiveDefinite)
