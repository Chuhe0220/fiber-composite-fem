# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 10:18:33 2021

@author: Hannes Grimm-Strele
"""

import numpy as np
import matplotlib.pyplot as plt

tolerance = 0.01

# discretization of triangle
lambda1 = np.array([1./3., 1./3.])
lambda5 = np.array([0.5 - 0.5*tolerance, 0.5 - 0.5*tolerance])
lambda15 = [1.0 - 2*tolerance, tolerance]

lambda2 = lambda1 + 0.25 * (lambda5 - lambda1)
lambda3 = lambda1 + 0.50 * (lambda5 - lambda1)
lambda4 = lambda1 + 0.75 * (lambda5 - lambda1)
lambda6 = lambda1 + 0.25 * (lambda15 - lambda1)
lambda9 = lambda5 + 0.25 * (lambda15 - lambda5)
lambda10 = lambda1 + 0.50 * (lambda15 - lambda1)
lambda12 = lambda5 + 0.50 * (lambda15 - lambda5)
lambda13 = lambda1 + 0.75 * (lambda15 - lambda1)
lambda14 = lambda5 + 0.75 * (lambda15 - lambda5)

plt.plot([lambda1[0], lambda15[0]], [lambda1[1], lambda15[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda2[0], lambda14[0]], [lambda2[1], lambda14[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda3[0], lambda12[0]], [lambda3[1], lambda12[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda4[0], lambda9[0]], [lambda4[1], lambda9[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda1[0], lambda5[0]], [lambda1[1], lambda5[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda6[0], lambda9[0]], [lambda6[1], lambda9[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda10[0], lambda12[0]], [lambda10[1], lambda12[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda13[0], lambda14[0]], [lambda13[1], lambda14[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda2[0], lambda6[0]], [lambda2[1], lambda6[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda3[0], lambda10[0]], [lambda3[1], lambda10[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda4[0], lambda13[0]], [lambda4[1], lambda13[1]], color = 'k', linewidth = 1, linestyle = 'solid')
plt.plot([lambda5[0], lambda15[0]], [lambda5[1], lambda15[1]], color = 'k', linewidth = 1, linestyle = 'solid')

plt.savefig('OrientationTriangle.png', dpi=300, bbox_inches="tight")
plt.close()