# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'wong2'
startp = np.array([2.0,3.0,5.0,5.0,1.0,2.0,7.0,3.0,6.0,10.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 5
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    x = x.reshape(-1,1)
    F = np.zeros(9)
    F[0] = x[0]**2 + x[1]**2 + x[0]*x[1] - 14.0*x[0]- 16.0*x[1] + (x[2]-10)**2 + 4.0*(x[3]-5.0)**2+(x[4]-3.0)**2 + 2.0*(x[5]-1.0)**2 + 5.0*x[6]**2 + 7.0*(x[7]-11.0)**2+2.0*(x[8]-10)**2 + (x[9]-7.0)**2 + 45.0
    F[1] = F[0] + 10*(3.0*(x[0]-2.0)**2 + 4.0*(x[1]-3.0)**2 + 2.0*x[2]**2 -7.0*x[3]-120.0)
    F[2] = F[0] + 10*(5.0*x[0]**2 + 8.0*x[1] + (x[2]-6.0)**2 - 2.0*x[3] -40.0)
    F[3] = F[0] + 10*(0.5*(x[0]-8.0)**2 + 2.0*(x[1]-4.0)**2 + 3.0*x[4]**2 -x[5]-30.0)
    F[4] = F[0] + 10*(x[0]**2 + 2.0*(x[1]-2.0)**2 - 2.0*x[0]*x[1] + 14.0*x[4] - 6.0*x[5])
    F[5] = F[0] + 10*(-3.0*x[0] + 6.0*x[1] + 12.0*(x[8]-8.0)**2 - 7.0*x[9])
    F[6] = F[0] + 10*(4.0*x[0]+5.0*x[1]-3.0*x[6] +9.0*x[7]-105.0)
    F[7] = F[0] + 10*(10*x[0]-8.0*x[1]-17.0*x[6]+2.0*x[7])
    F[8] = F[0] + 10*(-8.0*x[0]+2.0*x[1]+5.0*x[8] -2.0*x[9]-12.0)
    y = np.max(F)
    return y
