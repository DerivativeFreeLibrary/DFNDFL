# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'wong3'
startp = np.array([2.0,3.0,5.0,5.0,1.0,2.0,7.0,3.0,6.0,10.0,2.0,2.0,6.0,15.0,1.0,2.0,1.0,2.0,1.0,3.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 10
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    x = x.reshape(-1,1)
    f = np.zeros(18)
    f[0] = x[0]**2 + x[1]**2 + x[0]*x[1] - 14.0*x[0] - 16.0*x[1] + (x[2] - 10.0)**2  \
         + 4.0*(x[3] - 5.0)**2 + (x[4] - 3.0)**2 + 2.0*(x[5] - 1.0)**2 + 5.0*x[6]**2 \
         + 7.0*(x[7] - 11.0)**2 + 2.0*(x[8] - 10.0)**2 + (x[9] - 7.0)**2 + (x[10] - 9.0)**2 \
         + 10.0*(x[11] - 1.0)**2 + 5.0*(x[12] - 7.0)**2 + 4.0*(x[13] - 14.0)**2  \
         + 27.0*(x[14] - 1.0)**2 + x[15]**4 + (x[16] - 2.0)**2 + 13.0*(x[17] - 2.0)**2 \
         + (x[18] - 3.0)**2 + x[19]**2 + 95.0
    f[1] = f[0] + 10.0*(3.0*(x[0] - 2.0)**2 + 4.0*(x[1] - 3.0)**2 + 2.0*x[2]**2 - 7.0*x[3] - 120.0)
    f[2] = f[0] + 10.0*(5.0*x[0]**2 + 8.0*x[1] + (x[2] - 6.0)**2 - 2.0*x[3] - 40.0)
    f[3] = f[0] + 10.0*(0.5*(x[0] - 8.0)**2 + 2.0*(x[1] - 4.0)**2 + 3.0*x[4]**2 - x[5] - 30.0)
    f[4] = f[0] + 10.0*(x[0]**2 + 2.0*(x[1] - 2.0)**2 - 2.0*x[0]*x[1] + 14.0*x[4] - 6.0*x[5])
    f[5] = f[0] + 10.0*(4.0*x[0] + 5.0*x[1] - 3.0*x[6] + 9.0*x[7] - 105.0)
    f[6] = f[0] + 10.0*(10.0*x[0] - 8.0*x[1] - 17.0*x[6] + 2.0*x[7])
    f[7] = f[0] + 10.0*(-3.0*x[0] + 6.0*x[1] + 12.0*(x[8] - 8.0)**2 - 7.0*x[9])
    f[8] = f[0] + 10.0*(-8.0*x[0] + 2.0*x[1] + 5.0*x[8] - 2.0*x[9] - 12.0)
    f[9]= f[0] + 10.0*(x[0] + x[1] + 4.0*x[10] - 21.0*x[11])
    f[10]= f[0] + 10.0*(x[0]**2 + 5.0*x[10] - 8.0*x[11] - 28.0)
    f[11]= f[0] + 10.0*(4.0*x[0] + 9.0*x[1] + 5.0*x[12]**2 - 9.0*x[13] - 87.0)
    f[12]= f[0] + 10.0*(3.0*x[0] + 4.0*x[1] + 3.0*(x[12] - 6.0)**2 - 14.0*x[13] - 10.0)
    f[13]= f[0] + 10.0*(14.0*x[0]**2 + 35.0*x[14] - 79.0*x[15] - 92.0)
    f[14]= f[0] + 10.0*(15.0*x[1]**2 + 11.0*x[14] - 61.0*x[15] - 54.0)
    f[15]= f[0] + 10.0*(5.0*x[0]**2 + 2.0*x[1] + 9.0*x[16]**4 - x[17] - 68.0)
    f[16]= f[0] + 10.0*(x[0]**2 - x[1] + 19.0*x[18] - 20.0*x[19] + 19.0)
    f[17]= f[0] + 10.0*(7.0*x[0]**2 + 5.0*x[1]**2 + x[18]**2 - 30.0*x[19])
    y = np.max(f);
    return y
