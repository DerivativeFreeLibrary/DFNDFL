############################################################################################
#
#description  (problems for mixed-integer problems with AT LEAST
#              2 discrete variables and 2 continuous variables)
#
#This python module exports the following objects:
#- problem : a python class
#- prob_collection : A dictionary of entries "probname" => problem object
#
# a problem object is a structured type that has the following attributes:
# name   : string - name of the problem
# startp : numpy array - the starting point for the continuous problem
# lb     : numpy array - the lower bounds of the continuous problem
# ub     : numpy array - the upper bounds of the continuous problem
# n      : int - the total number of variables (>= 4)
# nint   : int - the number of discrete variables (>= 2)
# ncont  : int - the number of continuous variables (>= 2), N.B. n = nint+ncont
#        : BEWARE the variables are so intended
#        :        x[0] ... x[ncont-1] are the continuous variables
#        :    x[ncont] ... x[n-1]     are the discrete variables
# lbmix  : numpy array - the actual lower bounds of the mixed integer problem
# ubmix  : numpy array - the actual upper bounds of the mixed integer problem
# x_initial : numpy array - the actual initial point of the mixed integer problem
# xmix   : numpy array - a temporary array used for computation
# feval  : function handle - function to compute the objective function value
#        : N.B. the point must be reconstructed through the use of reconstruct_xmix
#        :      before calling feval!
# m      : dictionary with entries char => number of constraints
#        : it is equal to {'a': n-2, 'b': n-2, 'c': n-1, 'd': n-1, 'e': n-2, 'f': 1, 'z': 0}
#        : it is used to record the number of constraints for the given problem and for
#        : each of the six families of constraints a,b,c,d,e,f. 'z' means no constraints
#
############################################################################################
import numpy as np
import importlib

prob_collection = {}
list_prob_names = ["prob210", "kowalik"]

class problem:
    name = ""
    startp = None
    lb = None
    ub = None
    n = 0
    nint = 0
    ncont = 0
    lbmix = None
    ubmix = None
    x_initial = None
    xmix = None
    feval = None
    m = {}

    def reconstruct_xmix(self,x):
        self.xmix[self.ncont:] = self.lb[self.ncont:] + ((self.ub[self.ncont:] - self.lb[self.ncont:])/(self.ubmix[self.ncont:] - self.lbmix[self.ncont:]))*(x[self.ncont:]-self.lbmix[self.ncont:])
        self.xmix[:self.ncont] = x[:self.ncont]

    def func_f(self, x):
        self.reconstruct_xmix(x)
        return self.feval(self.xmix)

    def fconstr_a(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return  (3-2*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] - 2*self.xmix[J+2] + 1

    def fconstr_b(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return  (3-2*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] - 2*self.xmix[J+2] + 2.5

    def fconstr_c(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-1)
        return self.xmix[J]**2 + self.xmix[J+1]**2 + self.xmix[J]*self.xmix[J+1] - 2*self.xmix[J] - 2*self.xmix[J+1] + 1

    def fconstr_d(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-1)
        return self.xmix[J]**2 + self.xmix[J+1]**2 + self.xmix[J]*self.xmix[J+1] - 1

    def fconstr_e(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return (3-0.5*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] -2*self.xmix[J+2] +1

    def fconstr_f(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return np.array([np.sum((3-0.5*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] -2*self.xmix[J+2] +1)])

    def fconstr_z(self,x):
        return np.array([-1.0])

    def __init__(self,name,startp,lb,ub,n,nint,m,ncont,lbmix,ubmix,x_initial,xmix,feval):
        self.name = name
        self.startp = startp
        self.lb = lb
        self.ub = ub
        self.n = n
        self.nint = nint
        self.m = m
        self.ncont = ncont
        self.lbmix = lbmix
        self.ubmix = ubmix
        self.x_initial = x_initial
        self.xmix = xmix
        self.feval = feval

for pname in list_prob_names:
    mname = importlib.import_module(pname)
    if mname.n >= 3:
        prob_collection[pname] = problem(
            name   = mname.name,
            startp = mname.startp,
            lb     = mname.lb,
            ub     = mname.ub,
            n      = mname.n,
            nint   = mname.nint,
            m      = {'a': mname.n-2,
                      'b': mname.n-2,
                      'c': mname.n-1,
                      'd': mname.n-1,
                      'e': mname.n-2,
                      'f': 1,
                      'z': 0},
            ncont  = mname.ncont,
            lbmix  = mname.lbmix,
            ubmix  = mname.ubmix,
            x_initial = mname.x_initial,
            xmix   = mname.xmix,
            feval  = mname.feval
        )
    if mname.n <= 2:
        prob_collection[pname] = problem(
            name   = mname.name,
            startp = mname.startp,
            lb     = mname.lb,
            ub     = mname.ub,
            n      = mname.n,
            nint   = mname.nint,
            m      = {'c': mname.n-1,
                      'd': mname.n-1,
                      'f': 1,
                      'z': 0},
            ncont  = mname.ncont,
            lbmix  = mname.lbmix,
            ubmix  = mname.ubmix,
            x_initial = mname.x_initial,
            xmix   = mname.xmix,
            feval  = mname.feval
        )
