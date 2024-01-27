# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:02:10 2020

@author: tomma
"""
import sys
sys.path.append('../PY_PROBLEMS')  #add path to probs.py module
import argparse
import probs
import numpy as np
import mixed_DFN_DFLINT as dfn_dflint

class cache:
    # nnf = number of entries used in numpy array xf
    # xf  = numpy array where computed values are stored
    # hits= number of times item found in cache
    # maxf= maximum number of function values that can be stored in F
    # F   = numpy array of function values used for performance profiles
    nnf  = 0
    hits = 0
    xf   = []
    maxf = 0
    F    = []
    def __init__(self,max_fun,dim):
        self.nnf  = 0
        self.hits = 0
        self.F    = []
        self.xf   = np.inf*np.ones((10*max_fun,dim)) # dim = n+ncont+m+1
        self.maxf = max_fun

parser = argparse.ArgumentParser(prog='run_all',description='Run DFO codes')
parser.add_argument('-v','--version',action='store_true',help='Print version number')
parser.add_argument('--constrained',action='store_true',help='Solve constrained problems')
parser.add_argument('--alg',nargs=1,type=str,
                             default=['DFN_DFL'],
                             choices=["DFN_DFL", "DFL"],
                             help='Name of algo to be run')
parser.add_argument('--max_fun',nargs=1,type=int,
                             default=[5000],
                             help='Maximum number of function avaluations')
parser.add_argument('--outlev',nargs=1,type=int,
                             default=[1],
                             help='Output level')
parser.add_argument('-M','--NM_memory',nargs=1,type=int,
                             default=[3],
                             help='History size for nonmonotone linesearch')
args = parser.parse_args(sys.argv[1:])

#print(args.__dict__)
#for key in args.__dict__.keys():
#    print(key.ljust(25), args.__dict__[key])

###################################################################
# PRINT PARAMETERS
###################################################################
print()
print('------------------------------------------------------------')
print('Parameter           : value')
print('------------------------------------------------------------')
print('CONSTRAINED         : %s' % args.constrained)
print('ALGORITHM           : %s' % args.alg[0])
print('MAX_FUN             : %d' % args.max_fun[0])
print('OUTLEV              : %d' % args.outlev[0])
print('NM_MEMORY           : %d' % args.NM_memory[0])
print('------------------------------------------------------------')
print()

if args.version:
    print('\nrun_all.py version 0.1\n')
input()

prb = probs.prob_collection
keys= list(prb.keys())

if __name__ == '__main__':

    #set constrained
    constrained = args.constrained

    min_dim     = 2

    #set number of problems
    nprob = 0

    #choose which algorithm to run
    alg = [args.alg[0]]

    nsolver  = len(alg)

    #set nonmonotone memory size M (=1 for monotone)
    M = args.NM_memory[0]

    #set the maximum number of function evaluations
    max_fun = args.max_fun[0]

    #set the output level
    outlev = args.outlev[0]

    ii = 0
    for a in alg:
        #pick cyclically the problem to be solved
        for i,V in prb.items():
            if V.n >= min_dim:

                if a in ['DFN_DFL', 'DFL']:
                    print("alg: ",a," | prob: ",i) #,V.m.items())

                    for J,m in V.m.items():
                        if (constrained and m > 0) or (not constrained and m == 0):
                            print("alg: ",a," | prob: ",i," contype:",J,' n:',V.n,' m:',V.m[J])
                            CACHE = cache(max_fun,V.n+m+1)
                            x, f, stopfl, Dused, nf = dfn_dflint.box_DFN_DFL(a, M, J, V, CACHE, max_fun, outlev)
                            print("f: %15.6e | sotfl: %1d | len(CACHE.F): %6d | nf: %6d " % (f,stopfl,len(CACHE.F),nf))
                            print(' x=',x)
                            ii += 1
                            if outlev >= 0:
                                print('Hit return to continue ...',end='')
                                input()

        ii = 0
        if outlev > 2:
            print('Hit return to continue ...',end='')
            input()
