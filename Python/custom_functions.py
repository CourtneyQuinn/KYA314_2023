#!/usr/bin/env python
# coding: utf-8
# module with custom functions for dynamical systems analysis

import numpy as np
import scipy as scipy
import scipy.stats as stats
import scipy.linalg as linalg
from pprint import pprint
import sys
import os
import copy
import string
import glob
import xarray as xr

import warnings

#################### Functions ######################
def MyJacobian(f,x,h):
    '''Jacobian of a function with arbitrary dimensions
    2nd Order Central Difference
    Parameters
    ----------
    f : function handle
        Takes input x
    x : array (n, N)
        Points at which Jacobian is to be calculated
        If N = 1, x is one point
    h : float
        Step size for finite differences
    Returns
    -------
    df : array, shape (m, n, N)
        Array containing Jacobian matrix in first 2 dimensions
        3rd dimension denotes which point in original x array
    '''

    n = x.shape[0]
    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)

    N = x.shape[1]
    ftemp = f(x)
    m = ftemp.shape[0]
    
    df = np.empty((m,n,N))
    df[:] = np.nan

    for i in np.arange(n):
        xi1 = x.copy()
        xi2 = x.copy()
        xi1[i, :] = x[i, :] + h
        xi2[i, :] = x[i, :] - h
        dfi = (f(xi1) - f(xi2)) / (2 * h)
        df[:, i, :] = dfi
        
    return df

def MySolve(f,x0,df,tol,maxit):
    '''Newton iteration to find zeros of a nonlinear system
    of equations, 0 = f(x)
    2nd Order Central Difference
    Input
    ----------
    f  : function handle
         Takes input x
    x0 : array (n,)
         Initial guess
    df : function handle 
         Jacobian
    tol: float
         tolerance for convergence
    maxit : integer
            number of mximum iterations for
            Newton solver
    Returns
    -------
    x : array (n,)
        solution of f(x)=0
    converged : logical
                1 - converged
                0 - not converged
    J : array (n,n)
        Jacobian at solution
    '''
    
    if x0.ndim == 1:
        x0 = np.expand_dims(x0, axis=1)
    
    J = df(x0)
    x1 = x0 - np.matmul(np.linalg.inv(J),f(x0))
    
    ek1 = linalg.norm(x1 - x0)
    rk = linalg.norm(f(x0))
    print('Error is %g, Residual is %g' % (ek1,rk))
    x0 = x1
    
    for i in np.arange(0,maxit):
        J = df(x0)
        if i == maxit-1:
            converged = False;
            x = x0
            break
        elif ek1 < tol and rk < tol:
            converged = True
            J = df(x0)
            x = x0
            break
        else:
            x1 = x0 - np.matmul(np.linalg.inv(J),f(x0))
            ek1 = np.linalg.norm(x1 - x0)
            rk = np.linalg.norm(f(x0))
            x0 = x1
            print('Error is %g, Residual is %g' % (ek1,rk))

    if ek1 < tol and rk < tol:
        converged = True
        J = df(x0)
        x = x0
    
    return x, converged, J