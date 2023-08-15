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
