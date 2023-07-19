#!/usr/bin/env python
# coding: utf-8
# module with example functions for KYA314

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
def LogisticMap(x,lam):
    '''1D map exhibiting complex behaviour
    Input
    ----------
    x : array (1, N)
        Points to iterate
        If N = 1, x is one point
    lam : float
          parameter value
    Returns
    -------
    xout : array, shape (n, N)
           Array containing 1-step interate of map
           2nd dimension denotes iterates for different starting points
    '''
    n = x.shape[0]
    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    xout = lam*x[0,:]*(1-x[0,:]);
    
    return xout