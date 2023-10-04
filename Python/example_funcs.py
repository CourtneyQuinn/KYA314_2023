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
    xout : array, shape (1, N)
           Array containing 1-step interate of map
           2nd dimension denotes iterates for different starting points
    '''

    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    xout = lam*x[0,:]*(1-x[0,:]);
    
    return xout

def HenonMap(x,p):
    '''2D map exhibiting complex behaviour
    Input
    ----------
    x : array (2, N)
        Points to iterate
        If N = 1, x is one point
    p : array (2,)
        parameter values
    Returns
    -------
    xout : array, shape (2, N)
           Array containing 1-step interate of map
           2nd dimension denotes iterates for different starting points
    '''

    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    n = x.shape
    xout = np.empty(n)
    xout[:] = np.nan 
    
    alpha = p[0]
    beta = p[1]
    
    xout[0,:] = 1 - alpha*x[0,:]**2+x[1,:];
    xout[1,:] = beta*x[0,:]
    
    return xout

def PredatorPrey(t,x,p):
    '''2D ODE system exhibiting limit cycles
    Input
    ----------
    x : array (2, N)
        state space values
    p : array (4,)
        parameter values
    Returns
    -------
    xout : array, shape (2, N)
           Array containing derivative at x
           2nd dimension denotes derivatives for different points
    '''

    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    n = x.shape
    xout = np.empty(n)
    xout[:] = np.nan 
    
    r = p[0]
    m = p[1]
    a = p[2]
    b = p[3]
    
    xout[0,:] = r*x[0,:]-m*x[0,:]*x[1,:];
    xout[1,:] = -a*x[1,:]+b*x[0,:]*x[1,:];
    
    xout = xout.squeeze()
    
    return xout

def Stommel(t,x,p):
    '''2D ODE system for the AMOC
    Input
    ----------
    t : float
        time value
    x : array (2, N)
        state space values
    p : array (3,)
        parameter values
    Returns
    -------
    xout : array, shape (2, N)
           Array containing derivative at x
           2nd dimension denotes derivatives for different points
    '''

    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    n = x.shape
    xout = np.empty(n)
    xout[:] = np.nan 
    
    F = p[0,]
    alpha = p[1,]
    mu = p[2,]
    
    xout[0,:] = -alpha*(x[0,:]-1)-x[0,:]*(1+mu**2*(x[0,:]-x[1,:])**2);
    xout[1,:] = F-x[1,:]*(1+mu**2*(x[0,:]-x[1,:])**2);
    
    return xout

def Van_Duff(t,x,p):
    '''2D ODE system exhibiting relaxation oscillations
    Input
    ----------
    t : float
        time value
    x : array (3, N)
        state space values
    p : array (5,)
        parameter values
    Returns
    -------
    xout : array, shape (3, N)
           Array containing derivative at x
           2nd dimension denotes derivatives for different points
    '''

    if x.ndim == 1:
        x = np.expand_dims(x, axis=1)
    
    n = x.shape
    xout = np.empty(n)
    xout[:] = np.nan 
    
    mu = p[0,]
    alpha = p[1,]
    beta = p[2,]
    K = p[3,]
    omega = p[4,]
    
    xout[0,:] = x[1,:];
    xout[1,:] = mu*(1-x[0,:]**2)*x[1,:]-alpha*x[0,:]-beta*x[0,:]**3+K*np.cos(x[2,:]);
    xout[2,:] = omega*x[2,:]**0
    
    return xout
