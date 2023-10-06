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

def MyIVP(f,x0,tspan,h):
    '''Solving initial-value problems for ODEs of form 
    x*(t)=f(t,x(t))
    4th-order Runge-Kutta
    Input
    ----------
    f : function handle
        Takes input x
    x0 : array (n, m)
        initial conditions
        If m = 1, x is one point
    tspan : array (2,)
        start and end times of integration
    h : float
        Step size for numerical integration
    Returns
    -------
    xt : array (n, m, N+1)
        Array containing solution x(t) at times in t
        2nd dimension denotes trajectory corresponding to 
        x0 values
    t : array (N+1,)
        time steps for integration
    xend : array (n, m)
        final point in trajectory for each x0 value
    '''

    n = x0.shape[0]
    if x0.ndim == 1:
        x0 = np.expand_dims(x0, axis=1)

    N = int((tspan[1]-tspan[0])/h)
    m = x0.shape[1]
    
    xt = np.empty((n,m,N+1))
    xt[:] = np.nan
    xt[:,:,0]= x0

    t = np.empty((N+1))
    t[:] = np.nan
    t[0] = tspan[0]
    
    for k in np.arange(N):
        tk = t[k]
        xk = xt[:,:,k]
        j1 = f(tk,xk)
        j2 = f(tk+h/2,xk+(h/2)*j1)
        j3 = f(tk+h/2,xk+(h/2)*j2)
        j4 = f(tk+h,xk+h*j3)
        
        xt[:,:,k+1] = xk + (h/6)*(j1+2*j2+2*j3+j4)
        t[k+1] = tspan[0]+k*h
    
    xend = xt[:,:,-1]
    
    return xt, t, xend

def MyTrackCurve(userf,userdf,y0,ytan,**kwargs):
    '''Track a curve of n equations for n+1 variables
    0=f(y), f: R(n+1) to R(n)
    Input
    ----------
    userf : function handle
            Takes input y
    userdf : function handle
            Takes input y
    y0 : array (n+1,)
         initial guess for first point on curve
    ytan : array (n+1,)
           approx. tangent to curve at y0
    Returns
    -------
    ylist : array (n+1, nmax)
        array containing points along curve
        y_k = ylist(:,k)
    '''
    options={'stepsize':0.01,'nmax':100,'tol':1e-5,'maxit':10};
    options.update(kwargs)  
    s = options.get('stepsize')
    nmax = options.get('nmax')
    tol = options.get('tol')
    maxit = options.get('maxit')
    
    m = y0.shape[0]
    
    if y0.ndim == 1:
        y0 = np.expand_dims(y0, axis=1)
    if ytan.ndim == 1:
        ytan = np.expand_dims(ytan, axis=1)
    
    n = userf(y0).shape[0]
    ylist = np.zeros((m,nmax+1))*np.nan
    ylist[:,0] = y0.squeeze(axis=1)
    for j in np.arange(0,nmax):
        yj = y0 + s*ytan
        fj = lambda y: np.matmul(np.transpose(ytan),(y - yj))
        F = lambda y: np.concatenate((userf(y),fj(y)),axis=0)
        def df(y):
            J = MyJacobian(fj,y,1e-5)
            if J.ndim > 2:
                J = J.squeeze(axis=2)
            dfout = np.concatenate((userdf(y),J),axis=0)
            return dfout
        
        yk,converged,J = MySolve(F,yj,df,tol,maxit)
        while converged == 0:
            s = np.max([s/2,1e-8])
            yj = y0 + s*ytan
            fj = lambda y: np.matmul(np.transpose(ytan),(y - yj))
            F = lambda y: np.concatenate((userf(y),fj(y)),axis=0)
            yk,converged,J = MySolve(F,yj,df,tol,maxit)
        
        s = np.min([s*2,options.get('stepsize')])
        dfk = userdf(yk)
        zeros0 = np.zeros((n,))
        zeros1 = np.append(zeros0,np.array([1.]),axis=0)
        fk = np.concatenate((dfk,np.transpose(ytan)),axis=0)
        z = np.linalg.lstsq(fk,zeros1,rcond=None)[0]
        mult = np.sign(np.matmul(np.transpose(z),ytan))
        ytan = np.divide(z,linalg.norm(z))*mult
        ytan = np.expand_dims(ytan, axis=1)
        ylist[:,j+1] = np.squeeze(yk,axis=1)
        y0 = yk
    
    return ylist

def LyapQR(M,xini,N,**kwargs):
    '''Compute Lyapunov exponents for either a map 
    or an array of tangent linear propagators
    Input
    ----------
    M : function handle or array (n,n,N+1)
        Takes input x
    xini : array (n,)
           initial point
    N : integer
        number of steps over which to compute exponents
    Returns
    -------
    lams : array (n,)
           final Lyapunov exponents
    Rdiag : array (n,N)
            diagonal pf R at each step
    Lambda : array (n,N)
             Lyapunov exponents at increasing windows
    x : array (n,N+1)
        state space trajectory
    '''
    options={'dM':lambda x:MyJacobian(M,x,1e-6),'h':1};
    options.update(kwargs)  
    dM = options.get('dM')
    h = options.get('h')

    n = xini.shape[0]
    Q = np.eye(n,dtype = float)
    Rdiag = np.empty((n,N))*np.nan
    Lambda = np.empty((n,N))*np.nan
    x = np.empty((n,N+1))*np.nan
    x[:,0] = xini
    
    if xini.ndim == 1:
        xini = np.expand_dims(xini, axis=1)
    
    if type(M)=='function':
        for j in np.arange(0,N):
            x[:,j+1] = M(x[:,j])
            Aj = dM(x[:,j])
            B = np.matmul(Aj,Q)
            [Q,R] = linalg.qr(B)
            Rj = np.diag(R)
            l = np.where(Rj<0)[0]
            if len(l) != 0:
                Q[:,l] = -1*Q[:,l]
            Ri = abs(Rj)
            Rdiag[:,j] = Rj
            Lambda[:,j] = 1/j*np.sum(np.log(Rdiag[:,:j+1]),axis=1)
    else:
        for j in np.arange(0,N):
            B = np.matmul(M[:,:,j],Q)
            [Q,R] = linalg.qr(B)
            Rj = np.diag(R)
            l = np.where(Rj<0)[0]
            if len(l) != 0:
                Q[:,l] = -1*Q[:,l]
            Rj = abs(Rj)
            Rdiag[:,j] = Rj
            Lambda[:,j] = 1/(j*h)*np.sum(np.log(Rdiag[:,:j+1]),axis=1)
            
    lams = Lambda[:,-1]
    
    return lams, Rdiag, Lambda, x
