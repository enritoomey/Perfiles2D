# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 00:12:28 2015

@author: Enriquito
"""
import numpy as np

def atmosfera_estandar(calculo,*param):
    #datos
    b_suth = 1.458*10**(-6) #[Kg/(m*s*sqrt(K))]
    s_suth = 110.4 #[K]
    g  =9.81 #m/s**2
    R = 287.0 #J/(kg*K)
    k = 1.4
    Z = [0.0,11000.0,20100.0,32200.0,52400.0,61600.0,80000.0,95000.0]#[m]
    T = [288.0,216.5,215.6,228.5,270.5,252.5,180.5,180.5]#[K]
    P = [101325.0,22628.3619,5386.81,840.76,52.62,15.822,0.8455,0.04951]#Pa
    if calculo == 'altura':
        h = param[0]
        deltaT = param[1]   
        n = 1
        c = 0
        while c == 0 | n <7:
            if h < Z[n]:
                c =1
            else: n=n+1
        
        if T[n]==T[n-1]:
            t = T[n]+deltaT
            p = P[n-1]*np.e**(g*(h-Z[n-1])/R/T[n])
            rho = p/t/R
        else:
            m = (T[n]-T[n-1])/(Z[n]-Z[n-1])
            temp = m*(h-Z[n-1])+T[n-1]
            t = temp+deltaT
            p = P[n-1]*(temp/T[n-1])**(-g/m/R)
            rho = p/t/R
    
    if calculo == 'presion':
        p = param[0]
        deltaT = param[1]
        n = 1
        c = 0
        while c == 0 | n <7:
            if p > P[n]:
                c =1
            else: n=n+1
            
        if T[n]==T[n-1]:
            t = T[n]+deltaT
            h = R*T[n]/g*np.log(p/P[n-1])+Z[n-1]
            rho = p/t/R
        else:
            m = (T[n]-T[n-1])/(Z[n]-Z[n-1])
            h = ((p/P[n-1])**(-m*R/g)*T[n-1]-T[n-1])/m+Z[n-1]
            t = m*(h-Z[n-1])+T[n-1]+deltaT
            rho = p/t/R                        
        
        
        
    mu = b_suth*np.sqrt(t)/(1+s_suth/t)
    vson = np.sqrt(k*R*t)
    
    
    
    return h,deltaT,p,t,rho,mu,vson
    
    
    
 