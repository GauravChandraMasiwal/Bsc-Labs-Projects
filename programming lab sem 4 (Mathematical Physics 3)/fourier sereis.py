# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 10:40:15 2022

@author: hp
"""

#Name=VINAYAK JOSHI
#College Roll No. = 2020PHY1251
#University Roll no. = 20068567065

from MyIntegration import MySimp
from MyIntegration import MyTrap
from MyIntegration import MyLegQuadrature,MyLegQuadrature_tol
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

from scipy import integrate

#from legendre import Lege
from scipy.special import legendre
from sympy import *
from sympy import simplify,lambdify
import scipy.integrate as integrate

#  function1 defined for a range [-l,l]
def func1(x):
    if (x<0 ):
        return 0
    elif (x>=0):
        return 1
    
def fourier(li, lf, n, f):
    l = (lf-li)/2
    # Constant term
    a0=1/l*integrate.quad(lambda x: f(x), li, lf)[0]
    # Cosine coefficents
    A = 0
    # Sine coefficents
    B = 0
     
    for i in range(1,n+1):
        A=1/l*MyLegQuadrature(lambda x: f(x)*np.cos(i*np.pi*x/l), li, lf,40,1)[0]
        B=1/l* integrate.quad(lambda x: f(x)*np.cos(i*np.pi*x/l), li, lf)[0]
 
    return [a0/2.0, A, B]
 
if __name__ == "__main__":
 
 
    # Limits for the functions
    li = -2.5
    lf = 2.5
 
    # Number of harmonic terms
    nt = [3]
 
    # Fourier coeffficients for various functions
    for i in range (len(nt)):
        coeffs = fourier(li,lf,nt[i],func1)
        print('Fourier coefficients for function 1 are\n')
        print("a0 for n=",nt[i], "is="+str(coeffs[0]))
        print("an for n=",nt[i],"is=" +str(coeffs[1]))
        print("bn for n=",nt[i],"is="+str(coeffs[2]))
        print('-----------------------\n\n')