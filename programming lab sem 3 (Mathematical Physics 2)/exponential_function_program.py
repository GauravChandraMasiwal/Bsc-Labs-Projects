# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 13:40:24 2022

@author: chand
"""

import numpy as np
import matplotlib.pyplot as plt

def factorial(n):
    prod = 1
    for i in range(2,n+1):
        prod = i*prod
    return prod
def exp(x,n):
    sum1 = 0
    for j in range(n+1):
        sum1 += x**j/factorial(j) 
    return sum1

def for_plot(n):
    X = np.linspace(-5,5,50)
    
    for i in range(X):
        exp_f = []
        exp_f.append(exp(i,n))
    
def main_func():
    
    n = 2
    tol = 0.5 * 10**(-5)
    list_1 = []
    nlist = []
    
    array = np.array([for_plot(n),])
    
    x = 
    sum1 = exp(X,n)
    list_1.append(sum1)
    nlist.append(n)
    while True:
        n = n+1
    
        sum2 = exp(X,n)
    
        if abs(sum2 - sum1) > tol:
            n = n+1
            sum1 = sum2
            sum2 = exp(X,n)
            list_1.append(sum1)
            nlist.append(n)
        else :
            break
      
    plt.plot
    
        