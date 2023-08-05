import numpy as np
from scipy.special import legendre
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd
from sympy import *
from sympy import simplify

# 1(a) Gamma Function
def gamma_(s):               #This works only for positive integers and positive half integers
    n= float(s)
    if n<0 or n%0.5!=0 :
       return ("Value must be either a positive integer or positive half integer.")
    else:
        if n==0:
           return "INFINITE"
        elif n==0.5:
           return np.sqrt(np.pi)
        elif n==1:
           return n
        else:
           return (n-1)*gamma_(n-1)    

# 1(b) Legendre 
def Lege(n):
    x=symbols('x')
    polyn=0
    if n%2 == 0:       #even
        m=int(n/2)
    else :             #odd
        m = int((n-1)/2)   
    s = np.arange(0,m+1,1)
    for i in s:
        polyn += ( (-1)**i*gamma_(2*n-2*i+1)*(x**(n-2*i)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
    polyn=simplify(polyn)
    fx=lambdify(x,polyn,"math")
    return polyn
from scipy.special import legendre
x = symbols('x')
print((legendre(3)[6]))
print(legendre(2)[2]*x**2+legendre(2)[0])
#print(Lege(2))
