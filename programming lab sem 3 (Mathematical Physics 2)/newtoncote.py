# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:05:49 2021

@author: chand
"""

import matplotlib.pyplot as plt
from scipy.integrate import quad  #inbuilt function for integtration

strt = int(input("enter the start limit : "))
end = int(input("enter the end limit : "))
n = int(input("enter the number n : "))  #h = (end - strt)/n

def func(x):
    return x**2

def trapezoid_method(a,b,n):
    sum1 = 0
    h = (b-a)/n
    for i in range(1,n):
        sum1 += func(a +(i*h))
    INT = (h/2) * ((func(a)+func(b)) + 2 * sum1)   
    return INT

def simpson_method(a,b,n):
    sum_odd = 0
    sum_even = 0
    h = (b-a)/(2*n)
    for i in range(1,2*n,2):
        sum_odd +=func(a + (i*h))
    for j in range (2,2*n,2):
        sum_even += func(a + (j*h))
    INT2 = (h/3) * ((func(a)+func(b)) + 4*sum_odd + 2*sum_even)
    return INT2

I1,I2,I3 = trapezoid_method(strt, end, n),simpson_method(strt, end, n),quad(func,strt,end)

Integration,Integration2 ,H= [],[],[]

n2 = eval(input("enter the list of number n in list format : "))  #input from user for x axis values on plot

for i in n2:
    H.append((end-strt)/i)     

for i in range (len(n2)):
    Integration.append(trapezoid_method(strt, end, n2[i]))
    Integration2.append(simpson_method(strt, end, n2[i]))

plt.title("LOG GRAPH 0F h VS I(h)")
plt.plot(H,Integration)
plt.plot(H,Integration2)
plt.legend(['trapezoidal','simpson'],loc = 'best')
plt.xscale('log')  #both x and y scale here are log scale
plt.yscale("log")
plt.xlabel("h in log")
plt.ylabel('integration for h  in log')
plt.grid()
plt.scatter(H,Integration)
plt.scatter(H,Integration2)
plt.savefig('integration.png')
plt.show()

print("-:" ,'\t', "INTEGRATION OF FUNCTION WITH LIMITS FROM ",strt,' to ',end,' :-')
print("numerical value of integration is : ",I3)
print("integration result using trapezoidal method is : ",I1)
print("integration result using simpson method is : ",I2)
print("error with trapezoidal method is : ",I1 - I3[0])
print("error with simpson method is : ",I2 - I3[0])

#TO INTEGRATE USING TRAPEZOID AND SIMPSON METHOD FOR GIVEN DATA

X1 = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]    #current
Y1 = [0,0.5,2.0,4.05,8,12.5,18,24.5,32,40.5,50]    #voltage

#trapezoid method

n = len(X1) -1      #for n terms ,there will be n-1 subintervals
h = (X1[-1] - X1[0])/(n)

SUM = 0
for i in range (1,n):
    SUM += Y1[i]

INT = (h/2)*((Y1[0] + Y1[-1]) + 2 *SUM)
print(SUM,h)
print("the result of integration using trapezoid method is : ",INT)


#simpson method

n = len(X1) - 1      

h = (X1[-1] - X1[0])/(n)

sum_even = 0
sum_odd = 0
for i in range (1,n,2):
    sum_odd  += Y1[i]
for i in range (2,n,2):
    sum_even += Y1[i]
    
INT2 = (h/3) * ((Y1[0] + Y1[-1]) + 4 * sum_odd + 2 * sum_even)

print("the result of integration using simpson method is : ",INT2)