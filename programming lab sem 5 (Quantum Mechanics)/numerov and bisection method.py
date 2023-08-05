import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
from scipy.special import hermite
from scipy.stats import linregress

u_odd = [0]
u_even=[1]

# POTENTIAL FUNCTION
def V2(x):
    V = (x**2)/2
    return V

# NUMEROV
def numerov(a,b,e,i_cs, n, V):
    x = np.linspace(a, b, n)
    h = x[1] - x[0]
    i_cs.append(i_cs[0] + ((h ** 2) / 2) + ((3 * (h ** 4)) / 24))
    u0, u1 = i_cs[0], i_cs[1]

    Vx = []
    for i in x:
        Vx.append(V(i))
    Vx = np.array(Vx)
    alpha = 2*(-Vx + e + 0.5)
    ci = 1 + (((h ** 2) / 12) * alpha)
    p = 0
    u = [u0, u1]
    for i in range(1,len(x)-1):
        p+=1
        u2 = (((12 - 10*ci[i]) * u[i]) - (ci[i - 1] * u[i-1])) / (ci[i + 1])
        u.append(u2)
        if p==len(x)-2:
            break
    #x.append(t[-1])
    return x,np.array(u),ci

# FINDING EIGENVALUE FOR REQ NODES

def energy(u,ics,emax,nodes,tol):
    emin=0
    p=0
    while abs(emax-emin)>=tol:
        p+=1
        u2=u[:-1]
        u3=u[1:]
        count=0
        #print("emin",emin)
        for i,j in zip(u2,u3):
            if i*j<0:
                count+=1

        #print(int(nodes/2))
        if count<=int(nodes/2):
            emin=(emax+emin)/2
        else:
            emax=(emax+emin)/2
        e_n=(emax+emin)/2
        u=numerov(0,4,e_n,ics,400,V2)[1]
    return (emax+emin)/2 + 0.5,p,emin,emax

#PARITY
def parity(n,u):
    if n%2==0:
        p_u=u[::-1]
        p_u1=p_u[1:]
    else:
        p_u=-1*u[::-1]
        p_u1 =p_u[1:]
    return p_u

def u_total(n,u,x_total):
    u_con= np.concatenate((parity(n,u),u))
    return u_con

i_cs_back=[0]
i_cs_back1=[0]

# EXAMPLE:

x,u,c=numerov(6,0,5,i_cs_back, 200, V2)  # ----- BACKWARD
x_par=np.linspace(0,-4,200)
y_par=parity(5,u)
X=np.concatenate((x,x_par))            # ----- ALL X
Y=np.concatenate((u,y_par))            # ----- ALL Y
plt.plot(X,Y)                          # ----- PLOT
plt.show()

# MULTIPLYING BY FACTOR

def factor(a,b,n,i_cs,N,V):
    i_cs_odd=[0]
    i_cs_even=[1]
    if n%2==0:
        x,u,c=numerov(a,b,n,i_cs_even,N,V)
    else:
        x, u, c = numerov(a, b, n, i_cs_odd, N, V)
    #x,u,c=numerov(a,b,n,i_cs,N,V)
    xcl1 = round(np.sqrt(2 * n + 1), 2)
    p = 0
    for i in x:
        if (round(i, 2)) == xcl1:
            break
        else:
            p += 1
    c=x[p-1]
    #print(c)
    h = (b - a) / N
    N2= int((c+h - a)/h)
    #N = N / 2
    if n%2==0:
        #h1 = (c - a) / N
        x1,u1,c1=numerov(a,c+h,n,i_cs_even,N2,V)
    else:
        #h1 = (c - a) / N
        x1, u1, c1 = numerov(a, c+h, n, i_cs_odd, N2, V)
    N3 = int((b - (c-h)) / h)
    x2,u2,c2=numerov(b,c-h,n,i_cs,N3,V)

    h=x2[1]-x2[0]
    k= u1[-1]/u2[-1]
    u2_new= k*u2
    u2=u2_new
    #print(u1)
    #print(u1[-1],u2[-1])
    return u2[-3], u2[-2], u1[-3], c,u1[-2], u2[-1], u1[-1], h

print(factor(0,6,2,i_cs_back,500,V2))

def derivative(factor):
    num= factor[2] + factor[0] - (12*factor[3] - 10)*factor[1]
    d_der= num/factor[4]
    return d_der

#print(derivative(factor(0,6,3,i_cs_back,400,V2)))

def bisection(a,b,nmax,i_cs,N,V,derivative,tol,factor):
    nmin=0
    d_der=derivative(factor(a,b,nmax,i_cs,N,V))
    while abs(d_der)>=tol:
        if d_der<=0:
            nmin=(nmin+nmax)/2
        else:
            nmax=(nmin+nmax)/2
        n_new= (nmin+nmax)/2
        print(n_new)
        d_der = derivative(factor(a, b, n_new, i_cs, N, V))
        #print(d_der)
    return (nmax+nmin/2 + 0.5), d_der

print(bisection(0,8,5,i_cs_back,1000,V2,derivative,10**(-10),factor)[0])

#2ND APPROACH
def alternate(factor):
    num1=factor[6]-factor[2] - factor[0] + factor[5]
    #den1=2*factor[7]*factor[-2]
    #d_diff=num1/den1
    return num1

print(bisection(0,8,5,i_cs_back,10000,V2,alternate,10**(-5),factor)[0])
