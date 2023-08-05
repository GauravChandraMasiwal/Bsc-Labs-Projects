import numpy as np
import matplotlib.pyplot as plt
import math
from hermite_quad1 import MyHerQuad1
import pandas as pd
from MyIntegration2 import MySimp,MySimp_tol
from prettytable import PrettyTable
#from hermite_quad import limitsimp
eps_0=0.4
def f(x,n):
    eps=eps_0/2**n
    g=(1/2*math.sqrt(np.pi*eps))*np.exp(-x**2/2*eps)
    h=eps/(np.pi*(x**2+(eps)**2))
    return g,h
n=[1,2,3,4,5]
x=np.linspace(-7,7,50)
g1,g2,g3,g4,g5,h1,h2,h3,h4,h5=[],[],[],[],[],[],[],[],[],[]   
def plot():
    for i in x:
        g1.append(f(i,n[0])[0])
        g2.append(f(i,n[1])[0])
        g3.append(f(i,n[2])[0])
        g4.append(f(i,n[3])[0])
        g5.append(f(i,n[4])[0])
        h1.append(f(i,n[0])[1])
        h2.append(f(i,n[1])[1])
        h3.append(f(i,n[2])[1])
        h4.append(f(i,n[3])[1])
        h5.append(f(i,n[4])[1])
    fig,axs=plt.subplots(1,2,figsize=(12, 6))
    axs[0].plot(x,g1,c='red',label='n=1')
    axs[0].plot(x,g2,c='green',label='n=2')
    axs[0].plot(x,g3,c='orange',label='n=3')
    axs[0].plot(x,g4,c='blue',label='n=4')
    axs[0].plot(x,g5,c='pink',label='n=5')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('g(x)')
    axs[0].set_title("for $\sin$(x/$\epsilon$)/$\pi \epsilon$")
    axs[0].grid('true')
    axs[0].set_xlim([-20,20])
    axs[0].legend()
   
    axs[1].plot(x,h1,c='red',label='n=1')
    axs[1].plot(x,h2,c='green',label='n=2')
    axs[1].plot(x,h3,c='orange',label='n=3')
    axs[1].plot(x,h4,c='blue',label='n=4')
    axs[1].plot(x,h5,c='pink',label='n=5') 
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('h(x)')
    axs[1].set_title("for $\epsilon$/$\pi(x^2+\epsilon^2)$")
    axs[1].grid('true')
    axs[1].legend()
    plt.show()    
plot() 
b = 10
a = -10
n1=200
N=[1,2,3,4,5]
e=[]
q1=[]
q2=[]
def func1():
    for i in range(0,5):
         eps=eps_0/2**N[i]
         e.append(eps)
         z=lambda x:(1/2*math.sqrt(np.pi*e[i]))*np.exp(-x**2/2*e[i])
         z1=lambda x: e[i]/(np.pi*(x**2+(e[i])**2))
         r=MyHerQuad1(n1,z,z1)
         q1.append(r[0])
         q2.append(r[1])
    data = {"Value of N":N,"Integral I1":q1,"Integral I2":q2}          
    print(pd.DataFrame(data))
func1()   


    
def limitsimp(f,a,b,N = 2, N_max= 10000,tol = 0.5e-3):
    #b_b = [b]
    a_a = [a]
    trial = 0
    val = [MySimp_tol(f,a,b,N,N_max,tol)[0]]
    b0 = [b]
    while True:
        val.append(MySimp_tol(f,a,10*b,N,N_max,tol)[0])
        T = val[-1]-val[-2]
        b0.append(10*b)
        a_a.append(10*a)
        if abs(T) <= abs(tol*val[-1]) or round(abs(T*val[-1]),10) == 0:
            break
        
        
        elif abs(T) > abs(tol*val[-1]):
            b = 10*b
            a = 10*a
                
        else:
            trial = 1
            
            break
        

    if trial == 0:
        return val[-1],b0,a_a,val
    else:
        msg = "Error! Tolerance is not defined"
        print(msg)
b = 10
a = -10
print("FOR N=1")
def delta(x):
    s=lambda n: (1/np.pi)*(eps_0/2**n)/(x**2+(eps_0/2**n)**2)
    return s(1)
data = {"Value of a": limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[2],"Value of b":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[1],"Value of Integral 1":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[3]}
print(pd.DataFrame(data))
print("FOR N=2")    
def delta(x):
    s=lambda n: (1/np.pi)*(eps_0/2**n)/(x**2+(eps_0/2**n)**2)
    return s(2)
data = {"Value of a": limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[2],"Value of b":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[1],"Value of Integral 1":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[3]}
print(pd.DataFrame(data))    
print("FOR N=3")    
def delta(x):
    s=lambda n: (1/np.pi)*(eps_0/2**n)/(x**2+(eps_0/2**n)**2)
    return s(3)
data = {"Value of a": limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[2],"Value of b":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[1],"Value of Integral 1":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[3]}
print(pd.DataFrame(data))      
    
print("FOR N=4")    
def delta(x):
    s=lambda n: (1/np.pi)*(eps_0/2**n)/(x**2+(eps_0/2**n)**2)
    return s(4)
data = {"Value of a": limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[2],"Value of b":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[1],"Value of Integral 1":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[3]}
print(pd.DataFrame(data)) 
     
print("FOR N=5")    
def delta(x):
    s=lambda n: (1/np.pi)*(eps_0/2**n)/(x**2+(eps_0/2**n)**2)
    return s(5)
data = {"Value of a": limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[2],"Value of b":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[1],"Value of Integral 1":limitsimp(delta,a,b, N_max= 10**8,tol = 0.5e-3)[3]}
print(pd.DataFrame(data))    
    