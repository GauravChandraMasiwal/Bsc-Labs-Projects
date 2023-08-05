import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def trap(quad,V,T,a,b):
    n=1000
    h=float((b-a)/n)
    result = 0.5*(quad(a,V,T))+0.5*(quad(b,V,T))
   
    for i in range(1,n):
        result+=(quad(a+i*h,V,T))
    result *= h
    return result

def simp(quad,V,T,a,b):
    n=1000
    h = float((b-a)/n)
    result = (1/3)*(quad(a,V,T)+quad(b,V,T))
    for i in range(1,n,2):
        result+= 4*(quad(a+i*h,V,T))
    for j in range(2,n-1,2):
        result+=2*(quad(a+j*h,V,T))
    result*=h/3
    return result

def integrand(n_j,V,T):
    power=(-(h*n_j)**2)/(8*mo*k*T*V**(2/3))
    I=(n_j**2)*np.exp(power)
    return I

def Partition(Vlist,Tlist,a,b):
    Z=[]
    for V in Vlist:
        Z0=[]
        for T in Tlist:
            z=(np.pi/2)*simp(integrand,V,T, a, b)
            Z0.append(z)
        Z.append(Z0)
    return np.log(np.array(Z))

def Z_analytical(Vlist,Tlist):
    Z=[]
    for V in Vlist:
        Z0=[]
        for T in Tlist:
            z=V*(((2*np.pi*mo*k)/(h**2))*(3/2))*(T**(3/2))
            Z0.append(z)
        Z.append(Z0)
    return np.log(np.array(Z))

def forward_d(x,y):                  
    derive=[]
    for i in range(len(x)-1):
        dy=y[i+1]-y[i]
        dx=x[i+1]-x[i]
        derive.append(dy/dx)
    return np.array(derive)
    
def Pressure(Vlist,Tlist,a,b):
    Z=Partition(Vlist, Tlist, a, b)
    Plist=[]
    for i in range(len(Vlist)-1):
        d=forward_d(Vlist,Z.T[i])
        P=N*k*Tlist[i]*d
        Plist.append(P)
    return np.array(Plist)

def Internal_Energy(Vlist,Tlist,a,b):
    Z=Partition(Vlist, Tlist, a, b)
    Ulist=[]
    for i in range(len(Tlist)-1):
        d=forward_d(Tlist,Z[i])
        U=N*k*Tlist[i]*Tlist[i]*d.T[i]
        Ulist.append(U)
    return np.array(Ulist).T


def Entropy(Vlist,Tlist,a,b):
    U=Internal_Energy(Vlist, Tlist, a, b)
    Z=Partition(Vlist, Tlist, a, b)
    Slist=[]
    for i in range(len(Vlist)-1):
        
        S=(U[i]/Tlist[i])+N*k*(Z[i]-np.log(N)+1)
        Slist.append(S)
    return np.array(Slist)



h=6.626e-34      #Js
k=1.38e-23      #J/K
c=3e8           #m/s2
me=9.11e-31      #kg(mass of electron)
mp=1.67e-27         #kg(mass of proton)
mn=mp             #kg(mass of neutron)
mo=mp+mn
Na=6.022e23      # avagadro number
N=Na
ul=10e11
ll=0
n=20
q=4
deg=1
Vlist=np.linspace(20e-3,50e-3,n)
Tlist=np.linspace(150,450,n)
# print(Vlist,Tlist)


# a=simp(integrand,Vlist[0],Tlist[0],ll,ul)
Z=Partition(Vlist, Tlist, ll, ul)
Za=Z_analytical(Vlist, Tlist)
P=Pressure(Vlist, Tlist, ll, ul)
U=Internal_Energy(Vlist, Tlist, ll, ul)
S=Entropy(Vlist, Tlist, ll, ul)
slope,intercept=np.polyfit(Tlist[:n-1], U, deg)

print("The specific heat of capacity is :",slope/N,"J/K")

# print(P)
# print(S)
# print(len(Z),len(Tlist)-1)
# print(U)
# print(Z)
# print(z)

for i in range(len(Tlist)):
    if i%(n/q)==0:
        plt.plot(Tlist,Z[i],label=f' V = {round(Vlist[i],3)}$m^3$')
        plt.plot(Tlist,Za[i],color="black",linestyle="dotted")
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Temprature VS Partition function")
plt.xlabel("Temprature")
plt.ylabel("Partition function")
plt.legend()
plt.show()



for i in range(len(Vlist)):
    if i%(n/q)==0:
        plt.plot(Vlist,Z.T[i],label=f'T = {round(Tlist[i],3)}$\degree K$')
        plt.plot(Vlist,Za.T[i],color="black",linestyle="dotted")
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Volume VS Partition function")
plt.xlabel("Volume")
plt.ylabel("Partition function")
plt.legend()
plt.show()

for i in range(len(Vlist)):
    if i%(n/q)==0:
        plt.plot(Vlist[:n-1],P[i],label=f'T = {round(Tlist[i],3)}$\degree K$')
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Volume VS Pressure")
plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.legend()
plt.show()


for i in range(len(Tlist)):
    if i%(n/q)==0:
        plt.plot(Tlist[:n-1],P.T[i],label=f'V = {round(Vlist[i],3)}$m^3$')
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Temprature VS Pressure")
plt.xlabel("Temprature")
plt.ylabel("Pressure")
plt.legend()
plt.show()

for i in range(len(Vlist)):
    if i==0:
        plt.plot(Tlist[:n-1],U,label=f'V = {round(Vlist[i],3)}$m^3$')
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Temprature VS Internal Energy")
plt.xlabel("Temprature")
plt.ylabel("Internal Energy")
plt.legend()
plt.show()


for i in range(len(Vlist)):
    if i%(n/q)==0:
        plt.plot(Tlist,S[i],label=f'V = {round(Vlist[i],3)} $m^3$')
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Temprature VS Entropy")
plt.xlabel("Temprature")
plt.ylabel("Entropy")
plt.legend()
plt.show()

for i in range(len(Tlist)):
    if i%(n/q)==0:
        plt.plot(Vlist[:n-1],S.T[i],label=f'T = {round(Tlist[i],3)}$\degree K$')
plt.minorticks_on()
plt.grid(visible=True,axis="both",which="both")
plt.title("Volume VS Entropy")
plt.xlabel("Volume")
plt.ylabel("Entropy")
plt.legend()
plt.show()