# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 20:20:54 2023

@author: chand
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import integrate

c =3*10**(8)
h = 4.1357*10**(-15) #eV s
k = 8.617 * 10**(-5) # eV/K

    
#for bose-einstein
def f_be(X,alpha):
    
    Y = 1/(np.exp(alpha)*np.exp(X) - 1)
    return Y

#for fermi-dirac 
def f_fd(X,alpha):
    
    Y = 1/(np.exp(alpha)*np.exp(X) + 1)
    return Y

#distribution of particles 
def dN_de(e,T,case1,case2):  #case can be either R(relativistic) or NR(non relativistic)
    alpha=-U/(k*T)
    alpha2 =-ef/(k*T)
    X = e/(k*T)
    if case1 == 'R' and case2 =='B':
        return (4*V*np.pi/(h*c)**3)*f_be(X,alpha)*(e**2)
    
    elif case1 == 'NR' and case2 =='B':
        return (2*V*np.pi*(2*m)**(3/2)/h**3)*f_be(X,alpha)*(e**0.5)
        
    elif case1 == 'R' and case2 =='F':
        return (4*V*np.pi/(h*c)**3)*f_fd(X,alpha2)*(e**2)
    
    elif case1 == 'NR' and case2 =='F':
        return (2*V*np.pi*(2*m)**(3/2)/h**3)*f_fd(X,alpha2)*(e**0.5)
    
    else:
        print('ERROR? plz only enter the valid cases i.e (N,NR,B,F)for(relativistic,non-relativistic,bosons and fermions)')
    
    
    
def internal(T,case1,case2):
    internal_energy = []
    for i in T:
        alpha1=-U/(k*i)
        alpha2 =-ef/(k*i)
        #X = e/(k*i)
        
        if case1 == 'R' and case2 =='B':
            
            f=lambda e:(4*V*np.pi/(h*c)**3)*(1/(np.exp(alpha1)*np.exp(e/(k*i)) - 1))*(e**3)
            internal_energy.append(integrate.quad(f,U+0.0001,5)[0])
        
        elif case1 == 'NR' and case2 =='B':
            f=lambda e:(2*V*np.pi*(2*m)**(3/2)/h**3)*(1/(np.exp(alpha1)*np.exp(e/(k*i)) -1))*(e**1.5)
            internal_energy.append(integrate.quad(f,U+0.0001,5)[0])
        
        elif case1 == 'R' and case2 =='F':
            f=lambda e:(4*V*np.pi/(h*c)**3)*(1/(np.exp(alpha2)*np.exp(e/(k*i)) +1))*(e**3)
            internal_energy.append(integrate.quad(f,0.0001,10)[0])
            
        elif case1 == 'NR' and case2 =='F':
            f=lambda e:(2*V*np.pi*(2*m)**(3/2)/h**3)*(1/(np.exp(alpha2)*np.exp(e/(k*i)) +1))*(e**1.5)
            internal_energy.append(integrate.quad(f,0.0001,10)[0])        
        else:
            print('ERROR? plz only enter the valid cases i.e (N,NR,B,F)for(relativistic,non-relativistic,bosons and fermions)')
        
    return internal_energy




    
N=6.022 *10**(23)
U=1;ef=5
m =1
e = np.linspace(0,4,100)
V=1
#fermi_temp = h**2 /(8*m*k) * (3*N/np.pi /V)**(2/3)

#print(ef/k)
'''
#plot of dN/de vs e for fermions 
e = np.linspace(0,10,70)
T =[10,100,1000,10000,50000]

fig=plt.figure()
fig.set_figheight(5)
fig.set_figwidth(10)
plt.subplot(1,2,1)
plt.plot(e,dN_de(e, T[2], case1='NR', case2='F'),'o-',c='r',markersize=5,label='low T= '+str(T[2]))
plt.plot(e,dN_de(e, T[3], case1='NR', case2='F'),'o-',c='violet',markersize=5,label='high T= '+str(T[3]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.title("FOR NON-RELATIVISTIC FERMIONS")
#plt.ticklabel_format(useOffset=False,style='plain')

plt.subplot(1,2,2)
plt.plot(e,dN_de(e, T[2], case1='R', case2='F'),'o-',c='r',markersize=5,label='low T= '+str(T[2]))
plt.plot(e,dN_de(e, T[3], case1='R', case2='F'),'o-',c='violet',markersize=5,label='high T= '+str(T[3]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
#plt.ticklabel_format(useOffset=False,style='plain')
plt.title("FOR RELATIVISTIC FERMIONS")
plt.suptitle('PLOT OF DISTRIBUTION OF PARTICLES VS ENERGY FOR FERMIONS')

plt.show()
#fermi_level = h**2 /(8*m) * (3*N/np.pi /V)**(2/3)

#for internal energy plot

T=np.linspace(10,200000,200)

fig=plt.figure()
fig.set_figheight(7)
fig.set_figwidth(10)
plt.subplot(2,2,1)
plt.plot(T[:10]/1000,internal(T[:10], case1='NR', case2='F'),'o-',c='r',markersize=4,label='non-relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
plt.title("FOR LOW TEMPERATURE")
plt.legend()

plt.subplot(2,2,2)
plt.plot(T[10:]/1000,internal(T[10:], case1='NR', case2='F'),'o-',c='g',markersize=4,label='non-relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
plt.title("FOR HIGH TEMPERATURE")
plt.legend()
plt.suptitle('PLOT OF INTERNAL ENERGY VS TEMPERATURE FOR FERMIONS ')
#plt.show()
#print(internal(e, T=[100], case1='NR', case2='F'))

plt.subplot(2,2,3)
plt.plot(T[:10]/1000,internal(T[:10], case1='R', case2='F'),'o-',c='r',markersize=4,label='relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
#plt.title("FOR LOW TEMPERATURE")
plt.legend()

plt.subplot(2,2,4)
plt.plot(T[10:]/1000,internal(T[10:], case1='R', case2='F'),'o-',c='g',markersize=4,label='relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
#plt.title("FOR HIGH TEMPERATURE")
plt.suptitle('PLOT OF INTERNAL ENERGY VS TEMPERATURE FOR FERMIONS ')
plt.legend()
plt.show()

'''


#plot of dN/de vs e for bosons
U=0
e = np.linspace(U+0.1,1,100)

T =[10,100,1500,10000,50000]
T=1000
plt.plot(e,f_be(e/(k*T),alpha=-U/(k*T))*(e**0.5),'o-',c='r',markersize=5,label='low T= '+str(T))
#plt.plot(e,dN_de(e, T[3], case1='NR', case2='B'),'o-',c='violet',markersize=5,label='high T= '+str(T[3]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.title("FOR NON-RELATIVISTIC BOSONS")
#plt.ticklabel_format(useOffset=False,style='plain')
plt.show()
#fermi_level = h**2 /(8*m) * (3*N/np.pi /V)**(2/3)
'''
#for internal energy plot

T=np.linspace(10,200000,200)

fig=plt.figure()
fig.set_figheight(7)
fig.set_figwidth(10)
plt.subplot(2,2,1)
plt.plot(T[:10]/1000,internal(T[:10], case1='NR', case2='B'),'o-',c='r',markersize=4,label='non-relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
plt.title("FOR LOW TEMPERATURE")
plt.legend()

plt.subplot(2,2,2)
plt.plot(T[10:]/1000,internal(T[10:], case1='NR', case2='B'),'o-',c='g',markersize=4,label='non-relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
plt.title("FOR HIGH TEMPERATURE")
plt.legend()
plt.suptitle('PLOT OF INTERNAL ENERGY VS TEMPERATURE FOR BOSONS ')
#plt.show()
#print(internal(e, T=[100], case1='NR', case2='F'))

plt.subplot(2,2,3)
plt.plot(T[:10]/1000,internal(T[:10], case1='R', case2='B'),'o-',c='r',markersize=4,label='relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
#plt.title("FOR LOW TEMPERATURE")
plt.legend()

plt.subplot(2,2,4)
plt.plot(T[10:]/1000,internal(T[10:], case1='R', case2='B'),'o-',c='g',markersize=4,label='relativistic')
plt.grid()
plt.xlabel('Temperature (in K) (x10**3)')
plt.ylabel('internal energy(in eV) ')
#plt.title("FOR HIGH TEMPERATURE")
plt.suptitle('PLOT OF INTERNAL ENERGY VS TEMPERATURE FOR BOSONS ')
plt.legend()
plt.show()
'''