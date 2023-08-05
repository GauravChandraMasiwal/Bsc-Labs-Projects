#name : gaurav chandra
#rollno : 2020PHY1122

#to study the distribution of particles for energies for bosons and fermions and 
#calculate internal energy and specific heat using this.

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import integrate
import warnings

warnings.filterwarnings('ignore')
c =3*10**(8)
h = 4.1357*10**(-15) #eV s
k = 8.617 * 10**(-5) # eV/K
N=6.022 *10**(23)
m =1 ; V=1
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
        
        if case1 == 'R' and case2 =='B':
            
            f=lambda e:(4*V*np.pi/(h*c)**3)*(1/(np.exp(alpha1)*np.exp(e/(k*i)) - 1))*(e**3)
            internal_energy.append(integrate.quad(f,U+0.0001,10)[0])
        
        elif case1 == 'NR' and case2 =='B':
            f=lambda e:(2*V*np.pi*(2*m)**(3/2)/h**3)*(1/(np.exp(alpha1)*np.exp(e/(k*i)) -1))*(e**1.5)
            internal_energy.append(integrate.quad(f,U+0.0001,2)[0])
        
        elif case1 == 'R' and case2 =='F':
            f=lambda e:(4*V*np.pi/(h*c)**3)*(1/(np.exp(alpha2)*np.exp(e/(k*i)) +1))*(e**3)
            internal_energy.append(integrate.quad(f,0.0001,10)[0])
            
        elif case1 == 'NR' and case2 =='F':
            f=lambda e:(2*V*np.pi*(2*m)**(3/2)/h**3)*(1/(np.exp(alpha2)*np.exp(e/(k*i)) +1))*(e**1.5)
            internal_energy.append(integrate.quad(f,0.0001,10)[0])        
        else:
            print('ERROR? plz only enter the valid cases i.e (N,NR,B,F)for(relativistic,non-relativistic,bosons and fermions)')
        
    return internal_energy




#fermi_level = (h**2 /(8*m)) * (3*N/(np.pi *V))**(2/3)
#fermi_temp = fermi_level/k

#plot of dN/de vs e for fermions 
e = np.linspace(0,20,100)
T =[1000,20000,50000]
U=0;ef=5
print('The characteristic temperature used here is 1/k = ',1/k)
fig=plt.figure()
fig.set_figheight(7)
fig.set_figwidth(9)

plt.subplot(2,2,1)
plt.plot(e,dN_de(e, T[0], case1='NR', case2='F'),'o-',c='r',markersize=5,label='low T= '+str(T[0]))
plt.legend()
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.title("FOR NON-RELATIVISTIC FERMIONS")
plt.xlim(0,10)
plt.subplot(2,2,2)
plt.plot(e,dN_de(e, T[0], case1='R', case2='F'),'o-',c='r',markersize=5,label='low T= '+str(T[0]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.title("FOR RELATIVISTIC FERMIONS")
plt.xlim(0,10)
plt.subplot(2,2,3)
plt.plot(e,dN_de(e, T[1], case1='NR', case2='F'),'o-',c='violet',markersize=5,label='high T= '+str(T[1]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
#plt.title("FOR NON-RELATIVISTIC BOSONS")

plt.subplot(2,2,4)
plt.plot(e,dN_de(e, T[1], case1='R', case2='F'),'o-',c='violet',markersize=5,label='high T= '+str(T[1]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')

#plt.title("FOR RELATIVISTIC BOSONS")
plt.suptitle('PLOT OF DISTRIBUTION OF PARTICLES VS ENERGY FOR FERMIONS')

plt.show()

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




#plot of dN/de vs e for bosons
U=0
e = np.linspace(U+0.0001,7,100)

T =[1500,10000,50000]

fig=plt.figure()
fig.set_figheight(7)
fig.set_figwidth(9)

plt.subplot(2,2,1)
plt.plot(e,dN_de(e, T[0], case1='NR', case2='B'),'o-',c='r',markersize=5,label='low T= '+str(T[0]))
plt.legend()
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.xlim(0,1)
plt.ylim(0,1e44)
plt.title("FOR NON-RELATIVISTIC BOSONS")

plt.subplot(2,2,2)
plt.plot(e,dN_de(e, T[0], case1='R', case2='B'),'o-',c='r',markersize=5,label='low T= '+str(T[0]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.xlim(0,2)

plt.title("FOR RELATIVISTIC BOSONS")

plt.subplot(2,2,3)
plt.plot(e,dN_de(e, T[1], case1='NR', case2='B'),'o-',c='violet',markersize=5,label='high T= '+str(T[1]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')
plt.xlim(0,1)
plt.ylim([0,1e45])
#plt.title("FOR NON-RELATIVISTIC BOSONS")

plt.subplot(2,2,4)
plt.plot(e,dN_de(e, T[1], case1='R', case2='B'),'o-',c='violet',markersize=5,label='high T= '+str(T[1]))
plt.legend(loc='best')
plt.grid()
plt.xlabel('e (in eV)')
plt.ylabel('dN / de ')

#plt.title("FOR RELATIVISTIC BOSONS")
plt.suptitle('PLOT OF DISTRIBUTION OF PARTICLES VS ENERGY FOR BOSONS')

plt.show()
#fermi_level = h**2 /(8*m) * (3*N/np.pi /V)**(2/3)

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




slope1=stats.linregress(T[1:],internal(T[1:], case1='NR', case2='F'))[0]
slope2=stats.linregress(T[1:],internal(T[1:], case1='R', case2='F'))[0]
slope3=stats.linregress(T,internal(T, case1='NR', case2='B'))[0]
slope4=stats.linregress(T,internal(T, case1='R', case2='B'))[0]
specific_heat=[slope1,slope2,slope3,slope4]

print('The specific heat for fermi gas for non-relativistic fermions is ',slope1)
print('The specific heat for fermi gas for relativistic fermions is ',slope2)
print('The specific heat for boson gas for non-relativistic bosons is ',slope3)
print('The specific heat for boson gas for relativistic bosons is ',slope4)


#references :Thermal Physics by SC Garg ,R M Bansal & C K Ghosh book pageno:595 section 14.2,14.3,14.4  
