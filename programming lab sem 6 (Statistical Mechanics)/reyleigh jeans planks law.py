import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

#Reyleigh jeans and planks law and density of states plot


#plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')

#density of states

p = np.linspace(10,30,100)
#nu = 10**p
nu = np.logspace(10,30,21)

c = 2.9 * 10**(8)   # m/s
lo = 10**(-10) 
k = 8.61733 * 10**(-5)  # ev K^(-1)
h = 6.626 * 10**(-34)
e = 1.6 * 10**(-19)

g_nu = 8*np.pi*(nu**2)*(lo**3)/(c**3)

x = 2*lo*nu/c

g_x = np.pi*(x**2)

plt.figure(figsize = (10,7))

plt.subplot(1,2,1)
plt.plot(nu, g_nu)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('ν')
plt.ylabel('G (ν)')
plt.grid()

plt.subplot(1,2,2)
plt.plot(x, g_x)
plt.xlabel('x')
plt.ylabel('G (x)')
plt.xscale('log')
plt.yscale('log')
plt.grid()

plt.suptitle("Density of states: Rayleigh-Jeans and Planck's Law")

plt.savefig('04-1',dpi = 800)
plt.show()



#Rayleigh-Jeans: energy of states

x = np.linspace(0.000001,12,1000)

f_rj = x**2

#checking slope
print('slope (Rayleigh-Jeans):', linregress(np.log(x),np.log(f_rj))[0])



T_arr = [1200, 1500, 1800]

plt.figure(figsize = (10,7))

plt.subplot(1,2,1)
plt.plot(x,f_rj)
plt.xlabel('x')
plt.ylabel(r'$f_{RJ} (x)$')
plt.grid()


plt.subplot(1,2,2)
for T in T_arr:
    e_ = k*T
    nu_ = e_/h
    l_ = h*c/e_

    nu = nu_*x
    
    u_rj = 8*np.pi*e_ * f_rj / (l_**3 * nu_ * e) 
    
    
    plt.plot(nu,u_rj, label=f'T = {T}K')
    #plt.vlines(x= 4**14, ymin = min(u_rj), ymax = max(u_rj) )
    #plt.vlines(x= 8**14, ymin = min(u_rj), ymax = max(u_rj) )
      

plt.xlabel('ν')
plt.ylabel('U (ν)')
plt.legend()
plt.grid()

plt.suptitle('Rayleigh Jeans law: Energy of states')

plt.savefig('04-2',dpi=800)
plt.show()

#Plancks Radiation law: Energy of states

f_p = x**3 / (np.exp(x) - 1)

def fp(x):
    return x**3 / (np.exp(x) - 1)

plt.figure(figsize = (12,8))

plt.subplot(1,2,1)
plt.plot(x, f_p)
plt.xlabel('x')
plt.ylabel(r'$f_P (x)$')
plt.grid()

ii = list(f_p).index(max(f_p))

print(x[ii])



u_max = []
x_max = []

nu_vis = np.linspace(4*(10**14), 8*(10**14),100)
for T in T_arr:
    e_ = k*T
    nu_ = e_/h
    l_ = h*c/e_

    nu = nu_*x
    
    x_ = h*nu_vis/(k*T)
    nuu = nu_*x_
    
    u_p = 8*np.pi*e_ * fp(x) / (l_**3 * nu_ * e) 
    upp = 8*np.pi*e_ * fp(x_) / (l_**3 * nuu * e) 
    
    u_max.append(max(u_p))
    
    i = list(u_p).index(max(u_p))
    x_max.append(nu[i])
    
    
    
    
    plt.subplot(1,2,2)
    #plt.plot(x, u_p, label=f'T = {T}K')
    plt.scatter(nuu, upp)
    
plt.xlabel('ν')
plt.ylabel('U (v)')
plt.legend()
plt.grid()

plt.suptitle('Plancks radiation law: Energy of states')

plt.savefig('04-3',dpi=800)
plt.show()

print(u_max)
print(x_max)


