'''
Roll No.: 2020PHY1122
Name: Gaurav Chandra
'''

#This program is related to the particle distribution function 
#here i have made functions of maxwell distribution function, bose einstein and fermi dirac distribution function
#Also i have made plots for different temperatures and compared the plots of each of these functions

import numpy as np
import matplotlib.pyplot as plt

k = 8.6173 * 10**(-5)
X = np.linspace(-4,4,100)
alpha = 0
X1 = np.linspace(alpha+10**(-5),4,100)
def f_mb(X):
    
    Y = np.exp(-X)
    return Y

def f_be(X,alpha):
    
    Y = 1/(np.exp(alpha)*np.exp(X) - 1)
    return Y


def f_fd(X,alpha):
    
    Y = 1/(np.exp(alpha)*np.exp(X) + 1)
    return Y


#plot 1
plt.plot(X,f_mb(X),label = "Maxwell's Boltzmann")
plt.plot(X1,f_be(X1,alpha),label = "Bose Einstein")
plt.plot(X,f_fd(X,alpha),label = "Fermi Dirac")
plt.ylim([0,10])
plt.xlabel("E/KT")
plt.ylabel("F(E/KT)")
plt.title("PLOT OF PROBABILITY VS E/KT")
plt.grid()
plt.legend()
plt.show()

#plot 2
e_f = 1  #eV
T = [10,100,1000,5000]
e = np.linspace(-4,4,100)

for i in T:
    plt.plot(e,f_fd(e/(k*i) , alpha = -e_f/(k*i)),marker = 'o',label = 'Temp ='+str(i)+'  K')

plt.xlabel("E (in eV)")
plt.ylabel("F")
plt.legend()
plt.title("PROBABILITY PLOT OF FERMI DIRAC FUNCTION FOR DIFFERENT TEMPERATURE")
plt.grid()
plt.show()

#plot 3
U = 1  #eV
T = [10,100,1000,5000]
e = np.linspace(U+10**(-5),4,100)


for i in T:
    x = e/(k*i)
    alpha = -U/(k*i)
    
    plt.plot(e,f_be(x,alpha),marker = 'o',label = 'Temp ='+str(i)+'  K')

plt.xlabel("E (in eV)")
plt.ylabel("F")
plt.legend()
plt.ylim([0,10])
plt.xlim([U,2.5])
plt.title("PROBABILITY PLOT OF BOSE EINSTEIN FUNCTION FOR DIFFERENT TEMPERATURE")
plt.grid()
plt.show()


#plot 4

T = [500,1000,5000,10000]
e = np.linspace(-0.1,0.1,100)

for i in T:
    X = e/(k*i)
    plt.plot(e,f_mb(X),linewidth=4,label = 'Temp ='+str(i)+'  K')

plt.xlabel("E (in eV)")
plt.ylabel("F")
plt.legend()
plt.title("PROBABILITY PLOT OF MAXWELL BOLTZMANN FUNCTION FOR DIFFERENT TEMPERATURE")
plt.grid()
plt.show()
