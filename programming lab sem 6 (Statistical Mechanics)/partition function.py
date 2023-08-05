#name : gaurav chandra
#rollno : 2020phy1122
 
#In this program my main aim is to analyse plots of partition function,internal energy,
# entropy and population denstiy plots for high and low temperature.

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

k = 8.617 * 10**(-5) # eV/K

def Z(g,e,T):
    n = len(e)
    part = []
    
    for j in (T):
        
        z = 0

        for i in range(n):

            z =z+ g[i] * np.exp(-e[i]/(k*j))

        part.append(z)

    return part

def frac(g,e,T):
    n = len(e)
    FRAC = np.zeros(n*len(T) ).reshape(n, len(T))
   
    z = Z(g,e,T)
    
    for i in range(len(T)):

        for j in range(n):
            
            f1 = (g[j]*np.exp(-e[j]/(k*T[i]))) / z[i]
            FRAC[j][i] = f1
        
        
    return FRAC


T1 = np.linspace(0.0001,5000,1000)
T2 = np.linspace(5000,10**6,1000)

#2level
e = [0,1]
g=[1,1]
#3level
e_3 = [0,1,2]
g_3 = [1,1,1]

plt.subplot(2,2,1)
plt.plot(T1,Z(g,e,T1),c='r',label="2lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("Z")
plt.grid()
plt.legend()

plt.subplot(2,2,2)
plt.plot(T2,Z(g,e,T2),c='y',label="2lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("Z")
plt.grid()
plt.legend()

plt.subplot(2,2,3)
plt.plot(T1,Z(g_3,e_3,T1),c='b',label="3lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("Z")
plt.grid()
plt.legend()

plt.subplot(2,2,4)
plt.plot(T2,Z(g_3,e_3,T2),c='violet',label="3lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("Z")
plt.grid()

plt.suptitle("PLOT OF Z VS TEMP")
plt.savefig("ass6_1.png")
plt.legend()
plt.show()

#fractional plot

e_mid = 1/len(e)
e_mid_3 = 1/len(e_3)
E_mid ,E_mid_3= np.full( shape = len(T2) ,fill_value = e_mid),np.full( shape = len(T2) ,fill_value = e_mid_3)

plt.subplot(2,2,1)
for i in range(len(e)):
    plt.plot(T1,frac(g,e,T1)[i], label = "energy = "+str(e[i]) + " eV")

plt.xlabel("TEMPERATURE")
plt.ylabel("N_j / N")
plt.legend(loc=6)
plt.grid()


plt.subplot(2,2,2)
for i in range(len(e)):
    plt.plot(T2,frac(g,e,T2)[i], label = "energy = "+str(e[i]) + " eV")

plt.plot(T2,E_mid,'--')
plt.xlabel("TEMPERATURE")
plt.ylabel("N_j / N")
plt.legend(loc="best")
plt.grid()

plt.subplot(2,2,3)
for i in range(len(e_3)):
    plt.plot(T1,frac(g_3,e_3,T1)[i], label = "energy = "+str(e_3[i]) + " eV")

plt.xlabel("TEMPERATURE")
plt.ylabel("N_j / N")
plt.legend(loc=6)
plt.grid()


plt.subplot(2,2,4)
for i in range(len(e_3)):
    plt.plot(T2,frac(g_3,e_3,T2)[i], label = "energy = "+str(e_3[i]) + " eV")

plt.plot(T2,E_mid_3,'--')
plt.xlabel("TEMPERATURE")
plt.ylabel("N_j / N")
plt.legend(loc="best")
plt.grid()
plt.suptitle("PLOT N_j/N VS TEMPERATURE ")
plt.savefig("ass6_2.png")
plt.show()


#INTERNAL ENERGY

population_1 = frac(g,e,T1)
population_2 = frac(g,e,T2)
population_3 = frac(g_3,e_3,T1)
population_4 = frac(g_3,e_3,T2)
U_1,U_2,U_3,U_4 = 0,0,0,0
for i in range(len(population_1)) :
    U_1 += population_1[i]*e[i]

for i in range(len(population_2)) :
    U_2 += population_2[i]*e[i]

for i in range(len(population_3)) :
    U_3 += population_3[i]*e_3[i]

for i in range(len(population_4)) :
    U_4 += population_4[i]*e_3[i]
    
plt.subplot(2,2,1)
plt.plot(T1,U_1,c='r',label="2lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("U/N")
plt.grid()
plt.legend()

plt.subplot(2,2,2)
plt.plot(T2,U_2,c='y',label="2lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("U/N")
plt.grid()
plt.legend()

plt.subplot(2,2,3)
plt.plot(T1,U_3,c='b',label="3lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("U/N")
plt.grid()
plt.legend()

plt.subplot(2,2,4)
plt.plot(T2,U_4,c='g',label="3lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("U/N")
plt.grid()
plt.legend()
plt.savefig("ass6_3.png")
plt.suptitle("U/N VS TEMPERATURE")
plt.show()


#ENTROPY
z1 = Z(g,e,T1)
z2 = Z(g,e,T2)
z3 = Z(g_3,e_3,T1)
z4 = Z(g_3,e_3,T2)
N = 1

S1 = N*k*np.log(np.array(z1) / N) + U_1 / T1 
S2 = N*k*np.log(np.array(z2) / N) + U_2 / T2 
S3 = N*k*np.log(np.array(z3) / N) + U_3 / T1 
S4 = N*k*np.log(np.array(z4) / N) + U_4 / T2 

plt.subplot(2,2,1)
plt.plot(T1,S1,c='r',label="2lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("ENTROPY")
plt.legend()
plt.grid()

plt.subplot(2,2,2)
plt.plot(T2,S2,c='y',label="2lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("ENTROPY")
plt.legend()
plt.grid()

plt.subplot(2,2,3)
plt.plot(T1,S3,c='b',label="3lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("ENTROPY")
plt.legend()
plt.grid()

plt.subplot(2,2,4)
plt.plot(T2,S4,c='g',label="3lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("ENTROPY")
plt.grid()
plt.suptitle("PLOT ENTROPY VS TEMPERATURE")
plt.legend()
plt.savefig("ass6_4.png")
plt.show()


#HELMHOLTZ
F1 = -N*k*np.array(T1) * np.log(np.array(z1))
F2 = -N*k*np.array(T2) * np.log(np.array(z2))
F3 = -N*k*np.array(T1) * np.log(np.array(z3))
F4 = -N*k*np.array(T2) * np.log(np.array(z4))

plt.subplot(2,2,1)
plt.plot(T1,F1,c='r',label="2lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("HELMHOLTZ FNC")
plt.legend()
plt.grid()

plt.subplot(2,2,2)
plt.plot(T2,F2,c='y',label="2lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("HELMHOLTZ FNC")
plt.grid()
plt.legend()

plt.subplot(2,2,3)
plt.plot(T1,F3,c='b',label="3lvl and low T")
plt.xlabel("TEMPERATURE")
plt.ylabel("HELMHOLTZ FNC")
plt.legend()
plt.grid()

plt.subplot(2,2,4)
plt.plot(T2,F4,c='g',label="3lvl and high T")
plt.xlabel("TEMPERATURE")
plt.ylabel("HELMHOLTZ FNC")
plt.grid()
plt.legend()
plt.suptitle("PLOT F VS TEMP FOR HIGH TEMP")
plt.savefig("ass6_5.png")
plt.show()

result_1,result_2 = stats.linregress(T2,F2),stats.linregress(T2,F4)

print("The slope of the plot of F vs T for high temp for 2lvl system is :",result_1.slope)
print("The value obtained for entropy at high temperature :",np.max(S2))
print('')
print("The slope of the plot of F vs T for high temp for 3lvl system is :",result_2.slope)
print("The value obtained for entropy at high temperature :",np.max(S4))
