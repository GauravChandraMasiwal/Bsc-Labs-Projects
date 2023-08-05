#program for canonical ensemble maxwell boltzmann function




import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import stats

def forward_d(x,y):                  
    derive=[]
    for i in range(len(x)-1):
        dy=y[i+1]-y[i]
        dx=x[i+1]-x[i]
        derive.append(dy/dx)
    return np.array(derive)

k = 1.38*10**(-23)
h = 6.626*10**(-34)
N_a = 6.022*10**(23)
m = 1.6*10**(-27)
V = np.linspace(20*10**(-3),50*10**(-3),50)
T = np.linspace(150,450,50)

matrix = np.zeros(len(V)*len(T)).reshape(len(T),len(V))

for i in range(len(V)):
  for j in range(len(T)):
    v = V[i] 
    t = T[j]
    
    def z(n):
      z = (np.pi/2) * (n**2) * np.exp(-h**2 * n**2/(8*m*v**(2/3)*k*t))
      return z

    I = integrate.quad(z,0,10**(11))[0]
  
    matrix[i][j] = I
    
log_z = np.log(matrix)

#print(matrix)

# plot of log z vs temp
fig=plt.figure()
fig.set_figheight(6)
fig.set_figwidth(10)
plt.subplot(1,2,1)
plt.title("plot of log(z) vs T ")
plt.scatter(T,log_z[:,0],label = "for V= "+str(np.round(V[0],3)))
plt.scatter(T,log_z[:,4],label = "for V= "+str(np.round(V[4],3)))
plt.scatter(T,log_z[:,8],label = "for V= "+str(np.round(V[8],3)))
plt.xlabel("T") 
plt.ylabel("log(Z)")
plt.grid()
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.title("plot of log(z) vs log(T) ")
plt.scatter(np.log(T),log_z[:,0],label = "for log(V)= "+str(np.round(np.log(V[0]),3)))
plt.scatter(np.log(T),log_z[:,4],label = "for log(V)= "+str(np.round(np.log(V[4]),3)))
plt.scatter(np.log(T),log_z[:,8],label = "for log(V)= "+str(np.round(np.log(V[8]),3)))
plt.xlabel("log(T)") 
plt.ylabel("log(Z)")
plt.grid()
plt.legend(loc='best')
plt.show()
plt.savefig('fig7_1.jpeg',dpi=800)

# plot of log z vs volume
fig=plt.figure()
fig.set_figheight(6)
fig.set_figwidth(11)
plt.subplot(1,2,1)
plt.title("plot of log(z) vs V ")
plt.scatter(V,log_z[0],label = "for T= "+str(np.round(T[0],3)))
plt.scatter(V,log_z[4],label = "for T= "+str(np.round(T[4],3)))
plt.scatter(V,log_z[8],label = "for T= "+str(np.round(T[8],3)))
plt.xlabel("V") 
plt.ylabel("log(Z)")
plt.grid()
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.title("plot of log(z) vs log(V) ")
plt.scatter(np.log(V),log_z[0],label = "for T= "+str(np.round(T[0],3)))
plt.scatter(np.log(V),log_z[4],label = "for T= "+str(np.round(T[4],3)))
plt.scatter(np.log(V),log_z[8],label = "for T= "+str(np.round(T[8],3)))
plt.xlabel("log(V)") 
plt.ylabel("log(Z)")
plt.grid()
plt.legend(loc='best')
plt.show()
plt.savefig('fig7_2.jpeg',dpi=800)

# pressure matrix
pressure = []
for i in range(len(T)):
    der = forward_d(V, log_z[i])
    P = N_a*k*T[i] * der
    pressure.append(P)
pressure = np.array(pressure).reshape(len(T),len(V)-1)
#print(pressure)

fig=plt.figure()
fig.set_figheight(6)
fig.set_figwidth(12.5)
plt.subplot(1,2,1)
plt.title("plot of Pressure vs Volume ")
plt.scatter(V[:len(V)-1],pressure[0],label = "for T= "+str(np.round(T[0],3)))
plt.scatter(V[:len(V)-1],pressure[4],label = "for T= "+str(np.round(T[4],3)))
plt.scatter(V[:len(V)-1],pressure[8],label = "for T= "+str(np.round(T[8],3)))
plt.xlabel("V") 
plt.ylabel("Pressure")
plt.grid()
plt.legend(loc='best')

plt.subplot(1,2,2)
plt.title("plot of Pressure vs Temperature ")
plt.scatter(T,pressure[:,0],label = "for V= "+str(np.round(V[0],3)))
plt.scatter(T,pressure[:,4],label = "for V= "+str(np.round(V[4],3)))
plt.scatter(T,pressure[:,8],label = "for V= "+str(np.round(V[8],3)))
plt.xlabel("Temperature") 
plt.ylabel("Pressure")
plt.grid()
plt.legend(loc='best')
plt.show()
plt.savefig('fig7_3.jpeg',dpi=800)
# energy matrix
cv=[]
for i in range(3):
    energy = []
    der = forward_d(T, log_z[:,i])
    
    #energy.append(der)
    for j in range(len(T)-1):
        energy.append(k*T[j]**2 *der[j])
    plt.scatter(T[:len(T)-1],energy,label = "for V= "+str(np.round(V[i],4)))
    cv.append(stats.linregress(T[:len(T)-1],energy)[0])
    
plt.title("plot of Energy vs Temperature ")
plt.xlabel("Temperature") 
plt.ylabel("energy")
plt.grid()
plt.legend(loc='best')
plt.show()
plt.savefig('fig7_4.jpeg',dpi=800)


print('the specific heat of this ideal gas obtained is ',np.average(cv))


