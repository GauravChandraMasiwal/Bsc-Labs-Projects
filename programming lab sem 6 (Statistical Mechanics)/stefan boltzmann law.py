#name:gaurav chandra
#name : 2020phy1122

#Program of Stefan boltzmann law

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy import stats


h = 6.626*10**(-34)
c = 3*10**(8)
k = 8.61733 *10**(-5)* 1.6 *10**(-19)


x = np.linspace(0.01,12,100)

def f_p(a):
    return (a**3)/(np.exp(a)-1)

plt.plot(x,f_p(x))
plt.xlabel("x")
plt.ylabel("f_p(x) = x**3/(e**x -1)")
plt.grid()
plt.title("PLOT OF F_p VS X FOR PEAK")
plt.savefig("fig_1_a5")
plt.show()

i = list(f_p(x)).index(max(f_p(x)))
xp = x[i]
print("the value of peak(dimentionless) is : ",xp)

b = (h*c)/(k*xp)
print("the value of b is : ",b)


#part b

inte = integrate.quad(f_p,0.1,100)[0]
print("The value of integration obtained using python is ",inte)

inte_cal = (np.pi**4 /15)

print("The value of integration obtained using the numerical method is ",inte_cal)

def C(T) :
    l = h*c/(k*T)
    C = 8*np.pi*(k*T)/l**3
    return C

Temp = np.arange(100,10000,500)
C_t = C(Temp)
#print(C_t)

U = inte*(C_t)
F = c*U/4



plt.plot(Temp,F,'o-')
plt.xlabel("TEMPERATURE")
plt.ylabel("RADIANT FLUX")
plt.title("PLOT OF RADIANT FLUX VS TEMPERATURE")
plt.grid()
plt.savefig("fig_2_a5")
plt.show()


plt.plot(np.log(Temp),np.log(F),'o-')
plt.xlabel("TEMPERATURE")
plt.ylabel("RADIANT FLUX")
plt.title("LOG PLOT OF RADIANT FLUX VS TEMPERATURE")
plt.xlim([-1,10])
plt.ylim([-20,20])
plt.grid()
plt.savefig("fig_3_a5")
plt.show()

result = stats.linregress(np.log(Temp),np.log(F))
print("THE SLOPE IS :",result.slope,"AND THE INTERCEPT IS :",result.intercept)

sigma_cal = c*8*np.pi**5*k**4/(4*15*(c*h)**3)
sigma_num = np.exp(result.intercept)

print("THE VALUE OF SIGMA CALCULATED FROM STEFAN BLOTZMANN LAW IS :",sigma_cal," Wm**-2 * T**-4")
print("THE VALUE OF SIGMA CALCULATED FROM THE PLOT IS :",sigma_num," Wm**-2 * T**-4")
if np.round(sigma_num,10) == np.round(sigma_cal,10):
    print("AS THE SIGMA CALCULATED AND SIGMA FROM PLOT ARE EQUAL,SO THE STEFAN BOLTZMANN LAW IS PROVED")

#to calculate mean point


a = 0.001
b = 12
x = np.linspace(a,b,11)
for i in x:
    int_l = integrate.quad(f_p,a,i)[0]
    int_r = integrate.quad(f_p,i,b)[0]

    if abs(int_l - int_r)/int_r <= 0.1:
        x_mean = i
        print("THE X VALUE WHICH DIVIDES THE AREA IN EQUAL PARTS IS :",x_mean)
        break
    

    
b_mean = (h*c)/(k*x_mean)
print("THE MEAN B VALUE IS ;",b_mean)
