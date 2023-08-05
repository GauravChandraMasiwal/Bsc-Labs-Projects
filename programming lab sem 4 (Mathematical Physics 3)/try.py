import numpy as np
import matplotlib.pyplot as plt
from MyIVP import Eul
from MyIVP import RK2
from MyIVP import RK4
import pandas as pd

def func(x,Y):
    y1,y2,y3=Y
    f1=y2-y3+x
    f2=3*x*x
    f3=y2+np.exp(-x)
    return np.array([f1,f2,f3])


def main1(func,IC,tmin,tmax,N):
    
    rk4=RK4(func,IC,tmin,tmax,N)
    Z=rk4[0].T
    eu=Eul(func,IC,tmin,tmax,N)
    X=eu[0].T
    rk2=RK2(func,IC,tmin,tmax,N)
    Y=rk2[0].T
    t=rk4[1]
    
    return X,Y,Z,t

def plot(X,t,tit):
    for i in range(len(X)):
        plt.plot(t,X[i],label="$y_{0}$".format(i))  
    plt.title(tit)
    plt.xlabel("x")
    plt.ylabel("$y_i$")
    plt.legend()
    plt.grid()


def plot3(X,t,N):
    for i in range(len(X)):
        plt.plot(t,X[i],label="$N={0}$".format(N))  
    
def plot_ana(Ana,t):
    for g,i in zip(Ana,range(len(Ana))):
        plt.plot(t,g(t),label="$Analytic\_y_{0}$".format(i),linestyle="dashed") 
    plt.legend()
    

s1= lambda x : -0.05*x**5 + 0.25*x**4 + x +2 - np.exp(-x)
s2= lambda x: x**3 +1
s3 = lambda x :  0.25*x**4 + x - np.exp(-x)

IC=[1,1,-1];tmax=1;tmin=0;N=100
Ana=[s1,s2,s3]
X,Y,Z,t=main1(func,IC,tmin,tmax,N)
tit1="Euler Method for N={0}".format(N)
tit2="RK2 Method for N={0}".format(N)
tit3="RK4 Method for N={0}".format(N)

print("#-#-#-#-#-#-#-#  EULER'S METHOD FOR N=100 #-#-#-#-#-#-#-#-#")
data1={"x":t,"y1":X[0],"y2":X[1],"y3":X[2]}
#print(pd.DataFrame(data1))

print("#-#-#-#-#-#-#-#  RK2 METHOD FOR N=100 #-#-#-#-#-#-#-#-#")
data2={"x":t,"y1":Y[0],"y2":Y[1],"y3":Y[2]}
#print(pd.DataFrame(data2))

print("#-#-#-#-#-#-#-#  RK4 METHOD FOR N=100 #-#-#-#-#-#-#-#-#")
data3={"x":t,"y1":Z[0],"y2":Z[1],"y3":Z[2]}
#print(pd.DataFrame(data3))

plot(X,t,tit1)
plot_ana(Ana,t)
plt.show()
plot(Y,t,tit2)
plot_ana(Ana,t)
plt.show()
plot(Z,t,tit3)
plot_ana(Ana,t)
plt.show()


N_arr=[10**1,10**2,10**3,10**4,10**5]

E_err_1=[];E_err_2=[];E_err_3=[]
RK2_err_1=[];RK2_err_2=[];RK2_err_3=[]
RK4_err_1=[];RK4_err_2=[];RK4_err_3=[]

for N in N_arr:
    X,Y,Z,t=main1(func,IC,tmin,tmax,N)
    plot3(X,t,N)
    f1=s1(t)
    f2=s2(t)
    f3=s3(t)
    E_err_1.append(max(np.array(f1)-np.array(X[0])))
    E_err_2.append(max(np.array(f2)-np.array(X[1])))
    E_err_3.append(max(np.array(f3)-np.array(X[2])))

plot_ana(Ana,t)   
plt.title("Euler Method for different N")
plt.xlabel("x")
plt.ylabel("$y_i$")
plt.legend()
plt.grid()
plt.show()

for N in N_arr:
    X,Y,Z,t=main1(func,IC,tmin,tmax,N)
    plot3(Y,t,N)
    f1=s1(t)
    f2=s2(t)
    f3=s3(t)
    RK2_err_1.append(max(np.array(f1)-np.array(Y[0])))
    RK2_err_2.append(max(np.array(f2)-np.array(Y[1])))
    RK2_err_3.append(max(np.array(f3)-np.array(Y[2])))

plot_ana(Ana,t)
plt.title("RK2 Method for different N")
plt.xlabel("x")
plt.ylabel("$y_i$")
plt.legend()
plt.grid()
plt.show()


for N in N_arr:
    X,Y,Z,t=main1(func,IC,tmin,tmax,N)
    plot3(Z,t,N)
    f1=s1(t)
    f2=s2(t)
    f3=s3(t)
    RK4_err_1.append(max(np.array(f1)-np.array(Z[0])))
    RK4_err_2.append(max(np.array(f2)-np.array(Z[1])))
    RK4_err_3.append(max(np.array(f3)-np.array(Z[2])))

plot_ana(Ana,t)   
plt.title("RK4 Method for different N")
plt.xlabel("x")
plt.ylabel("$y_i$")
plt.legend()
plt.grid()
plt.show()
#print("#-#-#-#-#-#-# N and E= max(|y_ana -y_num|) values for y0,y1 and y2 for all three methods #-#-#-#-#-#-#")
print()
data={"N":N_arr,"E_y0(Euler)":E_err_1,"E_y0(RK2)":RK2_err_1,"E_y0(RK4)":RK4_err_1,"E_y1(Euler)":E_err_2,"E_y1(RK2)":RK2_err_2,"E_y1(RK4)":RK4_err_2,"E_y2(Euler)":E_err_3,"E_y2(RK2)":RK2_err_3,"E_y2(RK4)":RK4_err_3}
#print(pd.DataFrame(data))
fig, (ax1,ax2,ax3) = plt.subplots(3,sharex=True)
fig.suptitle("log(E) vs log(N) plot",fontsize=15,c="r")
ax1.plot(N_arr,E_err_1,label="$y_{0}(Euler)$")
ax1.plot(N_arr,RK2_err_1,label="$y_{0}(RK2)$")
ax1.plot(N_arr,RK4_err_1,label="$y_{0}(RK4)$")
ax1.set(ylabel="log(E)= log(max(|y_ana -y_num|))",yscale="log",xscale="log")
ax1.grid()
ax1.legend()

ax2.plot(N_arr,E_err_2,label="$y_{1}(Euler)$")
ax2.plot(N_arr,RK2_err_2,label="$y_{1}(RK2)$")
ax2.plot(N_arr,RK4_err_2,label="$y_{1}(RK4)$")
ax2.set(ylabel="log(E)= log(max(|y_ana -y_num|))",xscale="log",yscale="log")
ax2.grid()
ax2.legend()

ax3.plot(N_arr,E_err_3,label="$y_{2}(Euler)$")
ax3.plot(N_arr,RK2_err_3,label="$y_{2}(RK2)$")
ax3.plot(N_arr,RK4_err_3,label="$y_{2}(RK4)$")
ax3.set(xlabel="log(N)",ylabel="log(E)= log(max(|y_ana -y_num|))",xscale="log",yscale="log")
ax3.grid()
ax3.legend()
plt.show()

print(E_err_1)