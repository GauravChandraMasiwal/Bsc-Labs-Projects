#NAME : GAURAV CHANDRA
#ROLLNO : 2020PHY1122
#PARTNER NAME : KHUSHI 
#PARTNER ROLLNO : 2020PHY1155

import numpy as np
import matplotlib.pyplot as plt
from MyIntegration import MySimp,MyHermiteQuad
import pandas as pd    #for making tables
from scipy.integrate import quad


e_0 = 0.4

def seq1(x,a,n):  #sequence 1  phi = e^(-x^2/2e0)/(2 pi e0)^1/2
    e = e_0/(2**n)
    x = np.array(x)
    gauss = np.exp(-1*(x-a)**2/(2*e))/np.sqrt(2*np.pi*e)
    
    return gauss

def seq2(x,a,n):   #sequence 2  phi = e^(-|x-a|/e0)/(2 e0)
    e = e_0/(2**n)
    x = np.array(x)
    seq2 = np.exp(-1*abs(x-a)/e)/(2*e)
    
    return seq2

a_1 ,a_2 = 2,-2
x = np.linspace(-5,5,1000)

n= np.array([1,2,3,4,5])
y_1,y_2,y_3,y_4 = [],[],[],[]

for i in range(len(n)):
    y_1.append(seq1(x, a_1, n[i]))
    y_2.append(seq2(x, a_1, n[i]))
    y_3.append(seq1(x, a_2, n[i]))
    y_4.append(seq2(x, a_2, n[i]))
    
#PART A. i
    
fig, axs = plt.subplots(1,2)
width = [5,4,3,2,1]
for i in range(len(n)):
    axs[0].plot(x,y_1[i],linewidth = width[i],label = 'n = '+str(n[i]))
    axs[1].plot(x,y_3[i],linewidth = width[i],label = 'n = '+str(n[i]))

axs[0].set_xlabel("x")
axs[0].set_ylabel("value of dirac delta fucntion")
axs[1].set_xlabel("x")
axs[1].set_ylabel("value of dirac delta function")
fig.suptitle("DIREC DELTA FUNCTION USING SEQUENCES OF FUNCTION 1")
axs[0].set_title("a = "+str(a_1))
axs[1].set_title("a = "+str(a_2))
axs[0].grid()
axs[1].grid()
axs[0].legend(loc='best')
axs[1].legend(loc='best')
plt.savefig("a7_1.png")
plt.tight_layout()
plt.show()

fig, axs = plt.subplots(1,2)
width = [5,4,3,2,1]
for i in range(len(n)):
    axs[0].plot(x,y_2[i],linewidth = width[i],label = 'n = '+str(n[i]))
    axs[1].plot(x,y_4[i],linewidth = width[i],label = 'n = '+str(n[i]))

axs[0].set_xlabel("x")
axs[0].set_ylabel("value of dirac delta fucntion")
axs[1].set_xlabel("x")
axs[1].set_ylabel("value of dirac delta function")
fig.suptitle("DIREC DELTA FUNCTION USING SEQUENCES OF FUNCTION 2")
axs[0].set_title("a = "+str(a_1))
axs[1].set_title("a = "+str(a_2))
axs[0].grid()
axs[1].grid()
axs[0].legend(loc='best')
axs[1].legend(loc='best')
plt.savefig("a7_2.png")
plt.tight_layout()
plt.show()

#PART A. ii

epsilon = e_0/(2**n)

def simp(g,seq,a,n,check = 1,A=1,N=1000,tol = 0.5*10**-6) :  #this is to limit the simpson function by varing the range till tolerance is achieved
    if seq == 1:  #check the type of sequence entered
        f = lambda x:seq1(x, a, n)*g(x)
    elif seq == 2:
        f = lambda x:seq2(x, a, n)*g(x)
    else:
        print("invalid sequence name")
    
    e = e_0/(2**n)
    if check == 1:      #if check is 1,then value of A is obtained for first integral whose value is known
        
        val = 2
        T = abs(val - 1)
        while T > tol:
            val = MySimp(f,-1*e*A,A*e,N)[0]
            
            T=abs(val -1)
            
            A+=1
    
    
    val = MySimp(f,-1*e*A,A*e,N)[0]

    return [val,A]
        
#sequence 1
def hermite_seq(g,seq,a,n,N):
    if seq == 1:
        f = lambda x:seq1(x, a, n)*g(x)*np.exp(x**2)
    elif seq == 2:
        f = lambda x:seq2(x, a, n)*g(x)*np.exp(x**2)
    else:
        print("invalid sequence name")

    return MyHermiteQuad(f, N)

def true_val(f,F,a,n):
    if f==1 and F==0:
        h = lambda x: seq1(x, a, n)*g(x)
    elif f==1 and F==1:
        h = lambda x: seq1(x, a, n)*g1(x)
    elif f==1 and F==2:
        h = lambda x: seq1(x, a, n)*g2(x)
    elif f==2 and F==0:
        h = lambda x: seq2(x, a, n)*g(x)
    elif f==2 and F==1:
        h = lambda x: seq2(x, a, n)*g1(x)
    else:
        h = lambda x: seq2(x, a, n)*g2(x)
    
    return quad(h,-np.inf,np.inf)[0]


g = lambda x:1
g1 = lambda x:(x+1)**2
g2 = lambda x:3*((x-1)/3)**2   #in this case integral is transformed to x_new = 3*x+1 ,such that dx_new = 3 dx,and so x = (x_new -1)/3
 

I11,I12,I13,I21,I22,I23 = [],[],[],[],[],[]     #for storing integration values using simpson
i11,i12,i13,i21,i22,i23 = [],[],[],[],[],[]     #for storing integration values using hermite gauss
T11,T12,T13,T21,T22,T23 = [],[],[],[],[],[]     #for storing true value of integration using quad

vals1 = []
vals2 = []

for i in range(len(n)):
    vals1.append(simp(g, 1, 0, n[i])[1])
    vals2.append(simp(g, 2, 0, n[i])[1])


for i in range(len(n)):
    I11.append(simp(g , 1, 0, n[i],check=2,A = vals1[i])[0])
    I12.append(simp(g1, 1, 0, n[i],check=2,A = vals1[i])[0])
    I13.append(simp(g2, 1, 0, n[i],check=2,A = vals1[i])[0])
    I21.append(simp(g , 2, 0, n[i],check=2,A = vals2[i])[0])
    I22.append(simp(g1, 2, 0, n[i],check=2,A = vals2[i])[0])
    I23.append(simp(g2, 2, 0, n[i],check=2,A = vals2[i])[0])
    i11.append(hermite_seq(g, 1, 0, n[i],300))
    i12.append(hermite_seq(g1, 1, 0, n[i],300))
    i13.append(hermite_seq(g2, 1, 0, n[i],300))
    i21.append(hermite_seq(g, 2, 0, n[i],350))
    i22.append(hermite_seq(g1, 2, 0, n[i],340))
    i23.append(hermite_seq(g2, 2, 0, n[i],350))
    
    T11.append(true_val(1, 0, 0, n[i]))
    T12.append(true_val(1, 1, 0, n[i]))
    T13.append(true_val(1, 2, 0, n[i]))
    T21.append(true_val(2, 0, 0, n[i]))
    T22.append(true_val(2, 1, 0, n[i]))
    T23.append(true_val(2, 2, 0, n[i]))



pd.set_option("display.max_columns",None)
print("DATA FOR SEQUENCE 1")
print()
data1 = {'e':epsilon,'I1_simp':I11,'I2_simp':I12,'I3_simp':I13,'I1_herm':i11,'I2_herm':i12,'I3_herm':i13,'I1_true':T11,'I2_true':T12,'I3_true':T13}
print(pd.DataFrame(data1))   #for creating tables

print()


print("DATA FOR SEQUENCE 2")
print()
data1 = {'e':epsilon,'I1_simp':I21,'I2_simp':I22,'I3_simp':I23,'I1_herm':i21,'I2_herm':i22,'I3_herm':i23,'I1_true':T21,'I2_true':T22,'I3_true':T23}
print(pd.DataFrame(data1))