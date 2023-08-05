#Name=Gaurav
#College Roll No. = 2020PHY1122
#University Roll no. = 20068567021


from MyIntegration import MyLegQuadrature,MyLegQuadrature_tol
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sympy import *

def gamma_(s):            
    n= float(s)
    if n<0 or n%0.5!=0 :
       return ("Value must be either a positive integer or positive half integer.")
    else:
        if n==0.5:
           return np.sqrt(np.pi)
        elif n==1:
           return n
        else:
           return (n-1)*gamma_(n-1) 
def Lege(n):
    x=symbols('x')
    poly=0
    if n%2 == 0:       #even
        m=int(n/2)
    else :             #odd
        m = int((n-1)/2)   
    list1 = np.arange(0,m+1,1)
    for i in list1:
        poly += ( (-1)**i*gamma_(2*n-2*i+1)*(x**(n-2*i)))/(2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1))
    polyn=simplify(poly)
    fx=lambdify(x,polyn,"math")
    return fx

#(a)
#calculation of coefficient

def coeff_cal(func,n,n0=4,m=1):
  x=symbols('x')
  func2=Lege(n)
  def func1(*args):
    return func(*args) * func2(*args)
  Int=((2*n+1)/2)*MyLegQuadrature(func1,-1,1,n0,m)[0]
  return Int

#Expansion of function
def Expansion(func,n):
    l=[]
    x=symbols('x')
    for i in range(0,n):
        func2=Lege(i)
        def term(*args):
            return coeff_cal(func,i,n0=10,m=100) * func2(*args)
        l.append(term(x))
    k=sum(l)
    p_x=simplify(k)
    fx=lambdify(x,p_x,"math")
    return fx


#(b)        
f1=lambda x : 2 + 3*x + 2*x**4 
 
f2=lambda x : np.sin(x)*np.cos(x)

l1=[];l2=[];l3=[]
for i in range(0,5):
     
     g=coeff_cal(f1,i,n0=10,m=10)
     l3.append(g)
     if g!=0:
        pass
     else:
        l1.append(g)    
          
for i in range(0,10):
     g=coeff_cal(f2,i,n0=10,m=10)
     l2.append(g)
      
n=1;tol=0.1e-6
X=np.linspace(-np.pi,np.pi,100)
old=[]
s=Expansion(f2,n)

for x in X:
    old.append(s(x))
new=[]

while n<100:
      n=n+1
      u=Expansion(f2,n)
      for x in X:
          new.append(u(x))
      e=[]
      for a,b in zip(old,new):
           if b< 6e-15:
              err=abs(b-a)
              e.append(err)
           else:
              e.append((b-a)/b)
      maxv=max(e)
      if maxv<=tol:
         break
      else:
         old=new
         new=[]
n_tol=n
z1=["C0","C1","C2","C3","C4"]
z2=["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]
data1={"Coefficient corresponding to nth Legendre Polynomial":z1,"Value of Coefficient":l3}
print(pd.DataFrame(data1))
print()
print("Non-zero coefficients in the expansion of f (x) = 2 + 3x + 2x**4 :",l1)
print()
data2={"Coefficient corresponding to nth Legendre Polynomial":z2,"Value of Coefficient":l2}
print(pd.DataFrame(data2))
print("Number of terms required in the expansion of sin(x)*cos(x) which result in accuracy of 6 significant digits = ",n_tol)

#(c)
x_a=np.linspace(-2,2,100)
n_a=[1,2,3,4,5]
d1=[]
for n in n_a :
    d1.append(Expansion(f1,n))
e1=[];e2=[];e3=[];e4=[];e5=[];e6=[]
for x in x_a:
    e1.append(d1[0](x))
    e2.append(d1[1](x))
    e3.append(d1[2](x))
    e4.append(d1[3](x))
    e5.append(d1[4](x))
    e6.append(f1(x))
    
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Ques 2(c)')
  
ax1.plot(x_a,e1,linestyle='--',label="n=1",color = 'orange')
ax1.plot(x_a,e2,linestyle='--',label="n=2",color = 'red')
ax1.plot(x_a,e3,linestyle='--',label="n=3",color = 'blue')
ax1.plot(x_a,e4,linestyle='--',label="n=4",color = 'green')
ax1.plot(x_a,e5,linestyle='--',label="n=5",color = 'yellow')
ax1.plot(x_a,e6,label="exact",color = 'k')
ax1.grid()
ax1.legend()
ax1.set(xlabel="x",ylabel="f(x)",title="Series calculated for func1")



n_a=[2,4,6,8,10]
d2=[]
for n in n_a :
    d2.append(Expansion(f2,n))
E1=[];E2=[];E3=[];E4=[];E5=[];E6=[]
for x in x_a:
    E1.append(d2[0](x))
    E2.append(d2[1](x))
    E3.append(d2[2](x))
    E4.append(d2[3](x))
    E5.append(d2[4](x))
    E6.append(f2(x))
    
ax2.plot(x_a,E1,linestyle='--',label="n=2",color = 'orange')
ax2.plot(x_a,E2,linestyle='--',label="n=4",color = 'red')
ax2.plot(x_a,E3,linestyle='--',label="n=6",color = 'blue')
ax2.plot(x_a,E4,linestyle='--',label="n=8",color = 'green')
ax2.plot(x_a,E5,linestyle='--',label="n=10",color = 'yellow')
ax2.plot(x_a,E6,label="exact",color='k')
ax2.grid()
ax2.legend()
ax2.set(xlabel="x",ylabel="f(x)",title="Series calculated for func2")
plt.show()

