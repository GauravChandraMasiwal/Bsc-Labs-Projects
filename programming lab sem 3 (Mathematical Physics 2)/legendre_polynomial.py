#In this program ,I made Gamma , legendre and its derivative function and proving the legendre relations 

#1a GAMMA FACTORIAL
import numpy as np
import matplotlib.pyplot as plt
n = 2
x = 1
def gamma_(n):
    if n==0.5:
        return np.sqrt(np.pi)
    elif n==1:
        return n
    elif n > 0:
        return (n-1)*gamma_(n-1)
    else:
        return "ENTER A POSITIVE VALUE" 
print("The value Gamma Factorial is", gamma_(6))
#1b Legendre Polynomial
def legen(n,x):
    pn=0
    if n%2 == 0:
        m=(n/2)
    else :
        m = ((n-1)/2)
    s = np.arange(0,m+1,1)
    for i in s:
        a1 = (-1)**i*gamma_(2*n-2*i+1)
        a2 = x**(n-2*i)
        a3 = 2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1)
        pn += (a1*a2)/a3
    return pn
print("The value of Legendre Polynomial P",n,'(',x,')', " :  " ,legen(n,x))
print()
#1c Derivative of Legendre Polynomial
def de_legen(n,x):
    pn_d=0
    if n%2 == 0:
        m=(n/2)
    else :
        m = ((n-1)/2)
    s = np.arange(0,m+1,1)
    for i in s:
        a1 = (-1)**i*gamma_(2*n-2*i+1)
        #a2 = x**(n-2*i)
        a3 = 2**n*gamma_(i+1)*gamma_(n-i+1)*gamma_(n-2*i+1)
        a4 = (n-2*i)*x**(n-(2*i)-1)
        pn_d += (a1*a4)/a3
    return pn_d
print("The Differential of Legendre Polynomial P'",n,'(',x,')', " :  " ,de_legen(n,x))
#1d Inbuilt LEGENDRE POLYNOMIALS and the Differentials
from scipy.special import legendre
print()
p_x = legendre(n)
print("The value of P",n,"(x):")
print()
print(p_x) 
print()   #Polynomial type array
#For coefficient type 

#FOR DERIVATIVES
poly = legendre(n)  # coefficients of n^th degree Legendre polynomial 
polyd= poly.deriv() # coefficients of derivative of n^th degree Legendre Polynomial
print("Derivative of Legendre Polynomial P'",n,'(',x,')', ": ", polyd)

#2a
x1=np.linspace(-1,1,100)
poly_val =[]
for i in x1:
           z1=legen(n,i)
           poly_val.append(z1)
data= np.array([x1,legen(n=0,x=x1),legen(n=1,x=x1),legen(n=2,x=x1)],dtype="double")
np.savetxt('eleg00.dat',data.T,delimiter=',',fmt='%.12e')
ldata=np.loadtxt('eleg00.dat',delimiter=',',dtype='double').T

#2b
x2=np.linspace(-1,1,100)
poly_val_d =[]
for i in x2:
            z2=de_legen(n,i)
            poly_val_d.append(z2)
data1= np.array([x2,de_legen(n=0,x=x2),de_legen(n=1,x=x2),de_legen(n=2,x=x2),de_legen(n=3,x=x2)],dtype="double")
np.savetxt('eleg01.txt',data1.T,delimiter=',',fmt='%.12e')
ldata1=np.loadtxt('eleg01.txt',delimiter=',',dtype='double').T

def plot(plt,title=''):
    plt.spines['left'].set_position('zero')
    plt.spines['bottom'].set_position('zero')
    plt.spines['right'].set_color('none')
    plt.spines['top'].set_color('none')
    plt.xaxis.set_ticks_position('bottom')
    plt.yaxis.set_ticks_position('left')
    plt.legend(loc="lower right")
    plt.set_title(title)
    plt.grid()
fig,(plt1,plt2) = plt.subplots(1,2)
for i in range(1,4):
        plt1.plot(ldata[0],ldata[i],label="$P_{0}$".format(i-1))    
for i in [1,3]:
        plt2.plot(ldata1[0],ldata1[i],label="$P'_{0}$".format(i-1))
plt2.plot(ldata[0],ldata[2],label="$P_1$")
plot(plt1,'(a)'),plot(plt2,"(b)")
fig.suptitle("Legendre series", size=16)

#2c(i)
#Relation 1     nPn(x) = xP0n(x) − P0n−1(x)
t=[]
t1=[]
t2=[]
for w in range(100):
    t.append(float(2))
for w in range(100):
    t1.append(float(1))   
for w in range(100):
    t2.append(3) 
p=np.array(t) 
q=np.array(t1)
r=np.array(t2)
for i in range(-1,1,100):

    LHS = n*legen(n=2,x=x1)
    RHS = (x*de_legen(n=2,x=x2)) - de_legen(n=1,x=x2)
    
    if LHS[i] == RHS[i]:
      print("The Relation 1 is hence verified as LHS=RHS")
    else:
      print("The relation is not verified")

np.savetxt('1leg02.dat',(x2,p,q,legen(n=2,x=x1),de_legen(n=2,x=x2),de_legen(n=1,x=x2)),fmt='%.12e')

#2c(ii)
#Relation 2    (2n + 1)xPn(x) = (n + 1)Pn+1 + nPn−1(x)
for i in range(-1,1,100):
    LHS = (2*n+1)*x*legen(n=2,x=x1)
    RHS = ((n+1)*legen(n=3,x=x1) + n*legen(n=1,x=x1))  #n+1
    
    if LHS[i] == RHS[i]:
      print("The Relation 2 is hence verified as LHS=RHS")
    else:
      print("The relation is not verified")  
np.savetxt('2leg03.dat',(x2,p,q,r,legen(n=2,x=x1),legen(n=3,x=x1),legen(n=1,x=x1)),fmt='%.12e')

#2c(iii)     
#Relation 3    nPn(x) = (2n − 1)xPn−1(x) − (n − 1)Pn−2(x)
for i in range(-1,1,100):
    LHS = n*legen(n=3,x=x1)
    RHS = (2*n-1)*x*legen(n=2,x=x1) -(n-1)*legen(n=1,x=x1)
    
    if LHS[i] == RHS[i]:
      print("The Relation 3 is hence verified as LHS=RHS")
    else:
      print("The relation is not verified")
np.savetxt('3leg02.dat',(x2,p,q,r,legen(n=3,x=x2),legen(n=2,x=x2),legen(n=1,x=x2)),fmt='%.12e')


from scipy.integrate import quad
A=[] ; B =[]
for n in range(3):
    for m in range(3):
        if n == m:
           A.append(2/(2*n+1))
        else:
           A.append(0)
        f=legendre(n)*legendre(m)
        inte , err = quad(f, -1, 1)
        B.append(inte)
RHS = np.array(B).reshape(3,3)
LHS = np.array(A).reshape(3,3)
print("RHS = ",RHS)
print("LHS = ",LHS)
if np.allclose(LHS,RHS):
   print("Orthogonality verified")
else:
    print("Orthogonality not verified")
plt.show()
