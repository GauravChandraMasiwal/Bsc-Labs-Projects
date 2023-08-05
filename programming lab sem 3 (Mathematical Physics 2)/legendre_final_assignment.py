import matplotlib.pyplot as plt
import numpy as np
from scipy.special import eval_legendre,legendre       #inbuilt func for legendre
from scipy.integrate import quad     #inbuilt func for integration
 
def gamma_func(n):     #gamma function
    gamma_one = 1
    gamma_half = np.sqrt(np.pi)
    if n == 1:
        return gamma_one
    
    elif n == 0.5:
        return gamma_half
    
    else :
        return (n-1)*gamma_func(n-1)

#legendre function

def legen_poly(n,x):
    m = 0
    if (n %2) == 0 :     #check for even or odd
        m = int(n/2)
    else : 
        m = int((n-1)/2)
    p = 0
    for i in range(m+1):
        p += (((-1)**i) * (gamma_func(2*n-2*i +1))* (x**(n-2*i)))/((2**n) * gamma_func(i+1) * gamma_func(n-i+1) * gamma_func(n-2*i +1))    
    return p

#legendre differentiation function

def legen_diff(n,x):
    m = 0
    if (n %2) == 0 :
        m = int(n/2)
    else : 
        m = int((n-1)/2)
    p = 0
    for i in range(m+1):
        p += ((n-2*i)*((-1)**i) * (gamma_func(2*n-2*i +1))* (x**(n-2*i-1)))/((2**n) * gamma_func(i+1) * gamma_func(n-i+1) * gamma_func(n-2*i +1))
    return p
n = float(input("enter the positive integer n : "))
x = float(input("enter the value of x : "))

P = legen_poly(n, x)
P_derv = legen_diff(n, x)

inbuilt_P = eval_legendre(n,x)  #comparison with inbuilt function
p2 = legendre(n)
p2_derv =p2.deriv()
evaluated_derv = np.polyval(p2_derv,x)

error = inbuilt_P - P    #errors
error2  = evaluated_derv - P_derv

print('the value of legendre function for n = ',n,'x =',x,'is : ',P)  
print('the value of differentiation legendre function for n = ',n,'x =',x,'is : ',P_derv)  
print('the evaluated value of legendre function  for n = ',n,'x =',x,'is : ',inbuilt_P)   
print('the evaluated value of differentiation legendre function for n = ',n,'x =',x,'is : ',evaluated_derv)
print('the error in legendre function for n = ',n,'x =',x,'is :', error)
print("the error in differentiation of legendre function for n = ",n,'x =',x,'is :' ,error2)

X = np.linspace(-0.999999999,0.999999999,100)      #for x in (-1,1)

data= np.array([X,legen_poly(0,X),legen_poly(1,X),legen_poly(2,X)],dtype="double")
header = "       X       ,      P(0)      ,     P(1)    ,     P(2)  "       #spaces so that to kept headings in the same column as data
data2= np.array([X,legen_diff(0,X),legen_diff(1,X),legen_diff(2,X),legen_diff(3,X)],dtype="double")
header2= "    X       ,      P'(0)     ,     P'(1)   ,     P'(2)  ,   P'(3)"
np.savetxt('D:\python work\progg class\\leg00.dat',data.T,header=header,delimiter=' , ',fmt='%1.10f')
np.savetxt('D:\python work\progg class\\leg01.dat',data2.T,header=header2,delimiter=' ,',fmt='%1.10f')
retdata=np.loadtxt('D:\python work\progg class\\leg00.dat',delimiter=' , ',dtype='double').T
retdata1=np.loadtxt('D:\python work\progg class\\leg01.dat',delimiter=',',dtype='double').T

#plotting graph
plt.title("PLOT OF X VS P(n) ",c= 'm')                                   
plt.plot(retdata[0],retdata[1])
plt.plot(retdata[0],retdata[2])
plt.plot(retdata[0],retdata[3])
plt.legend(['P(0)','P(1)','P(2)'],loc = 'best')
plt.xlabel('X')
plt.ylabel('P(n)')
plt.savefig('legen_plot.png')
plt.grid()
plt.show()

plt.title("PLOT OF X VS DIFF P(n)", c = 'r')
plt.plot(retdata[0],retdata[2])
plt.plot(retdata[0],retdata1[1])
plt.plot(retdata[0],retdata1[3])
plt.legend(['P1',"P'(0)","P'(2)"],loc = 'best')
plt.xlabel('X')
plt.ylabel("P'(n)")
plt.savefig('legen_diff.png')
plt.grid()
plt.show()

#2.c     
n = []
n_1 = []
N = []
N_1 = []      #for n=2,n-1=1,n=3,n+1=4   lists to store in data files

for i in range(len(retdata[0])):
    n.append(2)
    n_1.append(1)
    N.append(3)
    N_1.append(4)
    
def check(l,r):
    if np.allclose(l,r,10):
        print("R.H.S IS EQUAL TO L.H.S  \n hence verified")
    else :
        print("not verified \n ")  
    
# relation 1
lhs,rhs = [],[]
print("TO PROVE : n*Pn(x) = x*P'n(x) - P'n-1(x) ",'\n','FOR n=2')
print('L.H.S')
for i in range(len(retdata[0])):
    lhs.append(n[i]*retdata[3][i])
    rhs.append(retdata[0][i]*retdata1[3][i] - retdata1[2][i])
#print(lhs)
print(' ')
print("R.H.S")
#print(rhs)


check(lhs,rhs)
data3= np.array([X,n,n_1,retdata[3],retdata1[3],retdata1[2],lhs,rhs],dtype="double")
header3 = "      X        ,        n     ,       n-1    ,      P(n)   ,      P'(n)  ,     P'(n-1)  ,      L.H.S  ,  R.H.S  "
np.savetxt('D:\python work\progg class\\leg02.dat' ,data3.T,header=header3,delimiter=' , ',fmt='%1.10f')

# relation 2
lhs2,rhs2 = [],[]
print("TO PROVE : (2n+1)*x*Pn(x) = (n+1)*Pn+1(x) + n*Pn-1(x) ",'\n','FOR n=2')
print("L.H.S")
for i in range(len(retdata[0])):
    lhs2.append((2*n[i]+1)*X[i]*retdata[3][i])
    rhs2.append(legen_poly(3,retdata[0][i]) * N[i] + n[i]*retdata[2][i])
#print(lhs2)
print(' ')
print("R.H.S")
#print(rhs2)
check(lhs2,rhs2)

data4= np.array([X,n,n_1,N,retdata[3],retdata[2],legen_poly(3,X),lhs2,rhs2],dtype="double")
header4 = "   X    ,         n    ,      n-1      ,       n+1     ,       P(n)   ,      P(n-1)    ,     P(n+1)   ,      L.H.S    ,     R.H.S  "
np.savetxt('D:\python work\progg class\\leg03.dat' ,data4.T,header=header4,delimiter=' , ',fmt='%1.10f')

# relation 3
lhs3,rhs3 = [],[]
print("to prove : n*Pn(x) = (2n-1)*x*Pn-1(x) - (n-1)*Pn-2(x) ",'\n','FOR n=3')
print("L.H.S")
for i in range(len(retdata[0])):
    lhs3.append(N[i]*legen_poly(3,retdata[0][i]))
    rhs3.append((2*N[i]-1)*retdata[0][i]*retdata[3][i] - n[i]*retdata[2][i])
#print(lhs3)
print(' ')
print("R.H.S")
#print(rhs3)
check(lhs3,rhs3)

data5= np.array([X,N,n,N_1,legen_poly(3,X),retdata[3],retdata[2],legen_poly(4,X),lhs3,rhs3],dtype="double")
header5 = "   X   ,          n    ,      n-1      ,      n+1     ,      P(n)   ,      P(n-1)    ,     P(n-2)  ,       P(n+1)      ,     L.H.S    ss ,   R.H.S  "
np.savetxt('D:\python work\progg class\\leg04.dat' ,data5.T,header=header5,delimiter=' , ',fmt='%1.10f')


#ORTHOGONALITY PROPERTY


FIRST, SECOND =[] , []
for n in range(3):
    for m in range(3):
        if n == m:
           FIRST.append(2/(2*n+1))
        else:
           FIRST.append(0)
        fnc=legendre(n)*legendre(m)
        inte , err = quad(fnc, -1, 1)
        SECOND.append(inte)
RHS = np.array(SECOND).reshape(3,3)
LHS = np.array(FIRST).reshape(3,3)
print("RHS : \n ",RHS)
print("LHS : \n ",LHS)
if np.allclose(LHS,RHS):
   print("Orthogonality verified")
else:
    print("Orthogonality not verified")
plt.show()