#In this program , I have made Trapeziodal and Simpson function to integrate a given function

import matplotlib.pyplot as plt
from scipy.integrate import quad

initial = int(input("enter the start limit : "))
end = int(input("enter the end limit : "))
N =  int(input("enter the number N : "))

def func(x):
    return x**2

def trapezoid_method(A,B,N):
    sum1 = 0
    t = (B-A)/N
    for i in range(1,N):
        sum1 += func(A +(i*t))
    INT = ((t/2) * ((func(A)+func(B)) + 2 * sum1))
    return INT
    
def simpson_method(A,B,N):
    sum_odd = 0
    sum_even = 0
    t = (B-A)/(2*N)
    for i in range(1,2*N,2):
        sum_odd += func(A + (i*t))
    for j in range(2,2*N,2):
        sum_even += func(A + (j*t))
    INT2 = (t/3) * ((func(A)+func(B)) + 4*sum_odd + 2*sum_even)
    return INT2

L1,L2,L3 = trapezoid_method(initial, end, N), simpson_method(initial, end, N), quad(func, initial,end)

Integration, Integration2, T= [],[],[]
N2 = eval(input("enter the list of number N in list format : "))

for i in N2 :
    T.append((end-initial)/i)
    
for i in range (len(N2)):
    Integration.append(trapezoid_method(initial, end,N2[i]))
    Integration2.append(simpson_method(initial,end, N2[i]))
   
plt.title("Log Graph of t VS I(t)")
plt.plot(T,Integration)
plt.plot(T,Integration2)
plt.legend(['trapezoidal', 'simpson'], loc = 'best')
plt.xscale('log')
plt.yscale("log")
plt.xlabel("t in log")
plt.ylabel('integration for t in log')
plt.grid()
plt.scatter(T,Integration)
plt.scatter(T,Integration2)
plt.savefig('integration.png')
plt.show()

print("-:",'\t', "INTEGRATION OF FUNCTION WIT LIMITS FROM", initial,'to',end,':-')
print("numerical value of integration is : ", L3)
print("integration result using trapezoidal method is : ",L1)
print("integration reult using simpson method is : ",L2)
print("error with trapezoidal method is : ",L1-L3[0])
print("error with simpson method is : ",L2-L3[0])

            

    



