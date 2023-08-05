# Trapezoidal Method And Simpson When Function is given

"""Trapezoidal rule 
I ≈h/2{f(x0) + f(xn) + 2 [f(x1) + f(x2) + ... + f(xn−1)]}


Simpson Rule
 I ≈h/3{f(x0) + f(x2n) + 4 [f(x1) + f(x3) + ... + f(x2n−1)] + 
        2 [f(x2) + f(x4) + ... + f(x2n−2)]}
  where xi = a+ i*h 
"""
from scipy.integrate import quad
import matplotlib.pyplot as plt
# we define a function here to integrate
def f(x):
    return x**2

"""defining recursive function for trapezoidal with 
upper limit b lower limit a and n intervals"""

def trpz(a,b,n):
    # calculating step size
    h = (b - a) / n
    
    trpzint = f(a) + f(b)
    
    for i in range(1,n):
        k = a + i*h
        trpzint = trpzint + 2 * f(k)
    
    # multiply h/2 with the obtained integration to get trapezoidal integration 
    trpzint =trpzint * h/2
    
    return trpzint
"""recursive function for Simpson with 
upper limit b lower limit a and n intervals"""
def simpson(a,b,n):
    # calculating step size
    h = (b - a) / n
    
    simpint = f(a) + f(b)
    
    for i in range(1,n):
        k = a + i*h
        
        if i%2 == 0:
            simpint = simpint + 2 * f(k)
        else:
            simpint = simpint + 4 * f(k)          
    
    # multiply h/2 with the obtained integration to get Simpson integration 
    simpint =simpint * h/3
    
    return simpint
          
# command to get input
l_l = float(input("Enter lower limit of integration: "))
u_l = float(input("Enter upper limit of integration: "))
s_i = int(input("Enter number of sub intervals: "))

# calling the function and get the value
value1 = trpz(l_l, u_l, s_i)
value2 = simpson(l_l,u_l,s_i)
value3 = quad(f,l_l,u_l)

print("Integration result by Trapezoidal method is: %0.6f" % (value1) )
print("Integration result by Simpson method is: %0.6f" % (value2) )
print("Numerical value is: ",(value3) )

n2 = eval(input("enter the list of number n : "))
H = []
for i in n2:
    H.append((u_l-l_l)/i)
    
    
trapz = []

for i in range (len(n2)):
    trapz.append(trpz(l_l, u_l, n2[i]))
    

plt.plot(H,trapz)

simps = []

for i in range (len(n2)):
    simps.append(simpson(l_l, u_l, n2[i]))
plt.scatter(H,simps)

plt.plot(H,simps)
plt.yscale("log")
plt.xlabel("h")
plt.scatter(H,trapz)
plt.scatter(H,simps)