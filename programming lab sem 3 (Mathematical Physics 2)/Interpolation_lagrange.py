import numpy as np
from scipy.interpolate import lagrange
import pandas as pd
Beta = np.array([0.00, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]) 
J0 = np.array([1.0, 0.99, 0.96, 0.91, 0.85, 0.76, 0.67, 0.57, 0.46, 0.34, 0.22, 0.11, 0.00, -0.10, -0.18, -0.26])   
#Lagrange Interpolation Basis Function
def intrplt(x,y):  
    xp=float(input("Enter the point to be calculated(xp)- "))
    #xp=2.3
    p=0
    for i in range(len(x)):
        Lag=1                       #x beta, j0 bessel
        for j in range(len(x)):
            if j!=i:
                Lag = Lag * (xp -x[j]) / (x[i] - x[j])
        p = p + y[i]*Lag                                      #Lagrange interpolating polynomial
    return p
#print("The value of function at the point is- %.4f"  % intrplt(Beta,J)) 

#using inbuilt function
xp=2.3
Lag1 = lagrange(Beta, J0) 
#print("The value of function using inbuilt function is- %.4f"  % Lag)
 
def invrs_intrplt(x,y):             #Inverse Lagrange Interpolation
    yq=float(input("Enter the point to be calculated(yq)- "))
    #yp=0.5
    s1=0
    for i in range(len(x)):
        Lag=1
        for j in range(len(x)):
            
            if j!=i:
                Lag = Lag * (yq -y[j]) / (y[i] - y[j])
        s1 = s1 + x[i]*Lag                                       
    return s1                                               #Inverse Lagrange interpolating polynomial
#print("The value of function at the point is- %.4f"  % invrs_intrplt(Beta,J))

yq=0.5           #inbuilt function
Lag =lagrange(J0,Beta)
#print("The value of function using inbuilt function is- %.4f"  % Lag(yq))

#Calculating intensity
I=np.array([2.81, 3.24, 3.80, 4.30, 4.37, 5.29, 6.03])
V=np.array([0.5, 1.2, 2.1, 2.9, 3.6, 4.5, 5.7])
print("The Laser Intensity at V= 2.4 is- %.4f"  % intrplt(V,I)) 
            


L1=[intrplt(Beta,J0),invrs_intrplt(Beta,J0)]
L2= float(Lag1(xp)),float (Lag(yq))
data={'Value of function':L1,'Value by inbuilt func':L2}
print(pd.DataFrame(data))

