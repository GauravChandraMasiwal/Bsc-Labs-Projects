
#NAME : GAURAV CHANDRA
#ROLLNO : 2020PHY1122
#COLLG ROLL NO : 20068567021

import numpy as np
import matplotlib.pyplot as plt
from MyIntegration import MyHermiteQuad,MySimp_tol
import pandas as pd
from scipy.integrate import quad

#B PART

def validate(y,n):      #this to verify the gauss hermite for mth order for poly. if n = 2*m-1 
    
    f = eval("lambda x:"+y)
    g = lambda x: f(x)*np.exp(-x**2)
    cal = quad(g, -np.inf,np.inf )[0]
    ans1 = MyHermiteQuad(f, n)
    print("The evaluated integration value of polynomial is : ",ans1)
    print("the absolute error between values evaluated from method and inbuilt func is : ",round(abs(ans1-cal),10))
    if round(abs(ans1-cal),10)==0 :
        print("The computational method gives exact result for given func for no. of interval=",n)
    else :
        print("The method doesn't give exact result")
    
    print("")
    

def limitsimp(f,N=2,N_max = 10**6,tol=0.5e-3) :  #this is to limit the simpson function by varing the range till tolerance is achieved
    b = 10
    trial = 0 
    val = [MySimp_tol(f, -b, b,N,N_max,tol)[0]]
    while True:
        val.append(MySimp_tol(f, -10*b, 10*b,N,N_max,tol)[0])
        
        T = val[-1] - val[-2]
        if abs(T) <= abs(tol*val[-1]) or round(abs(T*val[-1]),10) == 0:     #break the loop if tolerence is achieved
                       
            break
        elif abs(T) > abs(tol*val[-1]) : #if tolerence isn't achieved continue with b = 2*b
            
            b = 10*b
        else :
            trial =1              
            break
    
    
    if trial == 0 :
        return val[-1]
    else :
        msg = "ERROR !!! tolerance is not obtained" 
        print(msg)




if __name__ == "__main__":
    
    #b - i
    
    print("TO VALIDATE COMPUTATIONAL METHOD FOR POLYNOMIAL OF DIFFERENT ORDERS")
    print("for validation of hermite gauss quadrature for n = 2 and n = 4,use order upto 8")
    print("")
    
    y_3 = input("Enter a polynomial in terms of x of order 3:")
    
    y_7 = input("Enter a polynomial in terms of x of order 7:")
    
    print("---:VALIDATION FOR GAUSS HERMITE QUADRATURE METHOD FOR ORDER :3 AND N = 2 :---")
    print("")
    validate(y_3, 2)
    print("")
    print("---:VALIDATION FOR GAUSS HERMITE QUADRATURE METHOD FOR ORDER:7 AND N = 4 :---")
    print("")
    validate(y_7, 4)
     
    #b - ii
    
    n = []
    a = 2
    while a <=128:
        n.append(a)
        a = a*2
    
    f_1 = lambda x:1/(1+x**2)   #for integrand 1
    f_2 = lambda x:np.exp(x**2)/(1+x**2)    #for integrand 2
    
    I_1,I_2 = [],[]
    
    for i in n:
        I_1.append(MyHermiteQuad(f_1, i))
        I_2.append(MyHermiteQuad(f_2, i))
        
    
    data = {'n':n,'I1':I_1,'I2':I_2}
    print(pd.DataFrame(data))   #for creating tables
    
    dat1= np.array([n,I_1,I_2],dtype="double")
    np.savetxt('quad-herm-1122.txt',dat1.T,delimiter=',',fmt='%.12e')   #for saving the table as txt file
    
    #C PART
    
    g_1 = lambda x: np.exp(-x**2)*f_1(x)    #integrand 1 for simpson
        
    I_3 = limitsimp(g_1)
    I_4 = limitsimp(f_1)
    print("")
    print("THE VALUE OF INTEGRATION OF FUNC 1 OBTAINED USING SIMPSON'S RULE IS : ",I_3)
    print("THE VALUE OF INTEGRATION OF FUNC 1 OBTAINED USING SIMPSON'S RULE IS : ",I_4)

    #D PART
    
    fig, axs = plt.subplots(1,2)
    axs[0].plot(n,I_1,'o--',label = 'hermite gauss')
    axs[0].plot(n,[I_3]*len(n),label = 'simpson method')
 
    axs[1].plot(n,I_2,'o--',label = 'hermite gauss')
    axs[1].plot(n,[I_4]*len(n),label = 'simpson method')
    
    axs[0].set_xlabel("n")
    axs[0].set_ylabel("Integration value")
    axs[1].set_xlabel("n")
    axs[1].set_ylabel("Integration value")
    axs[0].set_title("FUNCTION 1")
    axs[1].set_title("FUNCTION 2")
    axs[0].grid()
    axs[1].grid()
    axs[0].legend(loc='best')
    axs[1].legend(loc = 'best')
    plt.savefig("a5_1.pdf")
    plt.tight_layout()
    plt.show()
    
    
    