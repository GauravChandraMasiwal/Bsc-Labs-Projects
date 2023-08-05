
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def factorial(n):      #function to calculate factorial for given n
    prod = 1
    for i in range(2,n+1):
        prod = i*prod
    return prod


def MySinSeries(x,n):   #function for taylor series expansion of sin series
    X = np.array(x)     #arguments are x (a list of values of x ) and n (no.of terms in the series)
    
    if min(abs(X)) > 0 :  #check if the list has zero or not
        factor = 1        #if not ,then i have used the method of calculating the second-last terms of series by multiplying the last terms by a factor
                          #this reduces the run time period of program for higher vals of n,as it uses the factorial func only one time
                          
        nth = (-1)**(n-1)*X**(2*n-1)/factorial(2*n-1)
        sum_sin = np.zeros(len(X))
        while n>0:
            sum_sin =sum_sin + nth*factor
            nth = nth*factor
            n = n-1
            factor = (-1)*2*n*(2*n +1)/X**2
        
        return sum_sin
                            #if there is zero in the array,then we cant take the factor,so i have used the iterative method to expand the series
    else:
        sum_sin = np.zeros(len(X))
        for i in range(1,n+1):
            sum_sin = sum_sin + (-1)**(i-1)*X**(2*i-1)/factorial(2*i-1)
        return sum_sin
        
def MyCosSeries(x,n):
    X = np.array(x)
    
    if min(abs(X)) > 0 :
        factor = 1
    
        nth = (-1)**(n-1)*X**(2*n-2)/factorial(2*n-2)
        sum_cos = np.zeros(len(X))
        while n>0:
            sum_cos =sum_cos + nth*factor
            nth = nth*factor
            n = n-1
            factor = (-1)*2*n*(2*n -1)/X**2
        
        return sum_cos
    else:
        sum_cos = np.zeros(len(X))
        for i in range(1,n+1):
            sum_cos = sum_cos + (-1)**(i-1)*X**(2*i-2)/factorial(2*i-2)
        return sum_cos
#que 1

def graph_skech(x0,X1,Y1,X2,Y2,labels,savef,Type):
    if Type == "SIN":
        fig, axs = plt.subplots(1,2)
        axs[0].plot(X1,np.sin(X1),color = "magenta",label = labels[0])
        axs[1].plot(X2,np.array([np.sin(x0)]*len(X2)),'-ro',label = "sin(\u03C0 /4)")
        axs[1].plot(X2,Y2,'--bo',label = "calculated points")
        for i in range(len(Y1)):
            axs[0].scatter(X1,Y1[i],s=20,label = labels[i+1])
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("f(x)")
        axs[1].set_xlabel("n")
        axs[1].set_ylabel("sin(\u03C0 /4)")
        axs[0].set_title("TAYLOR EXPANSION OF SIN(X)")
        axs[1].set_title("VARIATION OF SIN(\u03C0 /4) WITH n")
        axs[0].grid()
        axs[0].legend()
        axs[1].grid()
        axs[1].legend()
        plt.savefig(savef)
        plt.tight_layout()
        plt.show()
    
    else:
        fig, axs = plt.subplots(1,2)
        axs[0].plot(X1,np.cos(X1),c = "magenta",label = labels[0])
        axs[1].plot(X2,np.array([np.cos(x0)]*len(X2)),'-ro',label = "cos(\u03C0 /4)")
        axs[1].plot(X2,Y2,'--bo',label = "calculated points")
        for i in range(len(Y1)):
            axs[0].scatter(X1,Y1[i],s=20,label = labels[i+1])
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("f(x)")
        axs[1].set_xlabel("n")
        axs[1].set_ylabel("cos(\u03C0 /4)")
        axs[0].set_title("TAYLOR EXPANSION OF SIN(X)")
        axs[1].set_title("VARIATION OF COS(\u03C0 /4) WITH n")
        axs[0].grid()
        axs[0].legend()
        axs[1].grid()
        axs[1].legend()
        plt.savefig(savef)
        plt.tight_layout()
        plt.show()

def part_a():
    m = [1,2,5,10,20]
    ang = np.linspace(-2*np.pi,2*np.pi,50)
    vals1,vals2,vals3,vals4 = [],[],[],[]  
    #here vals1 and vals2 are lists calculated for sin and cos for part a of que 1 respectively
    #and vals3 and vals4 are lists for sin and cos for part b and c of que 1 respectively
    
    for i in m:
        vals1.append(MySinSeries(ang, i))
        vals2.append(MyCosSeries(ang, i))
        
    x_b = [np.pi/4]
    m_b = np.arange(2,21,2)
    for i in m_b:
        vals3.append(MySinSeries(x_b, i)[0])
        vals4.append(MyCosSeries(x_b, i)[0])
    
    label = ["inbuilt","n = 1","n = 2","n = 5","n = 10","n = 20"]
    savefig = ["mp1_plo1.png","mp1_plot2.png"]
    graph_skech(x_b[0], ang, vals1, m_b, vals3, label, savefig[0], 'SIN')    
    graph_skech(x_b[0], ang, vals2, m_b, vals4, label, savefig[1], 'COS')    
  
    return None

#que 2

def modi_sin(x,tol = 0.5e-3):   
    calc = []
    no_terms = []
    
    for i in x:
        n = 2
        while True:
            y1 = MySinSeries([float(i)], n)
            y2 = MySinSeries([float(i)], 2*n)
            
            err = abs(y2[0] - y1[0])
            
            if err > tol:
                n=2*n
            else :
                break
            
        calc.append(y1[0])
        no_terms.append(n)
        
    inbuilt = np.sin(x)
    
    plt.title("PLOT SIN(X) VS X")
    plt.plot(x,inbuilt,c = 'r')
    plt.scatter(x,calc )
    plt.legend(["inbuilt sin func","calculated points"])
    plt.xlabel("x")
    plt.ylabel("sin(x)")
    plt.savefig("mp1_plot3.png")
    plt.grid()
    plt.show()
    print(pd.DataFrame({"x":x,"sin(x)_calc":calc,"n":no_terms,"sin(x)_inbuilt":inbuilt}))
    
    return None

arg = np.arange(0,np.pi+np.pi/8,np.pi/8)
part_a()
modi_sin(arg)
