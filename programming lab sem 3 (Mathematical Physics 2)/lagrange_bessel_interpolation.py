import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import lagrange as lag   #inbuilt lagrange function

#A) BESSEL FUNCTION

beta = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3]
j = [1,0.99,0.96,0.91,0.85,0.76,0.67,0.57,0.46,0.34,0.22,0.11,0,-0.10,-0.18,-0.26]

#B) LINEAR INTERPOLATION

I = [2.81,3.24,3.8,4.3,4.37,5.29,6.03]
V = [0.5,1.2,2.1,2.9,3.6,4.5,5.7]


I2 = [2.81,6.03]    #FOR LINEAR INTERPOLATION,WE SHOULD HAVE ORDER 1
V2 = [0.5,5.7]


#METHOD/CODE

def lagrange(X,Y,x):         #(X,Y) IS THE GIVEN DATA SET AND x IS THE VALUE OF X AT WHICH WE WANT TO FIND Y
    
    y = 0
    
    for i in range(len(X)):
        l = 1
        for k in range(len(X)):
            if k != i :
                l = l*(x - X[k])/(X[i] - X[k])
            else : 
                continue
        y += l*Y[i]   #interpolated value
    f = lag(X,Y)      #inbuilt interpolation
    
    return [y,f(x),f(x)-y]  #[interpolated value,inbuilt interpolated value,relative error]

#APPLICATIONS

print('\t','\t',"-:  a) BESSEL FUNCTION   :-")

#i)
a = 2.3     # a stores the given value for which we want the interpolated value
first = lagrange(beta, j,a)
print('\n',"THE VALUE OF BESSEL FUNCTION AT  GIVEN beta VALUE IS : ",first[0])
print("THE EVALUATED VALUE IS : ",first[1],'\n',"ERROR IS : ",first[2])

x= np.linspace(beta[0],beta[-1],50)  #plotting
y,Y = [],[]
for i in range(len(x)):
    y.append(lagrange(beta,j,x[i])[0])
    Y.append(lagrange(beta,j,x[i])[1])
plt.plot(x,y,'orange',linewidth = 3)  #width is such that both orange and blue curves appear
plt.plot(x,Y,'b--')
plt.scatter(beta,j,c='red')
plt.scatter(a,first[0],c='g')
ann = '('+str(a)+','+str(round(first[0],2))+')'   #to annotate in coordinate form
plt.annotate(ann,xy=(a,first[0]), xytext = (a,first[0]+0.3), arrowprops = dict(facecolor = 'green',shrink = 0.01),)
plt.legend(['interpolated lagrange function','inbuilt lagrange function','given data points','interpolated point'])
plt.title("BESSEL FUNCTION PLOT")
plt.xlabel("BETA")
plt.ylabel("J(BETA)")
plt.grid()
plt.savefig("bessel_plot.png")
plt.show()

#ii)
a = 0.5
second = lagrange(j, beta,a)
print('\n',"THE VALUE OF beta FOR GIVEN VALUE OF BESSEL FUNCTION IS : ",second[0])
print("THE EVALUATED VALUE IS : ",second[1],'\n',"ERROR IS : ",second[2])

y1 =np.linspace(j[0],j[-1],50)   
x1,X1= [],[]
for i in range(len(y1)):
    x1.append(lagrange(j,beta,y1[i])[0])
    X1.append(lagrange(j,beta,y1[i])[1])
plt.plot(y1,x1,c='orange',linewidth=3)
plt.plot(y1,X1,'b--')
plt.scatter(j,beta,c = 'r')
plt.scatter(a,second[0],c='g')
ann = '('+str(a)+','+str(round(second[0],2))+')'
plt.annotate(ann,xy=(a,second[0]), xytext = (a,second[0]+0.5), arrowprops = dict(facecolor = 'green',shrink = 0.01),)
plt.legend(['interpolated inverse lagrange','inbuilt inverse lagrange','given data points','interpolated point'],loc = (0.1,0))
plt.title("INVERSE BESSEL FUNCTION PLOT")
plt.ylabel("BETA")
plt.xlabel("J(BETA)")
plt.grid()
plt.savefig("inverse_plot.png")
plt.show()


#iii)using given data points
print('\n','\t'," -: b) INTERPOLATION FOR GIVEN DATA POINTS OF V AND I :-")
a = 2.4
third = lagrange(V,I,a)
print('\n',"VALUE OF INCIDENT LAZER INTENSITY FOR GIVEN VOLTAGE IS : ",third[0])
print("THE EVALUATED VALUE IS : ",third[1],'\n',"ERROR IS : ",third[2])

x2= np.linspace(V[0],V[-1],50)
y2,Y2 = [],[]
for k in range(len(x2)):
    y2.append(lagrange(V,I,x2[k])[0])
    Y2.append(lagrange(V,I,x2[k])[1])
plt.plot(x2,y2,'orange',linewidth = 3)
plt.plot(x2,Y2,'b--')
plt.scatter(V,I,c='red')
plt.scatter(a,third[0],c='g')
ann = '('+str(a)+','+str(round(third[0],2))+')'
plt.annotate(ann,xy=(a,third[0]), xytext = (a,third[0]-1), arrowprops = dict(facecolor = 'green',shrink = 0.01),)
plt.legend(['interpolation fnc','inbuilt interpolation','given data points','interpolated point'])
plt.title("PLOT OF VOLTAGE VS CURRENT")
plt.ylabel("VOLTAGE")
plt.xlabel("CURRENT")
plt.grid()
plt.savefig("linear_interp.png")
plt.show()


#iv)using data points of order 1
a = 2.4
fourth = lagrange(V2,I2,a)
print('\n','\t'," -: LINEAR INTERPOLATION FOR DATA SET OF ORDER 1 :-")

print('\n',"VALUE OF INCIDENT LAZER INTENSITY FOR GIVEN VOLTAGE IS : ",fourth[0])
print("THE EVALUATED VALUE IS : ",fourth[1],'\n',"ERROR IS : ",fourth[2])

x3= np.linspace(V2[0],V2[-1],50)
y3,Y3 = [],[]
for k in range(len(x3)):
    y3.append(lagrange(V2,I2,x3[k])[0])
    Y3.append(lagrange(V2,I2,x3[k])[1])
plt.plot(x3,y3,'orange',linewidth = 3)
plt.plot(x3,Y3,'b--')
plt.scatter(V2,I2,c='red')
plt.scatter(a,fourth[0],c='g')
ann = '('+str(a)+','+str(round(fourth[0],2))+')'
plt.annotate(ann,xy=(a,fourth[0]), xytext = (a,fourth[0]-1), arrowprops = dict(facecolor = 'green',shrink = 0.01),)
plt.legend(['interpolation fnc','inbuilt interpolation','given data points','interpolated point'])
plt.title("LINEAR INTERPOLATION")
plt.ylabel("VOLTAGE")
plt.xlabel("CURRENT")
plt.grid()
plt.savefig("linear_interp2.png")
plt.show()

print('\n',"THE DIFFERENCE BETWEEN INTERPOLATED VALUE FOR GIVEN DATA SET AND DATA OF ORDER 1 IS :",third[0]-fourth[0])