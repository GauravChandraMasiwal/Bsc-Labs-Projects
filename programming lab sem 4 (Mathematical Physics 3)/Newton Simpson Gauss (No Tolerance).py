from sympy.abc import x
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
import inspect
from sympy import lambdify, simplify, poly, degree

def getval(expr, **kwargs):
    # technique_name = globals()[inspect.stack()[1].function]  # Integration Technique
    func = eval("lambda x:"+ expr)      # Lambdifying the function
    c = 0       # Checks if the user passed an array or defined the limits. 
    
    for key in kwargs.keys():
        if ((key == "x") or (key == "xarr") or (key == "xint")):
            xarr = list(kwargs.get(key))
            a = xarr[0]
            b = xarr[-1]
            n = len(xarr)
            if (("tol" in kwargs.keys()) and (kwargs.get("tol") != None)):
                tol = kwargs.get("tol")
            else:
                tol = 0.5*10**(-8)  # Default tol
            c = 1
            break
        
        elif ((key == "a") or (key == "b")):
            a = kwargs.get("a")
            b = kwargs.get("b")
            break
        
        else:
            raise(KeyError("Bounds not found! Provide Lower and Upper Bounds or the Complete Interval to proceed"))
    
    if(c != 1):
        n = 100             # Default n
        tol = 0.5*10**(-8)  # Default tol
        
        if (("n" in kwargs.keys()) and (kwargs.get("n") != None)):
            n = kwargs.get("n")

        if (("tol" in kwargs.keys()) and (kwargs.get("tol") != None)):
            # n = Tol(technique_name, func, kwargs.get("tol"), a, b)    # Use with "Tol" function. 
            tol = kwargs.get("tol")

        else:
            raise(KeyError("Number of Nodes not found!"))
    
    return(func, a, b, n, tol)

# def Tol(technique, func, tol, a, b): # Unified tolerance: single function for all the integration techniques.
#     n = 100
#     if (tol >= np.abs( (technique(func, a = a, b = b, n = 2*n) - technique(func, a = a, b = b, n = n) ) / ( technique(func, a = a, b = b, n = n))) ):
#         return(n)
#     else:
#         while(tol < np.abs( (technique(func, a = a, b = b, n = 2*n) - technique(func, a = a, b = b, n = n) ) / (technique(func, a = a, b = b, n = n))) ):
#             n *= 2
#     return(n)

def trapz(expr, **kwargs):
    func, a, b, n, tol = getval(expr, **kwargs)
    h = (b-a)/n
    y=[(func(a)+func(b))/2]            
    
    for i in range(1, n):
        y.append(func(a+i*h))       # Value of f(x) at the nodal points
    trap = h*sum(y)
    
    return(trap)

def simps(expr,**kwargs):
    func, a, b, n, tol = getval(expr, **kwargs)
    h = (b-a)/(2*n)
    y = [(func(a)+func(b))/3]
    
    for i in range(1,2*n): 
        if(i%2 == 0):
            y.append(2*func(a+i*h)/3)  # y at Even Nodes (at even indexes) 
            print(2*func(a+i*h)/3)
        
        elif(i%2 == 1):
            y.append(4*func(a+i*h)/3)  # y at Odd Nodes (at odd indexes)
    
    simp = h*sum(y)
    
    return(simp)

# def gausstol(expr, tol, n, m, m_max, xx):
#     if (tol <= np.abs((gaussquad(expr, n, 2*m, xarr = xx) - gaussquad(expr, n, m, xarr = xx))/(gaussquad(expr, n, m, xarr = xx))) ):
#         return(m)
#     elif(m <= m_max):
#         while (tol > np.abs((gaussquad(expr, n, 2*m, xarr = xx) - gaussquad(expr, n, m, xarr = xx))/(gaussquad(expr, n, m, xarr = xx))) ):
#             m *= 2
#         return(m)
#     else:
#         return(ValueError("Tolerance could not be reached for m_max"))

def gaussquad(expression, n = None, m = None, **kwargs):
    func, a, b, *other_params = getval(expression, **kwargs)
    print(other_params)
    a = 0
    b = 2
    if (n == None):
        if (simplify(expression).is_polynomial() == False):
            raise(KeyError("Expression is not a Polynomial. To Approximate the integral, enter n."))
        else:
            expr = poly(expression, x, domain = 'ZZ')
            xvals, weights = leggauss(degree(expr, gen = x))
            xvals = ((b-a)*xvals + (b+a))/2
            w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
        result = (b-a)*sum(w_fx)/2

    elif (n >= 1):
        xvals, weights = leggauss(n)
        tol = other_params[1]
        m_max = 2**10
        if (m != None):
            h = (b - a)/m
            result = 0
            for i in range(m):
                b = a + h
                new_xvals = (h*xvals + (b+a))/2
                w_fx = [weights[j]*func(new_xvals[j]) for j in range(len(weights))]
                result += h*sum(w_fx)/2
                a = b
            # tol_m = gausstol(expression, tol, n, m, m_max, xx = np.linspace(a, b, other_params[0]))
            # print("Tolerance reached at m = ", tol_m)
        else:
            xvals = ((b-a)*xvals + (b+a))/2
            w_fx = [weights[i]*func(xvals[i]) for i in range(len(weights))]
            result = (b-a)*sum(w_fx)/2

    return(result)


if __name__ == '__main__':
    function = str(input("Enter the Function f(x): "))
    xx = np.linspace(0, 1, 50)
    t = gaussquad(expression = function, n = 3, m = 1, xarr = xx, tol = 10e-5)
    print(quad(eval("lambda x:"+'x**3+5'),0,2))
    # t = trapz(func, a = 0, b = 2, n=200)
    # t = simps(func, a = 0, b = 2, tol = 10**(-10))
    print("final = ", t)