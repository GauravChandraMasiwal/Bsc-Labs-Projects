
# EULER COMPUTATIONAL METHOD TO SOLVE A SYSTEM OF n FIRST ORDER DIFFERENTIAL EQUAITIONS

#y_vector(n+1) = y_vector(n) + h * f(t(n),y_vector(n))

import numpy as np

def euler(ini_cond,inde_f,func,Nf = 10E5,tol = 0.5*10E-3):        #ini_condi is an arrray consisting all the initial conditions
                                        #inde_f is the final value of independent variable and func is the vector function
    time_vect = np.array([ini_cond[0]],dtype=np.float64)
    a = []
    for i in range(1,len(ini_cond)):    #this is to make the vector for dependent variables
        a.append([ini_cond[i]])
    y_vect = np.array(a,dtype=np.float64)
    N = 30                           #initial value of N,we'll then change it such that tolerance is achieved
    trial = 0                       #if somehow tolerence couldn't achieved then we'll chnge its value to 1 
    
    while N <= Nf : 
        nlist = np.array([])       #we will append the final dependent values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):  #i have only 2 possible values i.e N and 2*N
            y = y_vect     #so that we can have different arrays for N and 2*N,otherwise it'd get appended in the original array
            t = time_vect
            for j in range(i):
                H = abs((inde_f - ini_cond[0])/(i))    #H 'd be different for N and 2*N
                m1_vect = H*func(t,y)
                t = t + H
                y = y + m1_vect 
                          
            nlist = np.append(nlist,y)               #we will append the final values for N and 2*N in this array
        nlist = nlist.reshape(2,len(ini_cond)-1)         #nlist  = [[x(n)   y(n)   z(n)  ]
                                                         #          [x(2*n) y(2*n) z(2*n)]]
       
        if abs(max((nlist[1]-nlist[0])/nlist[1])) <= tol:
            
            break
        elif abs(max((nlist[1]-nlist[0])/nlist[1])) > tol and 2*N <=Nf:         #if tolerence isn't achieved but 2*N <= Nf then continue with N = 2*N
            
            N = 2*N
        else :
            trial = 1
            
            break
#trial is just to show if the tolerance is achieved or not,if trial is zero then program will go on to give output otherwise ,there is an error
    if trial == 0 :
        h = abs((inde_f - ini_cond[0])/(N))  
        for i in range(N):
            
            m1 = h*func(time_vect[i],y_vect[:,i])   #m1 would be a vector too

            T = time_vect[i] + h
            t_vect = np.append(time_vect,T)
            time_vect = t_vect
        
            y_next = y_vect[:,i] + m1
            Y = []
            for j in range(len(y_vect)):
                y = np.append(y_vect[j],y_next[j])
                Y.append(y)
        
            y_vect = np.array(Y,dtype=np.float64) 
        return [time_vect,y_vect]
        
    else : 
        print("ERROR!!!  tolerance value is not reached within given grid point limit", '\n',"CAUSE : either tolerance or grid limit is too low")
