
# RK4 COMPUTATIONAL METHOD TO SOLVE A SYSTEM OF n FIRST ORDER DIFFERENTIAL EQUAITIONS

import numpy as np


#y_vector(n+1) = y_vector(n) + h * mrk4_vector
#mrk4_vector = (1/6)*(m1_vector + 2*m2_vector + 2*m3_vector + m4_vector)
#m1_vector = h*func(t(n),y_vector(n))
#m2_vector = h*func(t(n)+(h/2),y_vector + (m1_vector/2)
#m3_vector = h*func(t(n)+(h/2),y_vector + (m2_vector/2)
#m4_vector = h*func(t(n)+h,y_vector + m3)

def rk4(ini_cond,inde_f,func,Nf = 10E5,tol = 0.5*10E-3):        #ini_condi is an arrray consisting all the initial conditions
                                                        #inde_f is the final value of independent variable and func is the vector function
    
    time_vect = np.array([ini_cond[0]],dtype=np.float64)       #time_vect makes an array out of initial value of independent variable,we will append next values in this array
    a = []                      
    for i in range(1,len(ini_cond)):                    #this is to make the vector for dependent variables
        a.append([ini_cond[i]])             #a = [[x0][y0][z0]..]
    
    y_vect = np.array(a,dtype=np.float64)             
    N = 30                              #initial value of N,we'll then change it such that tolerance is achieved
       
    trial = 0       #if somehow tolerence couldn't achieved then we'll chnge its value to 1 
    while N <= Nf :     
        nlist = np.array([])     #we will append the final dependent values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):    #i have only 2 possible values i.e N and 2*N
            y = y_vect    #so that we can have different arrays for N and 2*N,otherwise it'd get appended in the original array
            t = time_vect
            for j in range(i):
                H = abs((inde_f - ini_cond[0])/(i))     #H 'd be different for N and 2*N
                m1 = H*func(t,y)
                m2 = H*func(t+(H/2),y + (m1/2))
                m3 = H*func(t+(H/2),y + (m2/2))
                m4 = H*func(t+H,y+m3)
                mrk4 = (1/6)*(m1 + 2*m2 + 2*m3 + m4)
                t = t + H
                y = y + mrk4 
                          
            nlist = np.append(nlist,y)         #we will append the final values for N and 2*N in this array
        nlist = nlist.reshape(2,len(ini_cond)-1)  
                                                    #nlist  = [[x(n)   y(n)   z(n)  ]
                                                    #          [x(2*n) y(2*n) z(2*n)]]
       
        if abs(max((nlist[1]-nlist[0])/nlist[1])) <= tol:     #break the loop if tolerence is achieved
           
            break
        elif abs(max((nlist[1]-nlist[0])/nlist[1])) > tol and 2*N <=Nf:    #if tolerence isn't achieved but 2*N <= Nf then continue with N = 2*N
            
            N = 2*N
        else :
            trial =1              #if tolerence is not achieved within Nf,chnge trial value to 1,and show error! 
            
            break

#trial is just to show if the tolerance is achieved or not,if trial is zero then program will go on to give output otherwise ,there is an error
   
    if trial == 0 :
        h = abs((inde_f - ini_cond[0])/(N))
        for a in range(N):
            m1_vect = h*func(time_vect[a],y_vect[:,a])
            m2_vect = h*func(time_vect[a] + (h/2),y_vect[:,a] + (m1_vect/2))
            m3_vect = h*func(time_vect[a] + (h/2),y_vect[:,a] + (m2_vect/2))
            m4_vect = h*func(time_vect[a] + h,y_vect[:,a]+m3_vect)
            mrk4 = (1/6)*(m1_vect + 2*m2_vect + 2*m3_vect + m4_vect)
            
            t = time_vect[a] + h
            t_vect = np.append(time_vect,t)
            time_vect = t_vect
        
            y_next = y_vect[:,a] + mrk4
            Y = []
            for b in range(len(y_vect)):
                y = np.append(y_vect[b],y_next[b])
                Y.append(y)
        
            y_vect = np.array(Y,dtype=np.float64)
        
        
        return [time_vect,y_vect]
    else :
        print("ERROR!!!  tolerance value is not reached within given grid point limit", '\n',"CAUSE : either tolerance or grid limit is too low")
 