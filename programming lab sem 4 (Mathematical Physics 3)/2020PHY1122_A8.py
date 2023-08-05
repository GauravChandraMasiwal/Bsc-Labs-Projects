
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def func_vector(inde,dep):  #the given function
    
    dy1_dx = dep[1] - dep[2] + inde

    dy2_dx = 3*inde**2

    dy3_dx = dep[1] + np.exp(-inde)
    
    f_vector = np.array([dy1_dx,dy2_dx,dy3_dx])
    
    return f_vector

def analytic(inde):     #the given analytic solution for the given function
    x = np.array(inde)
    y1 = -0.05*(x**5) + 0.25*(x**4) + x + 2 - np.exp(-x)
    
    y2 = x**3 + 1
    
    y3 = 0.25*(x**4) + x - np.exp(-x)
        
    return np.array([y1,y2,y3])



def euler(ini_cond,inde_f,func,N):  #ini_condi is an arrray consisting all the initial conditions
                                    #inde_f is the final value of independent variable and func is the vector function
    
    h = (inde_f - ini_cond[0])/(N)
    time_vect = np.array([ini_cond[0]])  
    a = []
    for i in range(1,len(ini_cond)):        #this is to make the vector for dependent variables
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)  
    for i in range(N):
        m_vect = h*func(time_vect[i],y_vect[:,i])  
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + m_vect
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)       
        y_vect = np.array(Y)
    return [time_vect,y_vect]


def rk2(ini_cond,inde_f,func,N):
    h = (inde_f - ini_cond[0])/(N)        
    time_vect = np.array([ini_cond[0]])
    
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i])
        m2_vect = h*func(time_vect[i] + h,y_vect[:,i]+m1_vect)
        mrk2 = (1/2)*(m1_vect + m2_vect)
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + mrk2
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]

def rk4(ini_cond,inde_f,func,N):
    h = (inde_f - ini_cond[0])/(N)       
    time_vect = np.array([ini_cond[0]])
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    
    for i in range(N):
        m1_vect = h*func(time_vect[i],y_vect[:,i])
        m2_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m1_vect/2))
        m3_vect = h*func(time_vect[i] + (h/2),y_vect[:,i] + (m2_vect/2))
        m4_vect = h*func(time_vect[i] + h,y_vect[:,i]+m3_vect)
        mrk4 = (1/6)*(m1_vect + 2*m2_vect + 2*m3_vect + m4_vect)
        t = time_vect[i] + h
        t_vect = np.append(time_vect,t)
        time_vect = t_vect
        y_next = y_vect[:,i] + mrk4
        Y = []
        for j in range(len(y_vect)):
            y = np.append(y_vect[j],y_next[j])
            Y.append(y)
        y_vect = np.array(Y)
    return [time_vect,y_vect]

#tolerance functions

def euler_tol(ini_cond,inde_f,func,Nf = 10E5,tol = 0.5*10E-3):        #ini_condi is an arrray consisting all the initial conditions
                                        #inde_f is the final value of independent variable and func is the vector function
    
    N = 2                           #initial value of N,we'll then change it such that tolerance is achieved
    trial = 0                       #if somehow tolerence couldn't achieved then we'll chnge its value to 1 
    
    while N <= Nf : 
        nlist = np.array([])       #we will append the final dependent values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):  #i have only 2 possible values i.e N and 2*N
                 #so that we can have different arrays for N and 2*N,otherwise it'd get appended in the original array
            y = euler(ini_cond, inde_f, func, i)
            nlist = np.append(nlist,y[-1][-1])               #we will append the final values for N and 2*N in this array
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
        
        return [y[0],y[1]]
        
    else : 
        print("ERROR!!!  tolerance value is not reached within given grid point limit", '\n',"CAUSE : either tolerance or grid limit is too low")

def rk2_tol(ini_cond,inde_f,func,Nf = 10E5,tol = 0.5*10E-3):        #ini_condi is an arrray consisting all the initial conditions
                                        #inde_f is the final value of independent variable and func is the vector function
    
    N = 2                           #initial value of N,we'll then change it such that tolerance is achieved
    trial = 0                       #if somehow tolerence couldn't achieved then we'll chnge its value to 1 
    
    while N <= Nf : 
        nlist = np.array([])       #we will append the final dependent values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):  #i have only 2 possible values i.e N and 2*N
                 #so that we can have different arrays for N and 2*N,otherwise it'd get appended in the original array
            y = rk2(ini_cond, inde_f, func, i)
            nlist = np.append(nlist,y[-1][-1])               #we will append the final values for N and 2*N in this array
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
        
        return [y[0],y[1]]
        
    else : 
        print("ERROR!!!  tolerance value is not reached within given grid point limit", '\n',"CAUSE : either tolerance or grid limit is too low")

def rk4_tol(ini_cond,inde_f,func,Nf = 10E5,tol = 0.5*10E-3):        #ini_condi is an arrray consisting all the initial conditions
                                        #inde_f is the final value of independent variable and func is the vector function
    
    N = 2                           #initial value of N,we'll then change it such that tolerance is achieved
    trial = 0                       #if somehow tolerence couldn't achieved then we'll chnge its value to 1 
    
    while N <= Nf : 
        nlist = np.array([])       #we will append the final dependent values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):  #i have only 2 possible values i.e N and 2*N
                 #so that we can have different arrays for N and 2*N,otherwise it'd get appended in the original array
            y = rk4(ini_cond, inde_f, func, i)
            nlist = np.append(nlist,y[-1][-1])               #we will append the final values for N and 2*N in this array
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
        
        return [y[0],y[1]]
        
    
    
    else : 
        print("ERROR!!!  tolerance value is not reached within given grid point limit", '\n',"CAUSE : either tolerance or grid limit is too low")




def Graph_Sketch(X,Y1,Y2,Y3,title,savef,label):  #FOR PLOTTING FOR DIFFERENT N
    for i in range(len(X)):
        plt.plot(X[i],Y1[i],label=label[i])
        plt.plot(X[i],Y2[i],label =label[i])
        plt.plot(X[i],Y3[i],label=label[i])
    plt.title(title);plt.legend(loc = 'best')
    plt.savefig(savef)
    plt.grid()
    plt.show()

def Graph(X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,title,savef):    #FOR PLOTTING FOR DIFFERENT X_F
    plt.plot(X,Y1,label="y1(xf = 2.5)") 
    plt.plot(X,Y2,label="y2(xf = 2.5)")
    plt.plot(X,Y3,label="y3(xf = 2.5)")
    plt.plot(X,Y4,label="y1(xf = 5)")
    plt.plot(X,Y5, label="y2(xf = 5)")
    plt.plot(X,Y6, label="y3(xf = 5)")
    plt.plot(X,Y7,label="y1(xf = 7.5)")
    plt.plot(X,Y8, label="y2(xf = 7.5)")
    plt.plot(X,Y9, label="y3(xf = 7.5)")
    plt.plot(X,Y10 , label="y1(xf = 10)")
    plt.plot(X,Y11, label="y2(xf = 10)")
    plt.plot(X,Y12 , label="y3(xf = 10)")
    plt.grid()
    plt.legend()
    plt.xlabel("$y_i$")
    plt.ylabel("$x$")
    plt.title(title)
    plt.savefig(savef)
    plt.show()
    
    
if __name__ == '__main__':
    
    y = [0,1,1,-1]
    x_f = 1
    
    E1 = euler(y, x_f, func_vector, N=50)
    R21 = rk2(y,x_f , func_vector , N=50)
    R41 = rk4(y,x_f , func_vector , N=50)
    
    data1 = {"x":E1[0],"y1 EULER":E1[1][0],"y2 (EULER)":E1[1][1],"y3 (EULER)":E1[1][2]}
    
    data2 = {"x":R21[0],"y1 (RK2)":R21[1][0],"y2 (RK2)":R21[1][1],"y3 (RK2)":R21[1][2]}
    
    data3 = {"x":R41[0],"y1 (RK4)":R41[1][0],"y2 (RK4)":R41[1][1],"y3 (RK4)":R41[1][2]}
    
    print(pd.DataFrame(data1))
    print(pd.DataFrame(data2))
    print(pd.DataFrame(data3))
    
    
    N_list = np.array([10,100,1000,10000])
    
    E2_X,E2_y1,E2_y2,E2_y3 = [],[],[],[]
    R22_y1,R22_y2,R22_y3 = [],[],[]
    R42_y1,R42_y2,R42_y3 = [],[],[]
    
    E_err_1=[];E_err_2=[];E_err_3=[]
    RK2_err_1=[];RK2_err_2=[];RK2_err_3=[]
    RK4_err_1=[];RK4_err_2=[];RK4_err_3=[]
    
    for i in N_list:
        X = euler(y, x_f, func_vector, i)[0]
        EY1 = euler(y, x_f, func_vector, i)[1][0]
        EY2 = euler(y, x_f, func_vector, i)[1][1]
        EY3 = euler(y, x_f, func_vector, i)[1][2]
        R2Y1 = rk2(y, x_f, func_vector, i)[1][0]
        R2Y2 = rk2(y, x_f, func_vector, i)[1][1]
        R2Y3 = rk2(y, x_f, func_vector, i)[1][2]
        R4Y1 = rk4(y, x_f, func_vector, i)[1][0]
        R4Y2 = rk4(y, x_f, func_vector, i)[1][1]
        R4Y3 = rk4(y, x_f, func_vector, i)[1][2]
        
        E2_X.append(X)
        E2_y1.append(EY1)
        E2_y2.append(EY2)
        E2_y3.append(EY3)
        R22_y1.append(R2Y1)
        R22_y2.append(R2Y2)
        R22_y3.append(R2Y3)
        R42_y1.append(R4Y1)
        R42_y2.append(R4Y2)
        R42_y3.append(R4Y3)
    
        E_err_1.append(max(analytic(X)[0]-EY1))
        E_err_2.append(max(analytic(X)[1]-EY2))
        E_err_3.append(max(analytic(X)[2]-EY3))
        RK2_err_1.append(max(analytic(X)[0]-R2Y1))
        RK2_err_2.append(max(analytic(X)[1]-R2Y2))
        RK2_err_3.append(max(analytic(X)[2]-R2Y3))
        RK4_err_1.append(max(analytic(X)[0]-R4Y1))
        RK4_err_2.append(max(analytic(X)[1]-R4Y2))
        RK4_err_3.append(max(analytic(X)[2]-R4Y3))
        

    N_leg = ['N=10','N=100','N=1000','N=10000']
    Graph_Sketch(E2_X, E2_y1, E2_y2, E2_y3, "euler method for different N ", "A8_1.png",N_leg)
    Graph_Sketch(E2_X, R22_y1, R22_y2, R22_y3, "rk2 method for different N", "A8_2.png",N_leg)
    Graph_Sketch(E2_X, R42_y1, R42_y2, R42_y3, "rk4 method for different N", "A8_3.png",N_leg)
    
    fig,(ax1,ax2,ax3) = plt.subplots(3,sharex=True)
    fig.suptitle("LOG(E) vs LOG(N) PLOT",fontsize=10)
    ax1.plot(N_list,E_err_1,label="$y_{1}(Euler)$")
    ax1.plot(N_list,RK2_err_1,label="$y_{1}(RK2)$")
    ax1.plot(N_list,RK4_err_1,label="$y_{1}(RK4)$")
    ax1.set(yscale="log",xscale="log")
    ax1.grid()
    ax1.legend(loc = 'best')

    ax2.plot(N_list,E_err_2,label="$y_{2}(Euler)$")
    ax2.plot(N_list,RK2_err_2,label="$y_{2}(RK2)$")
    ax2.plot(N_list,RK4_err_2,label="$y_{2}(RK4)$")
    ax2.set(ylabel="LOG(E)= LOG(max(|y_ana -y_num|))",xscale="log",yscale="log")
    ax2.grid()
    ax2.legend(loc = 'best')

    ax3.plot(N_list,E_err_3,label="$y_{3}(Euler)$")
    ax3.plot(N_list,RK2_err_3,label="$y_{3}(RK2)$")
    ax3.plot(N_list,RK4_err_3,label="$y_{3}(RK4)$")
    ax3.set(xlabel="LOG(N)",xscale="log",yscale="log")
    ax3.grid()
    ax3.legend(loc = 'best')
    fig.savefig("A8_4.png")
    plt.show()
    
    
    euler_1=euler(y,1,func_vector,50)
    euler_11=euler(y,2.5,func_vector,50)
    euler_12=euler(y,5,func_vector,50)
    euler_13=euler(y,7.5,func_vector,50)
    euler_14=euler(y,10,func_vector,50)
    
    Graph(euler_1[0],euler_11[1][0],euler_11[1][1],euler_11[1][2],euler_12[1][0],euler_12[1][1],euler_12[1][2],euler_13[1][0],euler_13[1][1],euler_13[1][2],euler_14[1][0],euler_14[1][1],euler_14[1][2],"$y_i$ vs x Plot EULER","A8_5.png")

    rk_2 = rk2(y,1,func_vector,50)
    rk2_11=rk2(y,2.5,func_vector,50)
    rk2_12=rk2(y,5,func_vector,50)
    rk2_13=rk2(y,7.5,func_vector,50)
    rk2_14=rk2(y,10,func_vector,50)
    
    Graph(rk_2[0],rk2_11[1][0],rk2_11[1][1],rk2_11[1][2],rk2_12[1][0],rk2_12[1][1],rk2_12[1][2],rk2_13[1][0],rk2_13[1][1],rk2_13[1][2],rk2_14[1][0],rk2_14[1][1],rk2_14[1][2],"$y_i$ vs x Plot RK2","A8_6.png")

    rk_4 = rk4(y,1,func_vector,50)
    rk4_11=rk4(y,2.5,func_vector,50)
    rk4_12=rk4(y,5,func_vector,50)
    rk4_13=rk4(y,7.5,func_vector,50)
    rk4_14=rk4(y,10,func_vector,50)
    
    Graph(rk_4[0],rk4_11[1][0],rk4_11[1][1],rk4_11[1][2],rk4_12[1][0],rk4_12[1][1],rk4_12[1][2],rk4_13[1][0],rk4_13[1][1],rk4_13[1][2],rk4_14[1][0],rk4_14[1][1],rk4_14[1][2],"$y_i$ vs x Plot RK4","A8_7.png")
    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.expand_frame_repr',False)
    
    
    
    