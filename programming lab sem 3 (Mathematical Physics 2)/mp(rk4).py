import numpy as np
from tabulate import tabulate as tab
import matplotlib.pyplot as plt
def euler(ini_cond,inde_f,func,h):
    time_vect = np.array([ini_cond[0]])  
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)  
    N = int((inde_f - ini_cond[0])/h)  
    s_no = [1,]
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
        s_no.append(i+2)
    return [s_no,time_vect,y_vect]
def rk2(ini_cond,inde_f,func,h):        
    time_vect = np.array([ini_cond[0]])
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    N = int((inde_f - ini_cond[0])/h)  
    s_no = [1,]
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
        s_no.append(i+2)
    return [s_no,time_vect,y_vect]
def rk4(ini_cond,inde_f,func,h):       
    time_vect = np.array([ini_cond[0]])
    a = []
    for i in range(1,len(ini_cond)):
        a.append(ini_cond[i])
    y_vect = np.array(a).reshape(len(ini_cond)-1,1)
    N = int((inde_f - ini_cond[0])/h)  
    s_no = [1,]
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
        s_no.append(i+2)
    return [s_no,time_vect,y_vect]
def func_vector(inde,dep):
    dx_dt = dep[1] + dep[0] -dep[0]**3
    dy_dt = -dep[0]
    return np.array([dx_dt,dy_dt])
def func_vector2(inde,dep):
    di1_dt = -4*dep[0] + 4*dep[1] +12
    di2_dt = -1.6*dep[0]+1.2*dep[1]+4.8                               
    return np.array([di1_dt,di2_dt])
def graph_sketch(indep,dep,title,axis,labels,savefigs,queno):  
    if queno == 1 : 
        fig,axs = plt.subplots(2, 2)
        fig.tight_layout()
        ax = [[0,0],[0,1],[1,0],[1,1]]
        scat = ["*","^","s"]
        for i in range(len(dep)):
            a = ax[i]
            for j in range(len(dep[i])):
                axs[a[0],a[1]].scatter(indep,dep[i][j][0],marker=scat[j],label=labels[0][j][0],s=6) 
                axs[a[0],a[1]].scatter(indep,dep[i][j][1],marker=scat[j],label=labels[0][j][1],s =6) 
            axs[a[0],a[1]].set_xlabel(axis[0][0])
            axs[a[0],a[1]].set_ylabel(axis[0][1])
            axs[a[0],a[1]].set_title(title[i],c="k")
            axs[a[0],a[1]].grid()
            axs[a[0],a[1]].legend(loc = "upper right")
        plt.savefig(savefigs[0])
        plt.show()
        fig,axs = plt.subplots(2, 2)
        fig.tight_layout()
        for i in range(len(dep)):
            a = ax[i]
            for j in range(len(dep[i])):
                axs[a[0],a[1]].scatter(dep[i][j][0],dep[i][j][1],marker=scat[j],label=labels[1][j],s =6)
                axs[a[0],a[1]].set_xlabel(axis[1][0])
            axs[a[0],a[1]].set_ylabel(axis[1][1])
            axs[a[0],a[1]].set_title(title[i],c="k")
            axs[a[0],a[1]].grid()
            axs[a[0],a[1]].legend(loc = "upper right")
        plt.savefig(savefigs[1])
        plt.show()
    else : 
        scat = ["*","^","s"]
        for i in range(len(dep)):
            for j in range(len(dep[i])):
                plt.plot(indep,dep[i][j][0],marker=scat[j],label=labels[0][j][0]) 
                plt.plot(indep,dep[i][j][1],marker=scat[j],label=labels[0][j][1]) 
            plt.xlabel(axis[0][0])
            plt.ylabel(axis[0][1])
            plt.title(title[i],c="k")
            plt.grid()
        plt.legend(loc = "upper right")
        plt.savefig(savefigs[0])
        plt.show()
        for i in range(len(dep)):
            for j in range(len(dep[i])):
                plt.plot(dep[i][j][0],dep[i][j][1],marker=scat[j],label=labels[1][j])
            plt.xlabel(axis[1][0])
            plt.ylabel(axis[1][1])
            plt.title(title[i],c="k")
            plt.grid()
        plt.legend(loc = "upper right")
        plt.savefig(savefigs[1])
        plt.show()        
def q1():
    cond1,cond2,cond3,cond4,cond = [0,0,-1],[0,0,-2],[0,0,-3],[0,0,-4],[0,0,0]
    euler1 = euler(cond1,15, func_vector,0.1);euler2 = euler(cond2,15, func_vector,0.1);euler3 = euler(cond3,15, func_vector,0.1);euler4 = euler(cond4,15, func_vector,0.1)
    rk2_1 = rk2(cond1,15, func_vector,0.1);rk2_2 = rk2(cond2,15, func_vector,0.1);rk2_3 = rk2(cond3,15, func_vector,0.1);rk2_4 = rk2(cond4,15, func_vector,0.1)
    rk4_1 = rk4(cond1,15, func_vector,0.1);rk4_2 = rk4(cond2,15, func_vector,0.1);rk4_3 = rk4(cond3,15, func_vector,0.1);rk4_4 = rk4(cond4,15, func_vector,0.1)
    Euler1 = euler(cond,6, func_vector2,0.1);Rk2_1 = rk2(cond,6, func_vector2,0.1);Rk4_1 = rk4(cond,6, func_vector2,0.1)
    print('\n',' FOR 1st CONDITION [t,x,y] = [0,0,-1]','\n')
    print(tab({"S.NO":euler1[0],"t":euler1[1],"x (EULER)":euler1[2][0],"y (EULER)":euler1[2][1],"x (RK2)":rk2_1[2][0],"y(RK2)":rk2_1[2][1],"x (RK4)":rk4_1[2][0],"y (RK4)":rk4_1[2][1]},headers = 'keys'))
    print('\n',' FOR 2nd CONDITION [t,x,y] = [0,0,-2]','\n')
    print(tab({"S.NO":euler2[0],"t":euler2[1],"x (EULER)":euler2[2][0],"y (EULER)":euler2[2][1],"x (RK2)":rk2_2[2][0],"y(RK2)":rk2_2[2][1],"x (RK4)":rk4_2[2][0],"y (RK4)":rk4_2[2][1]},headers = 'keys'))
    print('\n',' FOR 3rd CONDITION [t,x,y] = [0,0,-3]','\n')
    print(tab({"S.NO":euler3[0],"t":euler3[1],"x (EULER)":euler3[2][0],"y (EULER)":euler3[2][1],"x (RK2)":rk2_3[2][0],"y(RK2)":rk2_3[2][1],"x (RK4)":rk4_3[2][0],"y (RK4)":rk4_3[2][1]},headers = 'keys'))
    print('\n',' FOR 4th CONDITION [t,x,y] = [0,0,-4]','\n')
    print(tab({"S.NO":euler4[0],"t":euler4[1],"x (EULER)":euler4[2][0],"y (EULER)":euler4[2][1],"x (RK2)":rk2_4[2][0],"y(RK2)":rk2_4[2][1],"x (RK4)":rk4_4[2][0],"y (RK4)":rk4_4[2][1]},headers = 'keys'))
    print('\n',' FOR LINEAR ELECTRIC CIRCUIT','\n')
    print(tab({"S.NO":Euler1[0],"t":Euler1[1],"I1 (EULER)":Euler1[2][0],"I2 (EULER)":Euler1[2][1],"I1 (RK2)":Rk2_1[2][0],"I2(RK2)":Rk2_1[2][1],"I1 (RK4)":Rk4_1[2][0],"I2 (RK4)":Rk4_1[2][1]},headers = 'keys'))
    inde,Inde = euler1[1],Euler1[1]
    dep1 = np.array([euler1[2],rk2_1[2],rk4_1[2]]);dep2 = np.array([euler2[2],rk2_2[2],rk4_2[2]]);dep3 =np.array([euler3[2],rk2_3[2],rk4_3[2]])
    dep4 = np.array([euler4[2],rk2_4[2],rk4_4[2]]);dep = np.array([dep1,dep2,dep3,dep4])
    Dep = np.array([[Euler1[2],Rk2_1[2],Rk4_1[2]]])
    TITLE = [["x(0) = 0 and y(0) = -1"],["x(0) = 0 and y(0) = -2"],["x(0) = 0 and y(0) = -3"],["x(0) = 0 and y(0) = -4"]]
    AXIS = [["Time","Values"],["x","y"]]
    LABELS = [[["x(EULER)","y(EULER)"],["x(RK2)","y(RK2)"],["x(RK4)","y(RK4)"]],["EULER","RK2","RK4"]]
    FIG_NAME = ["MP_lab1.png","MP_lab2.png"]
    graph_sketch(inde, dep, TITLE, AXIS, LABELS, FIG_NAME,1)
    TITLE1 = ["t VS I(t) plot","I2(t) VS I1(t) plot"]
    AXIS1 = [["Time","I(t)"],["I1(t)","I2(t)"]]
    LABELS1 = [[["I1(EULER)","I2(EULER)"],["I1(RK2)","I2(RK2)"],["I1(RK4)","I2(RK4)"]],["EULER","RK2","RK4"]]
    FIG_NAME1 = ["MP_lab3.png","MP_lab4.png"]
    graph_sketch(Inde, Dep, TITLE1, AXIS1, LABELS1, FIG_NAME1,2)
q1()