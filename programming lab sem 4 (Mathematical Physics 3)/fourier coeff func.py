def FourierCoeff(f,int_name,L,N,args,d=3,var = 0):
    if int_name == "TRAP":  #checking the type of the integration method used
        met = MyTrap_tol
    elif int_name == "SIMP":
        met = MySimp_tol
    elif int_name == "GAUSS":
        met = MyLegQuadrature_tol
    else : 
        print("ENTER VALID METHOD TYPE")
    
    a_n,b_n = np.zeros(N+1),np.zeros(N)  #a_n includes a_0 too,
    
    if var == 1:   #even func
        for i in range(0,N+1):
            const = (i*np.pi)/L
            g = lambda x:f(x) * np.cos(const * x)  #function of two multiplicated functions
            a_n[i] = (2/L)*met(g,0,L,args,tol = 0.5*10**(-d))[0]
        return a_n,b_n
    
    elif var == -1:   #odd func
        for i in range(0,N):
            const = ((i+1)*np.pi)/L
            g = lambda x:f(x) * np.sin(const * x)
            b_n[i] = (2/L)*met(g,0,L,args,tol = 0.5*10**(-d))[0]
        return a_n,b_n
    
    else :     
        for i in range(0,N+1):
            const = (i*np.pi)/L
            g = lambda x:f(x) * np.cos(const * x)  
            a_n[i] = (1/L)*(met(g,-L,L,args,tol = 0.5*10**(-d))[0])
            
        for j in range(0,N):
            const = ((j+1)*np.pi)/L
            g = lambda x:f(x) * np.sin(const * x)
            b_n[j] = (1/L)*met(g,-L,L,args,tol = 0.5*10**(-d))[0]
        
        return a_n,b_n