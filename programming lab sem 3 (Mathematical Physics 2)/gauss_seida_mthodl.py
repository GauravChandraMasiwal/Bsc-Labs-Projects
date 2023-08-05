

import numpy as np
import sys


#gauss-seidal method

def GAUSS_SEIDAL(A,B,error = 0.0005):
    arr = np.array(A)
    arr2 = np.array(B)
    n = len(arr);m = len(arr[0])  #N,M are no.of rows and columns of coefficient matrix
    
    Diag = np.diag(np.abs(arr))
    
    # Find row sum without diagonal
    Diag_new = np.sum(arr, axis=1) - Diag
    if np.all(Diag > Diag_new):
        print("MATRIX IS DIAGONALLY DOMINANT")
        print("----------------------------------")
        
        X = np.array([ ])
        count = 0
        for i in range(m):
            sol = arr2[i]/arr[i][i]     #x_i(0) = b_i/a_ii 
            X = np.append(X,sol)
    
        while True :
            SOL = np.array(X)
            for i in range(m):
                A_TIMES_X = 0 
                for j in range(m) :
                    if i != j :    
                        A_TIMES_X += arr[i][j] * X[j]
                    else : 
                        continue
            
                X[i] = (arr2[i] - A_TIMES_X)/arr[i][i]
        
            DIFF = abs(X - SOL)
        
            count += 1
            if abs(max(DIFF)) <= error:
                break
        return X,count
    else:
        sys.exit("MATRIX IS NOT DIAGONALLY DOMINANT")
  
if __name__ == "__main__":
    
    coeff_mat = [[45,2,3],[-3,22,2],[5,1,20]]
    result_mat = [58,47,67]
    
    
    print("\n GAUSS-SEIDAL METHOD TO SOLVE A SYSTEM OF LINEAR EQUATIONS : ")
    X1,count = GAUSS_SEIDAL(coeff_mat,result_mat)
    for i in range(len(X1)):
        print("x"+str(i+1)+" = ",np.round(X1[i],6))
    print("no.of iterations is : ",count)
    print(" ")
    
    print("THE SOLUTION FROM INBUILT PYTHON PACKAGE")
    X2 = np.linalg.solve(coeff_mat,result_mat)
    for i in range(len(X2)):
        print("x"+str(i+1)+" = ",X2[i])
    print("\n THE RELATIVE PERCENTAGE ERROR IN BOTH OF THE METHODS IS : ")

    for i in range(len(X1)):
        print("x_err"+str(i+1)+" = ",round(abs((X1[i] - X2[i])/X2[i])*100,5),"%")