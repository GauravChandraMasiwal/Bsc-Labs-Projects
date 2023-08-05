# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 20:54:37 2021

@author: chand
"""
#this is a program for gauss elimination method to solve the set of matrix equations
import numpy as np

def row_swap(a,row):
    n=len(Aug_mat)
    a = np.array(a)
    count=1
    while (a[row][row]==0):
        if (row==n-1):
            a[[row, 0]]=a[[0, row]]
        else:
            a[[row+count, row]]=a[[row, row+count]]
        count=count+1
        if count==n:
            return a,0
    return a,count

def rank(a):
    n=len(Aug_mat)
    m=int(np.size(Aug_mat)/n)
    a = np.array(a)
    zero_row=len(np.where(~a.any(axis=1))[0])
    coeff_mat=np.zeros((n,m-1))
    for i in range(n):
        for j in range(m-1):
            coeff_mat[i][j]=Aug_mat[i][j]
    zero_rows2=len(np.where(~coeff_mat.any(axis=1))[0])
    rank1=n-zero_row
    rank2=n-zero_rows2
    return rank1,rank2

def remove_column(a,r_col):
    n=len(Aug_mat)
    m=int(np.size(Aug_mat)/n)
    new_a=np.zeros((n,m-1))
    for i in range(n):
        for j in range(1,m,1):
            new_a[i][j-1]=a[i][j]
    return new_a

def reduced_form(a):
    n=len(a) 
    m=int(np.size(a)/len(a))
    for i in range(n):
        if a[i][i]==0:
            print("row swap")
            a,c=row_swap(a,i)
            print(a)
            if c==0:
                print("removing redundant column")
                a=remove_column(a,i)
                m=m-1       
        for j in range(i+1,n):
            ratio=a[j][i]/a[i][i]  
            for k in range(m):  
                a[j][k]=round((a[j][k]-(ratio*a[i][k])),3)
            print("intermediate steps")
            print(np.array(a))
            rank1,rank2=rank(a)
            if rank1<n:
                return a,rank1,rank2
    return a,rank1,rank2

def back_substitution(a):
    n=len(Aug_mat)
    x=np.zeros([n])
    x[n-1] =round(a[n-1][n]/a[n-1][n-1],4) 
    for i in range(n-2,-1,-1):     
        x[i] = a[i][n]             
        for j in range(i+1,n):     
            x[i] = x[i] - a[i][j]*x[j]  
        x[i] = round(x[i]/a[i][i],3)   
    return x
    


if __name__ == "__main__":
    Aug_mat=[[4,-1,0,-1],[-1,5,-1,2],[0,1,-3,1]]
    n=len(Aug_mat)
    m=int(np.size(Aug_mat)/n)
    coeff_mat=np.zeros((n,m-1))

    print("coefficient matrix:")
    for i in range(n):
        for j in range(m-1):
            coeff_mat[i][j]=Aug_mat[i][j]
    print(np.array(coeff_mat))

    print("AUGMENTED MATRIX")
    print(np.array(Aug_mat))

    ans,r1,r2=reduced_form(Aug_mat)

    print("REDUCED ECHELON FORM")
    print(np.array(ans))
    print("RANK OF AUGMENTED MATRIX: ",r1)
    print("RANK OF COEFFICIENT MATRIX: ",r2)
    if m-1 == r1 : 
        print("UNIQUE SOLN EXIST")
        x = back_substitution(ans)
        for i in range(len(x)):
            print("x"+str(i+1)+": ",x[i])
    elif m-1 < r1 :
        print("over-determined case")
    else :
        print("under-determined case")
    
    
    Aug_mat=[[2,7,4,2],[8,5,7,2],[7,4,0,1],[1,6,4,7]]   #overdetermined

    n=len(Aug_mat)
    m=int(np.size(Aug_mat)/n)
    coeff_mat=np.zeros((n,m-1))

    print("coefficient matrix:")
    for i in range(n):
        for j in range(m-1):
            coeff_mat[i][j]=Aug_mat[i][j]
    print(np.array(coeff_mat))

    print("AUGMENTED MATRIX")
    print(np.array(Aug_mat))

    ans,r1,r2=reduced_form(Aug_mat)

    print("REDUCED ECHELON FORM")
    print(np.array(ans))
    print("RANK OF AUGMENTED MATRIX: ",r1)
    print("RANK OF COEFFICIENT MATRIX: ",r2)
    if m-1 == r1 : 
        print("UNIQUE SOLN EXIST")
        x = back_substitution(ans)
        for i in range(len(x)):
            print("x"+str(i+1)+": ",x[i])
    elif m-1 < r1 :
        print("over-determined case")
    else :
        print("under-determined case")
    
    Aug_mat=[[7,6,4,4,5],[2,3,4,6,5],[1,2,3,7,5]]        #undetermined
    
    n=len(Aug_mat)
    m=int(np.size(Aug_mat)/n)
    coeff_mat=np.zeros((n,m-1))

    print("coefficient matrix:")
    for i in range(n):
        for j in range(m-1):
            coeff_mat[i][j]=Aug_mat[i][j]
    print(np.array(coeff_mat))

    print("AUGMENTED MATRIX")
    print(np.array(Aug_mat))

    ans,r1,r2=reduced_form(Aug_mat)

    print("REDUCED ECHELON FORM")
    print(np.array(ans))
    print("RANK OF AUGMENTED MATRIX: ",r1)
    print("RANK OF COEFFICIENT MATRIX: ",r2)
    if m-1 == r1 : 
        print("UNIQUE SOLN EXIST")
        x = back_substitution(ans)
        for i in range(len(x)):
            print("x"+str(i+1)+": ",x[i])
    elif m-1 < r1 :
        print("over-determined case")
    else :
        print("under-determined case")
    
