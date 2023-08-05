print("FIBONACCI SERIES")

num=int(input("ENTER THE TOTAL COUNT OF SERIES:"))

start=0
second=1
num1=num-2
l1=[0,1]

print("HERE IS THE FIBONACCI SERIES FOR TOTAL COUNT",num)

#print(start)
#print(second)
while(num1 >0):
    fib = start + second
    #print(fib)
    l1.append(fib)

    num1-=1

    start = second
    second = fib
    
for i in range(num-1):
    print(l1[i],end=",")
    
print(l1[-1],"........")
