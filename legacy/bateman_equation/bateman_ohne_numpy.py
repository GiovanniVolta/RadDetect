#!/usr/bin/python3

from math import*
import matplotlib.pyplot as plt
import numpy as np

def bateman_c(n,la):
    enum=1
    c=[]
     
    for j in range(n):
        enum = enum*la[j]
    
    for m in range(n):
        denom=1
        for i in range(n):
            if i != m:
                denom = denom*(la[i]-la[m])
            else:
                pass
        c.append(enum/denom)
    return c

def bateman(t,n,N0,c,la):
    bateman_sum=0
    for m in range(n):
        bateman_sum = bateman_sum + c[m]*exp(-la[m]*t)
    activity = N0 * bateman_sum    
    return activity



tH_222Rn=3.8235*24*3600
tH_218Po=3.1*60
tH_214Pb=26.8*60
tH_214Bi=19.9*60
tH_214Po=164.3E-6
tH_210Pb=22.3*365*24*3600
tH_210Bi=5.013*24*3600
tH_210Po=138.376*24*3600

rn222_chain=[log(2)/tH_222Rn, log(2)/tH_218Po, log(2)/tH_214Pb, log(2)/tH_214Bi, log(2)/tH_214Po, log(2)/tH_210Pb, log(2)/tH_210Bi, log(2)/tH_210Po]

timebin=int(330350*10)
daughters=8
A=[]
time=[]


for i in range(daughters):
    dummy=[]
    coefficients=bateman_c(i+1, rn222_chain)
    for j in range(timebin):
        #print(j,i+1)
        dummy.append(bateman(j ,i+1 , 1*tH_222Rn/log(2), coefficients, rn222_chain))
    A.append(dummy)

#        A.append(daughters*[bateman(j+1 ,i+1 , 1*tH_222Rn/log(2), rn222_chain)])
#        print(i, j+1, A[j][i])

for i in range(timebin):
    time.append(i)

        
fig = plt.figure(figsize=(8,8))
for i in range(daughters):
    plt.plot(time,A[i])
#plt.legend()
plt.show()
