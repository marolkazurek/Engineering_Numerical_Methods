import math
import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import scipy
import scipy.linalg

nx=40
ny=40
N=(nx+1)*(ny+1)
de=1.0
dt=1.0
TA=40
TB=0
TC=30
TD=0
kB=0.1
kD=0.6
IT_MAX=2000
A=[[0 for i in range(N)] for j in range(N)]
B=[[0 for i in range(N)] for j in range(N)]
c=[0 for i in range(N)]
d=[0 for i in range(N)]
T=[0 for i in range(N)]
#######################################################################
for i in range(1, nx):
    for j in range(1, ny):
        l=i+j*(nx+1)
        A[l][l-nx-1]=dt/(2*de**2)
        A[l][l-1]=dt/(2*de**2)
        A[l][l+1]=dt/(2*de**2)
        A[l][l+nx+1]=dt/(2*de**2)

        A[l][l]=(-2*dt)/(de**2)-1

        B[l][l-nx-1]=-dt/(2*de**2)
        B[l][l-1]=-dt/(2*de**2)
        B[l][l+1]=-dt/(2*de**2)
        B[l][l+nx+1]=-dt/(2*de**2)

        B[l][l]=(2*dt)/(de**2)-1

#######################################################################
for i in [0,nx]:
    for j in range(0, ny+1):
        l=i+j*(nx+1)
        A[l][l]=1
        B[l][l]=1
        c[l]=0

#######################################################################
for i in range(1, nx):
    l=i+ny*(nx+1)
    j=ny
    A[l][l-nx-1]=-1/(kB*de)
    A[l][l]=1+1/(kB*de)
    c[l]=TB
    for k in range(len(B)):
        B[l][k]=0
#######################################################################
for i in range(1, nx):
    l = i + 0*(nx+1)
    j=0
    A[l][l+nx+1]=-1/(kD*de)
    A[l][l]=1+1/(kD*de)
    c[l]=TD
    for k in range(len(B)):
        B[l][k]=0

#########################################################################
for j in range(ny+1):
    l=0+j*(nx+1)
    T[l]=TA

for j in range(ny+1):
    l=nx+j*(nx+1)
    T[l]=TC

for i in range(1, nx):
    for j in range(0, ny):
        l=i+j*(nx+1)
        T[l]=0

# for i in range(len(A)):
#     print(A[i])

lu,piv=scipy.linalg.lu_factor(A)


for it in range(IT_MAX):
    d=np.add(np.matmul(B,T),c)
    result=scipy.linalg.lu_solve((lu, piv), d)
    T=result
    if it==100 or it==200 or it==500 or it==1000 or it==2000:
        TT= np.reshape(T,(nx+1,ny+1))
        helper=np.zeros((nx+1,ny+1))
        for i in range(1,nx):
            for j in range(1,ny):
                helper[i][j] = (TT[i+1,j]+TT[i-1][j]-4*TT[i][j]+TT[i][j+1]+TT[i][j-1])/de**2
        plt.figure()
        plt.imshow(TT)
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(r'it%d' %it)
        
        plt.figure()
        plt.imshow(helper)
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(r'itt%d' %it)




