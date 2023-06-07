import math
import numpy as np
import matplotlib.pyplot as plt

b=0.001
N=500
y=0.1
t_max=100
dt=0.1
u0=1
TOL=10**(-6)
# mi <= 20
a=b*N-y
#############################
def f(u):
    return a * u - b * u**2

u=[u0]
u1=[u0]

T=[i for i in np.arange(0,t_max+dt,dt)]

for t in T[1:]:
    u_n1=u[-1]
    u_n=0.
    for _ in range(20):     
        if  math.fabs(u_n1-u_n)<TOL:
            break
        u_n=u_n1
        u_n1=u[-1]+(dt/2.)*(f(u[-1])+f(u_n))
    u.append(u_n1)
    u_n1=u1[-1]
    u_n=0.
    for _ in range(20):     
        if  math.fabs(u_n1-u_n)<TOL:
            break
        u_n=u_n1
        u_n1=u_n-(u_n-u1[-1]-(dt/2.)*(f(u1[-1])+f(u_n)))/(1-(dt/2.)*(a-2*b*u1[-1]))
    u1.append(u_n1)  

plt.clf() 
plt.plot(T,u,label='u(t)')
plt.plot(T,[N-i for i in u],label='z(t)=N-u(t)')
plt.legend()
plt.ylabel('u')
plt.xlabel('t')
plt.savefig('picard.png')

plt.clf() 
plt.plot(T,u1,label='u(t)')
plt.plot(T,[N-i for i in u1],label='z(t)=N-u(t)')
plt.legend()
plt.ylabel('u')
plt.xlabel('t')
plt.savefig('newton.png')


def F(U1,U2,u,dt,a1,a2):
    return U1-u-dt*(a1*(a*U1-b*U1**2)+a2*(a*U2-b*U2**2))
def m1(a_n,U,dt):
    return 1-dt*a_n*(a-2*b*U)
def m2(a_n,U,dt):
    return -dt*a_n*(a-2*b*U)

u=[u0]
U1=U2=u0

a11=0.25
a12=0.25-(3**(.5)/6)
a21=0.25+(3**(.5)/6)
a22=0.25
b1=b2=.5
for t in T[1:]:
    u_n1=u[-1]
    u1_n = 0 
    u2_n = 0
    for _ in range(20):     
        if  math.fabs(u1_n-U1)<=TOL or math.fabs(u2_n-U2)<=TOL:
            break
        u1_n = U1 
        u2_n = U2
        F1=F(U1,U2,u[-1],dt,a11,a12)
        F2=F(U1,U2,u[-1],dt,a21,a22)
        dU1=(F2*m2(a12,U2,dt)-F1*m1(a22,U2,dt))/(m1(a11,U1,dt)*m1(a22,U2,dt)-m2(a12,U2,dt)*m2(a21,U2,dt))
        dU2=(F1*m2(a21,U1,dt)-F2*m2(a11,U1,dt))/(m1(a11,U1,dt)*m1(a22,U2, dt)-m2(a12,U2,dt)*m2(a21,U2,dt))
        U1+=dU1
        U2+=dU2
    u_n1 = u[-1] + dt*(b1*f(U1) + b2*f(U2))
    u.append(u_n1)  
plt.clf() 
plt.plot(T,u,label='u(t)')
plt.plot(T,[N-i for i in u],label='z(t)=N-u(t)')
plt.legend()
plt.ylabel('u')
plt.xlabel('t')
plt.savefig('RK2.png')