import math
import matplotlib.pyplot as plt
import numpy as np

dt = 10**(-4)
R = 100
L = 0.1
C = 0.001
w0 = 1/((L*C)**(.5))
T0 = (2*math.pi)/w0 
y=[1]
lambd = -1
Dt=[0.01, 0.1, 1]
d=[0]


def f(x):
    return math.exp(lambd*x)

for dt in Dt:
    for t in np.arange(dt,5,dt):
        y.append(y[-1]+dt*lambd*y[-1])
        d.append(y[-1]-f(t))
    if dt==0.01:
        y1 = y.copy()
        d1 = d.copy()
    if dt==0.1:
        y2 = y.copy()
        d2 = d.copy()
    if dt==1:
        y3 = y.copy()
        d3 = d.copy()
    y=[1]
    d=[0]
T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,y1, 'ro', label='y(t), dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,y2, 'bo', label='y(t), dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,y3, 'go', label='y(t), dt=1')
plt.legend()
plt.ylabel('y(t)')
plt.xlabel('t')
plt.savefig('euler.png')

T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,d1,label='dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,d2,label='dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,d3,label='dt=1')
plt.legend()
plt.ylabel('y_num(t)-y_dok(t)')
plt.xlabel('t')
plt.savefig('eblad.png')

for dt in Dt:
    for t in np.arange(dt,5,dt):
        k1=lambd * y[-1] 
        k2=lambd * (y[-1] + dt*k1)
        y.append(y[-1]+(dt/2)*(k1+k2))
        d.append(y[-1]-f(t))
    if dt==0.01:
        y1 = y.copy()
        d1 = d.copy()
    if dt==0.1:
        y2 = y.copy()
        d2 = d.copy()
    if dt==1:
        y3 = y.copy()
        d3 = d.copy()
    y=[1]
    d=[0]
T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,y1, 'ro', label='y(t), dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,y2, 'bo', label='y(t), dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,y3, 'go', label='y(t), dt=1')
plt.legend()
plt.ylabel('y(t)')
plt.xlabel('t')
plt.savefig('rk2.png')

T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,d1,label='dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,d2,label='dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,d3,label='dt=1')
plt.legend()
plt.ylabel('y_num(t)-y_dok(t)')
plt.xlabel('t')
plt.savefig('rk2blad.png')

for dt in Dt:
    for t in np.arange(dt,5,dt):
        k1=lambd * y[-1]
        k2=lambd * (y[-1] + (dt/2)*k1)
        k3=lambd * (y[-1] + (dt/2)*k2)
        k4=lambd * (y[-1] + dt*k3)
        y.append(y[-1]+(dt/6)*(k1+2*k2+2*k3+k4))
        d.append(y[-1]-f(t))
    if dt==0.01:
        y1 = y.copy()
        d1 = d.copy()
    if dt==0.1:
        y2 = y.copy()
        d2 = d.copy()
    if dt==1:
        y3 = y.copy()
        d3 = d.copy()
    y=[1]
    d=[0]
T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,y1, 'ro', label='y(t), dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,y2, 'bo', label='y(t), dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,y3, 'go', label='y(t), dt=1')
plt.legend()
plt.ylabel('y(t)')
plt.xlabel('t')
plt.savefig('rk4.png')

T=[i for i in np.arange(0,5,Dt[0])]
plt.clf() 
plt.plot(T,d1,label='dt=0.01')
T=[i for i in np.arange(0,5,Dt[1])]
plt.plot(T,d2,label='dt=0.1')
T=[i for i in np.arange(0,5,Dt[2])]
plt.plot(T,d3,label='dt=1')
plt.legend()
plt.ylabel('y_num(t)-y_dok(t)')
plt.xlabel('t')
plt.savefig('rk4blad.png')

dt=10**(-4)

Q=[0]
I=[0]

def V(wv,t):
    return 10*math.sin(wv*t)

W=[0.5*w0, 0.8*w0, 1.0*w0, 1.2*w0]

for wV in W:
    for t in np.arange(dt,4*T0,dt):
        k1Q=I[-1]
        k1I=V(wV, t)/L-1/(L*C)*Q[-1]-(R/L)*I[-1]

        k2Q=I[-1]+(dt/2)*k1I
        k2I=V(wV, t+dt/2)/L-1/(L*C)*(Q[-1]+dt/2*k1Q)-(R/L)*(I[-1]+dt/2*k1I)

        k3Q=I[-1]+dt/2*k2I
        k3I=V(wV, t+dt/2)/L-1/(L*C)*(Q[-1]+dt/2*k2Q)-(R/L)*(I[-1]+dt/2*k2I)

        k4Q=I[-1]+dt*k3I
        k4I=V(wV, t+dt)/L-1/(L*C)*(Q[-1]+dt*k3Q)-(R/L)*(I[-1]+dt*k3I)

        Q.append(Q[-1]+dt/6*(k1Q+2*k2Q+2*k3Q+k4Q))
        I.append(I[-1]+dt/6*(k1I+2*k2I+2*k3I+k4I))
    if wV==0.5*w0:
        Q1 = Q.copy()
        I1 = I.copy()
    if wV==0.8*w0:
        Q2 = Q.copy()
        I2 = I.copy()
    if wV==1.0*w0:
        Q3 = Q.copy()
        I3 = I.copy()
    if wV==1.2*w0:
        Q4 = Q.copy()
        I4 = I.copy()
    Q=[0]
    I=[0]


T=[i for i in np.arange(0,4*T0,dt)]
plt.clf() 
plt.plot(T,Q1, label=r'$Q(t), \omega_v=0.5*\omega_0$')
plt.plot(T,Q2, label=r'$Q(t), \omega_v=0.8*\omega_0$')
plt.plot(T,Q3, label=r'$Q(t), \omega_v=1.0*\omega_0$')
plt.plot(T,Q4, label=r'$Q(t), \omega_v=1.2*\omega_0$')
plt.legend()
plt.ylabel('Q(t)')
plt.xlabel('t')
plt.savefig('qt.png')

T=[i for i in np.arange(0,4*T0,dt)]
plt.clf() 
plt.plot(T,I1,label=r'$I(t), \omega_v=0.5*\omega_0$')
plt.plot(T,I2,label=r'$I(t), \omega_v=0.8*\omega_0$')
plt.plot(T,I3,label=r'$I(t), \omega_v=1.0*\omega_0$')
plt.plot(T,I4,label=r'$I(t), \omega_v=1.2*\omega_0$')
plt.legend()
plt.ylabel('I(t)')
plt.xlabel('t')
plt.savefig('it.png')
