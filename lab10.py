import math
import numpy as np
#import numba as nb
import matplotlib.pyplot as plt

def delta_K(x,y):
    if abs(x - y) < 10**(-6): 
        return 1
    else:
	    return 0


nx = 150
nt = 1000
d  = 0.1
dt = 0.05
xA = 7.5
s  = 0.5
xF = 2.5

t  = [dt * n for n in range(nt+1)]
x  = [ d * i for i in range(nx+1)]
u  = [0 for i in range(nx+1)]
v  = [0 for i in range(nx+1)]
vp = [0 for i in range(nx+1)]

u[ 0] = 0
u[-1] = 0
v[ 0] = 0
v[-1] = 0

for i in range(1,nx):
    u[i] = math.exp(-(((d * i) - xA)**2 / (2 * (s)**2)))

u0=u

par = [(0,0),(0,0.1),(0,1)]

for alpha,beta in par:
	a  = [ 0 for i in range(nx+1)]
	for i in range(1,nx):
		a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / d**2 - (beta * (u[i] - u0[i]) / dt) + (alpha * math.cos((50 * 0.0) / (dt * nt)) * delta_K((d * i), xF))
	for n in range(1, nt+1):
		for i in range(1, nx):
			vp[i] = v[i] + (dt * a[i] / 2)
		u0 = u
		for i in range(1, nx):
			u[i] = u[i] + (dt * vp[i])
		for i in range(1, nx):
			a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / d**2 - (beta * (u[i] - u0[i]) / dt) + (alpha * math.cos(50.0 * n / (dt * nt)) * delta_K((d * i), xF))
		
		for i in range(1, nx):
			v[i] = vp[i] + (dt * a[i] / 2)
		
alpha = 1
beta = 1

E = (d / 4.) * ((u[1] - u[0]) / d)**2 + ((u[nx] - u[nx - 1]) / d)**2 + (d / 2)*sum(v[i]**2 + ((u[i + 1] - u[i - 1]) / (2 * d))**2 for i in range(1, nx))