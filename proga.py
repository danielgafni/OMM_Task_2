from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

x1 = 0.
x2 = -1. #h < 0
t1 = 0.
t2 = 0.2

Nx = 100
Mt = 100
h = (x2-x1) / Nx
tau = (t2-t1) / Mt

epsilon = 1e-6
DRAWCHAR = False #show charachteristics plot

def a(u):
	if u > 10:
		return 2.0
	if u < -10:
		return 0
	return 2.0 * exp(2.0*u) / (1.0 + exp(2.0*u))

def F(U, m, n):
	if U[m][n] > 10:
		return 2*U[m][n]
	if U[m][n] < -10:
		return 0
	return log(1.0 + exp(2.0*U[m][n]))

def cha1(x0, tm):
	return [x0 + a(sin(pi*x0)) * t for t in tm]

def cha2(t0, tm):
	return [a(0)*(t-t0) for t in tm]

tm = np.arange(t1,t2,tau)
xn = np.arange(x1,x2,h)

if DRAWCHAR:
	for x0 in np.arange(-2,0,0.05):
		plt.plot(cha1(x0, tm), tm)
	for t0 in np.arange(-1,1,0.1):
		plt.plot(cha2(t0, tm), tm)
	plt.xlabel('x')
	plt.ylabel('t')
	plt.title('Characteristics')
	plt.show()

U = [[0 for x in xn] for t in tm]
# starting conditions
for n in range(Nx):
	U[0][n] = sin(pi*xn[n])
# boundary conditions
for m in range(Mt):
	U[m][0] = 0.
# first approximation
for n in range(Nx)[1:Nx]:
	for m in range(Mt)[1:Mt]:
		U[m][n] = U[0][n]

def f(U, mp1, np1):
	n = np1-1
	m = mp1-1
	return (U[mp1][n]-U[m][n] + U[mp1][np1]-U[m][np1]) / (2.*tau) - (F(U, mp1, np1)-F(U, mp1,n) + F(U, m, np1)-F(U,m,n)) / (2.*h) # 4-point
	#return (U[mp1][np1]-U[m][np1]) / (1.*tau) - (F(U, mp1, np1)-F(U, mp1,n)) / (1.*h) # 3-point implicit

def df(U, mp1, np1):
	return 1. / (2.*tau) - a(U[mp1][np1]) / (2.*h) # 4-point
	#return 1. / (1.*tau) - a(U[mp1][np1]) / (1.*h) # 3-point implicit

eps = epsilon + 1; # > epsilon
while eps > epsilon:
	eps = 0
	for m in range(Mt)[0:Mt-1]:
		for n in range(Nx)[0:Nx-1]:
			#print m, n, U[m+1][n+1]
			ep = f(U, m+1, n+1) / df(U, m+1, n+1)
			U[m+1][n+1] = U[m+1][n+1] - ep
			if abs(ep) > eps:
				eps = abs(ep)
	print eps

X, T = np.meshgrid(xn, tm)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, T, U, cmap='YlOrBr')
plt.title('Solution')
plt.xlabel('x')
plt.ylabel('t')
plt.show()
