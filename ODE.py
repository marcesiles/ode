#importar bibliotecas y definir la función a resolver

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def func(x,t):
    return -x**3 + np.sin(t)

#método de Euler
def euler(f, x0, times):
    h=times[1]-times[0]
    xn=np.zeros(times.size)
    xn[0]=x0
    for i in range(times.size-1):
        xn[i+1]=xn[i]+(h*f(xn[i],times[i]))
    return xn

#definición de las variables
x_0 = 0
x0=np.copy(x_0)
t_0 = 0

N1=20
N2=100

#Ejecución
y1_eu=euler(func,x_0,np.linspace(t_0, 10.0, N1))
y2_eu=euler(func,x_0,np.linspace(t_0, 10.0, N2))


#gráfico
plt.plot(np.linspace(t_0, 10.0, N1), y1_eu, 'r')
plt.plot(np.linspace(t_0, 10.0, N2), y2_eu, 'b')

plt.show()


#RK2
def rk2(f, x0, times):
    h = times[1] - times[0]
    xn=np.zeros(times.size)
    xn[0]=x0

    for i in range(times.size-1):
        k1=h*f(xn[i],times[i])
        k2=h*f(xn[i]+k1/2,times[i]+h/2)
        xn[i+1]=xn[i]+k2
    return xn

#Ejecución
y1_rk2=rk2(func,0,np.linspace(0, 10.0, N1))
y2_rk2=rk2(func,0,np.linspace(0, 10.0, N2))

#Gráfica
plt.plot(np.linspace(t_0, 10.0, N1), y1_rk2, 'r')
plt.plot(np.linspace(t_0, 10.0, N2), y2_rk2, 'b')

plt.show()

#RK4
def rk4(f, x0, times):
    h = times[1] - times[0]
    xn=np.zeros(times.size)
    xn[0]=x0

    for i in range(times.size-1):
        k1=h*f(xn[i],times[i])
        k2=h*f(xn[i]+k1/2,times[i]+h/2)
        k3=h*f(xn[i]+k2/2,times[i]+h/2)
        k4=h*f(xn[i]+k3,times[i]+h)
        xn[i+1]=xn[i]+(k1+2*k2+2*k3+k4)/6
    return xn

#Ejecución
y1_rk4=rk4(func,0,np.linspace(0, 10.0, N1))
y2_rk4=rk4(func,0,np.linspace(0, 10.0, N2))

#Gráfico
plt.plot(np.linspace(t_0, 10.0, N1), y1_rk4, 'r')
plt.plot(np.linspace(t_0, 10.0, N2), y2_rk4, 'b')

plt.show()

