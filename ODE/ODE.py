#ODE/ODE.py

#importar bibliotecas y definir la función a resolver

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Derivada como función de x y t
def func(x,t):
    """Function of the derivative $-x^{3}+\sin{t}$ 

    Examples:
        Consider the function $-x^{3}+\sin{t}$
        >>> func(2.0,3.0)
        -7.86

    Args:
        x (float): First argument
        t (float): Second argument

    Returns:
        float: Returns the derivative evaluate in 'x' and 't'

    """
    return -x**3 + np.sin(t)

#método de Euler
def euler(f, x0, times):
    """Euler method $x_i=x_{i-1}+h\cdot f(x_i)$

    Examples:
        Consider the function $-x^{3}+\sin{t}$
        >>> euler(func,0.0,np.linspace(0.0,10.0,20.0)
        [ 0.          0.          0.26439534  0.71189381  1.04830651  0.89488942
        0.77464602  0.52141013  0.17502349 -0.28921312 -0.80263945 -0.97897518
        -0.73458315 -0.50879969 -0.16038526  0.30726691  0.81787725  0.97386719
        0.72957635  0.49945695]

    Args:
        f (function): First argument, derivative function
        x0 (float) : Second argument, initial condition
        times (array): Third argument, array of steps

    Returns:
        array: Returns the array of approximations from x0 to the size of times minus one
    """
    h=times[1]-times[0]
    xn=np.zeros(times.size)
    xn[0]=x0
    for i in range(times.size-1):
        xn[i+1]=xn[i]+(h*f(xn[i],times[i]))
    return xn


#RK2
def rk2(f, x0, times):
    """ Runge-Kutta 2 method $x(t+h)=x(t)+hf(x+{hf(x,t)\over 2},t+{h\over 2})$

    Examples:
        Consider the function $-x^{3}+\sin{t}$
        >>> rk2(func,0.0,np.linspace(0.0,10.0,20.0)
        [ 0.          0.13691106  0.50040599  0.83221858  0.89696773  0.83638552
        0.684369    0.42791827  0.0377284  -0.46988042 -0.7896376  -0.78706438
        -0.65422559 -0.40234264 -0.0089812   0.49847514  0.79691697  0.78634182
        0.64914835  0.39298319]

    Args:
        f (function): First argument, derivative function
        x0 (float): Second argument, initial condition
        times (array): Third argument, array of steps

    Returns:
        array: Returns the array of approximations from x0 to the size of times minus one
    """
    h = times[1] - times[0]
    xn=np.zeros(times.size)
    xn[0]=x0

    for i in range(times.size-1):
        k1=h*f(xn[i],times[i])
        k2=h*f(xn[i]+k1/2,times[i]+h/2)
        xn[i+1]=xn[i]+k2
    return xn



#RK4
def rk4(f, x0, times):
    """Runge-Kutta 4 method 
    $k_1 = h\cdot f(x,t)$

    $k_2 = h\cdot f(x+ {k_1 \over 2},t+{h \over 2})$

    $k_3 = h\cdot f(x+ {k_2 \over 2},t+{h \over 2})$

    $k_4 = h\cdot f(x+k_3, t+h)$
    
    $x(t+h)=x(t)+{1 \over 6} (k_1+2k_2+2k_3+k_4)$

    Examples:
        Consider the function $-x^{3}+\sin{t}$
        >>> rk4(func,0.0,np,linspace(0.0,10.0,20.0)
        [ 0.          0.13505937  0.48487586  0.81979626  0.93102462  0.88144889
        0.72370945  0.46026479  0.06848691 -0.42773679 -0.78084333 -0.83096351
        -0.6993985  -0.43991777 -0.04506765  0.45089026  0.79057875  0.83047637
        0.69376603  0.43014634]

    Args:
        f (function): First argument, derivative function
        x0 (float): Second argument, initial condition
        times (array): Third argument, array of steps

    Returns:
        array: Returns the array of approximations from x0 to the size of times minus one
    """
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
