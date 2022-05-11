import numpy as np
import math as math
import matplotlib.pyplot as plt

def MSE(X,T,U):
    error = 0
    for x in range(0,len(X)):
        for t in range(0,len(T)):
            error += (U[t][x] - Solution(X[x],T[t]))**2
    return (error/(len(X)*len(T)))/(len(X)*len(T))

def plotU(x,U):
    figure, ax = plt.subplots(1,1,figsize=(8,6))
    len_ = len(U)
    len_arr_double = np.linspace(0,len_-1,6)
    len_arr_int = [int(elem) for elem in len_arr_double]
    for u_ind in len_arr_int:
        ax.plot(x,U[u_ind])
    legends = ["t="+str(index) for index in len_arr_int]
    ax.legend(legends)
    ax.set_title('Решение при различных t')

def plotErrors(errors, start, delta, type):
    x = range(start,start + delta)
    plt.plot(x,errors,c='black')
    plt.title('ошибки метода с различным шагом ' + type)
    plt.xlabel('число разбиений')
    plt.ylabel('функция ошибки решения')

# ошибка по времени
def ErrorT(X,U,TIME):
    error = []
    for t in range(0,len(TIME)):
        error_tmp = 0
        for i in range(0,len(X)):
            error_tmp += (U[t][i] - Solution(X[i],TIME[t]))**2
        error.append(np.sqrt(error_tmp))
    return error

# ошибка по пространству
def ErrorH(X,U,TIME):
    error = []
    for x in range(0,len(X)):
        error_tmp = 0
        for t in range(0,len(TIME)):
            error_tmp += (U[t][x] - Solution(X[x],TIME[t]))**2
        error.append(np.sqrt(error_tmp))
    return error

def F(x,t):
    return np.exp(-t)*np.sin(x)

def Ux0(t):
    return np.exp(-t)

def Uxl(t):
    return -np.exp(-t)

def U(x):
    return np.cos(x)

def ksi2(x):
    return -np.cos(x)

def Solution(x, t):
    return np.exp(-t) * np.cos(x)

def progonka(a, b, c, d, s):
    P = np.zeros(s)
    Q = np.zeros(s)

    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    k = s - 1

    for i in range(1, s):
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1])
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])
    P[k] = 0
    Q[k] = (d[k] - a[k] * Q[k - 1]) / (b[k] + a[k] * P[k - 1])

    x = np.zeros(s)
    x[k] = Q[k]

    for i in range(s - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
    return x

def autofill(x0, space_step, m, n, param_a, time_step, aprox_f):
    Uarray = np.zeros([n, m])

    tmp_x = x0
    for j in range(m):
        Uarray[0][j] = U(tmp_x)
		
        if aprox_f == 1:
            Uarray[1][j] = U(tmp_x) + ksi2(tmp_x)*time_step
			
        if aprox_f == 2:
            Uarray[1][j] = U(tmp_x) + ksi2(tmp_x)*time_step + param_a**2*(-np.cos(tmp_x))*time_step**2/2
				
        tmp_x += space_step
    
    return Uarray

"""## Схема крест"""

def explicit(t, m, n, x0, xl, aprox_f):

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)


    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2
	
    for k in range(1, n - 1):
        for j in range(1, m - 1):

            Uarray[k + 1][j] = \
                Uarray[k][j + 1] *(sigma + time_step**2*param_b/(2*space_step))/(1.5*time_step+1) +\
                Uarray[k][j] * (2 - 2*sigma + param_c*time_step**2) /(1.5*time_step+1) + \
                Uarray[k][j - 1] *(sigma - time_step**2*param_b/(2*space_step))/(1.5*time_step+1) + \
                (Uarray[k - 1][j] * (-1)*(1 - 1.5*time_step))/(1.5*time_step+1) + \
                F(space_step*j,time_step*k)*time_step**2/(1.5*time_step+1)

        # граничные условия не содержат производных,
        # но аппроксимации, перечисленные в задании реализованы
        # в пятой лабораторной работе
        Uarray[k + 1][0] = Ux0(time_step*k)
        Uarray[k + 1][m - 1] = Uxl(time_step*k)
    return Uarray

x0 = 0
xl = math.pi
param_a = 1
param_b = 1
param_c = -1
aprox_f = 2
m = 50
n = 50
tend = 1
space_step = (xl - x0) / (m - 1)
time_step = tend / (n - 1)

x = np.arange(x0,xl+space_step-0.0001,space_step)
t = np.arange(0,tend+time_step-0.0001,time_step)

U_exp = explicit(tend, m, n, x0, xl, aprox_f)

plotU(x,U_exp)

MSE(x,t,U_exp)

errors_time = []
delta = 100
for n_ in range(n,n + delta):
    time_step_local = tend/(n_-1)
    t_local = np.arange(0,1+time_step_local-0.0001,time_step_local)
    Uarray_with_different_time = (explicit(tend, m, n_, x0, xl, aprox_f))
    errors_time.append(MSE(x,t_local,Uarray_with_different_time))

plotErrors(errors_time,n,delta,'t')

errors_h = []
for m_ in range(m,m+delta):
    space_step_local = (xl-x0)/(m_-1)
    x_local = np.arange(x0,xl+space_step_local-0.0001,space_step_local)
    Uarray_with_different_h = explicit(tend, m_, n, x0, xl, aprox_f)
    errors_h.append(MSE(x_local,t,Uarray_with_different_h))

plotErrors(errors_h,m,delta,'h')

"""## Неявная схема"""

def implicit(t, m, n, x0, xl, aprox_f):

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)


    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2

    alpha = 0
    betta = 1
    gamma = 0
    delta = 1
    tmp = param_b*time_step**2/(2*space_step)

    for k in range(1, n - 1):
        a = np.zeros(m)
        b = np.zeros(m)
        c = np.zeros(m)
        d = np.zeros(m)

        for j in range(1, m - 1):
            a[j] = -sigma + tmp
            b[j] = (1.5*time_step + 1 + 2*sigma - param_c*time_step**2)
            c[j] = -sigma - tmp
            d[j] = Uarray[k - 1][j]*(1.5*time_step - 1) + 2*Uarray[k][j]

        # граничные условия не содержат производных,
        # но аппроксимации, перечисленные в задании реализованы
        # в пятой лабораторной работе
        b[0] = betta - alpha / space_step
        c[0] = alpha / space_step
        d[0] = Ux0((k + 1) * time_step)

        a[m - 1] = - gamma / space_step
        b[m - 1] = delta + gamma / space_step
        d[m - 1] = Uxl((k + 1) * time_step)


        Y = progonka(a, b, c, d, m)
        Uarray[k + 1] = Y
    return Uarray

x0 = 0
xl = math.pi
param_a = 1
param_b = 1
param_c = -1
aprox_f = 2
m = 50
n = 50
tend = 1
space_step = (xl - x0) / (m - 1)
time_step = tend / (n - 1)

x = np.arange(x0,xl+space_step-0.0001,space_step)
t = np.arange(0,tend+time_step-0.0001,time_step)

U_im = implicit(tend, m, n, x0, xl, aprox_f)

plotU(x,U_im)

MSE(x,t,U_im)

errors_time = []
delta = 100
for n_ in range(n,n + delta):
    time_step_local = tend/(n_-1)
    t_local = np.arange(0,1+time_step_local-0.0001,time_step_local)
    Uarray_with_different_time = (implicit(tend, m, n_, x0, xl, aprox_f))
    errors_time.append(MSE(x,t_local,Uarray_with_different_time))

plotErrors(errors_time,n,delta,'t')

errors_h = []
for m_ in range(m,m+delta):
    space_step_local = (xl-x0)/(m_-1)
    x_local = np.arange(x0,xl+space_step_local-0.0001,space_step_local)
    Uarray_with_different_h = implicit(tend, m_, n, x0, xl, aprox_f)
    errors_h.append(MSE(x_local,t,Uarray_with_different_h))

plotErrors(errors_h,m,delta,'h')

"""Рассмотрим решения двух методов в одно время"""

plt.figure(figsize=(10,6))
plt.plot(x,U_exp[25])
plt.plot(x,U_im[25])
plt.legend(['Explicit','Implicit'])
plt.title('Решения методов в момент времени t = ' + str(round(25*time_step,4)))