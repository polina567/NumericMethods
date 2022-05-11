import numpy as np
import math as math
import matplotlib.pyplot as plt

def MSE(X,Y,U):
    error = 0
    for x in range(0,len(X)):
        for y in range(0,len(Y)):
            error += (U[y][x] - Solution(X[x],Y[y]))**2
    return (error/(len(X)*len(Y)))/(len(X)*len(Y))

def plotU(x,U):
    figure, ax = plt.subplots(1,1,figsize=(8,6))
    len_ = len(U)
    len_arr_double = np.linspace(0,len_-1,6)
    len_arr_int = [int(elem) for elem in len_arr_double]
    for u_ind in len_arr_int:
        ax.plot(x,U[u_ind])
    legends = ["y="+str(index) for index in len_arr_int]
    ax.legend(legends)
    ax.set_title('Решение при различных y')

def plotErrors(errors, start, delta, type):
    x = range(start,start + delta)
    plt.plot(x,errors,c='black')
    plt.title('ошибки метода с различным шагом ' + type)
    plt.xlabel('число разбиений')
    plt.ylabel('функция ошибки решения')

# ошибка по Y
def ErrorY(X,U,Y):
    error = []
    for j in range(0,len(Y)):
        error_tmp = 0
        for i in range(0,len(X)):
            error_tmp += (U[j][i] - Solution(X[i],Y[j]))**2
        error.append(np.sqrt(error_tmp))
    return error

# ошибка по X
def ErrorX(X,U,Y):
    error = []
    for x in range(0,len(X)):
        error_tmp = 0
        for j in range(0,len(Y)):
            error_tmp += (U[j][x] - Solution(X[x],Y[j]))**2
        error.append(np.sqrt(error_tmp))
    return error

def Ux0(x):
    return np.cos(x)

def Uxl(x):
    return np.zeros_like(x)

def U0y(y):
    return np.exp(-y)*np.cos(y)

def Uly(y):
    return np.zeros_like(y)

def Solution(x, y):
    return np.exp(-y)*np.cos(x)*np.cos(y)

def autofill(m, n, lx):
    Uarray = np.zeros([n, m])
    x = np.linspace(0, lx, m)
    
    for j in range(n):
        Uarray[j,:] = Ux0(x) * (j/n) + Uxl(x) * (1 - (j/n))
    
    return Uarray

"""## Метод простых итераций (метод Либмана)"""

def iters(m, n, lx, ly, eps):
    x_step = lx / (m - 1)
    y_step = ly / (n - 1)
    x_step2 = x_step**2
    y_step2 = y_step**2
    
    A = (x_step2 - x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    B = (x_step2 + x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    C = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    D = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    
    Uarray = autofill(m, n, lx)
    
    iters = 5000
    
    while iters > 0:
        iters -= 1
        maxdiff = 0
        newU = np.copy(Uarray)
        
        for i in range(1, m - 1):
            for j in range(1, n - 1):
                newU[j][i] = Uarray[j - 1][i] * A + Uarray[j + 1][i] * B + Uarray[j][i - 1] * C + Uarray[j][i + 1] * D
                maxdiff = max(maxdiff, abs(newU[j][i] - Uarray[j][i]))
        
        Uarray = newU
        if maxdiff < eps:
            break
    
    return Uarray

lx = math.pi / 2
ly = math.pi / 2
m = 25
n = 25
x_step = lx / (m - 1)
y_step = ly / (n - 1)
eps = 0.0001

x = np.arange(0,lx+x_step-0.0001,x_step)
y = np.arange(0,ly+y_step-0.0001,y_step)

U_iter = iters(m, n, lx, ly, eps)

plotU(x,U_iter)

MSE(x,y,U_iter)

errors_y = []
delta = 15
for n_ in range(n,n + delta):
    y_step_local = ly/(n_-1)
    y_local = np.arange(0,ly+y_step_local-0.0001,y_step_local)
    Uarray_with_different_y = iters(m, n_, lx, ly, eps)
    errors_y.append(MSE(x,y_local,Uarray_with_different_y))

plotErrors(errors_y,n,delta,'y')

errors_x = []
for m_ in range(m,m+delta):
    x_step_local = lx/(m_-1)
    x_local = np.arange(0,lx+x_step_local-0.0001,x_step_local)
    Uarray_with_different_x = iters(m_, n, lx, ly, eps)
    errors_x.append(MSE(x_local,y,Uarray_with_different_x))

plotErrors(errors_x,m,delta,'x')

"""## Метод Зейделя"""

def zeidel(m, n, lx, ly, eps):
    x_step = lx / (m - 1)
    y_step = ly / (n - 1)
    x_step2 = x_step**2
    y_step2 = y_step**2
    
    A = (x_step2 - x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    B = (x_step2 + x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    C = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    D = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    
    Uarray = autofill(m, n, lx)
    
    iters = 5000
    
    while iters > 0:
        iters -= 1
        maxdiff = 0
        
        for i in range(1, m - 1):
            for j in range(1, n - 1):
                newval = Uarray[j - 1][i] * A + Uarray[j + 1][i] * B + Uarray[j][i - 1] * C + Uarray[j][i + 1] * D
                maxdiff = max(maxdiff, abs(newval - Uarray[j][i]))
                Uarray[j][i] = newval
        
        if maxdiff < eps:
            break
    
    return Uarray

lx = math.pi / 2
ly = math.pi / 2
m = 25
n = 25
x_step = lx / (m - 1)
y_step = ly / (n - 1)
eps = 0.0001

x = np.arange(0,lx+x_step-0.0001,x_step)
y = np.arange(0,ly+y_step-0.0001,y_step)

U_zeid = zeidel(m, n, lx, ly, eps)

plotU(x,U_zeid)

MSE(x,y,U_zeid)

errors_y = []
delta = 15
for n_ in range(n,n + delta):
    y_step_local = ly/(n_-1)
    y_local = np.arange(0,ly+y_step_local-0.0001,y_step_local)
    Uarray_with_different_y = zeidel(m, n_, lx, ly, eps)
    errors_y.append(MSE(x,y_local,Uarray_with_different_y))

plotErrors(errors_y,n,delta,'y')

errors_x = []
for m_ in range(m,m+delta):
    x_step_local = lx/(m_-1)
    x_local = np.arange(0,lx+x_step_local-0.0001,x_step_local)
    Uarray_with_different_x = zeidel(m_, n, lx, ly, eps)
    errors_x.append(MSE(x_local,y,Uarray_with_different_x))

plotErrors(errors_x,m,delta,'x')

"""## Метод простых итераций с верхней релаксацией"""

def relax(m, n, lx, ly, eps, omega):
    x_step = lx / (m - 1)
    y_step = ly / (n - 1)
    x_step2 = x_step**2
    y_step2 = y_step**2
    
    A = (x_step2 - x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    B = (x_step2 + x_step*y_step2) / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    C = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    D = y_step2 / (2 * x_step2 + 2 * y_step2 - 3 * x_step2 * y_step2)
    
    Uarray = autofill(m, n, lx)
    
    iters = 5000
    
    while iters > 0:
        iters -= 1
        maxdiff = 0
        newU = np.copy(Uarray)
        
        for i in range(1, m - 1):
            for j in range(1, n - 1):
                newU[j][i] = Uarray[j - 1][i] * A + Uarray[j + 1][i] * B + Uarray[j][i - 1] * C + Uarray[j][i + 1] * D
                newU[j][i] = omega * newU[j][i] + (1 - omega) * Uarray[j][i]
                maxdiff = max(maxdiff, abs(newU[j][i] - Uarray[j][i]))
        
        Uarray = newU
        if maxdiff < eps:
            break
    
    return Uarray

lx = math.pi / 2
ly = math.pi / 2
m = 25
n = 25
x_step = lx / (m - 1)
y_step = ly / (n - 1)
eps = 0.0001
omega = 1.0001

x = np.arange(0,lx+x_step-0.0001,x_step)
y = np.arange(0,ly+y_step-0.0001,y_step)

U_rel = relax(m, n, lx, ly, eps, omega)

plotU(x,U_rel)

MSE(x,y,U_rel)

errors_y = []
delta = 15
for n_ in range(n,n + delta):
    y_step_local = ly/(n_-1)
    y_local = np.arange(0,ly+y_step_local-0.0001,y_step_local)
    Uarray_with_different_y = relax(m, n_, lx, ly, eps, omega)
    errors_y.append(MSE(x,y_local,Uarray_with_different_y))

plotErrors(errors_y,n,delta,'y')

errors_x = []
for m_ in range(m,m+delta):
    x_step_local = lx/(m_-1)
    x_local = np.arange(0,lx+x_step_local-0.0001,x_step_local)
    Uarray_with_different_x = relax(m_, n, lx, ly, eps, omega)
    errors_x.append(MSE(x_local,y,Uarray_with_different_x))

plotErrors(errors_x,m,delta,'x')

"""Рассмотрим решения трёх методов одновременно"""

plt.figure(figsize=(10,6))
ydraw = 10
plt.plot(x,U_iter[ydraw])
plt.plot(x,U_zeid[ydraw])
plt.plot(x,U_rel[ydraw])
plt.legend(['Iterations','Zeidel','Relaxation'])
plt.title('Решения методов в y = ' + str(round(ydraw*y_step,4)))

