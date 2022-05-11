import numpy as np
import math as math
import matplotlib.pyplot as plt

def init(x, y, a, b, mu):
    return 0
def left(y, t, a, b, mu):
    return 0
def right(y, t, a, b, mu):
    return -math.sin(y) * math.sin(mu * t)
def up(x, t, a, b, mu):
    return -math.sin(x) * math.sin(mu * t)
def down(x, t, a, b, mu):
    return 0
def f(x, y, t, a, b, mu):
    return math.sin(x) * math.sin(y) * (mu * math.cos(mu * t) + (a + b) * math.sin(mu * t))
def solution(x, y, t, mu):
    return math.sin(x) * math.sin(y) * math.sin(mu * t)

"""# Решение трехдиагональной системы"""

def three_diag(a, b, c, d):
    P = [ -c[0] / b[0] ]
    Q = [ d[0] / b[0] ]
    x = [ 0 ]
    n = len(a)
    for i in range(1, n, 1):
        den = b[i] + a[i] * P[i - 1]
        P.append(-c[i] / den)
        Q.append((d[i] - a[i] * Q[i - 1]) / den)
        x.append(0)
    x[n - 1] = Q[n - 1]
    for i in range(n-2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
    return x

"""# Решатель задачи
method = 1 или 2  
approx = 1 или 2
"""

def solve(a, b, mu, hx, hy, tau, T, method, approx):
    u = []
    temp = []

    for x in np.arange(0.0, np.pi, hx):
        ux = []
        tempx = []
        for y in np.arange(0.0, np.pi, hy):
            uxy = []
            for t in np.arange(0, T, tau):
                uxy.append(init(x, y, a, b, mu) if tau == 0 else 0)
            ux.append(uxy)
            tempx.append(0)
        u.append(ux)
        temp.append(tempx)
    lx = len(u)
    ly = len(u[0])
    lz = len(u[0][0])

    if method == 1:
        for k in range(1, lz):
            for j in range(1, ly - 1):
                da = [0]
                db = [0]
                dc = [0]
                dd = [0]
                e2 = 0
                for i in range(1, lx - 1):
                    da.append(-a / hx / hx)
                    db.append(2 * a / hx / hx + 2 / tau)
                    dc.append(-a / hx / hx)
                    dd.append(2 * u[i][j][k - 1] / tau +
                        b / hy / hy * (u[i][j - 1][k - 1] - 2 * u[i][j][k - 1] + u[i][j + 1][k - 1]) +
                        f(i * hx, j * hy, (k + 0.5) * tau, a, b, mu))
                da[0] = 0
                db[0] = 1
                dc[0] = 0
                dd[0] = left(j * hy, (k + 0.5) * tau, a, b, mu)
                if approx == 1:
                    da.append(-1 / hx)
                    db.append(1 / hx)
                    dc.append(0)
                    dd.append(right(j * hy, (k + 0.5) * tau, a, b, mu))
                else:
                    e2 = 0.5 / hx
                    da.append(-2 / hx)
                    db.append(1.5 / hx)
                    dc.append(0)
                    dd.append(right(j * hy, (k + 0.5) * tau, a, b, mu))
                da[lx - 1] -= e2 * db[lx - 2] / da[lx - 2]
                db[lx - 1] -= e2 * dc[lx - 2] / da[lx - 2]
                dd[lx - 1] -= e2 * dd[lx - 2] / da[lx - 2]
                sol = three_diag(da, db, dc, dd)
                for i in range(0, lx):
                    temp[i][j] = sol[i]
            for i in range(0, lx):
                temp[i][0] = down(i * hx, (k + 0.5) * tau, a, b, mu)
                if approx == 1:
                    temp[i][ly - 1] = temp[i][ly - 2] + hy * up(i * hx, (k + 0.5) * tau, a, b, mu)
                else:
                    temp[i][ly - 1] = (-temp[i][ly - 3] + 4 * temp[i][ly - 2] + 2 * hy * up(i * hx, (k + 0.5) * tau, a, b, mu)) / 3
            for i in range(1, lx-1):
                da = [0]
                db = [0]
                dc = [0]
                dd = [0]
                e2 = 0
                for j in range(1, ly-1):
                    da.append(-b / hy / hy)
                    db.append(2 * b / hy / hy + 2 / tau)
                    dc.append(-b / hy / hy)
                    dd.append(2 * temp[i][j] / tau +
                        a / hx / hx * (temp[i - 1][j] - 2 * temp[i][j] + temp[i + 1][j]) +
                        f(i * hx, j * hy, (k + 0.5) * tau, a, b, mu))
                da[0] = 0
                db[0] = 1
                dc[0] = 0
                dd[0] = down(i * hx, (k + 0.5) * tau, a, b, mu)
                if approx == 1:
                    da.append(-1 / hy)
                    db.append(1 / hy)
                    dc.append(0)
                    dd.append(up(i * hx, (k + 0.5) * tau, a, b, mu))
                else:
                    e2 = 0.5 / hy
                    da.append(-2 / hy)
                    db.append(1.5 / hy)
                    dc.append(0)
                    dd.append(up(i * hx, (k + 0.5) * tau, a, b, mu))
                da[ly - 1] -= e2 * db[ly - 2] / da[ly - 2]
                db[ly - 1] -= e2 * dc[ly - 2] / da[ly - 2]
                dd[ly - 1] -= e2 * dd[ly - 2] / da[ly - 2]
                sol = three_diag(da, db, dc, dd)
                for j in range(0, ly):
                    u[i][j][k] = sol[j]
            for j in range(0, ly):
                u[0][j][k] = left(j * hy, (k + 0.5) * tau, a, b, mu)
                if approx == 1:
                    u[lx - 1][j][k] = u[lx - 2][j][k] + hx * right(j * hy, (k + 0.5) * tau, a, b, mu)
                else:
                    u[lx - 1][j][k] = (-u[lx - 3][j][k] + 4 * u[lx - 2][j][k] + 2 * hx * right(j * hy, (k + 0.5) * tau, a, b, mu)) / 3
    # fractional steps
    else:
        for k in range(1, lz, 1):
            for j in range(0, ly, 1):
                da = [0]
                db = [0]
                dc = [0]
                dd = [0]
                e2 = 0
                for i in range(1, lx-1, 1):
                    da.append(-a / hx / hx)
                    db.append(2 * a / hx / hx + 1 / tau)
                    dc.append(-a / hx / hx)
                    dd.append(u[i][j][k - 1] / tau +
                        f(i * hx, j * hy, k * tau, a, b, mu) / 2)
                da[0] = 0
                db[0] = 1
                dc[0] = 0
                dd[0] = left(j * hy, k * tau, a, b, mu)
                if approx == 1:
                    da.append(-1 / hx)
                    db.append(1 / hx)
                    dc.append(0)
                    dd.append(right(j * hy, k * tau, a, b, mu))
                else:
                    e2 = 0.5 / hx
                    da.append(-2 / hx)
                    db.append(1.5 / hx)
                    dc.append(0)
                    dd.append(right(j * hy, k * tau, a, b, mu))
                da[lx - 1] -= e2 * db[lx - 2] / da[lx - 2]
                db[lx - 1] -= e2 * dc[lx - 2] / da[lx - 2]
                dd[lx - 1] -= e2 * dd[lx - 2] / da[lx - 2]
                sol = three_diag(da, db, dc, dd)
                for i in range(0, lx, 1):
                    temp[i][j] = sol[i]
            for i in range(0, lx, 1):
                da = [0]
                db = [0]
                dc = [0]
                dd = [0]
                e2 = 0
                for j in range(1, ly-1, 1):
                    da.append(-b / hy / hy)
                    db.append(2 * b / hy / hy + 1 / tau)
                    dc.append(-b / hy / hy)
                    dd.append(temp[i][j] / tau +
                        f(i * hx, j * hy, (k + 1) * tau, a, b, mu) / 2)
                da[0] = 0
                db[0] = 1
                dc[0] = 0
                dd[0] = down(i * hx, k * tau, a, b, mu)
                if approx == 1:
                    da.append(-1 / hy)
                    db.append(1 / hy)
                    dc.append(0)
                    dd.append(up(i * hx, k * tau, a, b, mu))
                else:
                    e2 = 0.5 / hy
                    da.append(-2 / hy)
                    db.append(1.5 / hy)
                    dc.append(0)
                    dd.append(up(i * hx, k * tau, a, b, mu))
                da[ly - 1] -= e2 * db[ly - 2] / da[ly - 2]
                db[ly - 1] -= e2 * dc[ly - 2] / da[ly - 2]
                dd[ly - 1] -= e2 * dd[ly - 2] / da[ly - 2]
                sol = three_diag(da, db, dc, dd)
                for j in range(0, ly, 1):
                    u[i][j][k] = sol[j]

    return u

"""# Тестирование решателя"""

def MSE(u, a, b, mu, hx, hy, tau, T):
    mse = 0.0
    for i in range(0, len(u)):
        for j in range(0, len(u[0])):
            for k in range(0, len(u[0][0])):
                mse += (u[i][j][k] - solution(i*hx, j*hy, k*tau, mu))**2
    return (mse / (hx * hy))**0.5

def h_plot(a=1, b=1, mu=1, tau=0.1, T=1.0, method=1, approx=1):
    h = []
    e = []
    for N in range(5, 50):
        u = solve(a, b, mu, np.pi/N, np.pi/N, tau, T, method, approx)
        h.append(math.pi / N)
        e.append(MSE(u, a, b, mu, np.pi/N, np.pi/N, tau, T))
    plt.plot(h, e)
    plt.xlabel("шаг разбиения h")
    plt.ylabel("невязка")
    plt.grid()
    
def t_plot(a=1, b=1, mu=1, hx=0.1, hy=0.1, T=1.0, method=1, approx=1):
    t = []
    e = []
    for N in range(5, 50):
        u = solve(a, b, mu, hx, hy, T/N, T, method, approx)
        t.append(T / N)
        e.append(MSE(u, a, b, mu, hx, hy, T/N, T))
    plt.plot(t, e)
    plt.xlabel("шаг разбиения t")
    plt.ylabel("невязка")
    plt.grid()

"""## Ошибки метода 1 с различным шагом tau по t"""

t_plot(a=1, b=1, mu=1, method=1, approx=1)

"""## Ошибки метода 2 с различным шагом tau по t"""

t_plot(a=2, b=2, mu=2, method=2, approx=2)

"""## Ошибки метода 1 с различным шагом h по x и по y"""

h_plot(a=1, b=1, mu=1, method=1, approx=1)

"""## Ошибки метода 1 с различным шагом h по x и по y"""

h_plot(a=2, b=2, mu=2, method=2, approx=2)

"""# Визуальное тестирование на 3-х мерных графиках"""

def plot_3d(a=1, b=1, mu=1, N=20, tau=0.1, T=1.0, method=1, approx=1, t=0):
    u = solve(a, b, mu, np.pi/N, np.pi/N, tau, T, method, approx)
    X = [[i for i in np.arange(0, np.pi, np.pi/N)] for _ in np.arange(0, np.pi, np.pi/N)]
    Y = [i for i in np.arange(0, np.pi, np.pi/N) for _ in np.arange(0, np.pi, np.pi/N)]
    experience = []
    analit = []
    for i in range(0, len(u)):
        for j in range(0, len(u[0])):
            experience.append(u[i][j][t])
            analit.append(solution(i/N*np.pi, j/N*np.pi, t*tau, mu))
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, experience, label="not analytical")
    ax.scatter(X, Y, analit, label='analytical')
    plt.legend()
    plt.grid()

plot_3d(a=1, b=1, mu=1, N=30, tau=0.1, T=1.0, method=1, approx=1, t=9)

plot_3d(a=1, b=1, mu=1, N=30, tau=0.1, T=1.0, method=1, approx=2, t=9)

plot_3d(a=1, b=1, mu=1, N=30, tau=0.1, T=1.0, method=2, approx=1, t=9)

plot_3d(a=2, b=1, mu=1, N=30, tau=0.1, T=1.0, method=1, approx=1, t=9)

plot_3d(a=1, b=2, mu=1, N=30, tau=0.1, T=1.0, method=1, approx=1, t=9)

plot_3d(a=1, b=1, mu=2, N=30, tau=0.1, T=1.0, method=1, approx=1, t=9)

def plot_3d_compare_methods(a=1, b=1, mu=1, N=20, tau=0.1, T=1.0, approx=1, t=0):
    u1 = solve(a, b, mu, np.pi/N, np.pi/N, tau, T, 1, approx)
    u2 = solve(a, b, mu, np.pi/N, np.pi/N, tau, T, 2, approx)
    X = [[i for i in np.arange(0, np.pi, np.pi/N)] for _ in np.arange(0, np.pi, np.pi/N)]
    Y = [i for i in np.arange(0, np.pi, np.pi/N) for _ in np.arange(0, np.pi, np.pi/N)]
    experience_1 = []
    experience_2 = []
    analit = []
    for i in range(0, len(u1)):
        for j in range(0, len(u1[0])):
            experience_1.append(u1[i][j][t])
            experience_2.append(u2[i][j][t])
            analit.append(solution(i/N*np.pi, j/N*np.pi, t*tau, mu))
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, experience_1, label="Метод переменных направлений", s=1)
    ax.scatter(X, Y, experience_2, label="Метод дробных шагов", s=1)
    ax.scatter(X, Y, analit, label='analytical', s=1)
    plt.legend()
    plt.grid()

plot_3d_compare_methods(a=1, b=1, mu=1, N=30, tau=0.1, T=1.0, approx=1, t=9)

"""# Выводы

В ходе лабораторной работы были построены решатели эллиптических уравнений (метод дробных шагов, метод переменных направлений). Были построены графики погрешности. Погрешность при умеьшении шага по времени уменьшалась. Погрешность при уменьшении шага по по координатам увеличивалась, но очень медленно (сумма квадратов отклонений по всем t, x, y очень невилика). Были построены 3-д графики, по которым видно, что приближенные решения сходятся к аналитическому.
"""