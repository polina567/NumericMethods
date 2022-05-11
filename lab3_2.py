import numpy as np
import matplotlib.pyplot as plt

# Метод прогонки
def Lab1_2(a, b):
    p = np.zeros(len(b))
    q = np.zeros(len(b))

    p[0] = -a[0][1]/a[0][0]
    q[0] = b[0]/a[0][0]

    for i in range(1, len(p)-1):
        p[i] = -a[i][i+1]/(a[i][i] + a[i][i-1]*p[i-1])
        q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i] + a[i][i-1]*p[i-1])
    i = len(a)-1
    p[-1] = 0
    q[-1] = (b[-1] - a[-1][-2]*q[-2])/(a[-1][-1] + a[-1][-2]*p[-2])

    x = np.zeros(len(b))
    x[-1] = q[-1]
    for i in reversed(range(len(b)-1)):
        x[i] = p[i]*x[i+1] + q[i]
    return x

def S(x, X, F, bl):
    """Кубический сплайн"""

    h = np.array([])
    n = len(X)
    for i in range(1, n):
        h = np.append(h, X[i] - X[i - 1])

    # Система относительно Ci, i=2,...,n с трёхдиагональной матрицей:
    c_syst = np.array([])
    rows = np.array([])

    for i in range(n - 2):
        if i == 0:
            rows = np.append(rows, 2 * (h[i] + h[i + 1]))
        elif i == 1:
            rows = np.append(rows, h[i+1])
        else:
            rows = np.append(rows, 0)

    c_syst = np.append(c_syst, rows)

    for i in range(1, n - 3):
        rows = np.array([])
        for j in range(n - 2):
            if i - 1 == j:
                rows = np.append(rows, h[i-1])
            elif i == j:
                rows = np.append(rows, 2 * (h[i - 1] + h[i]))
            elif i + 1 == j:
                rows = np.append(rows, h[i])
            else:
                rows = np.append(rows, 0)
        c_syst = np.append(c_syst, rows)

    rows = np.array([])

    for i in range(2, n):
        if i == n - 2:
            rows = np.append(rows, h[i])
        elif i == n - 1:
            rows = np.append(rows, 2 * (h[i - 2] + h[i - 1]))
        else:
            rows = np.append(rows, 0)

    c_syst = np.append(c_syst, rows)

    b = np.array([])

    for i in range(2, n):
        b = np.append(b, 3 * ((F[i] - F[i - 1]) / h[i - 1] - (F[i - 1] - F[i - 2]) / h[i - 2]))
    c_matrix = np.zeros((3, 3))
    k = 0
    for i in range(3):
        for j in range(3):
            c_matrix[i][j] = c_syst[k]
            k += 1
    c = np.array([])
    c = np.append(c, 0)
    c = np.append(c, Lab1_2(c_matrix, b))

    a = np.array([])
    b = np.array([])
    d = np.array([])

    for i in range(n - 1):
        a = np.append(a, F[i])
        if i == n - 2:
            b = np.append(b, (F[i + 1] - F[i]) / h[i] - (2 / 3) * h[i] * c[i])
            d = np.append(d, - c[i] / (3 * h[i]))
        else:
            b = np.append(b, (F[i + 1] - F[i]) / h[i] - (1 / 3) * h[i] * (c[i + 1] + 2 * c[i]))
            d = np.append(d, (c[i + 1] - c[i]) / (3 * h[i]))

    if bl == 1:
        print("\nКоэффициенты сплайнов:")
        print("a = ", a)
        print("b = ", b)
        print("c = ", c)
        print("d = ", d)

    for i in range(n - 1):
        if (x >= X[i]) & (x <= X[i + 1]):
            res = a[i] + b[i] * (x - X[i]) + c[i] * (x - X[i]) ** 2 + d[i] * (x - X[i]) ** 3
            break
    return res

def main():
    
    print("Вариант 13:")
    print("-----------------------------------------------------")
    print("i |     0    |     1    |    2    |   3    |   4    |")
    print("xi|    0.0   |    1.0   |   2.0   |  3.0   |  4.0   |")
    print("fi|    1.0   |  1.5403  |  1.5839 |  2.01  | 3.3464 |")
    print("-----------------------------------------------------")
    print("x* = 1.5")
    x = 1.5
    X = np.array([ 0.0, 1.0, 2.0, 3.0, 4.0])
    F = np.array([1.0, 1.5403, 1.5839, 2.01, 3.3464])

    print('f(x*) = ', S(x, X, F, 1))

    # ГРАФИК ФУНКЦИИ
    xmin = 0.0
    xmax = 4.0
    dx = 0.05

    xarr = np.arange(xmin, xmax, dx)
    ylist = [S(x_, X, F, 0) for x_ in xarr]


    fig = plt.figure(figsize=(10, 8))
    grid = plt.grid(True)

    plt.title('S(x)')
    plt.plot(xarr, ylist)
    plt.plot(X, F, '*')
    plt.plot([-0.5, 1], [0, 0])
    plt.plot([0, 0], [-1, 2])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(['S(x)', 'Xi', 'Ox', 'Oy'])
    plt.show()
main()