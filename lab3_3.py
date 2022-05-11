import numpy as np
import matplotlib.pyplot as plt

# LU разложение:
def LU(A):
    # Получаем размерность матрицы
    n = len(A)

    # Созжаем две нулевые матрицы размером nxn
    U = np.matrix(np.zeros([len(A), len(A)]))
    L = np.matrix(np.zeros([len(A), len(A)]))

    # Получаем матрицы U и L

    # Получаем верхне-треугольную матрицу U и добавляем коэффициенты в матрицу L
    for k in range(n):
        for j in range(k, n):
            U[k, j] = A[k, j] - L[k, :k] * U[:k, j]

        for i in range(k + 1, n):
            if i == k:
                L[i, k] = 1
            L[i, k] = (A[i, k] - L[i, : k] * U[: k, k]) / U[k, k]

    # Заполняем матрицы L единицами по диагонали
    for i in range(n):
        for k in range(n):
            if i == k:
                L[i, k] = 1
    return L, U

def Solve(A, b):
    n = len(A[0])
    L, U = LU(A)
    x = np.zeros(n)
    y = np.zeros(n)
    print("Коэффициенты:", A[0], "\n")
    def sum_y(i):
        res = 0
        for k in range(n):
            res += L[i, k] * y[k]
        return res
    for i in range(n):
        y[i] = b[i] - sum_y(i)

    def sum_x(i):
        res = 0
        for j in range(i, n):
            res += U[i, j]*x[j]
        return res

    for i in reversed(range(n)):
        x[i] = (y[i] - sum_x(i))/U[i, i]

    return x

# Метод наименьших квадратов
def MNK(X, Y, n):
    if len(X) != len(Y):
        raise Exception("Несоответсвие размерностей")
    A = np.zeros((n+1, n+1))
    b = np.zeros(n+1)

    for i in range(n + 1):
        for j in range(n + 1):
            A[i][j] = sum([x ** (i + j) for x in X])
        b[i] = sum([Y[k] * X[k] ** i for k in range(len(Y))])
    a = Solve(A, b)

    return lambda x: sum([a[i] * x ** i for i in range(n + 1)])

# Сумма квадратов ошибок
def err(f, X, Y):
    if len(X) != len(Y):
        raise Exception("Несоответсвие размерностей")
    return sum([(f(X[i]) - Y[i]) ** 2 for i in range(len(X))])

def main():
    print("Вариант 13:")
    print("---------------------------------------------------------------")
    print("i |     0    |     1    |     2    |    3    |   4    |   5    |")
    print("xi|   -1.0   |    0.0   |    1.0   |   2.0   |  3.0   |  4.0   |")
    print("fi|  -0.4597 |    1.0   |  1.5403  |  1.5839 | 2.010  | 3.3464 |")
    print("---------------------------------------------------------------")
    X = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
    Y = np.array([-1.4754, 1.0, 1.5403, 1.5839, 2.010, 3.3464])

    f_1 = MNK(X, Y, 1)
    f_2 = MNK(X, Y, 2)

    xmin = -2
    xmax = 5

    dx = 0.001
    xarr = np.arange(xmin, xmax, dx)

    print("Сумма квадратов ошибок для 1-ой степени = ", err(f_1, X, Y))
    ylist = [f_1(x) for x in xarr]
    plt.figure(figsize=(10,7))
    plt.title("Многочлен 1 степени")
    plt.grid()
    plt.scatter(X, Y)
    plt.plot(xarr, ylist)
    plt.plot([-2, 5], [0, 0])
    plt.plot([0, 0], [-6, 6])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(['f_1(x)', 'Ox', 'Oy', 'Yi'])
    plt.show()


    print("Сумма квадратов ошибок для 2-ой степени = ", err(f_2, X, Y))
    ylist = [f_2(x) for x in xarr]
    plt.figure(figsize=(10,7))
    plt.title("Многочлен 2 степени")
    plt.grid()
    plt.scatter(X, Y)
    plt.plot(xarr, ylist)
    plt.plot([-2, 5], [0, 0])
    plt.plot([0, 0], [-6, 6])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(['f_2(x)', 'Ox', 'Oy', 'Yi'])
    plt.show()



main()