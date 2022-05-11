import numpy as np
import matplotlib.pyplot as plt

# Функция для решения нелинейного уравнения
# методом квадратур. Используется  формула трапеций с равноотстоящими узлами
# Входные данные: K - ядро уравнения, f - правая
# часть (задаются аналитически), a - начало  отрезка интегрирования, b - конец отрезка, h - шаг
# сетки. Результат - вектор y приближений к  решению в узлах сетки


def f(x):
    return np.power(np.exp, -x + 4 * (np.power(x, 2)))


def y(x):
    return np.power(np.exp, (x - np.exp))


def df(x):
    return np.exp ** (-x + 4 * (np.power(x, 2))) * (-1 + 8 * x)


def trapeze(a, b, f, h):
    X = np.arange(a, b, h)
    return h * ((f(X[0]) + f(b)) / 2 + sum([f(X[i]) for i in range(1, len(X))]))


def main():
    print("Решение нелинейного интегрального уравнения\n")

    a = np.double(input("Введите левую границу а\n"))
    b = np.double(input("Введите правую границу b\n"))
    h = np.double(input("Введите шаг h\n"))

    r = int((a + b) / h)
    X = np.array(range(r), np.longdouble)
    for i in range(0, r):
        if i != int(b) & i != 0:
            X[i] = X[i-1] + h
        if i != int(b) & i == 0:
            X[i] = a

    K = np.array(range(r), np.longdouble)

    Yp = np.array(range(r), np.longdouble)
    for i in range(1, r):
        if i == 1:
            s = np.longdouble(0).value
            K[i] = (f(X[i])) * y(X[i]) * f(s)
        else:
            for j in range(1, i - 1):
                s = s + h * K[i] * K[j]
                K[i] = K[i - 1] + f(X[i]) * y(X[i]) * f(s)
        Yp[i] = (1 + np.sqrt(1 - 2 * h * (f(X[i]) + s))) / h


    K = np.array(range(r), np.longdouble)

    Ym = np.array(range(r), np.longdouble)
    for i in range(1, r):
        if i == 1:
            s = np.longdouble(0).value
            K[i] = f(X[i]) * y(X[i])
        else:
            for j in range(1, i - 1):
                s = s + h * K[i] * K[j]
                K[i] = K[i - 1] + f(X[i]) * y(X[i]) * f(s)
        Ym[i] = (1 - np.sqrt(1 - 2 * h * (f(X[i]) + s))) / h
    return Ym[i]

    print("Найденные значения y[i]:\n")

    print("1) Положительные значения y[i]:")
    for i in range(1, r):
        print("Y[", i, "] = ", Yp[i], "\n")

    print("2) Отрицательные значения y[i]:")
    for i in range(1, r):
        print("Y[", i, "] = ", Ym[i], "\n")


main()
