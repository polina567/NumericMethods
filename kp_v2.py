import numpy as np
import matplotlib.pyplot as plt


def K(x, s):
    return np.exp(-x - s)


def f(x):
    return np.exp(-x)


def y(x):
    return np.exp(-x)


def main():
    print("Решение нелинейного интегрального уравнения\n")

    a = np.double(input("Введите левую границу а\n"))
    b = np.double(input("Введите правую границу b\n"))
    h = np.double(input("Введите шаг h\n"))

    r = int((a + b) / h)
    X = np.array(range(r), np.longdouble)
    for i in range(0, r):
        if i != int(b) and i != 0:
            X[i] = X[i - 1] + h
        if i != int(b) and i == 0:
            X[i] = a

    Yp = np.array(range(r), np.longdouble)
    for i in range(1, r):
        s = 0
        for j in range(1, i - 1):
            s += h * K(X[i], X[j]) * (y(j) * y(j))

        Yp[i] = (1 + np.sqrt(1 - 2 * h / (f(X[i]) + s))) / h

    Ym = np.array(range(r), np.longdouble)
    for i in range(1, r):
        s = 0
        for j in range(1, i - 1):
            s += h * K(X[i], X[j]) * (y(j) * y(j))
        Ym[i] = (1 - np.sqrt(1 - 2 * h / (f(X[i]) + s))) / h

    print("Найденные значения y[i]:\n")

    for i in range(1, r):
        if not np.isnan(Ym[i]):
            print("Y[", i, "] =", Ym[i], "\n")

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 9
    plt.rcParams["figure.figsize"] = fig_size

    plt.title("Графическое представление")
    plt.grid()
    plt.scatter(b, max(Ym), c="orange")
    plt.scatter(X, Ym, c="orange")
    plt.scatter(0, 0, c='black')
    plt.plot([-0.25, 300], [0, 0])
    plt.plot([0, 0], [-3, 300])
    plt.xlim(-0.25, (b + 0.1))
    plt.ylim(-3, (max(Ym)+3))
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(['Ox', 'Oy', 'Y_i'])
    plt.show()


main()
