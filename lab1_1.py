import numpy as np
import math
import matplotlib.pyplot as plt


def f(x):
    return np.arcsin(x) + x

def omega(x, X):
    w = 1
    for i in range(len(X)):
        w *= (x - X[i])
    return w

def omega_(x, X):
    w = 1
    for i in range(len(X)):
        if X[i] != x:
            w *= (x - X[i])
    return w

def L(x, X):
    """Интерполяционный многочлен Лагранжа"""
    l = 0.0
    for i in range(len(X)):
        l += (f(X[i]) * omega(x, X)) / ((x - X[i]) * omega_(X[i], X))
    return l

def X_new(i, k, X):
    return [X[j] for j in range(i, k)]

def separate(X):
    if len(X) == 2:
        return (f(X[0]) - f(X[1])) / (X[0] - X[1])
    else:
        return (separate(X_new(0, len(X) - 1, X)) - separate(X_new(1, len(X), X))) / (X[0] - X[len(X) - 1])


def xxx(x, i, X):
    res = 1
    for j in range(i):
        res *= (x - X[j])
    return res

def P(x, X):
    p = f(X[0])
    for i in range(1, len(X)):
        X_ = X_new(0, i + 1, X)
        p += xxx(x, i, X) * separate(X_)
    return p




def main():
    print("Вариант 13: y=arcsin(x)+x, X* = 0.1")
    print("a) Xi = -0.4, -0.1, 0.2, 0.5")
    print("b) Xi = -0.4, 0, 0.2, 0.5")
    X_a = np.array([-0.4, -0.1, 0.2, 0.5])
    X_b = np.array([-0.4, 0.0, 0.2, 0.5])
    X = 0.1
    Yi = np.array([])
    print("Выберети входные данные (a или b):")
    while True:
        bl = input()
        if bl == "a":
            print("a) Xi = -0.4, -0.1, 0.2, 0.5")
            print("Xi = ", X_a)
            for i in range(0, 4):
                 Yi = np.append(Yi, f(X_a[i]))
            print("Yi = ", Yi)
            print("\nМногочлен Лагранжа:")
            print("f(x*) = ", f(X))
            print("L(x*) = ", L(X, X_a))
            print("Абсолютная погрешность delta = ", f(X) - L(X, X_a))
            print("\nМногочлен Ньютона")
            print("f(x*) = ", f(X))
            print("P(x*) = ", P(X, X_a))
            print("Абсолютная погрешность delta = ", f(X) - P(X, X_a))
            # ГРАФИК ФУНКЦИИ
            xmin = -1
            xmax = 1
            dx = 0.01

            xarr = np.arange(xmin, xmax, dx)
            ylist = [f(x) for x in xarr]
            y_X_a = [f(x) for x in X_a]
            Larr = X_a
            Llist = [L(x + 0.001, X_a) for x in Larr]
            Parr = X_a
            Plist = [P(x, X_a) for x in Parr]

            fig = plt.figure(figsize=(8, 6))
            grid = plt.grid(True)

            plt.title('f(x)')
            plt.plot(xarr, ylist)
            plt.plot(Larr, Llist)
            plt.plot(Parr, Plist)
            plt.plot(X_a, y_X_a, '*')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend(['f(x)', 'L(x)', 'P(x)', 'Xa', 'Ox', 'Oy'])
            plt.show()
        elif bl == "b":
            print("b) Xi = -0.4, 0, 0.2, 0.5")
            print("Xi = ", X_b)
            for i in range(0, 4):
                Yi = np.append(Yi, f(X_b[i]))
            print("Yi = ", Yi)
            print("\nМногочлен Лагранжа:")
            print("f(x*) = ", f(X))
            print("L(x*) = ", L(X, X_b))
            print("Абсолютная погрешность delta = ", f(X) - L(X, X_b))
            print("\nМногочлен Ньютона")
            print("\nМногочлен Ньютона")
            print("f(x*) = ", f(X))
            print("P(x*) = ", P(X, X_b))
            print("Абсолютная погрешность delta = ", f(X) - P(X, X_b))
            # ГРАФИК ФУНКЦИИ
            xmin = -1
            xmax = 1
            dx = 0.01

            xarr = np.arange(xmin, xmax, dx)
            ylist = [f(x) for x in xarr]
            y_X_b = [f(x) for x in X_b]
            Larr = X_b
            Llist = [L(x + 0.001, X_b) for x in Larr]
            Parr = X_b
            Plist = [P(x, X_b) for x in Parr]
            fig = plt.figure(figsize=(8, 6))
            grid = plt.grid(True)

            plt.title('f(x)')
            plt.plot(xarr, ylist)
            plt.plot(Larr, Llist)
            plt.plot(Parr, Plist)
            plt.plot(X_b, y_X_b, '*')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend(['f(x)', 'L(x)', 'P(x)', 'Xb', 'Ox', 'Oy'])
            plt.show()
        else:
            print("Неверный ввод!")


main()
