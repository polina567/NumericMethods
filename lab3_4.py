import math
import numpy as np

# Формулы из методички (По 3 точкам)
def diff(x, X, Y):
    if len(X) != len(Y):
        raise Exception("Несоответсвие размерностей")

    for i in range(len(X) - 1):
        if (x >= X[i]) & (x <= X[i + 1]):
            if (x == X[i]) & (i != 0):
                left = (Y[i] - Y[i - 1]) / (X[i] - X[i - 1])
                right = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])

                mean = (left + right) / 2

                print('Левая производная: {}'.format(left))
                print('Правая производная: {}'.format(right))
                print('Среднее: {}'.format(mean))

            elif (x == X[i + 1]) & (i != len(X) - 1):
                left = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])
                right = (Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1])

                mean = (left + right) / 2

                print('Левая производная: {}'.format(left))
                print('Правая производная: {}'.format(right))
                print('Среднее: {}'.format(mean))
            else:
                df_1 = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])

                print('1-я производная от функции первой степени: {}'.format(df_1))
            df_2 = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) + (((Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])) / (X[i + 2] - X[i])) * (2 * x - X[i] - X[i + 1])
            ddf = 2 * (((Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])) / (X[i + 2] - X[i]))

            print('1-я производная: {}'.format(df_2))
            print('2-я производная: {}'.format(ddf))
            break


# Разделённая разность f(x0,x1,x2,....)
def SeparatedDiff(X,Y):
    n = len(X)
    if n < 2:
        if n == 1: return Y[0]
        return None
    if n > 2:
        return (SeparatedDiff(X[:n-1], Y[:n-1]) - SeparatedDiff(X[1:], Y[1:]))/(X[0] - X[n-1])
    return (Y[0] - Y[1])/(X[0] - X[1])

# Произведение всех (x-xi),
def BracketsMult(x, X):
    ans = 1
    for i in range(len(X)):
        ans *= (x - X[i])
    return ans

# Многочлен Ньютона
def P(x,X,Y):
    n = len(X)
    ans = 0
    for i in range(0, n):
        ans += SeparatedDiff(X[:i+1], Y[:i+1])*BracketsMult(x,X[:i])
    return ans

# ((x-x1)(x-x2) + (x-x0)(x-x2) + (x-x0)(x-x1)) Нужно для нахождения первой производной
def BracketsMult1(x,X):
    ans = 0
    for i in range(len(X)):
        ans += BracketsMult(x, X[:i]+X[i+1:])
    return ans

# Нахождение первой производной
def Der1(x,X,Y):
    n = len(X)
    ans = 0
    for i in range(1,n):
        ans += SeparatedDiff(X[:i+1], Y[:i+1])*BracketsMult1(x,X[:i])
    return ans


def DerSeveralBrackets(x, X):
    n = len(X)
    if n == 1: return 1#если скобка одна
    return BracketsMult(x, X[1:]) + BracketsMult(x, X[:1])*DerSeveralBrackets(x, X[1:])

def BracketsMult2(x,X):
    ans = 0
    for i in range(len(X)):
        ans += DerSeveralBrackets(x, X[:i]+X[i+1:])
    return ans

# Нахождение второй производной
def Der2(x, X, Y):
    n = len(X)
    ans = 0
    for i in range(2, n):
        ans += SeparatedDiff(X[:i+1], Y[:i+1])*BracketsMult2(x,X[:i])
    return ans


def main():
    print("Вариант 13:")
    print("--------------------------------------------")
    print("i |    0    |    1   |    2   |    3    |    4   |")
    print("xi|   0.2   |   0.5  |   0.8  |   1.1   |   1.4  |")
    print("yi|  12.906 | 5.5273 | 3.8777 |  3.2692 | 3.0319 |")
    print("--------------------------------------------")
    print("x* = 0.8\n")
    x_star = 0.8
    X = [0.2, 0.5, 0.8, 1.1, 1.4]
    Y = [12.906, 5.5273, 3.8777,  3.2692, 3.0319]
    print("По 3 точкам:")
    print("_______________________________________")
    diff(x_star, X, Y)
    print("_______________________________________")

    print("\nПо 5 точкам:")
    print("_______________________________________")
    d1 = Der1(x_star, X, Y)
    d2 = Der2(x_star, X, Y)
    print(
        "Интерполяционный полином = {0}\
        \n1-я производная = {1}\
        \n2-я производная = {2}".format(P(x_star, X, Y), d1, d2))
    print("_______________________________________")
main()

