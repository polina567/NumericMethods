import numpy as np
import math

def zero_estimate():
    a = 0.5
    b = 1
    x_0 = (a + b) / 2
    return a, b, x_0

def function(x):
    f = math.log(x + 1, math.e) - 2 * x + 0.5
    return f

def derivative(x):
    df = 1 / (x + 1) - 2
    return df

def function_l(x):
    f = x + function(x) * 0.5
    return f

def derivative_l(x):
    df = 1 + derivative(x) * 0.5
    return df


def Newton(eps):
    a, b, x_cur = zero_estimate()
    x_nex = x_cur - function(x_cur) / derivative(x_cur)
    i = 1
    while abs(x_nex - x_cur) >= eps:
        x_cur = x_nex
        x_nex = x_cur - function(x_cur) / derivative(x_cur)
        i += 1
    return x_nex, i

def Simple(eps):
    a, b, x_0 = zero_estimate()
    x_nex = function_l(x_0)
    q = max(abs(derivative_l(a)), abs(derivative_l(b)))
    x_cur = x_0
    i = 1
    while (q / (1 - q) * abs(x_nex - x_cur) > eps):
        x_cur = x_nex
        x_nex = function_l(x_cur)
        i += 1
    return x_nex, i




def main():
    #ln(x + 1) - 2*x^2 + 1 = 0
    eps = 0.00000000001
    print('Метод Ньютона: ', Newton(eps))
    print('Метод Простых Итераций: ', Simple(eps))




if __name__ == "__main__":
    main()