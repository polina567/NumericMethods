import numpy as np

def y(x):
    return x**2/(x**3-27)

# Метод прямоугольников
def rectangle(x0, xk, f, h):
    X = np.arange(x0, xk, h)
    return h * sum([f(X[i] + h / 2) for i in range(len(X))])

# Метод трапеций
def trapeze(x0, xk, f, h):
    X = np.arange(x0, xk, h)
    return h * ((f(X[0]) + f(xk)) / 2 + sum([f(X[i]) for i in range(1, len(X))]))

# Метод Симпсона
def Simpson(x0, xk, f, h):
    res = 0
    x = x0 + h
    while x < xk:
        res += f(x - h) + 4 * f(x) + f(x + h)
        x += 2 * h
    return (h / 3) * res

# Метод Рунге-Ромберга-Ричардсона
def RRR(F1, F2, h1, h2, p):
    if h1 < h2:
        return F1 + (F1 - F2) / ((h2 / h1) ** p - 1)
    return F2 + (F2 - F1) / ((h1 / h2) ** p - 1)

def main():
    print("Вариант 13:")
    print("y = x^2/(x^3-27)")
    print("X0 = -2, Xk = 2, h1 = 1.0, h2 = 0.5\n\n")
    x0 = -2; xk = 2; h1 = 1.0; h2 = 0.5
    # Метод второго порядка
    p = 2
    p1 = rectangle(x0, xk, y, h1)
    t1 = trapeze(x0, xk, y, h1)
    s1 = Simpson(x0, xk, y, h1)

    p2 = rectangle(x0, xk, y, h2)
    t2 = trapeze(x0, xk, y, h2)
    s2 = Simpson(x0, xk, y, h2)

    rp = RRR(p1, p2, h1, h2, p)
    rt = RRR(t1, t2, h1, h2, p)
    rs = RRR(s1, s2, h1, h2, 4)
    print('h1 = {}'.format(h1))
    print("Прямоугольник: {}".format(p1))
    print("Трапеция: {}".format(t1))
    print("Симпсон: {}".format(s1))

    print('\nh2 = {}'.format(h2))
    print("Прямоугольник: {}".format(p2))
    print("Трапеция: {}".format(t2))
    print("Симпсон: {}".format(s2))

    print("\nТочное значение по Рунге-Ромбергу-Ридчардсону для прямоугольника:", rp)
    print("Погрешность: {0} и {1}".format(abs(rp - p1), abs(rp - p2)))

    print("\nТочное значение по Рунге-Ромбергу-Ридчардсону для трапеции:", rt)
    print("Погрешность: {0} и {1}".format(abs(rt - t1), abs(rt - t2)))

    print("\nТочное значение по Рунге-Ромбергу-Ридчардсону для Симпсона:", rs)
    print("Погрешность: {0} и {1}".format(abs(rs - s1), abs(rs - s2)))
main()