import numpy as np
import math as m
import cmath
import numpy.linalg

def norm(x):
    return m.sqrt(sum([x_i**2 for x_i in x]))

def s(a, j):
    n = a.shape[0]
    res = 0
    for i in range(j, n):
        res += a[i, j] ** 2
    return res ** 0.5

def s1(a, j):
    n = a.shape[0]
    res = 0
    for i in range(j+1, n):
        res += a[i, j] ** 2
    return res ** 0.5

def find_nu(A, i, j, nu):
    if i == j:
        nu[j] = A[i][j, j] + np.sign(A[i][j, j]) * s(A[i], j)
    elif j < i:
        nu[j] = 0
    else:
        nu[j] = A[i][j, i]
    return nu

def householder(a):
    n = a.shape[0]
    H = [np.zeros((n, n))]*(n-1)
    A = [np.zeros((n, n))]*(n)
    A[0] = a
    for i in range(n-1):
        nu = np.zeros(n)
        for j in range(n):
            nu = find_nu(A, i, j, nu)
        nu = nu.reshape((1, 3))
        nu = nu.transpose()
        H[i] = np.eye(n) - 2*np.multiply(nu, nu.transpose()) / np.dot(nu.transpose(), nu)
        A[i+1] = np.matmul((H[i]), (A[i]))
    Q = H[0]
    for i in range(1, n-1):
        Q = Q.dot(H[i])
    return Q, A[-1]

def answer(A, eps, n):
    js = []
    for i in range(n):
        if (s1(A, i)) < eps:
            sol = A[i, i]
            print(sol)
        else:
            js.append(i)
            if i > 0:
                js.append(2)
                break
    if i == [0, 1]:
        sol = A[i, i]
    a = 1
    b = -(A[js[0], js[0]] + A[js[1], js[1]])
    c = A[js[0], js[0]] * A[js[1], js[1]] - A[js[0], js[1]] * A[js[1], js[0]]

    d = (b ** 2) - (4 * a * c)
    print(d, 'd')

    sol1 = (-b - cmath.sqrt(d)) / (2 * a)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    print(sol1, sol2)


def QR(a, eps):
    n = 3
    q, r = householder(a)
    A = np.dot(r, q)
    i = 0
    eps = 0.001
    while i < 100 :
        q, r = householder(A)
        A = np.dot(r, q)
        i += 1
        print(A)
    sol, sol1, sol2 = answer(A, eps, n)
    return i, A, sol, format(sol1), format(sol2)



def main():
    with open('1_5.txt', 'r') as f:
        n = [int(el) for el in f.readline().split()]
        n = n[0]
        a = np.zeros((n, n))
        for i in range(0, n):
            a[i, :] = [int(el) for el in f.readline().split()]
        eps = 0.001

    print(QR(a, eps))



if __name__ == "__main__":
    main()