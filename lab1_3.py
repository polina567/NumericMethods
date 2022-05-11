import numpy as np

def ab(a, b):
    n = a.shape[0]
    alpha = np.zeros((n, n))
    beta = np.zeros(n)
    for i in range(n):
        alpha[i] = -a[i] / a[i][i]
        beta[i] = b[i]/a[i][i]
        alpha[i][i] = 0
    return alpha, beta
    
def norm_1(a):
    res = 0
    for i in range(a.shape[0]):
        s = 0
        try:
            for j in range(a.shape[0]):
                s += abs(a[j][i])
            if s > res:
                res = s
        except:
            res += abs(a[i])
    return res
    
def norm_2(a):
    res = 0
    for i in range(a.shape[0]):
        try:
            for j in range(a.shape[0]):
                res += a[i][j]**2
        except:
            res += a[i]**2
    return res**0.5
    
def norm_c(a):
    res = 0
    for i in range(a.shape[0]):
        s = 0
        try:
            for j in range(a.shape[0]):
                s += abs(a[i][j])
            if s > res:
                res = s
        except:
            if abs(a[i]) > res:
                res = abs(a[i])
    return res
    

print('1 - решить пример из задания \n2 - ввести свой')
answer = int(input())
while answer != 1 and answer != 2:
    answer = int(input('Введите 1 или 2\n'))
    
if answer == 1:
    a = np.zeros((4, 4))
    a[0][0] = 24; a[0][1] = -7; a[0][2] = -4; a[0][3] = 4
    a[1][0] = -3; a[1][1] = -9; a[1][2] = -2; a[1][3] = -2
    a[2][0] = 3; a[2][1] = 7; a[2][2] = 24; a[2][3] = 9
    a[3][0] = 1; a[3][1] = -6; a[3][2] = -2; a[3][3] = -15

    b = np.zeros(4)
    b[0] = -190; b[1] = -12; b[2] = 155; b[3] = -17
    
else:
    n = 'a'
    n = int(input('Введите размерность\n'))
    a = np.zeros((n, n))
    b = np.zeros(n)
    
    for i in range(n):
        s = input('Введите {}-ю строку\n'.format(i+1))
        a[i] = [int(s_i) for s_i in s.split()]
        b[i] = int(input('Введите правую часть\n'))
    
print('Матрица A', a); print('Правая чаасть', b)
alpha, beta = ab(a, b)
print('Alpha:\n', alpha)
print('beta\n', beta)
print('Нормы\n', norm_1(alpha), norm_2(alpha), norm_c(alpha))
eps = float(input('Введите эпсилон:\n'))
print('Метод простых итераций.')
prev = beta
curr = beta + alpha.dot(beta)
i = 0
n_alpha = norm_c(alpha)/(1 - norm_c(alpha))
while n_alpha*norm_c(curr - prev) > eps:
    prev = curr
    curr = beta + alpha.dot(prev)
    print('Текущее приближение: ', curr)
    i += 1
print('Решение закончилось на {}-й итерации'.format(i))

print('Решение x:\n', curr)

print('Проверка, A*x', a.dot(curr))

print('Метод Зейделя')
k = 0
prev = beta
converge = False
m = alpha.shape[0]
while not converge:
    curr = np.copy(prev)
    for i in range(m):
        s = 0
        for j in range(m):
            s += alpha[i, j] * curr[j]
        curr[i] = beta[i] + s
    print('Текущее приближение: ', curr)
    converge =  n_alpha*norm_c(curr - prev) <= eps
    k += 1
    prev = np.copy(curr)

print('Решение закончилось на {}-й итерации'.format(k))

print('Решение x:\n', curr)

print('Проверка, A*x', a.dot(curr))
