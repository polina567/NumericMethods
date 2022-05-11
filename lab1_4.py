import numpy as np

def sum_kv(a):
    res = 0
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if i != j:
                res += a[i, j]**2
    return res**0.5

def find_max(a):
    res = 0
    for i in range(a.shape[0]):
        for j in range(a.shape[0]):
            if i != j and abs(a[i, j]) > res:
                res = abs(a[i, j])
                res_i = i
                res_j = j
    return res, res_i, res_j

def U(a, i, j):
    u = np.eye(a.shape[0])
    fi = np.pi/4 if a[i, i] == a[j, j] else 0.5*np.arctan((2*a[i, j])/(a[i, i] - a[j, j]))
    u[i, j] = -np.sin(fi)
    u[j, i] = np.sin(fi)
    u[i, i] = u[j, j] = np.cos(fi)
    return u

def mult(u):
    res = u[0]
    for i in range(1, len(u)):
        res = res.dot(u[i])
    return res

print('1 - решить пример из задания \n2 - ввести свой')
answer = int(input())
while answer != 1 and answer != 2:
    answer = int(input('Введите 1 или 2\n'))

if answer == 1:
    n = 3
    a = np.zeros((3, 3))
    a[0][0] = 8;
    a[0][1] = 0;
    a[0][2] = -2;
    a[1][0] = 0;
    a[1][1] = 5;
    a[1][2] = 4;
    a[2][0] = -2;
    a[2][1] = 4;
    a[2][2] = -6


else:
    n = 'a'
    n = int(input('Введите размерность\n'))
    a = np.zeros((n, n))

    for i in range(n):
        s = input('Введите {}-ю строку\n'.format(i + 1))
        a[i] = [int(s_i) for s_i in s.split()]

print('Матрица A', a);
eps = float(input('Введите эпсилон:\n'))
print('Метод вращения.')
u = []
while sum_kv(a) > eps:
    max_a, i, j = find_max(a)
    u.append(U(a, i, j))
    print('U: \n', u[-1])
    a = u[-1].T.dot(a).dot(u[-1])
    for i in range(n):
        for j in range(n):
            a[i, j] = round(a[i, j], 10)
    print('A на {}-й итерации:\n'.format(len(u)), a)

print('Собственные значения: ', [a[i, i] for i in range(n)])
eigen = [mult(u)[:, j] for j in range(n)]
print('Собственные векторы: \n', eigen)
print('Скалярные произведения собственных векторов друг на друга:')
for i in range(n):
    for j in range(i, n):
        if i != j:
            print(eigen[i].dot(eigen[j]))