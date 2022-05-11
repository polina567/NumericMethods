import numpy as np

def sweep(a, b):

    p = np.zeros(len(b)); q = np.zeros(len(b))
    p[0] = -a[0][1]/a[0][0]
    q[0] = b[0]/a[0][0]
    for i in range(1, len(p)-1):
        p[i] = -a[i][i+1]/(a[i][i] + a[i][i-1]*p[i-1])
        q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i] + a[i][i-1]*p[i-1])
    i = len(a)-1
    p[-1] = 0
    q[-1] = (b[-1] - a[-1][-2]*q[-2])/(a[-1][-1] + a[-1][-2]*p[-2])
    
    x = np.zeros(len(b))
    x[-1] = q[-1]
    for i in reversed(range(len(b)-1)):
        x[i] = p[i]*x[i+1] + q[i]
    
    return x
    
print('1 - решить пример из задания \n2 - ввести свой')
answer = int(input())
while answer != 1 and answer != 2:
    answer = int(input('Введите 1 или 2\n'))
    
if answer == 1:
    
    a = np.zeros((5, 5))
    b = np.zeros(5)

    a[0][0] = 14; a[0][1] = 9;b[0] = 125
    a[1][0] = -8; a[1][1] = 14; a[1][2] = 6; b[1] = -56
    a[2][1] = -5; a[2][2] = -17; a[2][3] = 8; b[2] = 144
    a[3][2] = 1; a[3][3] = 5; a[3][4] = -2; b[3] = 36
    a[4][3] = -4; a[4][4] = -10; b[4] = 70
    
    print('Матрица', a); print('Правая чаасть', b)

    
else:
    n = 'a'
    n = int(input('Введите размерность\n'))
    a = np.zeros((n, n))
    b = np.zeros(n)
    a[0][0] = int(input('Введите b\n'))
    a[0][1] = int(input('Введите c\n'))
    b[0] = int(input('Введите правую часть\n'))
    for i in range(1, n-1):
        a[i][i-1] = int(input('Введите а\n'))
        a[i][i] = int(input('Введите b\n'))
        a[i][i+1] = int(input('Введите c\n'))
        b[i] = int(input('Введите правую часть\n'))
    a[-1][-2] = int(input('Введите a\n'))
    a[-1][-1] = int(input('Введите b\n'))
    b[-1] = int(input('Введите правую часть\n'))

print('Решение: ', sweep(a, b))
        
