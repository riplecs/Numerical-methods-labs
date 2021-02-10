import scipy
from scipy import signal
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

func = lambda x: x**5+x**4-2*x**3-9*x**2-3*x-2
#plt.grid()
#plt.axis([-10, 10, -30, 30])
#plt.plot(x, func(x))
#plt.show()


print('Теорема про кільце:')
a_n=[-2, -3, -9, -2, 1, 1]
A=abs(max(a_n[1:], key=abs))
B=abs(max(a_n[:(len(a_n)-1)], key=abs))
print(f'A = max(|a_i|) (i=0, 1, ..., 4) = {A}\nB = max(|a_i|) (i=1, 2, ..., 5) = {B}')
left=abs(a_n[0])/(B+abs(a_n[0]))
rigth=(abs(a_n[len(a_n)-1])+A)/abs(a_n[len(a_n)-1])
print(f'Всі корені лежать у кільці: {left} =< |x*| =< {rigth}')


print('Теорема про верхню межу додатніх коренів:')
a_neg=[]
for i in range(len(a_n)):
    if a_n[i]<0: a_neg.append(a_n[i])
B=abs(max(a_neg, key=abs))
print(f'B = max(|a_i|) (a_i<0; i=0, ..., n) = {B}')
def find(a):
    m=None
    for i in reversed(range(len(a))):
        if a[i]<0:
            m=i
            print(f'm = max(i) (a_i<0; i=0, ..., n) = {i}')        
            return m
            break
    if m is None: print('Рівняння не має додатніх коренів.')

R=1+math.pow((B/a_n[len(a_n)-1]), 1/(len(a_n)-1-find(a_n)))
print(f'R = {R} - верхня межа додатніх коренів.')

print('Для знаходження нижньої межі додатніх коренів зробимо заміну x=1/y і отримаємо наступне рівняння:')
a_n=[-2, -3, -9, -2, 1, 1]
a_n.reverse()
a_n=[-a for a in a_n]
print(a_n)
R=1+math.pow((B/a_n[len(a_n)-1]), 1/(len(a_n)-1-find(a_n)))
print(f'R = {1/(R)} - нижня межа додатніх коренів.')


print("Для знаходження нижньої межі від'ємних коренів зробимо заміну x=-x і отримаємо наступне рівняння:")
a_n=[2, -3, 9, -2, -1, 1]
R=1+math.pow((B/a_n[len(a_n)-1]), 1/(len(a_n)-1-find(a_n)))
print(f"R = {-R} - нижня межа від'ємних коренів.")


print("Для знаходження верхньої межі від'ємних коренів зробимо заміну x=-1/y і отримаємо наступне рівняння:")
a_n.reverse()
R=1+math.pow((B/a_n[len(a_n)-1]), 1/(len(a_n)-1-find(a_n)))
print(f"R = {(-1/R)} - верхня межа від'ємних коренів.")

print('Теорема Гюа про наявність комплесних коренів:')
a_n=[-2, -3, -9, -2, 1, 1]
for i in range(1, len(a_n)-1):
    if a_n[i]**2<a_n[i-1]*a_n[i+1]:
        print(f'Існує таке k, що (a_k)^2<a_(k-1)*a_(k+1), k = {i}, отже рівняння має комплексні корені')
        break
    else: 
        print('Рівняння не має комплексних коренів.')    


print('Теорема Штурма:')
n=len(a_n)
f=np.poly1d([1, 1, -2, -9, -3, -2])
f0=f
f1=np.poly1d([5, 4, -6, -18, -3])
f2=np.poly1d(scipy.signal.deconvolve(f0, f1)[1])
f3=np.poly1d(scipy.signal.deconvolve(f1, f2)[1])
f4=np.poly1d(scipy.signal.deconvolve(f2, f3)[1])
f5=np.poly1d(scipy.signal.deconvolve(f3, f4)[1])
mas=[f2, f3, f4, f5]
for j in range(len(mas)):
    for i in range(len(mas[j])+1):
        mas[j][i]=-mas[j][i]
print('Загальна формула: f_(i+1)=-[f_(i-1) mod f_i]')
print(f0, ' = f0')
print(f1, ' = f1')
print(mas[0], ' = f2')
print(mas[1], ' = f3')
print(mas[2], ' = f4')
print(mas[3], ' = f5')

a=0.4070873392637155
b=4
c=-0.18181818181818182
d=-10
eps=0.00001

df=pd.DataFrame({'f': ['f0', 'f1', 'f2', 'f3', 'f4', 'f5'], f'f({a})':[f'{f0(a)}', f'{f1(a)}', f'{f2(a)}', f'{f3(a)}', f'{f4(a)}', f'{f5(a)}'], f'f({b})':[f'{f0(b)}', f'{f1(b)}', f'{f2(b)}', f'{f3(b)}', f'{f4(b)}', f'{f5(b)}'], f'f({c})':[f'{f0(c)}', f'{f1(c)}', f'{f2(c)}', f'{f3(c)}', f'{f4(c)}', f'{f5(c)}'], f'f({d})':[f'{f0(d)}', f'{f1(d)}', f'{f2(d)}', f'{f3(d)}', f'{f4(d)}', f'{f5(d)}']})
print(df)
print(f'Бачимо, що у точках {c} і {d} подслідовності змінюють знак два рази, отже на цьому проміжку рівняння не має коренів, але у точці {a} послідовність змінює знак 2 рази, а у точці {b} - 1 раз. Отже на відрізку [{a}, {b}] буде лише один дійсний корінь(що також видно з графіку функції).')

def bisection(x, y):
    z=(x+y)/2
    i=0
    while abs(x-y)>eps:
        i=i+1
        z=(x+y)/2
        if func(x)*func(z)<=0:
            y=z
        else: x=z
        print(f'{i}) a = {x}, b = {y}')
    return (x+y)/2

print('Метод бісекції:')
x_bisection=bisection(a, b)
print(f'Корінь рівняння з проміжку [{a}, {b}] дорівнює: ', x_bisection)

def horda(x, y):
    i=0
    z=(x*func(y)-y*func(x))/(func(y)-func(x))
    while abs(func(z))>eps:
        i=i+1
        z=(x*func(y)-y*func(x))/(func(y)-func(x))
        if func(x)*func(z)<=0:
            y=z
        else: x=z 
        print(f'{i}) a = {x}, b = {y}')
    return (x*func(y)-y*func(x))/(func(y)-func(x))

print('Метод хорд:')
x_horda=horda(a, b)
print(f'Корінь рівняння з проміжку [{a}, {b}] дорівнює: ', x_horda) 

d_func=lambda x: 5*x**4+4*x**3-6*x**2-18*x-3
def newton(x, y):
    i=0
    x0=(x+y)/2
    x1=x0-func(x0)/d_func(x0)
    while(abs(func(x1))>eps):
        i=i+1
        print(f'{i}) a = {x0}, b = {x1}')
        x0=(x1+x0)/2
        x1=x0-func(x0)/d_func(x0)
    return x1

print('Метод Ньютона:')
x_newton=newton(a, b)
print(f'Корінь рівняння з проміжку [{a}, {b}] дорівнює: ', x_newton)  


res_b=func(x_bisection)
res_h=func(x_horda)
res_n=func(x_newton)
print(f'Бісекція: f({x_bisection}) = {res_b}\nХорди: f({x_horda}) = {res_h}\nНьютон: f({x_newton}) = {res_n}')
if min([res_b, res_h, res_n], key=abs)==res_b:
    print(f'Найбільш точне значення отримано методом бісекції.')
elif min([res_b, res_h, res_n], key=abs)==res_h:
    print(f'Найбільш точне значення отримано методом хорд.')
else:
    print(f'Найбільш точне значення отримано методом Ньютона.')
