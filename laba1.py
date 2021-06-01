import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

func = lambda x: x**5+x**4-2*x**3-9*x**2-3*x-2
x=np.linspace(-10, 10, 100)
plt.grid()
plt.axis([-10, 10, -30, 30])
plt.plot(x, func(x))
plt.show()

print('Метод Лагранжа:')
a_n=[-2, -3, -9, -2, 1, 1]; n=len(a_n)
print(np.poly1d(a_n[::-1]), ' = P(x)')
def Lagrange(f, z):
    k=0
    F=[]; Ф=[]
    for i in reversed(f):
        if i>0:
            F.insert(0, i)
            k=k+1
        if i<0: break      
    for i in reversed(range(len(f)-k)):
        if f[i]<0:
            F.insert(0, f[i])
        else: F.insert(0, 0)
    for i in reversed(f):
        if i in F:
            Ф.insert(0, 0)
        else: Ф.insert(0, i)
    return np.poly1d(F[::-1], variable=z), np.poly1d(Ф[::-1], variable=z)
F=Lagrange(a_n, 'x')[0]
Ф=Lagrange(a_n, 'x')[1]
print(F,' = F(x)\n', Ф, ' = Ф(x)')
def find_alpha(f):
    alpha=1
    while f(alpha)<=0:
        alpha=alpha+1
    return alpha
alpha=find_alpha(F)
print(f'\nα={alpha}>0: F(α) = {F(alpha)} > 0, отже ∀ x* =< {alpha}')


print('Для нижньої межі: зробимо заміну x=1/y і отримаємо наступне рівняння:')
print(np.poly1d(a_n, variable='y'))
def cheking(mas, z):
    if mas[0]<0:
        res=[]
        for i in range(len(mas)):
            res.append((-1*mas[i]))
        print('Домножимо на -1:\n', np.poly1d(res, variable=z))
    else: res=mas
    return res
a_nn=cheking(a_n, 'y')
F=Lagrange(a_nn[::-1], 'y')[0]
Ф=Lagrange(a_nn[::-1], 'y')[1]
print(F,' = F(x)\n', Ф, ' = Ф(x)')
alpha1=1/find_alpha(F)
print(f'\nα={1/alpha1}>0: F(α) = {F(1/alpha1)} > 0, отже ∀ x* >= {alpha1}')

print("Уточнимо нижню межу від'ємних коренів: зробимо заміну x=-x і отримаємо наступне рівняння:")
def sign(mas):
    res=[]
    for i in mas: res.append(i)
    for i in range(len(res)):
        if i%2==1: res[i]=-1*res[i]
    return res
print(np.poly1d(sign(a_n)[::-1]))
a_nn=cheking(sign(a_n)[::-1], 'x')
F=Lagrange(a_nn[::-1], 'x')[0]
Ф=Lagrange(a_nn[::-1], 'x')[1]
print(F,' = F(x)\n', Ф, ' = Ф(x)')
alpha2=-find_alpha(F)
print(f'\nα={-alpha2}>0: F(α) = {F(-alpha2)} > 0, отже ∀ x- >= {alpha2}')

print("Уточнимо верхню межу від'ємних коренів: зробимо заміну x=-1/y і отримаємо наступне рівняння:")
print(np.poly1d(sign(a_n[::-1])[::-1], variable='y'))
a_nn=cheking(sign(a_n[::-1])[::-1], 'y')
F=Lagrange(a_nn[::-1], 'y')[0]
Ф=Lagrange(a_nn[::-1], 'y')[1]
print(F,' = F(x)\n', Ф, ' = Ф(x)')
alpha3=-1/find_alpha(F)
print(f'\nα={-1/alpha3}>0: F(α) = {F(-1/alpha3)} > 0, отже ∀ x- <= {alpha3}')

print('Теорема про кільце:')
A=abs(max(a_n[1:], key=abs))
B=abs(max(a_n[:(len(a_n)-1)], key=abs))
print(f'A = max(|a_i|) (i=0, 1, ..., 4) = {A}\nB = max(|a_i|) (i=1, 2, ..., 5) = {B}')
left=abs(a_n[0])/(B+abs(a_n[0])); rigth=(abs(a_n[len(a_n)-1])+A)/abs(a_n[len(a_n)-1])
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

R1=1+math.pow((B/a_n[len(a_n)-1]), 1/(len(a_n)-1-find(a_n)))
print(f'R = {R1} - верхня межа додатніх коренів.')

print('Для знаходження нижньої межі додатніх коренів зробимо заміну x=1/y і отримаємо наступне рівняння:')
print(np.poly1d(a_n, variable='y'))
a_nn=cheking(a_n, 'y')
R=1+math.pow((B/a_nn[::-1][len(a_nn)-1]), 1/(len(a_nn)-1-find(a_nn[::-1])))
R2=1/R
print(f'R = {R2} - нижня межа додатніх коренів.')

print("Для знаходження нижньої межі від'ємних коренів зробимо заміну x=-x і отримаємо наступне рівняння:")
print(np.poly1d(sign(a_n)[::-1]))
a_nn=cheking(sign(a_n)[::-1], 'x')
R=1+math.pow((B/a_nn[::-1][len(a_nn)-1]), 1/(len(a_nn)-1-find(a_nn[::-1])))
R3=-1*R
print(f"R = {R3} - нижня межа від'ємних коренів.")

print("Для знаходження верхньої межі від'ємних коренів зробимо заміну x=-1/y і отримаємо наступне рівняння:")
print(np.poly1d(sign(a_n[::-1])[::-1], variable='y'))
a_nn=cheking(sign(a_n[::-1])[::-1], 'y')
R=1+math.pow((B/a_nn[::-1][len(a_nn)-1]), 1/(len(a_nn)-1-find(a_nn[::-1])))
R4=-1/R
print(f"R = {R4} - верхня межа від'ємних коренів.")


print('Теорема Гюа про наявність комплесних коренів:')
a_n.reverse()
for i in range(1, len(a_n)-1):
    if a_n[i]**2<a_n[i-1]*a_n[i+1]:
        print(f'Існує таке k, що (a_k)^2<a_(k-1)*a_(k+1), k = {i}: {a_n[i]}^2={pow(a_n[i], 2)}<{a_n[i-1]}*{a_n[i+1]}={a_n[i-1]*a_n[i+1]}, отже рівняння має комплексні корені.')
        break
    if i==len(a_n)-2: 
        print('Рівняння не має комплексних коренів.')    
       


print('Теорема Штурма:')
f=np.poly1d(a_n)
print(f)
f0=f
f1=np.poly1d([5, 4, -6, -18, -3])
f2=np.polymul(np.polydiv(f0, f1)[1], -1)
f3=np.polymul(np.polydiv(f1, f2)[1], -1)
f4=np.polymul(np.polydiv(f2, f3)[1], -1)
f5=np.polymul(np.polydiv(f3, f4)[1], -1)
mas=[f2, f3, f4, f5]
print('Загальна формула: f_(i+1)=-[f_(i-1) mod f_i]')
print(f0, ' = f0')
print(f1, ' = f1')
print(mas[0], ' = f2')
print(mas[1], ' = f3')
print(mas[2], ' = f4')
print(mas[3], ' = f5')

a=max(left, R2, alpha1)
b=min(alpha, rigth, R1)
c=-min(abs(alpha3), abs(R4))
d=-min(abs(R3), abs(alpha2))
points=[d, c, a, b]
eps=0.00001


df=pd.DataFrame({'f': ['f0', 'f1', 'f2', 'f3', 'f4', 'f5'], 
                 f'f({d})':[f'{f0(d)}', f'{f1(d)}', f'{f2(d)}', f'{f3(d)}', f'{f4(d)}', f'{f5(d)}'], 
                 f'f({c})':[f'{f0(c)}', f'{f1(c)}', f'{f2(c)}', f'{f3(c)}', f'{f4(c)}', f'{f5(c)}'], 
                 f'f({a})':[f'{f0(a)}', f'{f1(a)}', f'{f2(a)}', f'{f3(a)}', f'{f4(a)}', f'{f5(a)}'], 
                 f'f({b})':[f'{f0(b)}', f'{f1(b)}', f'{f2(b)}', f'{f3(b)}', f'{f4(b)}', f'{f5(b)}']})
df.set_index('f')

def count(num):
    k=0
    for i in range(len(df[f'f({num})'])-1):
        if (float(df[f'f({num})'][i])>0 and float(df[f'f({num})'][i+1])<0) or (float(df[f'f({num})'][i])<0 and float(df[f'f({num})'][i+1])>0):
            k=k+1
    return k

for i in points:
    print(f'ККЗ у точці {i}: ', count(i))
for i in range(len(points)-1):
    print(f'Кількість коренів на проміжку [{points[i]}, {points[i+1]}] = ', count(points[i])-count(points[i+1]))
    if count(points[i])-count(points[i+1])>0:
        a=points[i]
        b=points[i+1]
print(f'a = {a}\nb = {b}')


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

x=np.linspace(1, 3, 100) 
plt.figure(figsize=(3, 3))
plt.grid()
plt.axis([1, 3, -25, 50])
plt.plot(x, func(x))

def horda(x, y):
    i=0
    z=(x*func(y)-y*func(x))/(func(y)-func(x))
    while abs(func(z))>eps:
        i=i+1
        z=(x*func(y)-y*func(x))/(func(y)-func(x))
        if func(a)>0: y=z
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
