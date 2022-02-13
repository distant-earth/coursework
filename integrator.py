import numpy as np
from scipy.integrate import solve_ivp as solver
from astropy.constants import G
from astropy import units as u
from tabulate import tabulate
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

Gau = G.to('au3 / (M_jupiter yr2)')                                 # грав.постоянная в нужной размерности
#===========================================================================================================
# Начальные условия:
t0 = 0                                                              # начальное время [в годах]
m1 = (0.798 * u.solMass).to(u.M_jupiter).value                      # масса первого тела [в массах Юпитера]
m2 = 30.0                                                           # масса второго тела [в массах Юпитера]
m0 = m1 + m2        
a0 = 0.006                                                          # большая полуось [в а.е.]
ecc0 = 0.1                                                          # эксцентриситет
omega0 = 0.2                                                        # аргумент перицентра [в рад.]
f0 = 2.7                                                            # истинная аномалия [в рад.]
E0 = np.arccos((ecc0 + np.cos(f0)) / (1 + ecc0 * np.cos(f0)))       # эксцентрическая аномалия [рад]
# Конечное время [в годах]:
time = 10000
# Коэффициенты в законе Эддингтона-Джинса потери массы:
alpha = 1e-8                # должен быть положительным, близким к 0
n = 2                       # в пределах от 1.4 до 4.4
#===========================================================================================================

# Выражение для суммарной массы в момент t (если t0 = 0), полученное из модели Джинса.
def mass(t):
    return (alpha * (n - 1) * (t - t0) + m0**(1-n))**(1/(1-n))
# Для производной:
def dmass(t):
    return -alpha * mass(t)**n

# Правые части уравнений через истинную аномалию:
def rp_f(t, y):
    a = y[0]
    ecc = y[1]
    omega = y[2]
    f = y[3]
    res = np.empty(4)
    res[0] = a * ((1 + ecc**2 + 2 * ecc * np.cos(f)) / (ecc**2 - 1)) * (dmass(t) / mass(t))
    res[1] = - (ecc + np.cos(f)) * (dmass(t) / mass(t))
    res[2] = - (np.sin(f) / ecc) * (dmass(t) / mass(t))
    param = np.sqrt(Gau.value * mass(t)) / a**(3/2)
    res[3] = - res[2] + param * (1 + ecc * np.cos(f))**2 / (1 - ecc**2)**(3/2)
    return res

# Правые части уравнений через эксцентрическую аномалию:
def rp_E(t, y):
    a = y[0]
    ecc = y[1]
    omega = y[2]
    E = y[3]
    res = np.empty(4)
    res[0] = - a * ((1 + ecc * np.cos(E)) / (1 - ecc * np.cos(E))) * (dmass(t) / mass(t))
    res[1] = ((1 - ecc**2) * np.cos(E) / (ecc * np.cos(E) - 1)) * (dmass(t) / mass(t))
    res[2] = (np.sqrt(1 - ecc**2) * np.sin(E) / (ecc * (ecc * np.cos(E) - 1))) * (dmass(t) / mass(t))
    param = np.sqrt(Gau.value * mass(t)) / a**(3/2)
    res[3] = - res[2] / np.sqrt(1 - ecc**2) + param / (1 - ecc * np.cos(E))
    return res

#===========================================================================================================

tol = [1e-4, 1e-4, 1e-6, 1e-10]                        # абсолютные погрешности a, e, omega, f или E

init_f = np.array([a0, ecc0, omega0, f0])
sol_f = solver(rp_f, (t0, time), init_f, atol=tol, method='RK45')        # решение уравнений через f
a_arr1 = sol_f.y[0, 0:]
ecc_arr1 = sol_f.y[1, 0:]
omega_arr1 = sol_f.y[2, 0:] % (2*np.pi)
f_arr1 = sol_f.y[3, 0:] % (2*np.pi)

init_E = np.array([a0, ecc0, omega0, E0])
sol_E = solver(rp_E, (t0, time), init_E, atol=tol, method='RK45')        # решение уравнений через E
a_arr2 = sol_E.y[0, 0:]
ecc_arr2 = sol_E.y[1, 0:]
omega_arr2 = sol_E.y[2, 0:] % (2*np.pi)
E_arr2 = sol_E.y[3, 0:] % (2*np.pi)

print()
print('Статус решения системы через f:', sol_f.message)
print('Кол-во шагов:', sol_f.y.shape[1])
print('Статус решения системы через E:', sol_E.message)
print('Кол-во шагов:', sol_E.y.shape[1])
print()
print(tabulate([['Начальные данные:', t0, a0, ecc0, omega0 % (2*np.pi), f0 % (2*np.pi), E0 % (2*np.pi)], 
['Решение через f:', time, *sol_f.y[0:2,-1], *sol_f.y[2:,-1] % (2*np.pi)], 
['Решение через E:', time, *sol_E.y[0:2,-1], sol_E.y[2,-1] % (2*np.pi), None, sol_E.y[3,-1] % (2*np.pi)]], 
headers=['t [годы]', 'a [а.е.]', 'e', 'omega [рад.]', 'f [рад.]', 'E [рад.]']))

with PdfPages('result.pdf') as pdf:

    step = 10                           # шаг, с которым из результата берутся точки для графика

    plt.title("Большая полуось")
    plt.xlabel('t [годы]')
    plt.ylabel('a [а.е.]')
    plt.plot(sol_f.t[0::step], a_arr1[0::step], 'r.', label = 'Через f')
    plt.plot(sol_E.t[0::step], a_arr2[0::step], 'b.', label = 'Через E')
    plt.legend(loc = 'upper left', framealpha = 1, title = fr'$a_0 = {a0}$')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.title("Эксцентриситет")
    plt.xlabel('t [годы]')
    plt.ylabel('e')
    plt.gca().set_ylim([0, 1])
    plt.plot(sol_f.t[0::step], ecc_arr1[0::step], 'r.', label = 'Через f')
    plt.plot(sol_E.t[0::step], ecc_arr2[0::step], 'b.', label = 'Через E')
    plt.legend(loc = 'upper left', framealpha = 1, title = fr'$e_0 = {ecc0}$')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.title("Аргумент перицентра")
    plt.xlabel('t [годы]')
    plt.ylabel(r'$\omega$ [рад]')
    plt.gca().set_ylim([0, 2 * np.pi])
    plt.plot(sol_f.t[0::step], omega_arr1[0::step], 'r.', label = 'Через f')
    plt.plot(sol_E.t[0::step], omega_arr2[0::step], 'b.', label = 'Через E')
    plt.legend(loc = 'upper left', framealpha = 1, title = fr'$\omega_0 = {round(omega0, 3)}$')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    #plt.title("Аномалии")
    #plt.xlabel('t [годы]')
    #plt.ylabel('f, E [рад]')
    #plt.plot(sol_f.t[0::step], f_arr1[0::step], 'r.', label = 'f')
    #plt.plot(sol_E.t[0::step], E_arr2[0::step], 'b.', label = 'E')
    #plt.legend(loc = 'upper left', framealpha = 1, title = fr'$f_0 = {f0}, E_0 = {round(E0, 3)}$')
    #plt.tight_layout()
    #pdf.savefig()
    #plt.close()

    pdf.infodict()['Title'] = 'Параметры орбиты'