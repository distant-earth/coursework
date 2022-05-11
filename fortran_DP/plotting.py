import numpy as np
from matplotlib import pyplot as plt

t1 = []
a1 = []
ecc1 = []
omega1 = []
f = []
mass1 = []
with open('f_RESULT') as table:
    for line in table:
        t1.append(float(line.split()[1]))
        a1.append(float(line.split()[2]))
        ecc1.append(float(line.split()[3]))
        omega1.append(float(line.split()[4]))
        f.append(float(line.split()[5]))
        mass1.append(float(line.split()[6]))
        
t2 = []
a2 = []
ecc2 = []
omega2 = []
mass2 = []
with open('E_RESULT') as table:
    for line in table:
        t2.append(float(line.split()[1]))
        a2.append(float(line.split()[2]))
        ecc2.append(float(line.split()[3]))
        omega2.append(float(line.split()[4]))
        mass2.append(float(line.split()[6]))

step = 100    
# Шаг, с которым из выведенных в файл точек берутся точки для графика.

plt.title("Большая полуось")
plt.xlabel('t [лет]')
plt.ylabel('a [а.е.]')
plt.plot(t1[::step], a1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], a2[::step], 'b--', label = 'Через E')

# Для вывода теоретической кривой:
#def func(x):
#	y = [0.0747175 / (0.6175 - 1e-9 * i) for i in x] 
#	return y
#plt.plot(t1[::step], func(t1[::step]), 'k-', label = r'$a(t) = \frac{M_0 \cdot a_0}{M(t)}$')

# Для изменения формата чисел вдоль оси ординат:
#current_values = plt.gca().get_yticks()
#plt.gca().set_yticks(current_values)
#plt.gca().set_yticklabels(['{:.11f}'.format(x) for x in current_values])

plt.legend(loc = 'upper left', framealpha = 1)
plt.tight_layout()

plt.figure()
plt.title("Эксцентриситет")
plt.xlabel('t [лет]')
plt.ylabel('e')
plt.gca().set_ylim([0.70, 0.72])
plt.plot(t1[::step], ecc1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], ecc2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1)
plt.tight_layout()

plt.figure()
plt.title("Аргумент перицентра")
plt.xlabel('t [лет]')
plt.ylabel(r'$\omega$ [рад]')
plt.gca().set_ylim([4.12, 4.18])
plt.plot(t1[::step], omega1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], omega2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1)
plt.tight_layout()

g1 = np.array(a1) * (1 - np.array(ecc1))
g2 = np.array(a2) * (1 - np.array(ecc2))

plt.figure()
plt.title("Перицентрич. расстояние q = a(1 - e)")
plt.xlabel('t [лет]')
plt.ylabel('q [а.е.]')
plt.plot(t1[::step], g1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], g2[::step], 'b--', label = 'Через E')
#current_values = plt.gca().get_yticks()
#plt.gca().set_yticks(current_values)
#plt.gca().set_yticklabels(['{:.11f}'.format(x) for x in current_values])
plt.legend(loc = 'upper left', framealpha = 1)
plt.tight_layout()


ma1 = np.array(mass1) * np.array(a1)
ma2 = np.array(mass2) * np.array(a2)

plt.figure()
plt.title(r"Ma $\approx$ Ma(1 - e$^2$) = const")
plt.xlabel('t [лет]')
plt.ylabel('M$\cdot$a [$M_{Солнца}\cdot$ a.e.]')
plt.gca().set_ylim([0.07, 0.08])
plt.plot(t1[::step], ma1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], ma2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1)
plt.tight_layout()

plt.figure()
plt.title("Изменение массы")
plt.xlabel('t [лет]')
plt.ylabel('M [$M_{Солнца}$]')
plt.gca().set_ylim([0.58, 0.63])
plt.plot(t1[::step], mass1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], mass2[::step], 'b--', label = 'Через E')
#current_values = plt.gca().get_yticks()
#plt.gca().set_yticks(current_values)
#plt.gca().set_yticklabels(['{:.10f}'.format(x) for x in current_values])
plt.legend(loc = 'lower left', framealpha = 1)
plt.tight_layout()

plt.show()
