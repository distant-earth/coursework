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

step = 100                        # шаг, с которым из результата берутся точки для графика

plt.title("Большая полуось (теоретически)")
plt.xlabel('t [лет]')
plt.ylabel('a [а.е.]')
#plt.plot(t1[::step], a1[::step], 'r-', label = 'Через f')
#plt.plot(t2[::step], a2[::step], 'b--', label = 'Через E')
def func(x):
	y = [0.006 / (1 - 0.00001 * i) for i in x] 
	return y
plt.plot(t1[::step], func(t1[::step]), 'k-', label = 'a(t) = C / m(t)')
plt.legend(loc = 'upper left', framealpha = 1) #title = fr'$a_0 = {}$')
plt.tight_layout()

plt.figure()
plt.title("Эксцентриситет")
plt.xlabel('t [лет]')
plt.ylabel('e')
plt.gca().set_ylim([-0.1, 0.1])
plt.plot(t1[::step], ecc1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], ecc2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1) # title = fr'$e_0 = {ecc0}$')
plt.tight_layout()

plt.figure()
plt.title("Аргумент перицентра")
plt.xlabel('t [лет]')
plt.ylabel(r'$\omega$ [рад]')
plt.gca().set_ylim([0.1, 0.3])
plt.plot(t1[::step], omega1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], omega2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1) #title = fr'$\omega_0 = {round(omega0, 3)}$')
plt.tight_layout()

g1 = np.array(a1) * (1 - np.array(ecc1))
g2 = np.array(a2) * (1 - np.array(ecc2))

plt.figure()
plt.title("Перицентрич. расстояние g = a(1 - e)")
plt.xlabel('t [лет]')
plt.ylabel('g [а.е.]')
plt.plot(t1[::step], g1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], g2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1) #title = fr'$\omega_0 = {round(omega0, 3)}$')
plt.tight_layout()

ma1 = np.array(mass1) * np.array(a1)
ma2 = np.array(mass2) * np.array(a2)

plt.figure()
plt.title("Ma = const")
plt.xlabel('t [лет]')
plt.ylabel('M$\cdot$a [$M_{системы}\cdot$ a.e.]')
plt.gca().set_ylim([0.004, 0.008])
plt.plot(t1[::step], ma1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], ma2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'upper left', framealpha = 1) #title = fr'$\omega_0 = {round(omega0, 3)}$')
plt.tight_layout()

plt.figure()
plt.title("Изменение массы")
plt.xlabel('t [лет]')
plt.ylabel('M [$M_{системы}$]')
#plt.gca().set_ylim([0.9, 1.1])
plt.plot(t1[::step], mass1[::step], 'r-', label = 'Через f')
plt.plot(t2[::step], mass2[::step], 'b--', label = 'Через E')
plt.legend(loc = 'lower left', framealpha = 1) #title = fr'$\omega_0 = {round(omega0, 3)}$')
plt.tight_layout()

plt.show()
