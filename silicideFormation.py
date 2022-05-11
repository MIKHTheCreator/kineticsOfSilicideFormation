import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

T_1 = 400
T_2 = 500
T_3 = 600
k = 1.38E-23
h_me = 100E-9
po = 8.9
A = 58.7
dt = 0.5
dT = 50

def t_func(t):
    return 10 * t + 5 * math.pow(t, 2)


def funct(x, E_a_i):
    return np.exp((-E_a_i * 1.602E-19) / (k * 5 * x))


def x_j_i(U_i, E_a_i, T):
    intg, err = integrate.quad(funct, 0, np.inf, args=(E_a_i))
    return U_i * intg


def main(U_i, E_a_i, T):
    return U_i * funct(t, E_a_i, T)


def v_t(integral):
    val = math.sqrt(integral)
    return val


def conditions(left_side):
    condition = (po / A) * h_me
    return left_side <= condition


def border_conditions(x_j_1, x_j_2, x_j_3):
    return (po / A) * x_j_1 + (po / A) * (x_j_2 - x_j_1) + (po / A) * (x_j_3 - x_j_2)


def integration(U_i, E_a_i, T, t):
    function = U_i * np.exp((-E_a_i * 1.602E-19) / (k * T)) * t
    return function


def formation(U_i, E_a_i, T, time):
    return U_i  * math.pow(math.e, (-E_a_i * 1.602E-19 / (1.38E-23 * T))) * time


t_list = []
h_list = []
h2_list = []
h3_list = []
t1 = 10

# Расчет точек для графика:
for t in np.arange(0, t1, dt):
    x_j_1_sq = integration(0.16, 1.65, T_1, t)
    x_j_2_sq = integration(0.85, 1.55, T_1, t)
    x_j_3_sq = integration(6, 1.5, T_1, t)

    x_j_1 = v_t(math.fabs(x_j_1_sq))
    x_j_2 = v_t(math.fabs(x_j_2_sq))
    x_j_3 = v_t(math.fabs(x_j_3_sq))

    x_j_1_sq_2 = integration(0.16, 1.65, T_2, t)
    x_j_2_sq_2 = integration(0.85, 1.55, T_2, t)
    x_j_3_sq_2 = integration(6, 1.5, T_2, t)

    x_j_1_2 = v_t(math.fabs(x_j_1_sq_2))
    x_j_2_2 = v_t(math.fabs(x_j_2_sq_2))
    x_j_3_2 = v_t(math.fabs(x_j_3_sq_2))

    x_j_1_sq_3 = integration(0.16, 1.65, T_3, t)
    x_j_2_sq_3 = integration(0.85, 1.55, T_3, t)
    x_j_3_sq_3 = integration(6, 1.5, T_3, t)

    x_j_1_3 = v_t(math.fabs(x_j_1_sq_3))
    x_j_2_3 = v_t(math.fabs(x_j_2_sq_3))
    x_j_3_3 = v_t(math.fabs(x_j_3_sq_3))

    if (conditions(border_conditions(x_j_1_sq, x_j_2_sq, x_j_3_sq))):
        t_list.append(t)
        h_list.append(x_j_1 * 1000)
        h2_list.append(x_j_1_2 * 100)
        h3_list.append(x_j_1_3)

# Объявим форму для рисования и оси:
fig, axes = plt.subplots()

plt.plot(t_list, h_list, color='blue', label='Ni2Si(400C)')
plt.plot(t_list, h2_list, color='red', label='Ni2Si(500C)')
plt.plot(t_list, h3_list, color='orange', label='Ni2Si(600C)')

axes.set_xlim(0, 10)
axes.set_ylim(0, 700E-9)
plt.title('Kinetics of Silicide Formation')
plt.xlabel('t, hours')
plt.ylabel('х, nm')
plt.legend(loc=5)

# Добавление дополнительной ссетки:
axes.grid(which='major', color='#666666')
axes.minorticks_on()
axes.grid(which='minor', color='gray', linestyle=':')

# # Выведем график на экран:
plt.show()
