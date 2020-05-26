from matplotlib.widgets import Slider

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

b = 1/2 # contacts per day of infected individual
k = 1/3 # fraction recovered per day
n = 200 # days

N = 7900000000 # total pop.
s_nau = 1
i_nau = 10 / N

def dsdt(s, i):
    return -b*s*i

def didt(s, i):
    return b*s*i - k*i

def drdt(i):
    return k*i

def RK4(n, s_nau, i_nau, r_nau = 0, dt = 1):

    s_arr = [s_nau] * n
    i_arr = [i_nau] * n
    r_arr = [r_nau] * n

    for i in range(n - 1):

        si = s_arr[i]
        ii = i_arr[i]
        ri = r_arr[i]

        sk1 = dsdt(si, ii)
        ik1 = didt(si, ii)
        rk1 = drdt(ii)

        sk2 = dsdt(si + 0.5 * dt * sk1, ii + 0.5 * ik1)
        ik2 = didt(si + 0.5 * dt * sk1, ii + 0.5 * ik1)
        rk2 = drdt(ii + 0.5 * dt * ik1)

        sk3 = dsdt(si + 0.5 * dt * sk2, ii + 0.5 * ik2)
        ik3 = didt(si + 0.5 * dt * sk2, ii + 0.5 * ik2)
        rk3 = drdt(ii + 0.5 * dt * ik2)

        sk4 = dsdt(si + dt * sk3, ii + dt * sk3)
        ik4 = didt(si + dt * sk3, ii + dt * ik3)
        rk4 = drdt(ii + dt * ik3)

        s_arr[i + 1] = si + dt / 6 * (sk1 + 2 * sk2 + 2 * sk3 + sk4)
        i_arr[i + 1] = ii + dt / 6 * (ik1 + 2 * ik2 + 2 * ik3 + ik4)
        r_arr[i + 1] = ri + dt / 6 * (rk1 + 2 * rk2 + 2 * rk3 + rk4)
    
    return [s_arr, i_arr, r_arr]

[s, i, r] = RK4(n, s_nau, i_nau)

plt.plot(s)
plt.plot(i)
plt.plot(r)
plt.legend(('Susceptible', 'Infected', 'Recovered'), loc = 'upper right', fontsize = 11)
plt.title('SIR Model')
plt.xlabel('Days', fontsize = 11)
plt.ylabel('Fraction of pop. infected', fontsize = 11)
plt.show()