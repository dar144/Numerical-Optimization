#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt
import numpy as np

lmbd = -1
y_0 = 1
delta_t = [0.01, 0.1, 1.0]
t_min = 0.0
t_max = 5.0
colors = ['bo', 'ro', 'go']
diff_colors = ['b', 'r', 'g']

# dt = 0.0001
# R = 100
# L = 0.1
# C = 0.001
# omega_0 = 1/math.sqrt(L*C)
# T_0 = 2*math.pi/omega_0
# time_min = 0
# time_max = 3*T_0
# omegas = [0.5*omega_0, 0.8*omega_0, 1.0*omega_0, 1.2*omega_0]

# Q_0 = 0
# I_0 = 0
# V_0 = 0

def n_t(step):
    n = math.floor((t_max-t_min)/step)
    t = [t_min + step*i for i in range(n+1)]
    return n, t

# def RK4_2rzed(omega_v):
#     n, t = n_t(0.05)
#     Q = [Q_0]
#     I = [I_0]
#     for i in range(n):
#         V_n = 10*math.sin(omega_v*t[i])
#         kQ_1 = I[i]
#         kI_1 = V_n/L - Q[i]/L/C - R*I[i]/L

#         kQ_2 = I[i]+dt/2*kI_1
#         kI_2 = (10*math.sin(omega_v*t[i]+math.pi/2))/L - (Q[i]+dt/2*kQ_1)/L/C - R*(I[i]+dt/2*kI_1)/L

#         kQ_3 = I[i]+dt/2*kI_2
#         kI_3 = (10*math.sin(omega_v*t[i]+math.pi/2))/L - (Q[i]+dt/2*kQ_2)/L/C - R*(I[i]+dt/2*kI_2)/L

#         kQ_4 = I[i]+dt/2*kI_3
#         kI_4 = (10*math.sin(omega_v*t[i]+math.pi))/L - (Q[i]+dt*kQ_3)/L/C - R*(I[i]+dt*kI_3)/L

#         Q.append(Q[i]+dt/6*(kQ_1 + 2*kQ_2 + 2*kQ_3 + kQ_4))
#         I.append(I[i]+dt/6*(kI_1 + 2*kI_2 + 2*kI_3 + kI_4))
#     return Q, I, t

def Euler(step):
    n, t = n_t(step)
    y = [y_0]
    for i in range(n):
        y.append(y[i] + step*lmbd*y[i])
    return t, y

def RK2(step):
    n, t = n_t(step)
    y = [y_0]
    for i in range(n):
        k_1 = y[i]*lmbd
        k_2 = lmbd*(y[i]+step*k_1)
        y.append(y[i]+step/2*(k_1+k_2))
    return t, y

def RK4(step):
    n, t = n_t(step)
    y = [y_0]
    for i in range(n):
        k_1 = y[i]*lmbd
        k_2 = lmbd*(y[i]+step*k_1/2)
        k_3 = lmbd*(y[i]+step*k_2/2)
        k_4 = lmbd*(y[i]+step*k_3)
        y.append(y[i]+step/6*(k_1+2*k_2+2*k_3+k_4))
    return t, y


n_0, t_0 = n_t(delta_t[0])


# EULERA
fig, (ax1, ax2) = plt.subplots(1, 2)
for i in range(len(delta_t)):
    n = math.floor((t_max-t_min)/delta_t[i])
    t, y = Euler(delta_t[i])
    y_analitical = [math.exp(lmbd*t[i]) for i in range(n+1)]
    diff_y = [y[i] - y_analitical[i] for i in range(n+1)]
    ax1.plot(t, y, colors[i])
    ax2.plot(t, diff_y, diff_colors[i])
ax1.plot(t_0, [math.exp(lmbd*t_0[i]) for i in range(n_0+1)],'k')
ax1.set_xlabel('t')
ax2.set_xlabel('t')
ax1.set_ylabel('y(t)')
ax2.set_ylabel('y_num(t)-y_dok(t)')
ax1.legend(['dt = 0.01', 'dt = 0.1', 'dt = 1','analitycznie'])
ax1.set_title('z.1 - Metoda Eulera - rozwiazanie')
ax2.set_title('z.1 - Metoda Eulera - blat globalny')


# RK2
fig, (ax1, ax2) = plt.subplots(1, 2)
ax3 = fig.add_axes([.8, .65, .09, .2])
for i in range(len(delta_t)):
    n = math.floor((t_max-t_min)/delta_t[i])
    t, y = RK2(delta_t[i])
    y_analitical = [math.exp(lmbd*t[i]) for i in range(n+1)]
    diff_y = [y[i] - y_analitical[i] for i in range(n+1)]
    ax1.plot(t, y, colors[i])
    ax2.plot(t, diff_y, diff_colors[i])
    ax3.plot(t, diff_y, diff_colors[i])
ax1.plot(t_0, [math.exp(lmbd*t_0[i]) for i in range(n_0+1)],'k')
ax3.axis([0, 5, 0, 0.001])
ax1.set_xlabel('t')
ax2.set_xlabel('t')
ax1.set_ylabel('y(t)')
ax2.set_ylabel('y_num(t)-y_dok(t)')
ax1.legend(['dt = 0.01', 'dt = 0.1', 'dt = 1','analitycznie'])
ax1.set_title('z.1 - Metoda RK2 - rozwiazanie')
ax2.set_title('z.1 - Metoda RK2 - blat globalny')

# RK4
fig, (ax1, ax2) = plt.subplots(1, 2)
ax3 = fig.add_axes([.8, .65, .09, .2])
for i in range(len(delta_t)):
    n = math.floor((t_max-t_min)/delta_t[i])
    t, y = RK4(delta_t[i])
    y_analitical = [math.exp(lmbd*t[i]) for i in range(n+1)]
    diff_y = [y[i] - y_analitical[i] for i in range(n+1)]
    ax1.plot(t, y, colors[i])
    ax2.plot(t, diff_y, diff_colors[i])
    ax3.plot(t, diff_y, diff_colors[i])
ax1.plot(t_0, [math.exp(lmbd*t_0[i]) for i in range(n_0+1)],'k')
ax3.axis([0, 5, 0, 0.0000005])
ax1.set_xlabel('t')
ax2.set_xlabel('t')
ax1.set_ylabel('y(t)')
ax2.set_ylabel('y_num(t)-y_dok(t)')
ax1.legend(['dt = 0.01', 'dt = 0.1', 'dt = 1','analitycznie'])
ax1.set_title('z.1 - Metoda RK2 - rozwiazanie')
ax2.set_title('z.1 - Metoda RK2 - blat globalny')



# # 2 RZED
# fig, (ax1, ax2) = plt.subplots(1, 2)
# for omega in omegas:
#     Q, I, t = RK4_2rzed(omega)
#     ax1.plot(t, Q)
#     ax2.plot(t, I)

plt.show()
