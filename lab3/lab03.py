import math
import matplotlib.pyplot as plt

x0 = 0.01
v0 = 0
dt0 = 1
S = 0.75
p = 2
t_max = 40
alpha = 5
f = lambda v : v
g = lambda x, v : alpha*(1-x**2)*v-x
delta = 1e-10


def timeControl(method, TOL):
    t = [0]
    dt = [dt0]
    x_n = [x0]
    v_n = [v0]

    while True:
        x1, v1 = method(x_n[-1], v_n[-1], dt[-1])
        x2, v2 = method(x1, v1, dt[-1])

        x, v = method(x_n[-1], v_n[-1], 2*dt[-1])

        Ex = (x2 - x)/(2**p-1)
        Ev = (v2 - v)/(2**p-1)

        if(max(abs(Ex), abs(Ev))<TOL):
            t.append(t[-1] + 2*dt[-1])
            x_n.append(x2)
            v_n.append(v2)
            dt.append(dt[-1])
        
        dt[-1] *= (S*TOL/(max(abs(Ex),abs(Ev))))**(1/(p+1))

        if t[-1] >= t_max:
            break

    return t, x_n, v_n, dt


# # # METODA TRAPEZOW # # #

def methodTrapez(x_n, v_n, dt):
    x = x_n
    v = v_n
    a11 = 1
    a12 = -dt/2

    while True:

        a21 = -dt/2*(-2*alpha*x*v-1)
        a22 = 1 - dt/2*alpha*(1-x**2)

        F = x - x_n - dt/2*(f(v_n)+f(v))
        G = v - v_n - dt/2*(g(x_n, v_n)+g(x, v))

        dx = (-F*a22+G*a12)/(a11*a22-a12*a21)
        dv = (-G*a11+F*a21)/(a11*a22-a12*a21)

        x += dx 
        v += dv

        if abs(dx)>=delta and abs(dv)>=delta:
            break

    return x, v


# # # METODA RK2 # # #

def methodRK2(x_n, v_n, dt):
    k1x = f(v_n)
    k1v = g(x_n, v_n)

    k2x = f(v_n+dt*k1v)
    k2v = g(x_n+dt*k1x, v_n+dt*k1v)

    x = x_n + dt/2*(k1x+k2x)
    v = v_n + dt/2*(k1v+k2v)

    return x, v


# # # RYSOWANIE SUBPLOTU # # #

def subPlot(axs, x1, y1, x2, y2, str, x_min, x_max, y_min, y_max, x_name, y_name, num):
    if num == 3:
        axs.plot(x1, y1, 'ro-', x2, y2, 'b*-')
    else:
        axs.plot(x1, y1, 'r-', x2, y2, 'b-')
    axs.set_title(str)
    axs.set_xlim([x_min, x_max])
    axs.set_ylim([y_min, y_max])
    axs.set_xlabel(x_name)
    axs.set_ylabel(y_name)
    axs.legend(['$TOL=10^{-2}$', '$TOL=10^{-5}$'])


# # # PYSOWANIE PLOTU # # #

def drawPlot(t1, t2, x1, x2, v1, v2, dt1, dt2, str):
    fig, axs = plt.subplots(2, 2)

    subPlot(axs[0, 0], t1, x1, t2, x2, str, 0, 40, -2.5, 2.5, 't', 'x(t)', 1)
    subPlot(axs[0, 1], t1, v1, t2, v2, str, 0, 40, -8, 8, 't', 'v(t)', 2)
    subPlot(axs[1, 0], t1, dt1, t2, dt2, str, 0, 40, 0, 1, 't', 'dt(t)', 3)
    subPlot(axs[1, 1], x1, v1, x2, v2, str, -2.5, 2.5, -8, 8, 'x', 'v', 4)



# # # MAIN # # #

t1, x1, v1, dt1 = timeControl(methodTrapez, 1e-2)
t2, x2, v2, dt2 = timeControl(methodTrapez, 1e-5)
drawPlot(t1, t2, x1, x2, v1, v2, dt1, dt2, "Metoda Trapezow")

t1, x1, v1, dt1 = timeControl(methodRK2, 1e-2)
t2, x2, v2, dt2 = timeControl(methodRK2, 1e-5)
drawPlot(t1, t2, x1, x2, v1, v2, dt1, dt2, "Metoda RK2")

plt.show()
