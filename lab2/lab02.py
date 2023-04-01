import math
import matplotlib.pyplot as plt

dt = 0.1
beta = 0.001
N = 500
gamma = 0.1
alpha = beta * N - gamma 
t_max = 100
t = [i * dt for i in range(int(t_max/dt))]
u_0 = 1
TOL = 10**(-6)
mu_max = 20

c = [1/2-math.sqrt(3)/6, 1/2+math.sqrt(3)/6]
b = [1/2, 1/2]
a = [[1/4, 1/4-math.sqrt(3)/6], [1/4+math.sqrt(3)/6, 1/4]]

# funckja rysująca wykresy dla metody trapezów
def drawPlot(t, u, v, str):
    plt.figure(str)
    plt.plot(t, u, 'r-', t, v, 'b-')
    plt.xlabel('t')
    plt.ylabel('u(t),v(t)')
    plt.legend(['u(t)', 'v(t)'])
    plt.xlim([0, 100])
    plt.ylim([0, 500])
    plt.title(str)



# # # METODA TRAPEZOW # # #

# funckje liczące u_(n+1)^(mu+1)
def nextPicard(u_i, u_n):
    return u_i+dt/2*((alpha*u_i-beta*(u_i**2))+(alpha*u_n-beta*(u_n**2)))

def nextNewton(u_i, u_n):
    return u_n-(u_n-nextPicard(u_i, u_n))/(1-dt/2*(alpha-2*beta*u_n))

# funckje liczące finalne wartosci u, v oraz mu
def resultsTrapez(nextVal, name):
    u = [u_0]
    v = [N-u_0]

    for i in range(len(t)-1):
        u_curr = u[i] 
        u_next = nextVal(u[i], u_curr) 
        mu = 1

        while abs(u_next-u_curr)>TOL and mu <=20:
            u_curr = u_next
            u_next = nextVal(u[i], u_curr)
            mu += 1

        u.append(u_next)
        v.append(N - u_next)

    drawPlot(t, u, v, name)


# # # METODA RK2 # # #

def nextF(U_first, U_second, u_n):
    return U_first - u_n - dt*(a[0][0]*(alpha*U_first-beta*(U_first**2))
        +a[0][1]*(alpha*U_second-beta*(U_second**2)))

def findFun(u):
    return (beta*N - gamma)*u - beta*(u**2)

def nextU(u_n, U_1, U_2):
    m_11 = 1 - dt*a[0][0]*(alpha-2*beta*U_1)
    m_12 = -dt*a[0][1]*(alpha-2*beta*U_2)
    m_21 = -dt*a[1][0]*(alpha-2*beta*U_1)
    m_22 = 1 - dt*a[1][1]*(alpha-2*beta*U_2)

    F_1 = nextF(U_1, U_2, u_n)
    F_2 = nextF(U_2, U_1, u_n)

    U_1 += (F_2*m_12 - F_1*m_22)/(m_11*m_22-m_12*m_21)
    U_2 += (F_1*m_21 - F_2*m_11)/(m_11*m_22-m_12*m_21)

    return U_1, U_2


def resultsRK(nextVal, name):
    u = [u_0]
    v = [N-u_0]

    for i in range(len(t)-1):
        u_curr = u[i] 
        U_1 = u[i]
        U_2 = u[i]
        U_1_new, U_2_new = nextVal(u_curr, U_1, U_2) 
        mu = 1

        while abs(U_1_new-U_1)>TOL and abs(U_2_new-U_2)>TOL and mu <=20:
            U_1 = U_1_new
            U_2 = U_2_new
            U_1_new, U_2_new = nextVal(u_curr, U_1, U_2)
            mu += 1

        u.append(u_curr + dt*(b[0]*findFun(U_1)+b[1]*findFun(U_2)))
        v.append(N - u[i+1])

    drawPlot(t, u, v, name)



# main
resultsTrapez(nextPicard, 'Metoda Picarda')
resultsTrapez(nextNewton, 'Metoda Newtona')
resultsRK(nextU, 'Metoda RK2')

plt.show()
