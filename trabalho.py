# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math as m

def exata(x,alpha,t,case):
    y = np.empty(len(x))
    if case == 1:
        for i in range(len(x)):
            y[i] = 0
            xAux = x[i] - alpha * t
            if xAux >=0 and xAux < 0.2:
                y[i] = m.exp(-m.log(50,m.e) * ((((xAux) - 0.15)/0.05) ** 2))
            elif xAux > 0.3 and xAux < 0.4:
                y[i] = 1
            elif xAux > 0.5 and xAux < 0.55:
                y[i] = 20 * (xAux) - 10
            elif xAux >= 0.55 and xAux < 0.66:
                y[i] = -20 * (xAux) + 12
            elif xAux > 0.7 and xAux < 0.8:
                y[i] = m.sqrt(1 - (((xAux) - 0.75) / 0.05) ** 2)

    elif case == 2:
        for i in range(len(x)):
            xAux = x[i] - alpha * t
            if xAux >= 0 and xAux <= 0.2:
                y[i] = 1
            elif xAux > 0.2 and xAux <= 0.4:
                y[i] = 4 * (xAux) - 0.6
            elif xAux > 0.4 and xAux <= 0.6:
                y[i] = -4 * (xAux) + 2.6
            elif xAux > 0.6 and xAux <= 0.8:
                y[i] = 1
            else:
                y[i] = 0

    return y

def fou(uu, Nx, Nt, h, dt, alpha):
    u = np.empty([Nx+1,Nt+1])
    u[:,:] = uu[:,:]

    for t in range(0,Nt):
        for i in range(1,Nx):
            u[i,t+1] = u[i,t] - (dt/h*alpha) * ( facef(u,i,t,alpha) - faceg(u,i,t,alpha))

    return u

def topus(uu, Nx, Nt, h, dt, alpha):
    u = np.empty([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(0,Nt):
        for i in range(1,Nx):
            ux = u[i,t] - (dt/h*alpha) * ( facef(u,i,t,alpha) - faceg(u,i,t,alpha))
            if ux >= 0 and ux <= 1:
                u[i,t+1] = (2*(ux**4) - 3*(ux**3) + 2*(ux))*(u[i + 1,t] - u[i - 1,t]) + u[i - 1,t]
            else:
                u[i,t+1] = ux

    return u

def facef(u,i,t,alpha): # i + 1/2
    if alpha > 0:
        return u[i,t]
    else:
        return u[i+1,t]

def faceg(u,i,t,alpha): # i - 1/2
    if alpha > 0:
        return u[i-1,t]
    else:
        return u[i,t]

############################################################

def condicaoInicial(Nx, Nt, h, x, case):
    u = np.zeros([Nx+1,Nt+1])
    if case == 1:
        for i in range(Nx+1):
            if x[i] < 0.2:
                u[i,0] = m.exp(-m.log(50,m.e) * (((x[i] - 0.15)/0.05) ** 2))
            elif x[i] > 0.3 and x[i] < 0.4:
                u[i,0] = 1
            elif x[i] > 0.5 and x[i] < 0.55:
                u[i,0] = 20 * x[i] - 10
            elif x[i] > 0.55 and x[i] < 0.66:
                u[i,0] = -20 * x[i] + 12
            elif x[i] > 0.7 and x[i] < 0.8:
                u[i,0] = m.sqrt(1 - ((x[i] - 0.75) / 0.05) ** 2)
            else:
                u[i,0] = 0
    elif case == 2:
        for i in range(len(x)):
            if x[i] >= 0 and x[i] <= 0.2:
                u[i,0] = 1
            elif x[i] > 0.2 and x[i] <= 0.4:
                u[i,0] = 4 * x[i] - 0.6
            elif x[i] > 0.4 and x[i] <= 0.6:
                u[i,0] = -4 * x[i] + 2.6
            elif x[i] > 0.6 and x[i] <= 0.8:
                u[i,0] = 1
            else:
                u[i,0] = 0

    return u


def condicoesContorno(u,Nx,Nt):
    for t in range(0, Nt+1):
        u[0,0] = 0.0
        u[Nx-1,0] = 0.0
    return u

#Entrada
case = 2 # case1 ou case2
alpha = 1.0 # coeficiente convectivo
if case == 1:
    T = 1 #Tempo final case1 = 1 | case2 = 0.25
    A = 0 #Limite inferior para x
    B = 2 #Limite superior para x
elif case == 2:
    T = 0.25 #Tempo final case1 = 1 | case2 = 0.25
    A = -1 #Limite inferior para x
    B = 1 #Limite superior para x


Nx = 400 #Quantidade de elementos
h = (B - A)/np.float(Nx) #Discretização no espaço
dt = 0.0025 #Discretização no tempo s δt = 2.5 × 10−4 e δt = 2.5 × 10-3
Nt = np.int(T/dt) #Quantidade de iterações no tempo
x = np.linspace(A,B,Nx+1) #Para plot das aproximações
tj = 0.1 #Instante de tempo desejado
t = int(tj/dt) #Índice correspondente ///

print(alpha*dt/h)
xe = np.linspace(A,B,1000) #Para plot da solução exata
u_e = exata(xe, alpha, tj, case)

u = condicaoInicial(Nx,Nt,h,x,case)
u = condicoesContorno(u,Nx,Nt)

u_fou = topus(u,Nx,Nt,h,dt,alpha)
#u_laxFri,teste = laxFriedrichs(u,Nx,Nt,h,dt)
#u_kt,teste2 = KT(u,Nx,Nt,h,dt,t)

#Plota o grafico
plt.title('Exata x Aproximada ('+ str(Nx) +  ' elementos , dt = ' + str(dt))
if case == 1:
    plt.xlim (A ,B)
#plt.ylim (0 ,1.1)
plt.grid()

plt.plot(xe, u_e, 'r-', label = 'Exata')
plt.plot(x, u_fou[:,t], 'bx', label = 'First Order UpWind')

plt.xlabel( 'x' )
plt.ylabel ( 'u' )
plt.legend(loc='best')
diretorio = ""
dtAux = str(dt/h)
dtAux = dtAux.replace('.','-')
nomefig = 'FOU_' + dtAux + '_case_' + str(case) + '.png'
plt.savefig(diretorio+nomefig, dpi=200)
plt.show()
