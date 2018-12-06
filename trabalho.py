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
            u[i,t+1] = u[i,t] - (dt/h*alpha) * ( facef(u,i,t,alpha,Nx,0) - faceg(u,i,t,alpha,Nx,0))

    return u

def topus(uu, Nx, Nt, h, dt, alpha):
    u = np.empty([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(0,Nt):
        for i in range(1,Nx):
            u[i,t+1] = u[i,t] - (dt/h*alpha) * ( facef(u,i,t,alpha,Nx,1) - faceg(u,i,t,alpha,Nx,1))

    return u

def facef(u,i,t,alpha,Nx,type): # i + 1/2
    if type == 0: # FOU
        if alpha > 0:
            return u[i,t]
        else:
            return u[i + 1,t]
    elif type == 1: # TOPUS
        uD = u[i + 1,t]
        uU = u[i,t]
        uR = u[i - 1, t]
        if(uD == uR):
            return uU
        uHat = (uU - uR)/(uD - uR)
        if(uHat < 0 or uHat > 1):
            return uU
        return uR + (uD - uR)*( 2*(uHat ** 4) - 3*(uHat ** 3) + 2*uHat)

def faceg(u,i,t,alpha,Nx,type): # i - 1/2
    if type == 0: # FOU
        if alpha > 0:
            return u[i - 1,t]
        else:
            return u[i,t]
    elif type == 1: # TOPUS
        uD = u[i,t]
        uU = u[i - 1,t]
        if(i > 1):
            uR = u[i - 2,t]
        else:
            uR = u[i,t] # 0
        if(uD == uR):
            return uU
        uHat = (uU - uR)/(uD - uR)
        if(uHat < 0 or uHat > 1):
            return uU
        return uR + (uD - uR)*( 2*(uHat ** 4) - 3*(uHat ** 3) + 2*uHat)


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

######################################################################

#Entrada
alpha = 1.0 # coeficiente convectivo
plotType = 1 # 0 - FOU | 1 - TOPUS
case = 1 # case1 ou case2
tj = 0.25 #Instante de tempo desejado
dt = 0.0025 #Discretização no tempo s δt = 2.5 × 10−4 e δt = 2.5 × 10-3
Nx = 400 #Quantidade de elementos

if case == 1:
    T = 1 #Tempo final case1 = 1 | case2 = 0.25
    A = 0 #Limite inferior para x
    B = 2 #Limite superior para x
elif case == 2:
    T = 0.25 #Tempo final case1 = 1 | case2 = 0.25
    A = -1 #Limite inferior para x
    B = 1 #Limite superior para x

h = (B - A)/np.float(Nx) #Discretização no espaço
Nt = np.int(T/dt) #Quantidade de iterações no tempo
x = np.linspace(A,B,Nx+1) #Para plot das aproximações

t = int(tj/dt) #Índice correspondente ///

print(alpha*dt/h)
xe = np.linspace(A,B,1000) #Para plot da solução exata
u_e = exata(xe, alpha, tj, case)

u = condicaoInicial(Nx,Nt,h,x,case)
u = condicoesContorno(u,Nx,Nt)

u_fou = fou(u,Nx,Nt,h,dt,alpha)
u_topus = topus(u,Nx,Nt,h,dt,alpha)

#Plota o grafico
plt.title('Exata x Aproximada ('+ str(Nx) +  ' elementos , dt = ' + str(dt))
if case == 1:
    plt.xlim (A ,B)
#plt.ylim (0 ,1.1)
plt.grid()

plt.plot(xe, u_e, 'r-', label = 'Exata')
methodName = ''
if plotType == 0:
    methodName = 'FOU_'
    plt.plot(x, u_fou[:,t], 'bx', label = 'FOU')
elif plotType == 1:
    methodName = 'TOPUS_'
    plt.plot(x, u_topus[:,t], 'gx', label = 'TOPUS')

plt.xlabel( 'x' )
plt.ylabel ( 'u' )
plt.legend(loc='best')
diretorio = ""
dtAux = str(dt/h)
dtAux = dtAux.replace('.','-')
nomefig = methodName + dtAux + '_case_' + str(case) + '.png'
plt.savefig(diretorio+nomefig, dpi=200)
plt.show()
