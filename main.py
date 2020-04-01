import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import threading
from mpl_toolkits import mplot3d

R = 4
T = 20
k_ = 0.04
C = 1.84
a2 = k_ / C

I = 4
K = 1

def Psi(r):
    return np.cos(np.pi*r/R)

def solve(r,t):
    U = np.zeros((len(t),len(r)),dtype=np.float64)

    U[0] = Psi(r)
    for k in range(0, K - 1):
        #i == 0
        U[k+1][0] = 6 * a2 * h_t * (U[k][1] - U[k][0]) / (h_r ** 2) + U[k][0]

        #  1 < i < I-1
        for i in range(1, I - 1):
            c1 = (U[k][i+1]-2*U[k][i]+U[k][i-1])/(h_r**2)
            c2 = (U[k][i+1]-U[k][i-1])/(h_r*r[i])
            U[k+1][i] = a2 * h_t * (c1 + c2) + U[k][i]
        I_ = I - 1
        c1 = (U[k][I_-1] - 2 * U[k][I_] + U[k][I_-1]) / (h_r ** 2)
        c2 = (U[k][I_-1] - U[k][I_-1]) / (h_r * r[i])
        U[k + 1][i] = a2 * h_t * (c1 + c2) + U[k][i]
        U[k+1][I_] = a2 * h_t * (c1 + c2) + U[k][I_]
    return U


def plot_layer(x,y):
    X = np.array(x,dtype=float)
    Y = np.array(y,dtype=float)
    plt.plot(X,Y)
    plt.grid()
    plt.show()
    return

def plot_3d(x,y,z):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X=x,Y=y,Z=z, cmap='plasma',vmin=np.min(z),vmax=np.max(z))
    ax.set_xlabel('t')
    ax.set_ylabel('r')
    ax.set_zlabel('u')
    plt.show()

def set_vars(isR):
    global K,I,h_t,h_r
    if(isR):
        I = int(input("I:"))
        K = int(round((I**2)*6*a2*T/(R**2)))
        h_r = R/I
        h_t = T/K
        print("K:"+str(K))
        print("h_r:"+str(h_r))
        print("h_t:"+str(h_t))
    else:
        K = int(input("K:"))
        I = int(round(np.sqrt(K*R**2/(6*a2*T))))
        h_r = R / I
        h_t = T / K
        print("I:" + str(K))
        print("h_t:" + str(h_t))
        print("h_r:" + str(h_r))

def main():
    set_vars(True)
    r = np.linspace(0, R, I)
    t = np.linspace(0, T, K)
    U = solve(r,t)
    #layer = int(input("desired layer:"))
    #plot_layer(r,U[layer])
    r,t = np.meshgrid(r,t)
    plot_3d(t,r,U)


if __name__ == '__main__':
    main()