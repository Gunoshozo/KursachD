import numpy as np
import matplotlib.pyplot as plt

R = 4
T = 20
k_ = 0.04
C = 1.84
a2 = k_ / C

I = 4
K = 1

lab = ""


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
        a1 = 2*(U[k][I_-1] -  U[k][I_]) / (h_r ** 2)
        U[k+1][I_] = a2 * h_t * a1  + U[k][I_]
    return U


def plot_layer(ax,x,y):
    X = np.array(x,dtype=float)
    Y = np.array(y,dtype=float)
    ax.plot(X,Y,label = lab)
    ax.legend()

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

def get_solution_I(i):
    global K,I,h_t,h_r
    I = i
    K = int(round((I ** 2) * 6 * a2 * T / (R ** 2)))
    h_r = R / I
    h_t = T / K
    r = np.linspace(0, R, I)
    t = np.linspace(0, T, K)
    return r, t,solve(r,t)

def densing_net():
    global  lab
    ax = plt.subplot()
    for i in range(1,6):
        r,t,sol = get_solution_I(5*i)
        lab = "I="+str(i*5)
        half = len(sol)//2
        plot_layer(ax,r, sol[-1])
    plt.xlabel("r")
    plt.ylabel("u")
    plt.grid()
    plt.show()

def main():
    densing_net()
    # global lab
    # set_vars(True)
    # r = np.linspace(0, R, I)
    # t = np.linspace(0, T, K)
    # U = solve(r,t)

    #r,t,sol = get_solution_I(30)
    #n_max = 9
    #plt.grid()
#
    #res = []
    #plt.figure(1)
    #ax = plt.subplot()
    #for i in range(n_max+1):
    #    lab = "K = " + str(int(i/n_max*(K-1)))
    #    plot_layer(ax,r,sol[int(i/n_max*(K-1))])
    #    print(int(i/n_max*(K-1)))
    #    res = [res,sol]
    #ax.set_xlabel("r")
    #ax.set_ylabel("u")
    #plt.savefig(str(I) + " " + str(K)+".png")

    #layer = int(input("desired layer:"))
    # plot_layer(r,U[layer])
    #
    # fig = plt.figure()

    #r,t = np.meshgrid(r,t)
    #res = v_analytical(r,t)
    #ax = plt.axes(projection="3d")
    #ax.plot_surface(X=t, Y=r, Z=U, cmap='magma', vmin=np.min(U), vmax=np.max(U))
    #ax.plot_surface(X=t, Y=r, Z=res, cmap='plasma', vmin=np.min(U), vmax=np.max(U))
    #ax.set_xlabel('t')
    #ax.set_ylabel('r')
    #ax.set_zlabel('u')
    #plt.show()


if __name__ == '__main__':
    main()