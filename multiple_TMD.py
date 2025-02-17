import numpy as np
import matplotlib.pyplot as plt


#constants
m = 1.83
TMD_ratio = 0.1
ma = m * TMD_ratio
L = 0.2 #length
b = 0.08 # width
E = 210E9 # Young's Modulus
d = 0.001 # thickness
I = b*d*d*d/12 # second moment of area
k = (24*E*I)/(L*L*L) # static stiffness for each floor

n = 100  # number of TMDs
ka = np.zeros(n)
for i in range(n):
    ka[i] = k * TMD_ratio * (1 + i/100)
w = np.linspace(0,100,1000)
y = []
y1 = []
def M(n):
    entries = [m]
    for i in range(n):
        entries.append(ma/n)
    matrix = np.diag(entries)
    return matrix

def K(n):
    matrix = np.zeros((n+1,n+1))
    for i in range(n+1):
        if i == 0:
            matrix[i][i] = k + ka.sum()
        else:
            matrix[i][i] = ka[i-1]
    for i in range(1,n+1):
        matrix[0][i] = -ka[i-1]
        matrix[i][0] = -ka[i-1]
    return matrix

def F(n):
    matrix = np.zeros((n+1,1))
    matrix[0] = 1
    return matrix

def plot_with_TMD(n):
    Mn = M(n)
    Kn = K(n)
    Fn = F(n)
    print(Fn.size)
    print(Mn.size)
    print(Kn.size)
    for omega in w:
        B = -omega*omega*Mn + Kn
        #print(F)
        disp = np.linalg.inv(B) @ Fn
        print(disp)
        y.append(abs(disp[0][0]))
        y1.append(abs(disp[1][0]))
    #print(y)
    plt.ylabel('y/f')
    plt.xlabel('frequency/Hz')
    plt.ylim(0,0.02)
    plt.title("Response with %i TMDs"%n)
    plt.plot(w/(2*np.pi),y, label = 'displacement')
    #plt.plot(w/(2*np.pi),y1, label = 'displacement')
    #plt.axvline(10, color = 'orange', linestyle = '--', label = 'f_max')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    #print(F(n))
    plot_with_TMD(n)
    print(K(n))
