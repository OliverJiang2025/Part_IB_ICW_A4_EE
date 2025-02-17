import numpy as np
import matplotlib.pyplot as plt

#change this bool var to decide if the simulation contains a single tmd
include_TMD = True

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
ka = k * 0.3
w = np.linspace(0,100,1000)
y = []


def plot_without_TMD():
    y = abs(1/(-w*w*m + k))
    plt.ylabel('y/f')
    plt.xlabel('frequency/Hz')
    plt.ylim(0,0.02)
    plt.title("Response with no TMD")
    plt.plot(w/(2*np.pi),y)
    plt.show()


def plot_with_TMD():
    M = np.array([[ma, 0],
                  [0,m]])
    K = np.array([[ka, -ka],
                  [-ka, k+ka]])
    F = np.array([[0],
                  [1]])
    for omega in w:
        B = -omega*omega*M + K
        #print(F)
        disp = np.linalg.inv(B) @ F
        print(disp)
        y.append(abs(disp[0][0]))
    #print(y)
    plt.ylabel('y/f')
    plt.xlabel('frequency/Hz')
    plt.ylim(0,0.02)
    plt.title("Response with one TMD")
    plt.plot(w/(2*np.pi),y)
    plt.show()


if __name__ == "__main__":
    if not include_TMD:
        plot_without_TMD()
    elif include_TMD:
        plot_with_TMD()

