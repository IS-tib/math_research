import numpy as np
import matplotlib.pyplot as plt

l = 10
Nx = 100
dx = l/(Nx-1)
x = np.linspace(0,l, Nx)

time = 20
Nt = 2000
dt = 0.01
t = np.linspace(0, time, Nt)

T = np.zeros((Nt, Nx))
E = np.zeros((Nt, Nx))
A = np.zeros((Nt, Nx))

A0 = 5.5
Da = 1e-9
lamda=0.086
kt = 0.05
ke = 0.03
k=0.3/1e5

A[0, :] = 5.5
T[0, :] = 1e5
E[0, :] = 1

for i in range(Nt-1):
    for j in range(1,Nx-1):

        A[i+1][j] = A[i][j] + dt * (Da * (A[i][j+1] - 2*A[i][j] + A[i][j-1]) / dx**2 - lamda*A[i][j])
        T[i+1][j] = T[i][j] + dt * (-1*kt*A[i][j]*T[i][j])
        E[i + 1][j] = E[i][j] + dt * (-1 * ke * A[i][j] * E[i][j] - k * ((T[i + 1][j] - T[i][j]) / dt))
        #E[i+1][j] = E[i][j] + dt * (-1*ke*A[i][j]*E[i][j]-k*max(0,((T[i+1][j]-T[i][j])/dt)))

    A[i+1][0] = A[i+1][1]
    A[i+1][-1] = A[i+1][-2]
    T[i + 1][0] = T[i + 1][1]
    T[i + 1][-1] = T[i + 1][-2]
    E[i + 1][0] = E[i + 1][1]
    E[i + 1][-1] = E[i + 1][-2]

center = Nx // 2

plt.figure(1)
plt.plot(t, A[:, center], label="Iodine Activity Concentration at Center")
plt.xlabel("Time (days)")
plt.ylabel("Activation Concentration (GBq/mL)")
plt.title("Iodine Activity Concentration at Center")
plt.grid(True)
plt.legend()

plt.figure(2)
plt.plot(t, T[:, center], label="Thyroid Cell Density")
plt.xlabel("Time (days)")
plt.ylabel("Cell Density")
plt.title("Cell Density Over Time at Center")
plt.legend()
plt.grid(True)

plt.figure(3)
plt.plot(t, E[:, center], label="Salivary Gland Ejection Fraction")
plt.xlabel("Time (days)")
plt.ylabel("Percentage (%)")
plt.title("Ejection Fraction Over Time at Center")
plt.legend()
plt.grid(True)

plt.show()














