import numpy as np
import matplotlib.pyplot as plt

l = 10
Nx = 100
dx = l/(Nx-1)
x = np.linspace(0,l, Nx)

time = 20#simulate for twenty days
Nt = 2000
dt = 0.01 #if this gets bigger, it starts oscilltating
t = np.linspace(0, time, Nt)

density_T = np.zeros((Nt, Nx))
eject_frac = np.zeros((Nt, Nx))

radiation_dose = np.zeros((Nt, Nx))

#placeholders

A0 = 5.5
De = 0.2
r=0.3
K=1e6
Dt = 0.1
k=1e-6
kt = 0.005
density_T[0, :] = 1e5
eject_frac[0, :] = 1

for i in range(Nt):
    for j in range(1, Nx-1):
        radiation_dose[i][j] = A0 * np.exp(-0.11 * (i * dt))

#finite difference method

for i in range(Nt-1):
    for j in range(1,Nx-1):

        eject_frac[i+1][j] = eject_frac[i][j] + dt * (De * (eject_frac[i][j+1] -
                            2*eject_frac[i][j] + eject_frac[i][j-1]) / dx**2 -
                            k*eject_frac[i][j]*radiation_dose[i][j])

        density_T[i+1][j] = density_T[i][j] + dt * (Dt * (density_T[i][j+1] -
                            2 * density_T[i][j] + density_T[i][j-1]) /
                            dx**2 - k*((eject_frac[i+1][j] - eject_frac[i][j]) / dt) -
                            kt*radiation_dose[i][j]*density_T[i][j])

    density_T[i+1][0] = density_T[i+1][1]
    density_T[i+1][-1] = density_T[i+1][-2]

center = Nx // 2

plt.figure(1)
plt.plot(t, radiation_dose[:, center], label="Radiation Dose at Center")
plt.xlabel("Time (days)")
plt.ylabel("Radiation Dose (becquerel/day)")
plt.title("Radiation Dose Over Time at Center")
plt.grid(True)
plt.legend()

plt.figure(2)
plt.plot(t, density_T[:, center], label="Thyroid Cell Density")
plt.xlabel("Time (days)")
plt.ylabel("Cell Density")
plt.title("Cell Density Over Time at Center")
plt.legend()
plt.grid(True)

plt.figure(3)
plt.plot(t, eject_frac[:, center], label="Salivary Gland Ejection Fraction")
plt.xlabel("Time (days)")
plt.ylabel("Percentage")
plt.title("Ejection Fraction Over Time at Center")
plt.legend()
plt.grid(True)

plt.show()














