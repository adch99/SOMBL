import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.special import spherical_jn

L = 100
# qx = 1
# qy = 1
# theta = 0
theta = np.pi/4

def I_S_x(qx, qy, theta):
    S_x_avg = np.zeros(qx.shape)
    for x in range(L):
        for y in range(L):
            # phi = qx*x + qy*y + theta - 1.4*(np.pi/100) * x**2 + (2.3*np.pi/100) * y**2
            phi = qx*x + qy*y + theta
            S_xi = np.sin(2*phi)
            # a = spherical_jn(3, phi)
            # b = np.sqrt(1 - a**2)
            # S_xi = 2 * a * b
            if (x + y) % 2 == 0:
                S_x_avg += S_xi
            else:
                S_x_avg -= S_xi
    return S_x_avg / (L*L/2)

qx = qy = np.linspace(0, 4*np.pi, 200)
# qx = qy = np.linspace(-2*np.pi/100, 2*np.pi/100, 200) + np.pi/2
QX, QY = np.meshgrid(qx, qy)
Z = I_S_x(QX, QY, theta)
high_regions = np.where(np.abs(Z) > 0.1)
print(QX[high_regions]/np.pi)
print(QY[high_regions]/np.pi)

# norm = mpl.colors.LogNorm(vmin=1e-20, vmax=Z.max())
mappable = plt.contourf(QX, QY, Z)
# mappable = plt.imshow(Z)
plt.colorbar(mappable)

plt.show()