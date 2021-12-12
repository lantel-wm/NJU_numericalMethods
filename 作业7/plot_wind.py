#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# initial gird
x = np.linspace(-2, 2, 17)
y = np.linspace(-2, 2, 17)
x, y = np.meshgrid(x, y)

#read U, V
dfU = pd.read_excel('u.xls', usecols=range(1, 18))
dfV = pd.read_excel('v.xls', usecols=range(1, 18))

U = dfU.values
V = dfV.values

#read u_prime, v_prime, phi
up = []
vp = []
phi = []
with open('up.txt', 'r') as f:
    for line in f:
        up = list(map(float, line.split()))
    f.close()

with open('vp.txt', 'r') as f:
    for line in f:
        vp = list(map(float, line.split()))
    f.close()

with open('phi.txt', 'r') as f:
    for line in f:
        phi = list(map(float, line.split()))
    f.close()

up = np.array(up).reshape(17, 17)
vp = np.array(vp).reshape(17, 17)
phi = np.array(phi).reshape(17, 17)

# plot original wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, U, V)
plt.title('Original Wind Field')
plt.savefig('wind.png')
plt.close()

# plot div wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

#  contourf = plt.contourf(x, y, phi, cmap='flag')
contour = plt.contour(x, y, phi, np.arange(-0.8, 0.601, 0.1), colors='gray', linestyles='-')
plt.quiver(x, y, up, vp)
plt.clabel(contour, fontsize=10, colors='gray')
#  plt.colorbar(contourf, drawedges=True, orientation='vertical',spacing='uniform')
plt.title('Divergence Wind Field')
plt.savefig('div_wind.png')
plt.close()

# plot vor wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, U - up, V - vp)
plt.title('Vortex Wind Field')
plt.savefig('vor_wind.png')
plt.close()


