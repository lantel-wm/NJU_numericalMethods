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

#read u_prime, v_prime
up = []
vp = []

with open('up.txt', 'r') as f:
    for line in f:
        print(line.split())
        up = list(map(float, line.split()))
f.close()

with open('vp.txt', 'r') as f:
    for line in f:
        print(line.split())
        vp = list(map(float, line.split()))
f.close()

up = np.array(up).reshape(17, 17)
vp = np.array(vp).reshape(17, 17)

# plot original wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, U, V)

plt.savefig('wind.png')
plt.close()

# plot div wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, up, vp)

plt.savefig('div_wind.png')
plt.close()

# plot vor wind field
plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, U - up, V - vp)

plt.savefig('vor_wind.png')
plt.close()


