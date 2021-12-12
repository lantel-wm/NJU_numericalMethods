import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

x = np.linspace(-2, 2, 17)
y = np.linspace(-2, 2, 17)
x, y = np.meshgrid(x, y)

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

plt.subplots(figsize=(12, 8))

plt.xlabel('X')
plt.ylabel('Y')

plt.quiver(x, y, up, vp)

plt.savefig('div_wind.png')
