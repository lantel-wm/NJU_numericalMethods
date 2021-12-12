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
        print(type(line))
        up.append(map(int, line.split(' ')))
    f.close()

print(up)
