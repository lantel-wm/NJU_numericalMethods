import pandas as pd
import numpy as np
import re

x = np.linspace(-2, 2, 17)
y = np.linspace(-2, 2, 17)
grid = np.meshgrid(x, y)

U = pd.read_excel('u.xls', usecols=range(1, 18))
V = pd.read_excel('v.xls', usecols=range(1, 18))

np.set_printoptions(linewidth=np.inf)
with open('u.txt', 'w') as f:
    f.write(' ')
    f.write(re.sub('[\[\]]', '', np.array_str(U.values.transpose())))
    f.close()

with open('v.txt', 'w') as f:
    f.write(' ')
    f.write(re.sub('[\[\]]', '', np.array_str(V.values.transpose())))
    f.close()

with open('gird.txt', 'w') as f:
    f.write(' ')
    f.write(re.sub('[\[\]]', '', np.array_str(grid.values.transpose())))
    f.close()
