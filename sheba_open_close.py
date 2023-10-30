import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs/input/Icepack_data/forcing/SHEBA"
filename = "open_clos_lindsay.dat"

df = pd.read_csv(os.path.join(data_path, filename), names=['days', 'open', 'clos'], sep='\s+')

df.set_index('days', inplace=True)
df['net'] = df.open + df.clos

print('Check if opening is strictly positive')
print(np.all(df.open >= 0))

print('Check if closing is strictly negative')
print(np.all(df.clos <= 0))

ax = df.plot(xlabel='Time (days)', ylabel='Deformation (1/s)')
plt.show()
