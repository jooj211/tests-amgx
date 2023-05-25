import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('../build/normaL2.csv')

stepSize = data['stepSize'].tolist()
l2_error = data['l2_error'].tolist()

stepSize = np.array(stepSize)
l2_error = np.array(l2_error)

plt.loglog(stepSize, l2_error)
#plt.plot(stepSize, l2_error, 'bo')
plt.xlabel('Número de pontos de discretização (log scale)')
plt.ylabel('Norma L2 do erro (log scale)')
plt.grid(True)

# ajuste de uma reta aos dados e cálculo da inclinação
#m, b = np.polyfit(np.log(stepSize), np.log(l2_error), 1)
#plt.plot(stepSize, np.exp(b) * stepSize**m, 'r--', label=f'Slope = {m:.2f}')

# ajuste de uma reta aos dados e cálculo da inclinação
#m, b = np.polyfit(stepSize, l2_error, 1)
#plt.plot(stepSize, b + m * stepSize, 'r--', label=f'Slope = {m:.2f}')

plt.legend()
#plt.show()
plt.savefig('../build/normaL2.png')