import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, np.pi * 2)
y = np.sin(x)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, y)
plt.show()