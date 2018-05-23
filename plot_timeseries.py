import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

data = np.loadtxt("TimeSeries.txt", delimiter=' ')
plt.plot(data[:,0], data[:,1], '-')
plt.savefig("TimeSeries.png")
plt.clf()

