import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

data = np.loadtxt("TimeSeries.txt", delimiter=' ', ndmin=2)
plt.plot(data[:,0], data[:,1], '-')
plt.savefig("TimeSeries.png")
plt.clf()

plt.plot(data[:,0], data[:,2], '-')
plt.savefig("DegreeSeries.png")
plt.clf()

