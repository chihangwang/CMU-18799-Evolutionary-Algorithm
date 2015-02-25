import numpy as np
import matplotlib.pyplot as plt
import re

stats = open('stats', 'r')

fits = []
sizes = []

for line in stats.readlines():
    line = line.strip()
    fitness = re.search('^\d*\.\d*', line).group()
    size    = re.search('\d*\.\d*$', line).group()
    # print "fitness: %s size: %s" %   (fitness, size)
    fits.append(float(fitness))
    sizes.append(float(size))

stats.close()

x = sizes
y = fits
z = np.arange(0, 91, 1)
max = [89]*91
plt.xlabel('generation')
plt.plot(z, x, 'bs')
plt.annotate('fitness', xy=(15, 40), xytext=(25, 40),
            arrowprops=dict(facecolor='green', shrink=0.05))
plt.plot(z, y, 'g^')
plt.annotate('size', xy=(40, 30), xytext=(50, 30),
            arrowprops=dict(facecolor='blue', shrink=0.05))
plt.plot(z, max, 'r--')
plt.show()