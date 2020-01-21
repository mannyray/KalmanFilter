import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

results = [];
for assumedProcessNoise in np.arange(0.05,0.25,0.01):
    localResults = []
    for assumedSensorNoise in [0.01]:#np.arange(0.1,0.6,0.01):
        measurementFile = open("assumedP_"+"{0:.2f}".format(assumedProcessNoise)+
            "assumedR"+"{0:.2f}".format(assumedSensorNoise)+"measurements.txt")
        estimatesBeforeFile = open("assumedP_"+"{0:.2f}".format(assumedProcessNoise)+
            "assumedR"+"{0:.2f}".format(assumedSensorNoise)+"estimatesBefore.txt")
        measurements =  np.array([float(i) for i in measurementFile.readlines()])
        #measurements = measurements[1:500]
        estimates = np.array([float(i) for i in estimatesBeforeFile.readlines()])
        #estimates = estimates[1:500]
        localResults.append(sum(measurements - estimates)/len(measurements)) 
        print sum(measurements - estimates)/len(measurements) 
    results.insert(0,localResults)

print results
print measurements - estimates
print sum(measurements - estimates)


for i in results:
    for j in i:
        print round(j,4),
    print

fig,ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
im = ax.imshow(results,extent=(0.4,0.6,0.4,0.6))
fig.colorbar(im,cax=cax,orientation='vertical')
plt.show()

plt.plot(estimates)
plt.show()



