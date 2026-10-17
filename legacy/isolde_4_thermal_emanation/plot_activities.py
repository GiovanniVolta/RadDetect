import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data_points_out.dat", skiprows=1, delimiter="\t")

plt.figure()
plt.errorbar(	data[:,0], data[:,2], xerr=data[:,1], yerr=data[:,3], fmt="o",
		color="cornflowerblue")
plt.ylabel(r"²²²Rn emanation rate (Bq)", fontsize=15)
plt.xlabel(r"Emanation temperature (°C)", fontsize=15)
plt.gca().tick_params(labelsize=13)
plt.show()
