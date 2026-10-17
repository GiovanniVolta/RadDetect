import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 


data = pd.read_fwf("RANGE_3D.txt", header=None, skiprows=21, names=["i", "d", "x", "y"])

print("Mean:", np.mean(data["d"]))
print("std:",np.std(data["d"]))


markersize=13
plt.figure(figsize=(12,12), dpi=50)

plt.hist(data["d"]*1E-1, bins=200, label=r'$\mathrm{E_{kin} = 30\,keV}$', color=tableau20[3])



plt.xlabel('Projected range [nm]', fontsize=18, horizontalalignment='right', x=1.0)
plt.locator_params(nbins=8)

handles, labels = plt.gca().get_legend_handles_labels()
order = [0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=17, numpoints=1, ncol=2, frameon=False)
ax=plt.gca()

#ax.set_ylabel(r'$^{226}$Ra density [A.U.]', fontsize=18, verticalalignment='top', y=0.6375, labelpad=40)
ax.set_ylabel(r'$^{226}$Ra density [A.U.]', fontsize=18, verticalalignment='top', y=0.77, labelpad=40)

ax.tick_params(which='major', direction='in', length=15, width=1)
ax.tick_params(which='minor', direction='in', length=7, width=1)
for item in [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels() + ax.get_legend().get_texts():
        item.set_fontsize(32)
plt.yscale('log')
#plt.xlim(-0.07, 1.45)
#plt.ylim(20, 2100)

plt.show()
