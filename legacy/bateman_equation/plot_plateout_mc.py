import numpy as np
import matplotlib.pyplot as plt

ISOTOPES = np.array([	("226Ra", 1E6, 1600*365*24*3600, 1),
                        ("222Rn", 0, 3.8*24*3600, 1),
			("218Po", 0, 3.1*60, 1),
			("214Pb", 0, 26.8*60, 1),
			("214Bi", 0, 19.9*60, 1)],

#			("214Po", 0, 164E-6)],
			dtype=[('name', '<U10'), ('N_0', '<f4'), ('t_1_2', '<f4'), ('alpha', "<i1")])

S1_ONLY_RATIOS = [1,0.956592512208356, 0.7596310363537709,0.4389582202930006]

LAMBDAS = np.log(2) / ISOTOPES["t_1_2"]

PLATE_HALF_LIFE = 0.04*24*3600
PLATE_CONST = np.log(2) / PLATE_HALF_LIFE

N_ISOS = len(ISOTOPES["name"])
N_TIMESTEPS = int(1E5)
DELTA_T = 1

# T_12 = log(2) / lambda

numbers = np.zeros((N_TIMESTEPS, N_ISOS))


numbers[0,:] = ISOTOPES["N_0"]

for idx in range(N_TIMESTEPS - 1):
    if idx == 0:
        continue
    for jdx in range(N_ISOS):
       decay = LAMBDAS[jdx] * numbers[idx-1, jdx]
       plate = 0
       feed  = 0
       if jdx > 0:
           feed = LAMBDAS[jdx-1] * numbers[idx-1, jdx-1]
           if ISOTOPES['alpha'][jdx-1]:
               plate = PLATE_CONST * numbers[idx-1, jdx]

       numbers[idx, jdx] = numbers[idx-1, jdx] + feed - decay - plate


print(numbers)

plt.figure()

times = np.linspace(0, N_TIMESTEPS*DELTA_T, N_TIMESTEPS)

for jdx in range(N_ISOS):
    plt.plot(times, LAMBDAS[jdx] * numbers[:, jdx] / (LAMBDAS[1] * numbers[:,1]), label=ISOTOPES["name"][jdx])

for this_y in S1_ONLY_RATIOS:
    plt.axhline(this_y)

plt.ylim(0,1)
plt.legend()
plt.show()
