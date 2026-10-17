import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import scipy.optimize as opt 
from matplotlib.ticker import FuncFormatter
formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))

import matplotlib as mpl
import cycler

n = 10
color = plt.cm.hsv(np.linspace(0, 1,n))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)


n_lines = 10


# Approach for analytical matrix solution of Batemans equation
# Taken from this paper Algebraic approach to the radioactive decay equations
# DOI: 10.1119/1.1571834

# To have the most flexibility lets try to implement the complete chain

# Thorium chain
# 228Th->224Ra->220Rn->216Po->212Pb->212Bi->212Po

# 225Ra->225Ac->221Fr->217At->213Bi->213Po


labels = ["228Th","224Ra","220Rn","216Po","212Pb","212Bi","212Po"]
t_1_2 = np.array((	1.9116*365*24*3600,
			3.6319*24*3600,
			55.6,
                	145e-3,
			10.64*3600,
                	60.55*60,
                	299e-9))


labels = ["225Ra","225Ac","221Fr","217At","213Bi","213Po"]
t_1_2 = np.array((	15*24*3600,
			10*24*3600,
			4.9*60,
                	32E-3,
			46*60,
                	4.2E-6))

lambdas = np.log(2)/t_1_2
n_isotopes = len(lambdas)


# Start with construction of a matrix [A], describing the system of PDEs
# For an non-branching chain this has the -lambda_i on the diagonal and
# +lambda_i on the lower diagonal
A_lower = np.eye(n_isotopes, k=-1)*lambdas
A_diag = np.diag(-lambdas)
A = A_lower + A_diag

# Calculate the eigenvalues and eigenvectors of the [A]
A_eig_vals, V = np.linalg.eig(A)

# We also need [V]^-1, which is the inverse matrix of [V]
V_inv = np.linalg.inv(V)


def get_n_t(t, a0):

	# calculate the initial number of particles for each isotope
	n_0 = a0 / lambdas
	n_t = np.zeros((len(t), n_isotopes))

	# Now the vector holding the number of atoms [N] at time t is given by:
	# [N] = [V] * [Lambda] * [V]^-1 * [N_0]

	for idx, time in enumerate(t):
		# Where [Lambda] is given as
		Lambda = np.diag(np.exp(A_eig_vals*time))
		n_t[idx] = np.dot(np.dot(np.dot(V, Lambda), V_inv), n_0)
	return n_t


def get_a_t_po212(t, *a0):

	return get_n_t(t, a0)[:,-1] * lambdas[-1]





# Just draw some functional evaluation of the Bateman equation
n_time_bins = 1000
time_bounds = [1e-1,14*7*24*3600]
#times = np.linspace(time_bounds[0],time_bounds[1], n_time_bins)
times = np.geomspace(time_bounds[0],time_bounds[1], n_time_bins)

# Set initial activities and choose 1Bq of 224Ra present
a_0 = np.zeros_like(lambdas)
a_0[0] = .05
a_0[1] = 1

numbers = get_n_t(times, a_0)
activities = numbers*lambdas

isotope_colors = ["black", "royalblue", "tab:red", "tab:cyan", "darkviolet", "green"]
fontsize_1 = 20

indicator_box_y = [7e-3, 1e-2]

plt.figure(figsize=(8,6), constrained_layout=True, dpi=150)

for idx, isotope in enumerate(labels):
     if idx == 0 or idx == len(labels)-1:
         pass 
         #continue
     # Prevent drawing the lines over the boxes indicating the times
     # To do so, cut the times, where the activity is below the box
     time_cut_arr = activities.T[idx] > indicator_box_y[1]*1.01
     plt.plot(times[time_cut_arr]/3600., activities.T[idx,time_cut_arr], label = labels[idx], 
              color = isotope_colors[idx], linewidth = 3)

# Draw ranges at the top, indicating what are seconds, minutes, hours and days

indicator_box_color = ["dimgray", "gray", "darkgray", "lightgray"]
indicator_box_color = ["gray", "darkgray", "lightgray", "whitesmoke"]
indicator_box_text = ["seconds", "minutes", "hours", "days"]
indicator_box_range = [(1./3600., 60/3600.),
                       (60/3600., 1.),
                       (1., 24.),
                       (24., time_bounds[1]/3600.)]

for idx, this_text in enumerate(indicator_box_text):
    plt.fill_between([indicator_box_range[idx][0],indicator_box_range[idx][1]], 
                     indicator_box_y[0], indicator_box_y[1], color=indicator_box_color[idx])
    box_center_x = np.sqrt(indicator_box_range[idx][0] * indicator_box_range[idx][1])
    box_center_y = np.sqrt(indicator_box_y[0] * indicator_box_y[1])
    plt.text(box_center_x, box_center_y, this_text, ha="center", va="center",
             fontsize=fontsize_1)

# Draw text labels onto the lines
label_positions = [[0, 0],
                   [4e-4, 0.7],
                   [4e-4, 0.035],
                   [6e-4, 0.015],
                   [0.13, 0.015],
                   [1.45, 0.014],
                   [0,0],
                   [0,0]]

label_rotations = [0, 0, 65, 65, 65, 75, 0, 0]

label_text = labels


for idx, this_color in enumerate(isotope_colors):
    if idx == 0 or idx == len(labels)-1:
         continue
    plt.text(label_positions[idx][0], label_positions[idx][1], label_text[idx],
             rotation=label_rotations[idx], fontsize=fontsize_1, color=isotope_colors[idx])

plt.gca().tick_params(axis='both', labelsize=fontsize_1)

# plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')

plt.gca().xaxis.set_major_formatter(formatter)
plt.gca().yaxis.set_major_formatter(formatter)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top') 

plt.ylim(7e-3,1.5)
plt.xlim(1./3600,time_bounds[1]/3600.)
plt.ylabel("Activity (Bq)", fontsize=fontsize_1)
plt.xlabel("Time (hours)", fontsize=fontsize_1, labelpad=10)

plt.savefig("ra224_decay_isotope_evolution_1.pdf", bbox_inches="tight", dpi="figure")
plt.show()
