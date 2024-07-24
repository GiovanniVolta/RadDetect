import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import scipy.optimize as opt 

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

labels = ["228Th","224Ra","220Rn","216Po","212Pb","212Bi","212Po"]
t_1_2 = np.array((	1.9116*365*24*3600,
			3.6319*24*3600,
			55.6,
                	145e-3,
			10.64*3600,
                	60.55*60,
                	299e-9))

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
n_time_bins = 10000
time_bounds = [1e-3,14*24*3600]
times = np.linspace(time_bounds[0],time_bounds[1], n_time_bins)
a_0 = np.zeros_like(lambdas)
a_0[1] = 1

numbers = get_n_t(times, a_0)
activities = numbers*lambdas

plt.figure()
plt.plot(times/3600./24., activities.T[1], label = labels[1], color="black")


# Now sweep the 212Pb initial value and draw the according 212Po activity
a_212pb_initials = np.linspace(0,2*a_0[1], 10)
for idx, this_a_212pb in enumerate(a_212pb_initials):
	a_0[-3:] = this_a_212pb
	numbers = get_n_t(times, a_0)
	activities = numbers*lambdas
	plt.plot(times/3600./24., activities.T[-1], label = r"$A_0(^{212}Pb) = $" + str(np.round(this_a_212pb, 2)))

#plt.figure()
#for idx,act in enumerate(activities.T):
#	plt.plot(times, act, label=labels[idx])

plt.legend(loc=0)
#plt.xscale('log')
#plt.yscale('log')
plt.ylim(1e-9,2)
plt.ylabel("Rate (a.u.)")
plt.xlabel("Time (days)")
plt.show()
