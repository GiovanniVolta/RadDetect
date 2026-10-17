import numpy as np
from matplotlib import pyplot as plt

# Approach for analytical matrix solution of Batemans equation
# Taken from this paper Algebraic approach to the radioactive decay equations
# DOI: 10.1119/1.1571834

# To have the most flexibility lets try to implement the complete chain

# defining the vector of lambdas for the uranium and thorium chain

# Uranium chain
# 226Ra->222Rn->218Po->214Pb->214Bi->214Po->210Pb->210Bi->210Po

labels = ["226Ra","222Rn","218Po","214Pb","214Bi","214Po","210Pb","210Bi","210Po"]
t_1_2_uranium = np.array((	1602*365*24*3600,
				3.8235*24*3600,
				3.1*60,
                        	26.8*60,
				19.9*60,
                        	164.3E-6,
                        	22.3*365*24*3600,
                        	5.013*24*3600,
                        	138.376*24*3600 ))

lambdas_uranium = np.log(2)/t_1_2_uranium
n_isotopes = len(lambdas_uranium)

# Start with construction of a matrix [A], describing the system of PDEs
# For an non-branching chain this has the -lambda_i on the diagonal and
# +lambda_i on the lower diagonal
A_lower = np.eye(n_isotopes, k=-1)*lambdas_uranium
A_diag = np.diag(-lambdas_uranium)
A = A_lower + A_diag

# Calculate the eigenvalues and eigenvectors of the [A]
A_eig_vals, V = np.linalg.eig(A)

# We also need [V]^-1, which is the inverse matrix of [V]
V_inv = np.linalg.inv(V)





def get_n_t(t, a0):

	# calculate the initial number of particles for each isotope
	n_0 = a_0 / lambdas_uranium
	n_t = np.zeros((len(t), n_isotopes))

	# Now the vector holding the number of atoms [N] at time t is given by:
	# [N] = [V] * [Lambda] * [V]^-1 * [N_0]


	for idx, time in enumerate(t):
		# Where [Lambda] is given as
		Lambda = np.diag(np.exp(A_eig_vals*time))
		n_t[idx] = np.dot(np.dot(np.dot(V, Lambda), V_inv), n_0)
	return n_t


# Defining the vector of initial conditions for each isotope in the chain
a_0 = np.zeros((len(lambdas_uranium)))
# Give an initial number [N_0] of atoms for the head of the chain
a_0[1] = 1
a_0[2] = 1
a_0[3] = 1



n_time_bins = 10000
time_bounds = [1e-3,2000*365*24*3600]
times = np.geomspace(time_bounds[0],time_bounds[1], n_time_bins)
numbers = get_n_t(times, a_0)
activities = numbers*lambdas_uranium

plt.figure()
for idx,act in enumerate(activities.T):
	plt.plot(times, act, label=labels[idx])
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.show()

