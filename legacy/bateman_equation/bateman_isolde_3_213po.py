import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import scipy.optimize as opt 


# Define where the .root files are stored on the LFS1
base_path = "/d12/lin/xenon/radon_monitor/database/backup/"
file_ending = ".root"

# Define runs and diferent intervals

run_alpha_pos = ["Rn22122017", "Rn04012018", "Rn14012018"]

time_since_implantation = 91*24*3600

# define time intervals to cut
time_cut = 1009.1 
time_cut = time_cut * 3600*24

file_names = run_alpha_pos


po212_selection = [400,430]
po212_selection = [250,300]

bin_width = int(20 * 60. * 10)

bin_centers = np.array(())
bin_counts = np.array(())
start_times = np.zeros((len(file_names)))
stop_times = np.zeros((len(file_names)))


for i in range(len(file_names)):
	
	this_file_path = base_path + file_names[i] + "/" + file_names[i] + file_ending
	print(this_file_path)
	this_tree = up.open(this_file_path)["t"]

	time_stamp = this_tree.array("timestamp")
	channel = this_tree.array("channel")

	start_times[i] = np.min(time_stamp)
	stop_times[i] = np.max(time_stamp)
	duration = stop_times[i] - start_times[i]

	event_mask = np.less(channel, po212_selection[1]) * np.greater(channel, po212_selection[0])
	selected_timestamps = time_stamp[event_mask]


	bins, bins_width = np.linspace(start_times[i], stop_times[i], int(duration/bin_width), retstep=True)
	bin_conts, bin_edges = np.histogram(selected_timestamps, bins=bins)
	bin_cent = (bin_edges[:-1] + bin_edges[1:])/2.

	bin_counts = np.append(bin_counts, bin_conts)
	bin_centers = np.append(bin_centers, bin_cent)

t_min = start_times[0]
t_max = stop_times[-1]

# Approach for analytical matrix solution of Batemans equation
# Taken from this paper Algebraic approach to the radioactive decay equations
# DOI: 10.1119/1.1571834

# To have the most flexibility lets try to implement the complete chain

# 225Ra->225Ac->221Fr->217At->213Bi->213Po

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

def get_full_po213_evolution(t, a_0_225ra, a_0_225ac):

	# Convert input time from seconds to days
	t = t*(24*3600) + time_since_implantation

	# This function will allow to fit the 212Po activity for different
	# Data taking conditions (i.e. coated, uncoated, de-coated)
	# in a combined fashion

	# Start with initializing the vector of initial conditions to be zero
	# for all isotopes.
	# Now fill the ones, which we like to fit and/or insert constraints

	# The initial vector a_0_1 will be initial conditions for first interval
	# We assume equilibrium of activity between 224Ra -> 216Po
	# As well as between 212Pb->212Po. So we have two degrees of freedom
	
	a_0_0 = np.zeros((n_isotopes))
	a_0_0[0] = a_0_225ra
	a_0_0[1] = a_0_225ac
	a_0_0[2:] = 0

	# Calculate the 212Pb activity at the first transition
	
	n_t_0 = get_n_t(t, a_0_0)

	a_t_212po = n_t_0[:,-1] * lambdas[-1]
	return a_t_212po



interval_rejection_mask = (bin_centers-t_min) < time_cut

print(interval_rejection_mask)
n_last_point = 300

bin_centers_cut = bin_centers[interval_rejection_mask]
bin_counts_cut = bin_counts[interval_rejection_mask]
bins_width_cut = bins_width

# define initial values
init_vals = np.array((1,1))
p_opt, p_cov = opt.curve_fit(get_full_po213_evolution, (bin_centers_cut-t_min)/3600./24, bin_counts_cut/bins_width_cut * 1e3,
				sigma = np.sqrt(bin_counts_cut)/bins_width_cut*1e3, p0 = init_vals, bounds=[0, np.inf],
				absolute_sigma = True)

print(p_opt)
print(p_cov)

# Get array of x-values and corresponding fit values
x_vals = np.linspace(0, (t_max-bin_centers[0])/3600/24, 10000)
a_fit = get_full_po213_evolution(x_vals, *p_opt)

# also calculate the residuals between fit and data
residuals = (bin_counts/bins_width*1e3 - get_full_po213_evolution((bin_centers-t_min)/3600./24, *p_opt))/(np.sqrt(bin_counts)/bins_width*1e3)


print("A_0(225Ra): ", p_opt[0], "+/-", np.sqrt(p_cov[0,0]))
print("A_0(225Ac): ", p_opt[1], "+/-", np.sqrt(p_cov[1,1]))



fig = plt.figure(figsize=(12,8), dpi=80) 
gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)


ax1.plot(x_vals, a_fit, color="tab:red", linewidth=3)
ax1.errorbar((bin_centers-t_min)/3600./24., bin_counts/bins_width*1e3, 
              yerr = np.sqrt(bin_counts)/bins_width*1e3, xerr=bins_width/3600./24. * 0.5,
		fmt='.', color='black', linewidth=3, markersize=1)


# Add residual plot


ax2.plot((bin_centers-t_min)/3600./24., residuals, '.', color='black', markersize=10)
# ax2.grid(axis='y')
ax2.fill_between([-1,15], 2*[-1], 2*[1], color="tab:red", alpha=0.5)
ax2.fill_between([-1,15], 2*[-2], 2*[2], color="tab:red", alpha=0.25)


ax1.set_ylabel("Detected $^{213}$Po Rate (mHz)")
ax1.get_xaxis().set_visible(False)
ax2.set_ylabel("Residual")
ax2.set_xlabel("Measurement time (days)")

ax1.set_ylim(0.1,15.9)
ax2.set_ylim(-3.9,3.9)
ax2.set_xlim(-0.5,12)

for item in [ax1.title, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
	item.set_fontsize(25)  
for item in [ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels():
	item.set_fontsize(25)  

fig.tight_layout()

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig("po213_time_dependence_1.pdf",dpi="figure")
plt.show()


# Just draw some functional evaluation of the Bateman equation
'''
n_time_bins = 10000
time_bounds = [1e-3,2000*365*24*3600]
times = np.geomspace(time_bounds[0],time_bounds[1], n_time_bins)
numbers = get_n_t(times, a_0)
activities = numbers*lambdas

plt.figure()
for idx,act in enumerate(activities.T):
	plt.plot(times, act, label=labels[idx])
plt.tegend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-9,1e1)
plt.show()
'''