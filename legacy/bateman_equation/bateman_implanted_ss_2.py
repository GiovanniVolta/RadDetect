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
#runs_uncoated = ['Rn14122020', 'Rn15122020']
runs_uncoated = ['Rn18012021','Rn23012021']
runs_coated = ['Rn25012021', 'Rn30012021']
#runs_decoated = ['Rn29122020', 'Rn04012021']

# define time intervals to cut
time_cut = 1009.1 
time_cut = time_cut * 3600*24

#file_names = runs_uncoated + runs_coated + runs_decoated
file_names = runs_uncoated + runs_coated

po212_selection = [220,320]
bin_width = int(20 * 60.)

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

t_0 = 0
t_1 = 607443 #509.4E3
t_2 = 1210.0E3

def get_full_po212_evolution(t,
				a0_212pb_0, a0_224ra_0,
				a0_212pb_1, a0_224ra_1,
				a0_212pb_2, a0_224ra_2):

	# Convert input time from seconds to days
	t = t*(24*3600)

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
	a_0_0[1:4] = a0_224ra_0
	a_0_0[-3:] = a0_212pb_0
	
	#a_0_0[-3]  = a0_212pb_0
	#a_0_0[-2:] = a0_212bi_0


	# Calculate the 212Pb activity at the first transition
	a_212Pb_t_1 = a0_212pb_1
	#a_212Pb_t_1 = get_n_t((t_1,), a_0_0)[0,-3] * lambdas[-3]

	a_0_1 = np.zeros((n_isotopes))
	a_0_1[1:4] = a0_224ra_1
	a_0_1[-3:] = a_212Pb_t_1


	# Calculate the 212Pb activity at the second transition
	a_212Pb_t_2 = a0_212pb_2
	#a_212Pb_t_2= get_n_t((t_2,), a_0_1)[0,-3] * lambdas[-3]

	a_0_2 = np.zeros((n_isotopes))
	a_0_2[1:4] = a0_224ra_2
	a_0_2[-3:] = a_212Pb_t_2

	
	t_0s = t[t<t_1] - t_0
	t_1s = t[np.logical_and(t>t_1, t<t_2)] - t_1
	t_2s = t[t>t_2] - t_2

	n_t_0 = get_n_t(t_0s, a_0_0)
	n_t_1 = get_n_t(t_1s, a_0_1)
	n_t_2 = get_n_t(t_2s, a_0_2)

	n_t = np.concatenate((n_t_0, n_t_1, n_t_2), axis = 0)
	a_t_212po = n_t[:,-1] * lambdas[-1]
	return a_t_212po



interval_rejection_mask = (bin_centers-t_min) < time_cut

print(interval_rejection_mask)
n_last_point = 300

bin_centers_cut = bin_centers[interval_rejection_mask]
bin_counts_cut = bin_counts[interval_rejection_mask]
bins_width_cut = bins_width

# define initial values
init_vals = np.array((0,5, 0, 1, 0, 2))
p_opt, p_cov = opt.curve_fit(get_full_po212_evolution, (bin_centers_cut-t_min)/3600./24, bin_counts_cut/bins_width_cut,
				sigma = np.sqrt(bin_counts_cut+1)/bins_width_cut, p0 = init_vals, bounds=[0, np.inf],
				absolute_sigma = True)

print(p_opt)
print(p_cov)

# Get array of x-values and corresponding fit values
x_vals = np.linspace(0, (t_max-bin_centers[0])/3600/24, 10000)
a_fit = get_full_po212_evolution(x_vals, *p_opt)

# also calculate the residuals between fit and data
residuals = (bin_counts/bins_width - get_full_po212_evolution((bin_centers-t_min)/3600./24, *p_opt))/(np.sqrt(bin_counts+1)/bins_width)


print("A_0(212Pb): ", p_opt[0], "+/-", np.sqrt(p_cov[0,0]))
print("A_0(224Ra): ", p_opt[1], "+/-", np.sqrt(p_cov[1,1]))
print("A_1(212Pb): ", p_opt[2], "+/-", np.sqrt(p_cov[2,2]))
print("A_1(224Ra): ", p_opt[3], "+/-", np.sqrt(p_cov[3,3]))
print("A_2(212Pb): ", p_opt[4], "+/-", np.sqrt(p_cov[4,4]))
print("A_2(224Ra): ", p_opt[5], "+/-", np.sqrt(p_cov[5,5]))


# from the initial parameters, calculate the 224Ra at the transition
a_ra224_t_1 = p_opt[1] * np.exp(-lambdas[1]*t_1)
a_ra224_t_1_err = np.sqrt(p_cov[1,1]) * np.exp(-lambdas[1]*t_1)

a_ra224_t_2 = p_opt[3] * np.exp(-lambdas[1]*(t_2-t_1))
a_ra224_t_2_err = np.sqrt(p_cov[3,3]) * np.exp(-lambdas[1]*(t_2-t_1))

print("A_t1(224Ra): ", a_ra224_t_1, "+/-", a_ra224_t_1_err)
print("A_t2(224Ra): ", a_ra224_t_2, "+/-", a_ra224_t_2_err)

# Calculate the reduction factor from interval 0 -> 1
f_red_0_1 = a_ra224_t_1 / p_opt[3]
f_red_0_1_err = f_red_0_1 * np.sqrt((np.sqrt(p_cov[3,3])/p_opt[3])**2 + (a_ra224_t_1_err/a_ra224_t_1)**2)
print("R(0->1): ", f_red_0_1, "+/-", f_red_0_1_err)

# Calculate the reduction factor from interval 1 -> 2
f_red_1_2 = p_opt[5] / a_ra224_t_2
f_red_1_2_err = f_red_1_2 * np.sqrt((np.sqrt(p_cov[5,5])/p_opt[5])**2 + (a_ra224_t_1_err/a_ra224_t_1)**2)
print("R(1->2): ", f_red_1_2, "+/-", f_red_1_2_err)

fig = plt.figure(figsize=(12,8), dpi=80) 
gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)


ax1.plot(x_vals, a_fit)
ax1.errorbar((bin_centers-t_min)/3600./24., bin_counts/bins_width, yerr = np.sqrt(bin_counts+1)/bins_width,
		fmt='.', color='gray', linewidth=1, markersize=1)

ax1.plot(np.array((t_0,t_0))/3600./24, (-1, 5), '--',label=r"$\mathrm{T_0}$(implantation stop)")
ax1.plot(np.array((t_1,t_1))/3600./24, (-1, 5), '--',label=r"$\mathrm{T_1}$(fake coating)")
ax1.plot(np.array((t_2,t_2))/3600./24, (-1, 5), '--',label=r"$\mathrm{T_2}$(de-coating)")

ax1.legend(loc=0)
# Add residual plot


ax2.plot((bin_centers-t_min)/3600./24., residuals, '.', color='gray', markersize=2)
ax2.grid(axis='y')

ax1.set_ylabel("Rate (1/s)")
ax2.set_ylabel("Residuals")
ax2.set_xlabel("Time (days)")
ax1.set_title("224Ra implanted stainless steel sample #2")

ax1.set_ylim(5e-4,15)
ax1.set_yscale('log')
ax2.set_ylim(-5,5)


for item in [ax1.title, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels() + ax1.get_legend().get_texts():
	item.set_fontsize(20)  
for item in [ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels():
	item.set_fontsize(20)  

plt.subplots_adjust(wspace=0, hspace=0)

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
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-9,1e1)
plt.show()
'''
