import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import scipy.optimize as opt 


# Define where the .root files are stored on the LFS1
base_path = "/d12/lin/xenon/radon_monitor/database/backup/"
file_ending = ".root"

t0 = 0
n_days = 70
file_names = ['Mn09052023']

x_axis_lim = np.array([-0.5, n_days])
y_axis_act_lim = np.array([1e-5, 1e-2])


# define time intervals to cut
time_cut = 1e9
time_cut = time_cut * 3600*24

po212_selection = [300, 400]
re_bin = 60 * 60 * 10

bin_centers = np.array(())
bin_counts = np.array(())
start_times = np.zeros((len(file_names)))
stop_times = np.zeros((len(file_names)))


for i in range(len(file_names)):
	
	this_file_path = base_path + file_names[i] + "/" + file_names[i].lower() + file_ending
	print(this_file_path)
	this_tree = up.open(this_file_path)["t"]

	time_stamp = this_tree['timestamp'].array(library='np')
	channel = this_tree['channel'].array(library='np')

	start_times[i] = np.min(time_stamp)
	stop_times[i] = np.max(time_stamp)
	duration = stop_times[i] - start_times[i]
	print(duration)

	event_mask = np.less(channel, po212_selection[1]) * np.greater(channel, po212_selection[0])
	selected_timestamps = time_stamp[event_mask]

	bins = np.arange(int(start_times[i]), int(stop_times[i]), re_bin)
	bins_width = np.mean(bins[1:] - bins[:-1])

	bin_conts, bin_edges = np.histogram(selected_timestamps, bins=bins)
	bin_cent = (bin_edges[:-1] + bin_edges[1:])/2.

	print(bin_edges[1:] - bin_edges[:-1])

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
n_isotopes = len(labels)

def calc_matrix():

	t_1_2 = np.array((	1.9116*365*24*3600,
						3.6319*24*3600,
						55.6,
	                	145e-3,
						10.64*3600,
	                	60.55*60,
	                	299e-9))

	lambdas = np.log(2)/t_1_2
	# Multiply the 212Bi Half-live with the branching fraction
	po212_branching = 0.6406
	lambdas[-2] *= po212_branching

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
	
	return lambdas, A_eig_vals, V, V_inv

lambdas, A_eig_vals, V, V_inv = calc_matrix()


def get_n_t(t, a0):

	lambdas, A_eig_vals, V, V_inv = calc_matrix()

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

def evolution_fit_function(t, a_0_228th, a_0_224ra, a_0_212pb):

	a0_0 = np.zeros((n_isotopes))
	a0_0[0] = a_0_228th
	a0_0[1:3] = a_0_224ra
	a0_0[4:] = a_0_212pb
	
	return get_n_t(t-t0*24*3600, a0_0)[:,-1] * lambdas[-1]

def exp_228_decay(t, a0, t_0=0):
    # Calcualte simple decay of 228Thorium with a shift on the time axis of t_0
    return a0*np.exp(-(t-t_0)*np.log(2)/(1.8*365))

def exp_224_decay(t, a0, t_0=0):
    # Calcualte simple decay of 224Radium with a shift on the time axis of t_0
    return a0*np.exp(-(t-t_0)*np.log(2)/3.6319)

def exp_212_decay(t, a0, t_0=0):
    # Calcualte simple decay of 212Lead with a shift on the time axis of t_0
    return a0*np.exp(-(t-t_0)*np.log(2)/(12./24.))

interval_rejection_mask = (bin_centers-t_min) < time_cut

print(interval_rejection_mask)
n_last_point = 300

bin_centers_cut = bin_centers[interval_rejection_mask]
bin_counts_cut = bin_counts[interval_rejection_mask]
bins_width_cut = bins_width

# define initial values
# A_228Th, A_224Ra, A_212Po
init_vals = np.array((1e-4, 1e-1, 1))
p_opt, p_cov = opt.curve_fit(evolution_fit_function, (bin_centers_cut-t_min), bin_counts_cut/bins_width_cut,
				sigma = np.sqrt(bin_counts_cut+1)/bins_width_cut, p0 = init_vals, bounds=[0, np.inf],
				absolute_sigma = True)

# export datapoints for oleg
exp_data = np.stack(((bin_centers_cut-t_min)/3600./24, bin_counts_cut), axis=1)
np.savetxt("implanted_sample_2_212po.dat", exp_data, header="time [days] po212 counts per 2 minutes")

print(p_opt)
print(p_cov)

# Get array of x-values and corresponding fit values
x_vals = np.linspace(0, (t_max-bin_centers[0])/3600/24, 10000)
a_fit = evolution_fit_function(x_vals*3600*24, *p_opt)

# also calculate the residuals between fit and data
residuals = (bin_counts/bins_width - evolution_fit_function((bin_centers-t_min), *p_opt))/(np.sqrt(bin_counts+1)/bins_width)

# Initial activities at the start of measurement
a228_0 = np.array([p_opt[0], np.sqrt(p_cov[0,0])])
a224_0 = np.array([p_opt[1], np.sqrt(p_cov[1,1])])
a212_0 = np.array([p_opt[2], np.sqrt(p_cov[2,2])])

print("A_0(228Th)(t=0):", a228_0[0], "+/-", a228_0[1])
print("A_0(224Ra)(t=0):", a224_0[0], "+/-", a224_0[1])
print("A_0(212Pb)(t=0):", a212_0[0], "+/-", a212_0[1])


fig = plt.figure(figsize=(16,8), dpi=80) 
gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)


ax1.errorbar((bin_centers-t_min)/3600./24., bin_counts/bins_width, 
			yerr = np.sqrt(bin_counts)/bins_width, xerr=bins_width/3600./24. * 0.5,
			fmt='.', color='black', linewidth=3, markersize=1)
ax1.plot(x_vals, a_fit, color="tab:red", linewidth=3)

# Add lines to indicate the equilibrium value of the 224Ra activity also draw
# shaded boxes to indicate the uncertainty on the activity
act_line_cols = 3*['tab:gray']
# act_line_cols = ['tab:pink', 'tab:green', 'tab:purple']
act_line_width = 3
act_line_style = 'solid'

ax1.plot(x_axis_lim, exp_228_decay(x_axis_lim, a228_0[0], t0), "--")
ax1.plot(x_axis_lim, exp_224_decay(x_axis_lim, a224_0[0], t0), "--")
ax1.plot(x_axis_lim, exp_212_decay(x_axis_lim, a212_0[0], t0), "--")


ax2.plot((bin_centers-t_min)/3600./24., residuals, '.', color='black', markersize=10)
# ax2.grid(axis='y')
ax2.fill_between([-1,n_days*1.1], 2*[-1], 2*[1], color="tab:red", alpha=0.5)
ax2.fill_between([-1,n_days*1.1], 2*[-2], 2*[2], color="tab:red", alpha=0.25)

ax1.set_ylabel("Detected $^{212}$Po Rate (Hz)")
ax1.get_xaxis().set_visible(False)
ax2.set_ylabel("Residual")
ax2.set_xlabel("Measurement time (days)")

ax1.set_yscale("log")
ax1.set_ylim(y_axis_act_lim)
ax2.set_ylim(-3.9,3.9)
ax2.set_xlim(x_axis_lim)

for item in [ax1.title, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
	item.set_fontsize(25)  
for item in [ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels():
	item.set_fontsize(25)  

fig.tight_layout()

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig("spiking_test_2",dpi="figure")
plt.show()

