import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import scipy.optimize as opt 


# Define where the .root files are stored on the LFS1
base_path = "/d12/lin/xenon/radon_monitor/database/backup/"
file_ending = ".root"

# Toggle between the plot showing true rate or the 
# rate with a correction that cancels out the overall
# 224Ra decay (only affects plotting, not the underlying
# analysis!)
plot_corrected_rate = True
# plot_corrected_rate = False

# Define runs and diferent intervals

t0 = 0
runs_uncoated = ['Rn18012021','Rn23012021']
t1 = 6.88
runs_coated = ['Rn25012021', 'Rn30012021']
t2 = 14.033
# runs_decoated = ['Rn29122020', 'Rn04012021']

# The constant towards which the ratio of the two isotopes converges
ra224_po212_ratio = 0.867618

x_axis_lim = np.array([-0.5, 14])

if plot_corrected_rate:
    y_axis_act_lim = np.array([0.2, 23])
else:
    y_axis_act_lim = np.array([3.5e-3, 8])

n_days = 21
time_since_implantation = 91*24*3600
time_since_implantation = 0

# define time intervals to cut
time_cut = 1009.1 
time_cut = time_cut * 3600*24

file_names = runs_uncoated+runs_coated


po212_selection = [220,320]
bin_width = int(1 * 60. *2)
re_bin = 60 * 60

bin_centers = np.array(())
bin_counts = np.array(())
start_times = np.zeros((len(file_names)))
stop_times = np.zeros((len(file_names)))


for i in range(len(file_names)):
	
	this_file_path = base_path + file_names[i] + "/" + file_names[i] + file_ending
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


def calc_matrix(t_1_2_mods):


	t_1_2 = np.array((	1.9116*365*24*3600,
				3.6319*24*3600,
				55.6,
	                	145e-3,
				10.64*3600,
	                	60.55*60,
	                	299e-9))
	t_1_2 *= t_1_2_mods

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

lambdas, A_eig_vals, V, V_inv = calc_matrix(np.ones((len(labels))))


def get_n_t(t, a0, mods):

	lambdas, A_eig_vals, V, V_inv = calc_matrix(mods)

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

	return get_n_t(t, a0, mods)[:,-1] * lambdas[-1]

# def get_full_po213_evolution(t, a_0_224ra, a_0_212pb, mod_224ra, mod_220rn, mod_212pb, mod_212bi):
# def get_full_po213_evolution(t, a_0_224ra, a_0_212pb,
#                                 a_1_224ra, a_1_212pb,
#                                 a_2_224ra, a_2_212pb,
# 				mod_224ra=1, mod_220rn=1, 
# 				mod_212pb=1, mod_212bi=1):

def get_full_po213_evolution(t, a_0_224ra, a_0_212pb,
                                a_1_224ra,
                                a_2_224ra,
				mod_224ra=1, mod_220rn=1, 
				mod_212pb=1, mod_212bi=1):

	# Convert input time from seconds to days
	t = t*(24*3600) + time_since_implantation

	# This function will allow to fit the 212Po activity for different
	# Data taking conditions (i.e. coated, uncoated, de-coated)
	# in a combined fashion

	modificators = np.ones((len(labels)))
	modificators[1] = mod_224ra
	modificators[2] = mod_220rn
	modificators[4] = mod_212pb
	modificators[5] = mod_212bi

	lambdas, A_eig_vals, V, V_inv = calc_matrix(modificators)

	# Start with initializing the vector of initial conditions to be zero
	# for all isotopes.
	# Now fill the ones, which we like to fit and/or insert constraints

	# The initial vector a_0_1 will be initial conditions for first interval
	# We assume equilibrium of activity between 224Ra -> 216Po
	# As well as between 212Pb->212Po. So we have two degrees of freedom


	# Interval 0	
	a0_0 = np.zeros((n_isotopes))
	a0_0[0] = 0
	a0_0[1:3] = a_0_224ra
	a0_0[4:] = a_0_212pb

	n_t_0 = get_n_t(t-t0*24*3600, a0_0, modificators)

	# Interval 1

	# compute the activity from interval zero at t1
	n_t_0_to_1 = get_n_t((t1*24*3600,), a0_0, modificators)
	# And use it as initial activity
	a1_0 = n_t_0_to_1[0,:] * lambdas
	# Now set the new Ra224 activity in this interval
	a1_0[1] = a_1_224ra

	n_t_1 = get_n_t(t-t1*24*3600, a1_0, modificators)

	# Interval 2

	# compute the activity from interval one at t2
	n_t_1_to_2 = get_n_t(((t2-t1)*24*3600,), a1_0, modificators)
	# And use it as initial activity
	a2_0 = n_t_1_to_2[0,:] * lambdas
	# Now set the new Ra224 activity in this interval
	a2_0[1] = a_2_224ra
	
	n_t_2 = get_n_t(t-t2*24*3600, a2_0, modificators)

	# Now stitch together the evolutions from the different intervals
	t_0_mask = t < t1*24*3600
	t_1_mask = np.logical_and(t >= t1*24*3600, t < t2*24*3600)
	t_2_mask = t >= t2*24*3600

	n_t = np.concatenate((n_t_0[t_0_mask], n_t_1[t_1_mask], n_t_2[t_2_mask]), axis=0)

	a_t_212po = n_t[:,6] * lambdas[6]
	return a_t_212po


def exp_224_decay(t, a0, t_0=0):
    # Calcualte simple decay of 224Radium with a shift on the time axis of t_0
    return a0*np.exp(-(t-t_0)*np.log(2)/3.6319)

interval_rejection_mask = (bin_centers-t_min) < time_cut

print(interval_rejection_mask)
n_last_point = 300

bin_centers_cut = bin_centers[interval_rejection_mask]
bin_counts_cut = bin_counts[interval_rejection_mask]
bins_width_cut = bins_width

# define initial values
init_vals = np.array((1,1e-5,1,1))
# init_vals = np.array((1,1))
p_opt, p_cov = opt.curve_fit(get_full_po213_evolution, (bin_centers_cut-t_min)/3600./24, bin_counts_cut/bins_width_cut,
				sigma = np.sqrt(bin_counts_cut+1)/bins_width_cut, p0 = init_vals, bounds=[0, np.inf],
				absolute_sigma = True)

# export datapoints for oleg
exp_data = np.stack(((bin_centers_cut-t_min)/3600./24, bin_counts_cut), axis=1)
np.savetxt("implanted_sample_2_212po.dat", exp_data, header="time [days] po212 counts per 2 minutes")

print(p_opt)
print(p_cov)

# Get array of x-values and corresponding fit values
x_vals = np.linspace(0, (t_max-bin_centers[0])/3600/24, 10000)
a_fit = get_full_po213_evolution(x_vals, *p_opt)

# also calculate the residuals between fit and data
residuals = (bin_counts/bins_width - get_full_po213_evolution((bin_centers-t_min)/3600./24, *p_opt))/(np.sqrt(bin_counts+1)/bins_width)

# Initial radium activities at the respective interval boundary
a224_0 = np.array([p_opt[0], np.sqrt(p_cov[0,0])]) / ra224_po212_ratio
a224_1 = np.array([p_opt[2], np.sqrt(p_cov[2,2])]) / ra224_po212_ratio
a224_2 = np.array([p_opt[3], np.sqrt(p_cov[3,3])]) / ra224_po212_ratio

# Calculate the fitted 224Ra activities at time zero
a224_0_zero = exp_224_decay(0, a224_0[0], t0)
a224_1_zero = exp_224_decay(0, a224_1[0], t1)
a224_2_zero = exp_224_decay(0, a224_2[0], t2)
a224_0_zero_err = exp_224_decay(0, a224_0[1], t0)
a224_1_zero_err = exp_224_decay(0, a224_1[1], t1)
a224_2_zero_err = exp_224_decay(0, a224_2[1], t2)


print("A_0(224Ra)(t=0):", a224_0_zero, "+/-", a224_0_zero_err)
print("A_0(212Pb)(t=0):", p_opt[1], "+/-", np.sqrt(p_cov[1,1])) # not correct, if t0 != 0
print("A_1(224Ra)(t=0):", a224_1_zero, "+/-", a224_1_zero_err)
print("A_2(224Ra)(t=0):", a224_2_zero, "+/-", a224_2_zero_err)

# Calculate reduction factors and their uncertainties (neglect correlation!)
r_1 = a224_0_zero/a224_1_zero
r_1_err = r_1 * np.sqrt((a224_0_zero_err/a224_0_zero)**2 + (a224_1_zero_err/a224_1_zero)**2)
r_2 = a224_2_zero/a224_1_zero
r_2_err = r_2 * np.sqrt((a224_2_zero_err/a224_2_zero)**2 + (a224_1_zero_err/a224_1_zero)**2)

print("R_1:", r_1, "+/-", r_1_err)
print("R_2:", r_2, "+/-", r_2_err)


fig = plt.figure(figsize=(12,8), dpi=80) 
gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[3, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# define a correction factor that cancels out the overall decay of 224Ra in the
# sample to make the plot easier to comprehend (suggested by Hardy)
# Also apply this factor to the rate and fit function (if desired by user)

dec_corr = 1/(exp_224_decay((bin_centers-t_min)/3600./24., 1))
dec_corr_rate = bin_counts/bins_width * dec_corr 
dec_corr_err = np.sqrt(bin_counts)/bins_width * dec_corr

if plot_corrected_rate:
    ax1.errorbar((bin_centers-t_min)/3600./24., dec_corr_rate, 
                  yerr = dec_corr_err, xerr=bins_width/3600./24. * 0.5,
                  fmt='.', color='black', linewidth=3, markersize=1, zorder=4)
    ax1.plot(x_vals, a_fit*1/exp_224_decay(x_vals, 1), 
             color="tab:red", linewidth=3, zorder=3)
else:
    ax1.errorbar((bin_centers-t_min)/3600./24., bin_counts/bins_width, 
                  yerr = np.sqrt(bin_counts)/bins_width, xerr=bins_width/3600./24. * 0.5,
                  fmt='.', color='black', linewidth=3, markersize=1)
    ax1.plot(x_vals, a_fit, color="tab:red", linewidth=3)




# Add indicator boxes with labels for the different phases
indicator_box_y = [y_axis_act_lim[0], 0.35]
indicator_box_color = ["gray", "darkgray", "lightgray"]
indicator_box_text = ["Before submersion", "After submersion"]  
indicator_box_range = [(t0, t1), (t1,t2), (t2,x_axis_lim[1])]
                                                                                
for idx, this_text in enumerate(indicator_box_text):                            
    ax1.fill_between([indicator_box_range[idx][0],indicator_box_range[idx][1]], 
                     indicator_box_y[0], indicator_box_y[1], color=indicator_box_color[idx])
    box_center_x = 0.5*(indicator_box_range[idx][0] + indicator_box_range[idx][1])
    box_center_y = np.sqrt(indicator_box_y[0] * indicator_box_y[1])             
    ax1.text(box_center_x, box_center_y, this_text, ha="center", va="center",   
             fontsize=25)


# Draw arrows indicating the nomenclature of the reduction factors
text_offset = 0.5
arrow_offset_1 = 2
arrow_offset_2 = 4

arrowcol="darkorange"
arrowwidth = 6

ax1.text(t1-arrow_offset_1-text_offset, np.sqrt(a224_0_zero*a224_1_zero), "$\mathrm{R}$", 
         ha='center', va='center', rotation=90, size=25+2, fontweight='semibold')
ax1.annotate("", xy=(t1-arrow_offset_1, a224_1_zero), xytext=(t1-arrow_offset_1,a224_0_zero),
            arrowprops=dict(facecolor=arrowcol, width=arrowwidth,linewidth=0), 
            ha='center', zorder=2)

# Add lines to indicate the equilibrium value of the 224Ra activity also draw
# shaded boxes to indicate the uncertainty on the activity
act_line_cols = 3*['tab:gray']
# act_line_cols = ['tab:pink', 'tab:green', 'tab:purple']
act_line_width = 3
act_line_style = 'solid'
if plot_corrected_rate:
    ax1.hlines(a224_0_zero, t0, t1, colors=act_line_cols[0], zorder=1, 
               linewidth=act_line_width, linestyle=act_line_style) 
    ax1.hlines(a224_1_zero, t1-arrow_offset_1-text_offset, t1,
               colors=act_line_cols[1], zorder=1,
               linewidth=act_line_width, linestyle="--") 
    ax1.hlines(a224_1_zero, t2, t2+arrow_offset_2+text_offset, 
               colors=act_line_cols[1], zorder=1,
               linewidth=act_line_width, linestyle="--") 
    ax1.hlines(a224_1_zero, t1, t2, colors=act_line_cols[1], zorder=1,
               linewidth=act_line_width, linestyle=act_line_style) 
    ax1.hlines(a224_2_zero, t2, x_axis_lim[1], colors=act_line_cols[2], zorder=1,
               linewidth=act_line_width, linestyle=act_line_style) 

#    ax1.hlines([a224_0_zero, a224_1_zero, a224_2_zero], x_axis_lim[0], x_axis_lim[1],
#               linestyles="--", zorder=1)
else:
    ax1.plot(x_axis_lim, exp_224_decay(x_axis_lim, a224_0[0], t0), "--")
    ax1.plot(x_axis_lim, exp_224_decay(x_axis_lim, a224_1[0], t1), "--")
    ax1.plot(x_axis_lim, exp_224_decay(x_axis_lim, a224_2[0], t2), "--")

# Also add vertical lines to indicate the different periods of the mesurement
time_line_cols = ["tab:blue", "tab:red", "tab:orange"]
ax1.vlines([t0,t1,t2], 1e-10, 1e10, colors=time_line_cols, linewidth=2, linestyle="dotted")

# Add residual plot


ax2.plot((bin_centers-t_min)/3600./24., residuals, '.', color='black', markersize=10)
# ax2.grid(axis='y')
ax2.fill_between([-1,n_days*1.1], 2*[-1], 2*[1], color="tab:red", alpha=0.5)
ax2.fill_between([-1,n_days*1.1], 2*[-2], 2*[2], color="tab:red", alpha=0.25)

if plot_corrected_rate:
    ax1.set_ylabel("Corrected $^{212}$Po Rate (Hz)")
else:
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

plt.savefig("implanted_ss_2_initial_ra224_1.pdf",dpi="figure")
plt.show()

