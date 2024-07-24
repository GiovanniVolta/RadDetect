import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

base_path = "/d12/lin/xenon/radon_monitor/database/backup/"
file_ending = ".root"

file_names = ['Ap21022020', 'Ap24022020_preheatin', 'Ap24022020', 'Ap02032020_cooling', 'Ap02032020']
temperatures = np.array((20., -1, 150., -1, 20.))
colors = ["white", "gray", "red", "gray", "white"] 

trees = []

po216_selection = [200,230]
bin_width = 500 * 60.

plt.figure(figsize=(6, 4.5))

for i in range(len(file_names)):
	if i == 0:
		pass
 		#continue
	this_file_path = base_path + file_names[i] + "/" + file_names[i] + file_ending
	print(this_file_path)
	this_tree = up.open(this_file_path)["t"]
	trees.append(this_tree)
	this_tree.show()

	time_stamp = this_tree.array("timestamp")
	channel = this_tree.array("channel")
	t_min = np.min(time_stamp)
	t_max = np.max(time_stamp)
	duration = t_max - t_min
	bins = 1
	bins_width = duration
	if np.round(duration / bin_width) > 0:
		bins, bins_width = np.linspace(t_min, t_max, np.round(duration / bin_width), retstep=True)
	event_mask = np.less(channel, po216_selection[1]) * np.greater(channel, po216_selection[0])

	selected_timestamps = time_stamp[event_mask]
	bin_conts, bin_edges = np.histogram(selected_timestamps, bins=bins)
	bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
	bin_edges_dt = [dt.datetime.fromtimestamp(b) for b in bin_edges]
	bin_centers_dt = [dt.datetime.fromtimestamp(b) for b in bin_centers]
	print(bin_conts, bin_edges, bin_centers)
	plt.errorbar(bin_centers_dt,bin_conts/bins_width * 1000., yerr=np.sqrt(bin_conts)/bins_width *1000 , fmt='.', color='navy')
	print(t_min, t_max)
	t_min_dt = dt.datetime.fromtimestamp(t_min)
	t_max_dt = dt.datetime.fromtimestamp(t_max)
	if temperatures[i] > 20:
		plt.axvspan(t_min_dt, t_max_dt, alpha=0.25, color=colors[i], label='{:3.0f}'.format(temperatures[i]) +r'$^\circ$C')
	else:
		plt.axvspan(t_min_dt, t_max_dt, alpha=0.5, color=colors[i])
	

	#break

plt.legend(title=r'Temperature', loc="upper right")
plt.xticks(rotation=25)
plt.title("Sample 2")
plt.ylabel(r"Detected $\mathrm{{}^{216}Po}$ Activity in mBq")
plt.tight_layout()
plt.show()
