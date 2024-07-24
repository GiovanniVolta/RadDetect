import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from datetime import datetime
import matplotlib.patches as mpatches



#216Po
base_path = "/d12/lin/xenon/radon_monitor/database/backup/"
file_ending = ".root"


file_names = ['Ap17082020', 'Ap20082020', 'Ap04092020_cooling', 'Ap04092020','Ap09092020_heating', 'Ap10092020','Ap15092020',
              'Ap24092020','Ap19102020','Ap26102020_heating','Ap26102020','Ap05112020','Ap06112020']
temperatures = np.array((20., 20., -100, -30., -100, -5, 20.,20.,20.,-100.,75.,-100.,20.))
colors = ["gold", "gold", "black", "blue", "black", "green", "gold","gold","gold","black","red","black","gold"] 

trees = []
time = []

po216_selection = [200,225]
bin_width = 200 * 60.

plt.figure(figsize=(11, 4.5))

for i in range(len(file_names)):
    this_file_path =  base_path + file_names[i] + "/" + file_names[i] + file_ending
    #print(this_file_path)
    this_tree = up.open(this_file_path)["t"]
    trees.append(this_tree)
    #this_tree.show()

    time_stamp = this_tree.array("timestamp")
    runtime = this_tree.array("runtime")
    channel = this_tree.array("channel")
    t_min = int(np.min(time_stamp))
    t_max = int(np.max(time_stamp))
    
        
    time.append(t_min)
    time.append(t_max)
    t_min_dt = dt.datetime.fromtimestamp(t_min)
    t_max_dt = dt.datetime.fromtimestamp(t_max)
    
    duration = t_max - t_min
    bins = 1
    bins_width = duration
    if np.round(duration / bin_width) > 0:
        bins, bins_width = np.linspace(int(t_min), int(t_max), int(np.round(int(duration) / int(bin_width))), retstep=True)
    event_mask = np.less(channel, po216_selection[1]) * np.greater(channel, po216_selection[0])


    selected_timestamps = time_stamp[event_mask]
    bin_conts, bin_edges = np.histogram(selected_timestamps, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    bin_edges_dt = [dt.datetime.fromtimestamp(b) for b in bin_edges]
    bin_centers_dt = [dt.datetime.fromtimestamp(b) for b in bin_centers]
    if (t_min < 1599753402  or t_max > 1600178898 ):
        plt.errorbar(bin_centers_dt,bin_conts/bins_width * 1000., yerr=np.sqrt(bin_conts)/bins_width *1000 , fmt='.', color='navy')
        plt.axvspan(t_min_dt, t_max_dt, alpha=0.25, color=colors[i])

    else:
        bin_centers_dt_corrected = []
        bin_centers_dt_heating_5_20 = []
        for l in range(len(bin_centers_dt)): 
            if bin_centers_dt[l]<datetime(2020, 9, 14, 15, 30, 4):
                bin_centers_dt_corrected.append(bin_centers_dt[l])
            else:
                bin_centers_dt_heating_5_20.append(bin_centers_dt[l])
        
        new = np.array(bin_centers_dt_corrected)
        new2 = bin_conts[:len(bin_centers_dt_corrected)]/bins_width * 1000.
        new3 = bin_conts[len(bin_conts)-7:]/bins_width * 1000.
      
        plt.errorbar(new,new2, yerr=np.sqrt(bin_conts[:len(bin_centers_dt_corrected)])/bins_width *1000 , fmt='.', color='navy')
        plt.errorbar(bin_centers_dt_heating_5_20,new3,yerr=np.sqrt(bin_conts[len(bin_conts)-7:])/bins_width *1000 , fmt='.', color='navy')
        plt.axvspan(min(new), max(new), alpha=0.25, color='green')
        plt.axvspan(min(bin_centers_dt_heating_5_20), max(bin_centers_dt_heating_5_20), alpha=0.25, color='black')
    
    
    
    
    

gold_patch = mpatches.Patch(color='gold', alpha=0.25,label='20'+r'$^\circ$C')
blue_patch = mpatches.Patch(color='blue', alpha=0.25,label='-30'+r'$^\circ$C')
green_patch = mpatches.Patch(color='green',alpha=0.25, label='-5'+r'$^\circ$C')
red_patch = mpatches.Patch(color='red',alpha=0.25, label='75'+r'$^\circ$C')


plt.xticks(rotation=25)
plt.ylim(0,4.5)
plt.title("$\mathrm{{}^{216}Po}$ Activity of DLC coated rods")
plt.ylabel(r"Detected $\mathrm{{}^{216}Po}$ Activity in mBq")
#plt.legend( bbox_to_anchor=(1.05, 1),loc="upper left")
#plt.legend(handles=[gold_patch,blue_patch,green_patch],loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(handles=[gold_patch,blue_patch,green_patch,red_patch],loc="upper right")
#plt.legend(handles=[p1, p2], title='title', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
#plt.legend()
plt.tight_layout()
plt.savefig('Po216_DLC_temp_2.png',dpi=900)
plt.show()


#212Po

file_ending = ".root"

file_names = ['Ap17082020', 'Ap20082020', 'Ap04092020_cooling', 'Ap04092020','Ap09092020_heating', 'Ap10092020','Ap15092020',
              'Ap24092020','Ap19102020','Ap26102020_heating','Ap26102020','Ap05112020','Ap06112020']
temperatures = np.array((20., 20., -100, -30., -100, -5, 20.,20.,20.,-100.,75.,-100.,20.))
colors = ["gold", "gold", "black", "blue", "black", "green", "gold","gold","gold","black","red","black","gold"] 

trees = []
time = []

po216_selection = [250,310]
bin_width = 200 * 60.

plt.figure(figsize=(11, 4.5))

for i in range(len(file_names)):
    this_file_path =  base_path + file_names[i] + "/" + file_names[i] + file_ending
    #print(this_file_path)
    this_tree = up.open(this_file_path)["t"]
    trees.append(this_tree)
    #this_tree.show()

    time_stamp = this_tree.array("timestamp")
    runtime = this_tree.array("runtime")
    channel = this_tree.array("channel")
    t_min = int(np.min(time_stamp))
    t_max = int(np.max(time_stamp))
    
        
    time.append(t_min)
    time.append(t_max)
    t_min_dt = dt.datetime.fromtimestamp(t_min)
    t_max_dt = dt.datetime.fromtimestamp(t_max)
    
    duration = t_max - t_min
    bins = 1
    bins_width = duration
    if np.round(duration / bin_width) > 0:
        bins, bins_width = np.linspace(int(t_min), int(t_max), int(np.round(int(duration) / int(bin_width))), retstep=True)
    event_mask = np.less(channel, po216_selection[1]) * np.greater(channel, po216_selection[0])


    selected_timestamps = time_stamp[event_mask]
    bin_conts, bin_edges = np.histogram(selected_timestamps, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    bin_edges_dt = [dt.datetime.fromtimestamp(b) for b in bin_edges]
    bin_centers_dt = [dt.datetime.fromtimestamp(b) for b in bin_centers]
    if (t_min < 1599753402  or t_max > 1600178898 ):
        plt.errorbar(bin_centers_dt,bin_conts/bins_width * 1000., yerr=np.sqrt(bin_conts)/bins_width *1000 , fmt='.', color='navy')
        plt.axvspan(t_min_dt, t_max_dt, alpha=0.25, color=colors[i])

    else:
        bin_centers_dt_corrected = []
        bin_centers_dt_heating_5_20 = []
        for l in range(len(bin_centers_dt)): 
            if bin_centers_dt[l]<datetime(2020, 9, 14, 15, 30, 4):
                bin_centers_dt_corrected.append(bin_centers_dt[l])
            else:
                bin_centers_dt_heating_5_20.append(bin_centers_dt[l])
        
        new = np.array(bin_centers_dt_corrected)
        new2 = bin_conts[:len(bin_centers_dt_corrected)]/bins_width * 1000.
        new3 = bin_conts[len(bin_conts)-7:]/bins_width * 1000.
      
        plt.errorbar(new,new2, yerr=np.sqrt(bin_conts[:len(bin_centers_dt_corrected)])/bins_width *1000 , fmt='.', color='navy')
        plt.errorbar(bin_centers_dt_heating_5_20,new3,yerr=np.sqrt(bin_conts[len(bin_conts)-7:])/bins_width *1000 , fmt='.', color='navy')
        plt.axvspan(min(new), max(new), alpha=0.25, color='green')
        plt.axvspan(min(bin_centers_dt_heating_5_20), max(bin_centers_dt_heating_5_20), alpha=0.25, color='black')
    
    
    
    
    

gold_patch = mpatches.Patch(color='gold', alpha=0.25,label='20'+r'$^\circ$C')
blue_patch = mpatches.Patch(color='blue', alpha=0.25,label='-30'+r'$^\circ$C')
green_patch = mpatches.Patch(color='green',alpha=0.25, label='-5'+r'$^\circ$C')
red_patch = mpatches.Patch(color='red',alpha=0.25, label='75'+r'$^\circ$C')



plt.xticks(rotation=25)
plt.ylim(0,10.5)
plt.title("$\mathrm{{}^{212}Po}$ Activity of DLC coated rods")
plt.ylabel(r"Detected $\mathrm{{}^{212}Po}$ Activity in mBq")
#plt.legend( bbox_to_anchor=(1.05, 1),loc="upper left")
#plt.legend(handles=[gold_patch,blue_patch,green_patch],loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(handles=[gold_patch,blue_patch,green_patch,red_patch],loc="upper right")
#plt.legend(handles=[p1, p2], title='title', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
#plt.legend()
plt.tight_layout()
plt.savefig('Po212_DLC_temp_2.png',dpi=900)
plt.show()
