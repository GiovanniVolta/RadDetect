import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt

file_names = [	"2020_02_21_15_33_done.txt",
		"2020_02_24_15_34_done.txt",
		"2020_03_02_15_15_done.txt",
		"2020_03_09_16_12_done.txt",
		"2020_03_16_15_07_done.txt"]

for this_file in file_names:
	this_data = pd.read_csv(this_file, header=0, sep=" ", names=["date", "time", "temp"])
	
	date = this_data["date"]
	time = this_data["time"]
	form_str = '%m/%d/%y %H:%M:%S'
	date_times = date + " " + time
	datetimes = [dt.datetime.strptime(b, form_str) for b in date_times]
	temp = this_data["temp"]

	print(np.mean(temp), np.std(temp))
	plt.plot(datetimes, temp)
	plt.show()


