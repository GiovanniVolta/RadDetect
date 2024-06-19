import requests
from bs4 import BeautifulSoup
import urllib.request
import uproot
import tempfile
import os
import inspect

import urllib
from iminuit import cost, Minuit

import numpy as np
import scipy as sp
from scipy.stats import crystalball

import matplotlib.pyplot as plt
import matplotlib as mpl

class BlumchenAnalysis:
    
    def __init__(self, url, filename):
        # Get branches
        self.mca_ch, self.timestamp, self.runtime = self.get_data(url, filename)
        
    def get_data(self, url, filename):
        # Get available ROOT file
        # It should be one
        _url = f'{url}/{filename}/'
        root_file = _scrape_radon_db(_url)
        
        if len(root_file) != 1:
            raise ValueError(f"Expected 1 elements in 'root_file', but got {len(root_file)}.")
    
        # Download, process, and delete the ROOT file
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            
            tmp_filename = tmp_file.name + '.root'
            urllib.request.urlretrieve(root_file[0], tmp_filename)
            
            with uproot.open(tmp_filename) as _file:
                tree_name = _file.keys()[0]
                tree = _file[tree_name]
                mca_ch = tree["channel"].array().to_numpy()
                timestamp = tree["timestamp"].array().to_numpy()
                runtime = tree["runtime"].array().to_numpy()
                # Process data as needed
                print(f"Retriving data from {root_file[0]}")

        os.remove(tmp_filename)
        return mca_ch, timestamp, runtime
    
    def get_mca_histogram(self, MCA_range=[910, 1060], time_range=[240, np.inf], n_channels=None):
        
        print('Selection:\n'
            f'\t-MCA range: {MCA_range}\n'
            f'\t-time range: {time_range} minutes')
        
        mask = (self.mca_ch >= MCA_range[0]) & (self.mca_ch <= MCA_range[1])
        mask &= (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_mca_ch = self.mca_ch[mask]
        
        # Energy/MCA channel histogram
        ch_min = np.min(selected_mca_ch)
        ch_max = np.max(selected_mca_ch)
        if n_channels == None:
            channel_bins = np.linspace(ch_min, ch_max, int(ch_max - ch_min))
        else:
            channel_bins = np.linspace(ch_min, ch_max, n_channels)        
        data, _ = np.histogram(selected_mca_ch, bins=channel_bins)
        channels = 0.5*(channel_bins[1:] + channel_bins[:-1])
        return data, channels
    
    def get_time_evolution(self, MCA_range=[910, 1060], time_range=[240, np.inf], n_timestamp=None):
        
        print('Selection:\n'
            f'\t-MCA range: {MCA_range}\n'
            f'\t-time range: {time_range} minutes')
        
        # Time evolution histogram
        # of defined MCA range and time selection
        # time in minutes
        mask = (self.mca_ch >= MCA_range[0]) & (self.mca_ch <= MCA_range[1])
        mask &= (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_runtime = self.runtime[mask]
        
        t_min = np.min(selected_runtime)
        t_max = np.max(selected_runtime)
        
        if n_timestamp == None:
            time_bins = np.linspace(t_min, t_max, int(t_max - t_min))
        else:
            time_bins = np.linspace(t_min, t_max, n_timestamp)
        
        data_time_evolution, _ = np.histogram(selected_runtime, bins=time_bins)
        dt = np.diff(time_bins)  # bin width
        times = 0.5 * (time_bins[:-1] + time_bins[1:])
        rate = data_time_evolution / dt
        rate_err = np.sqrt(data_time_evolution) / dt
        # retuns in s, Hzm Hz
        return times * 60, rate * 60, rate_err * 60
    
    def get_base_plot(self, MCA_range=[910, 1060], time_range=[0, np.inf], n_channels=None, n_timestamp=None):
        
        data, channels = self.get_mca_histogram(MCA_range=[0, 1300], time_range=time_range, n_channels=n_channels)
        times, rate, rate_err = self.get_time_evolution(MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp)
        
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
        axs = axs.flatten()
        
        axs[0].plot(channels, data, ds='steps', color='black')
        axs[0].fill_between(channels, data, step='mid', color='black', alpha=0.3)
        axs[0].set_yscale('log')
        axs[0].set_xlim(0, 1300)
        axs[0].set_xlabel('MCA channel')
        axs[0].set_ylabel('Counts')
        axs[0].grid()
        
        axs[1].scatter(self.runtime, self.mca_ch, s=3, color = 'black', alpha=0.3)
        axs[1].set_ylabel('MCA channel')
        axs[1].set_xlabel('Runtime [minutes]')
        axs[1].set_ylim(0, 1300)
        axs[1].grid()
        
        axs[2].errorbar(times / 60 / 60 / 24, rate, yerr=rate_err, lw=0, color = 'black', marker='o', ms=2, elinewidth=1, label=f'Time evolution in {MCA_range} MCA ch')
        axs[2].set_xlabel('Times [day]')
        axs[2].set_ylabel('Rate [Hz]')
        axs[2].grid()
        axs[2].legend()

        plt.tight_layout()
        plt.show()
        
    def get_mca_spectrum_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[910, 1060], time_range=[240, np.inf], MCA_counts_limit=5, n_channels=None):        
        
        data, channels = self.get_mca_histogram(MCA_range=MCA_range, time_range=time_range, n_channels=n_channels)
        
        mask = (data > MCA_counts_limit)
        _data = data[mask]
        _channels = channels[mask]
        _init = _prepare_init_for_iminuit(model, init)
        cost_function = cost.LeastSquares(_channels, _data, np.sqrt(_data), model)
        m = Minuit(cost_function, *_init)
        if limits != None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed != None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m
    
    def get_time_evolution_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[910, 1060], time_range=[240, np.inf], n_timestamp=None, rate_limit=5):
        
        times, rate, rate_err = self.get_time_evolution(MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp)
        
        mask = (rate > 5)
        _times = times[mask]
        _rate = rate[mask]
        _rate_err = rate_err[mask]
        
        _init = _prepare_init_for_iminuit(model, init)
        
        cost_function = cost.LeastSquares(_times, _rate, _rate_err, model)
        m = Minuit(cost_function, *_init)
        if limits != None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed != None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m
    
    
def _scrape_radon_db(url):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith('.root')]
    
def _prepare_init_for_iminuit(Model, init):
    parameter_names = inspect.signature(Model).parameters.keys()
    _init = []
    for _p in parameter_names:
        if _p == 'x':
            pass
        else:
            _init.append(init[_p])
    return _init