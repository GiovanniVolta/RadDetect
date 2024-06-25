import requests
from bs4 import BeautifulSoup
import urllib.request
import uproot
import tempfile
import os
import inspect
from iminuit import cost, Minuit
import numpy as np
import matplotlib.pyplot as plt

class MonalphaAnalysis:
    """
    A class for analyzing data from ROOT files, specifically designed for Monalpha analysis.
    
    This class provides methods to retrieve data from a specified URL, generate histograms of MCA channel data,
    plot time evolution of the data, and fit the data using specified models.
    
    Attributes:
        mca_ch (numpy.ndarray): MCA channel data.
        timestamp (numpy.ndarray): Timestamp data.
        runtime (numpy.ndarray): Runtime data.
    """
    
    def __init__(self, file_path):
        """
        Initializes the MonalphaAnalysis class by retrieving data from the specified file path.
        
        Args:
            file_path (str): The path to the ROOT file.
        """
        self.mca_ch, self.timestamp, self.runtime = self.get_data(file_path)

    def get_data(self, file_path):
        """
        Retrieves data from a ROOT file located at the specified file path.
        
        Args:
            file_path (str): The path to the ROOT file.
            
        Returns:
            tuple: A tuple containing arrays of MCA channel data, timestamp data, and runtime data.
        """
        if not os.path.exists(file_path):
            raise ValueError(f"File does not exist at the specified path: {file_path}")

        with uproot.open(file_path) as _file:
            tree = _file[_file.keys()[0]]
            mca_ch = tree["channel"].array().to_numpy()
            timestamp = tree["timestamp"].array().to_numpy()
            runtime = tree["runtime"].array().to_numpy()
            print(f"Retrieving data from {file_path}")

        return mca_ch, timestamp, runtime

    def get_mca_histogram(self, MCA_range=[650, 705], time_range=[240, np.inf], n_channels=None):
        """
        Generates a histogram of MCA channel data within a specified range and time range.
        
        Args:
            MCA_range (list, optional): The range of MCA channels to include. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            n_channels (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            
        Returns:
            tuple: A tuple containing the histogram data and the channel bins.
        """
        # Apply filters based on MCA range and time range
        mask = (self.mca_ch >= MCA_range[0]) & (self.mca_ch <= MCA_range[1]) & \
            (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_mca_ch = self.mca_ch[mask]

        # Generate histogram
        ch_min, ch_max = selected_mca_ch.min(), selected_mca_ch.max()
        channel_bins = np.linspace(ch_min, ch_max, n_channels or int(ch_max - ch_min))
        data, channels = np.histogram(selected_mca_ch, bins=channel_bins)
        channels = 0.5 * (channel_bins[1:] + channel_bins[:-1])
        return data, channels



    def get_base_plot(self, MCA_range=[910, 1060], time_range=[0, np.inf], n_channels=None, n_timestamp=None):
        """
        Generates and displays base plots for the MCA channel data, including a histogram, scatter plot, and error bar plot.
        
        Args:
            MCA_range (list, optional): The range of MCA channels for the time evolution plot. Defaults to [910, 1060].
            time_range (list, optional): The range of time for the plots, in minutes. Defaults to [0, np.inf].
            n_channels (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            n_timestamp (int, optional): The number of timestamps for the time evolution plot. If None, the number is determined automatically.
        """
        # Generate data for plots
        data, channels = self.get_mca_histogram(MCA_range=[400, 1300], time_range=time_range, n_channels=n_channels)
        
        
        # Plotting
        fig, axs = plt.subplots(1, 2, figsize=(18, 5), dpi=150)
        axs = axs.flatten()
        
        # MCA histogram
        axs[0].plot(channels, data, ds='steps', color='black')
        axs[0].fill_between(channels, data, step='mid', color='black', alpha=0.3)
        axs[0].set_yscale('log')
        axs[0].set_xlim(400, 1300)
        axs[0].set_xlabel('MCA channel')
        axs[0].set_ylabel('Counts')
        axs[0].grid()
        
        # Scatter plot of runtime vs MCA channel
        axs[1].scatter(self.runtime, self.mca_ch, s=3, color='black', alpha=0.3)
        axs[1].set_ylabel('MCA channel')
        axs[1].set_xlabel('Runtime [minutes]')
        axs[1].set_ylim(400, 1300)
        axs[1].grid()
        

        plt.tight_layout()
        plt.show()

    def get_mca_spectrum_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[650, 705], time_range=[420, np.inf], MCA_counts_limit=5, n_channels=None):
        """
        Prepares a fitting object for the MCA spectrum using the specified model and initial parameters.
        
        This method filters the MCA channel data based on the specified MCA range, time range, and counts limit,
        and then prepares a Minuit object for fitting the filtered data to the specified model.
        
        Args:
            model (callable): The model function to fit the data to. Must be compatible with iminuit.
            init (dict): Initial parameter values for the fitting model.
            limits (dict, optional): Parameter limits for the fitting model. Defaults to None.
            fixed (dict, optional): Parameters to be held fixed during the fit. Defaults to None.
            MCA_range (list, optional): The range of MCA channels to include in the fit. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            MCA_counts_limit (int, optional): The minimum number of counts required for a channel to be included in the fit. Defaults to 5.
            n_channels (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            
        Returns:
            iminuit.Minuit: A Minuit object configured for fitting the specified model to the MCA spectrum data.
        """
        data, channels = self.get_mca_histogram(MCA_range=MCA_range, time_range=time_range, n_channels=n_channels)
        
        mask = (data > MCA_counts_limit)
        _data = data[mask]
        _channels = channels[mask]
        _init = self._prepare_init_for_iminuit(model, init)
        cost_function = cost.ExtendedBinnedNLL(_data, _channels, model)
        m = Minuit(cost_function, *_init)
        if limits is not None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed is not None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m

    @staticmethod
    def _prepare_init_for_iminuit(model, init):
        parameter_names = inspect.signature(model).parameters.keys()
        return [init[p] for p in parameter_names if p != 'x']