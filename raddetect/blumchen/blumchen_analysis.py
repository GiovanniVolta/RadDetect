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

class BlumchenAnalysis:
    """
    A class for analyzing data from ROOT files, specifically designed for Blumchen analysis.
    
    This class provides methods to retrieve data from a specified URL, generate histograms of MCA channel data,
    plot time evolution of the data, and fit the data using specified models.
    
    Attributes:
        mca_ch (numpy.ndarray): MCA channel data.
        timestamp (numpy.ndarray): Timestamp data.
        runtime (numpy.ndarray): Runtime data.
    """
    
    def __init__(self, url, filename):
        """
        Initializes the BlumchenAnalysis class by retrieving data from the specified URL and filename.
        
        Args:
            url (str): The URL where the ROOT file is located.
            filename (str): The name of the ROOT file to retrieve.
        """
        self.mca_ch, self.timestamp, self.runtime = self.get_data(url, filename)

    def get_data(self, url, filename):
        """
        Retrieves data from a ROOT file located at the specified URL and filename.
        
        Args:
            url (str): The URL where the ROOT file is located.
            filename (str): The name of the ROOT file to retrieve.
            
        Returns:
            tuple: A tuple containing arrays of MCA channel data, timestamp data, and runtime data.
        """
        root_file = self._scrape_radon_db(f'{url}/{filename}/')
        if len(root_file) != 1:
            raise ValueError(f"Expected 1 element in 'root_file', got {len(root_file)}.")

        with tempfile.NamedTemporaryFile(suffix='.root', delete=False) as tmp_file:
            tmp_filename = tmp_file.name
            urllib.request.urlretrieve(root_file[0], tmp_filename)

            with uproot.open(tmp_filename) as _file:
                tree = _file[_file.keys()[0]]
                mca_ch = tree["channel"].array().to_numpy()
                timestamp = tree["timestamp"].array().to_numpy()
                runtime = tree["runtime"].array().to_numpy()
                print(f"Retrieving data from {root_file[0]}")

        os.remove(tmp_filename)
        return mca_ch, timestamp, runtime

    def get_mca_histogram(self, MCA_range=[910, 1060], time_range=[240, np.inf], n_channels=None):
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

    def get_time_evolution(self, MCA_range=[910, 1060], time_range=[240, np.inf], n_timestamp=None):
        """
        Generates data for the time evolution of the MCA channel data within a specified range and time range.
        
        Args:
            MCA_range (list, optional): The range of MCA channels to include. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            n_timestamp (int, optional): The number of timestamps for the histogram. If None, the number is determined automatically.
            
        Returns:
            tuple: A tuple containing the times, rate, and rate error.
        """
        # Apply filters based on MCA range and time range
        mask = (self.mca_ch >= MCA_range[0]) & (self.mca_ch <= MCA_range[1]) & \
            (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_runtime = self.runtime[mask]

        # Generate time evolution data
        t_min, t_max = selected_runtime.min(), selected_runtime.max()
        time_bins = np.linspace(t_min, t_max, n_timestamp or int(t_max - t_min))
        data_time_evolution, _ = np.histogram(selected_runtime, bins=time_bins)
        dt = np.diff(time_bins)
        times = 0.5 * (time_bins[:-1] + time_bins[1:])
        rate = data_time_evolution / dt
        rate_err = np.sqrt(data_time_evolution) / dt
        return times * 60, rate * 60, rate_err * 60

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
        data, channels = self.get_mca_histogram(MCA_range=[0, 1300], time_range=time_range, n_channels=n_channels)
        times, rate, rate_err = self.get_time_evolution(MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp)
        
        # Plotting
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
        axs = axs.flatten()
        
        # MCA histogram
        axs[0].plot(channels, data, ds='steps', color='black')
        axs[0].fill_between(channels, data, step='mid', color='black', alpha=0.3)
        axs[0].set_yscale('log')
        axs[0].set_xlim(0, 1300)
        axs[0].set_xlabel('MCA channel')
        axs[0].set_ylabel('Counts')
        axs[0].grid()
        
        # Scatter plot of runtime vs MCA channel
        axs[1].scatter(self.runtime, self.mca_ch, s=3, color='black', alpha=0.3)
        axs[1].set_ylabel('MCA channel')
        axs[1].set_xlabel('Runtime [minutes]')
        axs[1].set_ylim(0, 1300)
        axs[1].grid()
        
        # Time evolution error bar plot
        axs[2].errorbar(times / 60 / 60 / 24, rate, yerr=rate_err, lw=0, color='black', marker='o', ms=2, elinewidth=1, label=f'Time evolution in {MCA_range} MCA ch')
        axs[2].set_xlabel('Times [day]')
        axs[2].set_ylabel('Rate [Hz]')
        axs[2].grid()
        axs[2].legend()

        plt.tight_layout()
        plt.show()

    def get_mca_spectrum_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[910, 1060], time_range=[240, np.inf], MCA_counts_limit=5, n_channels=None):
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
        cost_function = cost.LeastSquares(_channels, _data, np.sqrt(_data), model)
        m = Minuit(cost_function, *_init)
        if limits is not None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed is not None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m

    def get_time_evolution_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[910, 1060], time_range=[240, np.inf], n_timestamp=None, rate_limit=5):
        """
        Prepares a fitting object for the time evolution data using the specified model and initial parameters.
        
        This method filters the time evolution data based on the specified MCA range, time range, and rate limit,
        and then prepares a Minuit object for fitting the filtered data to the specified model.
        
        Args:
            model (callable): The model function to fit the data to. Must be compatible with iminuit.
            init (dict): Initial parameter values for the fitting model.
            limits (dict, optional): Parameter limits for the fitting model. Defaults to None.
            fixed (dict, optional): Parameters to be held fixed during the fit. Defaults to None.
            MCA_range (list, optional): The range of MCA channels to include in the fit. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            n_timestamp (int, optional): The number of timestamps for the histogram. If None, the number is determined automatically.
            rate_limit (int, optional): The minimum rate required for a time point to be included in the fit. Defaults to 5.
            
        Returns:
            iminuit.Minuit: A Minuit object configured for fitting the specified model to the time evolution data.
        """
        times, rate, rate_err = self.get_time_evolution(MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp)
        
        mask = (rate > rate_limit)
        _times = times[mask]
        _rate = rate[mask]
        _rate_err = rate_err[mask]
        
        _init = self._prepare_init_for_iminuit(model, init)
        
        cost_function = cost.LeastSquares(_times, _rate, _rate_err, model)
        m = Minuit(cost_function, *_init)
        if limits is not None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed is not None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m

    @staticmethod
    def _scrape_radon_db(url):
        page = requests.get(url).text
        soup = BeautifulSoup(page, 'html.parser')
        return [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith('.root')]

    @staticmethod
    def _prepare_init_for_iminuit(Model, init):
        parameter_names = inspect.signature(Model).parameters.keys()
        return [init[p] for p in parameter_names if p != 'x']
    