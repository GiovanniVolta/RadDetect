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
from scipy.optimize import curve_fit

import warnings

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
    
    def __init__(self, filename, energy_calibration=None):
        """
        Initializes the MonalphaAnalysis class by retrieving data from the specified file path.
        
        Args:
            filename (str): The path or the name the ROOT file.
            energy_calibration (None or list): The energy calibration parameters. It can be either None or a list [q, m] 
            where q and m are the linear energy calibration parameters.
        """
        self.energy_calibration = energy_calibration
        if self.energy_calibration is not None:
            warnings.warn("energy_calibration is not None. this will affect the MCA range selection. Be caferul!")
                        
        self.mca, self.timestamp, self.runtime = self.get_data(filename)

    def get_data(self, filename):
        """
        Retrieves data from a ROOT file located at the specified URL or locally.

        Args:
            filename (str): The name of the ROOT file to retrieve.

        Returns:
            tuple: A tuple containing arrays of MCA channel data, timestamp data, and runtime data.
        """
        if not filename.endswith('.root'):
            local_path = filename + '.root'
        else:
            local_path = filename
        url = "https://radon-srv1.mpi-hd.mpg.de/coating_db/resultfiles"
        
        # Check if the file exists locally
        if os.path.exists(local_path):
            file_to_use = local_path
            print(f"Using local file: {local_path}")
        else:
            # If the file does not exist locally, download it
            root_file = self._scrape_radon_db(f'{url}/{filename}/')
            if len(root_file) != 1:
                raise ValueError(f"Expected 1 element in 'root_file', got {len(root_file)}.")

            with tempfile.NamedTemporaryFile(suffix='.root', delete=False) as tmp_file:
                tmp_filename = tmp_file.name
                urllib.request.urlretrieve(root_file[0], tmp_filename)
                file_to_use = tmp_filename
                print(f"Retrieving data from {root_file[0]}")

        # Open the file using uproot
        with uproot.open(file_to_use) as _file:
            tree = _file[_file.keys()[0]]
            mca = tree["channel"].array().to_numpy()
            timestamp = tree["timestamp"].array().to_numpy()
            runtime = tree["runtime"].array().to_numpy()

        # If a temporary file was used, remove it
        if not os.path.exists(local_path):
            os.remove(tmp_filename)

        if self.energy_calibration is not None:
            mca = (mca - self.energy_calibration[1])/self.energy_calibration[0]
        
        return mca, timestamp, runtime

    def get_mca_histogram(self, MCA_range=[600, 715], time_range=[0, np.inf], n_mca=None):
        """
        Generates a histogram of MCA channel data within a specified range and time range.
        
        Args:
            MCA_range (list, optional): The range of MCA channels to include. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            n_mca (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            
        Returns:
            tuple: A tuple containing the histogram data and the channel bins.
        """
        # Apply filters based on MCA range and time range
        mask = (self.mca >= MCA_range[0]) & (self.mca <= MCA_range[1]) & \
            (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_mca = self.mca[mask]

        # Generate histogram
        mca_min, mca_max = selected_mca.min(), selected_mca.max()
        mca_bins = np.linspace(mca_min, mca_max, n_mca or int(mca_max - mca_min))
        data, _ = np.histogram(selected_mca, bins=mca_bins)
        mcas = 0.5 * (mca_bins[1:] + mca_bins[:-1])
        return data, mcas

    def get_time_evolution(self, MCA_range=[600, 715], time_range=[0, np.inf], n_timestamp=None):
        """
        Generates data for the time evolution of the MCA channel data within a specified range and time range.
        
        Args:
            MCA_range (list, optional): The range of MCA channels to include. Defaults to [910, 1060].
            time_range (list, optional): The range of time to include, in minutes. Defaults to [240, np.inf].
            n_timestamp (int, optional): The number of timestamps for the histogram. If None, the number is determined automatically.
            
        Returns:
            tuple: A tuple containing the times, rate, and rate error in seconds and hertz
        """
        # Apply filters based on MCA range and time range
        mask = (self.mca >= MCA_range[0]) & (self.mca <= MCA_range[1]) & \
            (self.runtime > time_range[0]) & (self.runtime < time_range[1])
        selected_runtime = self.runtime[mask]

        # Generate time evolution data
        t_min, t_max = selected_runtime.min(), selected_runtime.max()
        time_bins = np.linspace(t_min, t_max, n_timestamp or int(t_max - t_min))
        data_time_evolution, _ = np.histogram(selected_runtime, bins=time_bins)
        # in second and Hz
        dt = np.diff(time_bins) * 60
        times = 0.5 * (time_bins[:-1] + time_bins[1:]) * 60
        rate = data_time_evolution / dt
        rate_err = np.sqrt(data_time_evolution) / dt
        return times, rate, rate_err

    def get_base_plot(self, MCA_range=[600, 715], time_range=[0, np.inf], n_mca=None, n_timestamp=None):
        """
        Generates and displays base plots for the MCA channel data, including a histogram, scatter plot, and error bar plot.
        
        Args:
            MCA_range (list, optional): The range of MCA channels for the time evolution plot. Defaults to [910, 1060].
            time_range (list, optional): The range of time for the plots, in minutes. Defaults to [0, np.inf].
            n_mca (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            n_timestamp (int, optional): The number of timestamps for the time evolution plot. If None, the number is determined automatically.
        """
        # Plotting
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
        axs = axs.flatten()
        
        # Plot property that depends on the energy_calibration
        if self.energy_calibration is not None:
            label = f'Time evolution in {MCA_range} keV'
            axs[1].set_ylabel('Energy [keV]')
            axs[0].set_xlabel('Energy [keV]')
            _MCA_range = [(0 - self.energy_calibration[1])/self.energy_calibration[0], 
                        (1300 - self.energy_calibration[1])/self.energy_calibration[0]
                        ]
        else:
            label = f'Time evolution in {MCA_range} MCA ch'
            axs[1].set_ylabel('MCA channel')
            axs[0].set_xlabel('MCA channel')
            _MCA_range = [0, 1300]
            
        data, mcas = self.get_mca_histogram(MCA_range=_MCA_range, time_range=time_range, n_mca=n_mca)
        times, rate, rate_err = self.get_time_evolution(MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp)
                    
        # MCA histogram
        axs[0].plot(mcas, data, ds='steps', color='black')
        axs[0].fill_between(mcas, data, step='mid', color='black', alpha=0.3)
        axs[0].axvspan(*MCA_range, color='pink', lw=0, alpha=0.5)
        axs[0].set_yscale('log')
        axs[0].set_xlim(_MCA_range)                        
        axs[0].set_ylabel('Counts')
        axs[0].grid()
        
        # Scatter plot of runtime vs MCA channel
        axs[1].scatter(self.runtime, self.mca, s=3, color='black', alpha=0.3)
        axs[1].axhspan(*MCA_range, color='pink', lw=0, alpha=0.5)            
        axs[1].set_xlabel('Runtime [minutes]')
        axs[1].set_ylim(_MCA_range)    
        axs[1].grid()
        
        # Time evolution error bar plot           
        axs[2].errorbar(times / 60 / 60 / 24, rate, yerr=rate_err, lw=0, color='black', marker='o', ms=2, elinewidth=1, label=label)
        axs[2].set_xlabel('Times [day]')
        axs[2].set_ylabel('Rate [Hz]')
        axs[2].grid()
        axs[2].legend()
        
        plt.tight_layout()
        plt.show()

    def get_mca_spectrum_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[600, 715], time_range=[0, np.inf], 
                                        MCA_counts_limit=5, n_mca=None, prefit=True):
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
            n_mca (int, optional): The number of channels for the histogram. If None, the number is determined automatically.
            
        Returns:
            iminuit.Minuit: A Minuit object configured for fitting the specified model to the MCA spectrum data.
        """
        data, mcas = self.get_mca_histogram(MCA_range=MCA_range, time_range=time_range, n_mca=n_mca)
        
        mask = (data > MCA_counts_limit)
        _data = data[mask]
        _mcas = mcas[mask]
        
        _parameter_names, _init = self._prepare_init_for_fit(model, init)
        
        if prefit:
            print('Prefit with scipy for deriving inital values')
            _bounds = self._prepare_bounds_for_fit(model, init, fixed, limits)
            _init, _init_cov = curve_fit(model.total_model, _mcas, _data, sigma=np.sqrt(_data), 
                                        absolute_sigma=True, p0=_init, maxfev = 500000, bounds=_bounds)
            _init_err = np.sqrt(np.diag(_init_cov))
            self.print_table(_parameter_names, _init, _init_err)
        
        cost_function = cost.LeastSquares(_mcas, _data, np.sqrt(_data), model.total_model)
        m = Minuit(cost_function, *_init)
        if limits is not None:
            for k in limits.keys():
                m.limits[k] = limits[k]
        if fixed is not None:
            for k in fixed.keys():
                m.fixed[k] = fixed[k]            
        return m

    def get_time_evolution_fitting_object(self, model, init, limits=None, fixed=None, 
                                        MCA_range=[600, 715], time_range=[0, np.inf], n_timestamp=None, rate_limit=0, prefit=True):
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
        
        _parameter_names, _init = self._prepare_init_for_fit(model, init)
        
        if prefit:
            print('Prefit with scipy for deriving inital values')
            _bounds = self._prepare_bounds_for_fit(model, init, fixed, limits)
            _init, _init_cov = curve_fit(model.total_model, _times, _rate, sigma=_rate_err, absolute_sigma=True, p0=_init, maxfev = 500000, bounds=_bounds)
            _init_err = np.sqrt(np.diag(_init_cov))
            self.print_table(_parameter_names, _init, _init_err)
                
        cost_function = cost.LeastSquares(_times, _rate, _rate_err, model.total_model)
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
        # Sends a GET request to the specified URL and parses the content using BeautifulSoup.
        page = requests.get(url).text
        soup = BeautifulSoup(page, 'html.parser')
        # Returns a list of URLs, appending the href attribute of each <a> tag that ends with '.root' to the base URL.
        return [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith('.root')]

    @staticmethod
    def _prepare_init_for_fit(Model, init):
        # Retrieves the names of the parameters (excluding 'self') of the Model's constructor.
        parameter_names = list(inspect.signature(Model).parameters.keys())[1:]
        # Returns a list of initial values for these parameters based on the provided 'init' dictionary.
        return parameter_names, [init[p] for p in parameter_names]
    
    @staticmethod
    def _prepare_bounds_for_fit(model, init, fixed, limits, lower_bound=-np.inf, upper_bound=np.inf, epsilon=1e-9):
        # Retrieves the names of the parameters (excluding 'x') of the model's function.
        parameter_names = list(inspect.signature(model).parameters.keys())[1:]
        bounds = []
        
        # prepare limits and fixed.
        if limits == None:
            _limits = {}
        else:
            _limits = limits
        
        if fixed == None:
            _fixed = {}
        else:
            _fixed = fixed
            
        # Creates bounds for each parameter based on whether it is fixed (using 'epsilon' for tight bounds) 
        # or if it is specified in limits 
        # or free (using provided 'lower_bound' and 'upper_bound')
        for p in parameter_names:
            if _fixed.get(p, False):
                bounds.append((init[p] - epsilon, init[p] + epsilon))
            elif p in _limits.keys():
                bounds.append((_limits[p][0], _limits[p][1]))
            else:
                bounds.append((lower_bound, upper_bound))
        # Separates the bounds into two lists: one for lower bounds and one for upper bounds.
        lower_bounds, upper_bounds = zip(*bounds)  # Unzip the pairs into two lists
        return lower_bounds, upper_bounds
    
    @staticmethod
    def print_table(_parameter_names, _init, _init_err):
        print(f'{"| Parameters":12} | {"Value":11} | {"Error":11} |')
        print('-' * 42)
        for (l, i, j) in zip(_parameter_names, _init, _init_err):
            print(f'| {l:10} | {i:11.3f} | {j:11.3f} |')
        print('-' * 42)