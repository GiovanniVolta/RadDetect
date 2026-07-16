import inspect
import os
import tempfile
import urllib.request
import warnings

import matplotlib.pyplot as plt
import numpy as np
import requests
import uproot
from bs4 import BeautifulSoup
from iminuit import Minuit, cost
from scipy.optimize import curve_fit


class RadonAnalysis:
    """A base class for analyzing data from ROOT files.

    This class provides methods to retrieve data from a specified URL or local file,
    generate histograms of MCA channel data,
    plot time evolution of the data, and fit the data using specified models.

    Attributes:
        mca (numpy.ndarray): MCA channel data.
        timestamp (numpy.ndarray): Timestamp data in seconds.
        runtime (numpy.ndarray): Runtime data in minutes.
    """

    # Used in the full spectrum plot
    DEFAULT_MCA_RANGE = [0, 1300]
    DEFAULT_TIME_RANGE = [0, np.inf]
    
    # To select a specific line to plot the temporal evolution
    SELECTED_MCA_RANGE = [0, 1300]
    SELECTED_TIME_RANGE = [0, np.inf]
    
    def __init__(
        self,
        filename,
        energy_calibration=None,
        compute_runtime_from_timestamp=False,
        timestamp_interval=60,
    ):
        """Initializes the RadonAnalysis class by retrieving data from the specified
        file path or URL.

        Args:
            filename (str): The path or the name the ROOT file.
            energy_calibration (None or list): The energy calibration parameters. It can be either
            None or a list [q, m]
            where q and m are the linear energy calibration parameters.
            compute_runtime_from_timestamp (bool, optional): If True, computes runtime from the timestamp.
            Defaults to False.
            timestamp_interval (int, optional): The list file logging interval in seconds. Defaults to 60.
        """
        self.energy_calibration = energy_calibration
        self.compute_runtime_from_timestamp = compute_runtime_from_timestamp
        self.timestamp_interval = timestamp_interval
        if self.energy_calibration is not None:
            warnings.warn(
                "energy_calibration is not None. this will affect the MCA range selection. "
                "Be careful!"
            )

        self.mca, self.timestamp, self.runtime = self.get_data(filename)

    def get_data(self, filename):
        """Retrieves data from a ROOT file located at the specified URL or locally.

        Args:
            filename (str): The name of the ROOT file to retrieve.

        Returns:
            tuple: A tuple containing arrays of MCA channel data, timestamp data (in seconds),
            and runtime data (in minutes).
        """
        if not filename.endswith(".root"):
            local_path = filename + ".root"
        else:
            local_path = filename
        url = "https://radon-srv1.mpi-hd.mpg.de/coating_db/resultfiles"

        # Check if the file exists locally
        if os.path.exists(local_path):
            file_to_use = local_path
            print(f"Using local file: {local_path}")
        else:
            print("Using radon_db")
            self._radon_db_is_reachable(url)
            # If the file does not exist locally, download it
            root_file = self._scrape_radon_db(f"{url}/{filename}/")
            if len(root_file) != 1:
                raise ValueError(
                    f"Expected 1 element in 'root_file', got {len(root_file)}."
                )

            with tempfile.NamedTemporaryFile(suffix=".root", delete=False) as tmp_file:
                tmp_filename = tmp_file.name
                urllib.request.urlretrieve(root_file[0], tmp_filename)
                file_to_use = tmp_filename
                print(f"Retrieving data from {root_file[0]}")

        # Open the file using uproot
        with uproot.open(file_to_use) as _file:
            tree = _file[_file.keys()[0]]
            mca = tree["channel"].array().to_numpy()
            timestamp = tree["timestamp"].array().to_numpy()
            if self.compute_runtime_from_timestamp:
                # timestamp is in UNIX time in seconds and runtime is in minutes, so we divide by 60.
                # timestamp_interval is added because the list file is generated every interval (e.g. 60 seconds).
                runtime = (timestamp - timestamp[0] + self.timestamp_interval) / 60
            else:
                runtime = tree["runtime"].array().to_numpy()

        # If a temporary file was used, remove it
        if not os.path.exists(local_path):
            os.remove(tmp_filename)

        if self.energy_calibration is not None:
            mca = (mca - self.energy_calibration[1]) / self.energy_calibration[0]

        return mca, timestamp, runtime

    def get_mca_histogram(self, MCA_range=None, time_range=None, n_mca=None):
        """Generates a histogram of MCA channel data within a specified range and time
        range. If nothing is specified, the spectrum will be done using the DEFUALT 
        intervals.

        Args:
            MCA_range (list, optional): The range of MCA channels to include. Defaults
                                        to class DEFAULT_MCA_RANGE.
            time_range (list, optional): The range of time to include, in minutes.
                                        Defaults to class DEFAULT_TIME_RANGE.
            n_mca (int, optional): The number of channels for the histogram.
                                If None, the number is determined automatically.

        Returns:
            tuple: A tuple containing the histogram data and the channel bins.
        """
        MCA_range = MCA_range if MCA_range is not None else self.DEFAULT_MCA_RANGE
        time_range = time_range if time_range is not None else self.DEFAULT_TIME_RANGE

        # print('get_mca_histogram')
        # print(MCA_range)
        # print(time_range)
        
        # Apply filters based on MCA range and time range
        mask = (
            (self.mca >= MCA_range[0])
            & (self.mca <= MCA_range[1])
            & (self.runtime > time_range[0])
            & (self.runtime < time_range[1])
        )
        selected_mca = self.mca[mask]

        # Generate histogram
        if n_mca is None:
            # Bins are exactly 1 integer wide. 
            # +2 ensures the upper bound is fully included as a bin edge.
            mca_bins = np.arange(MCA_range[0], MCA_range[1] + 2)
        else:
            # If a specific number of bins is requested, calculate edges evenly
            mca_bins = np.linspace(MCA_range[0], MCA_range[1], n_mca + 1)
            
        data, _ = np.histogram(selected_mca, bins=mca_bins)
        mcas = 0.5 * (mca_bins[1:] + mca_bins[:-1])
        return data, mcas

    def get_time_evolution(self, MCA_range=None, time_range=None, n_timestamp=None):
        """Generates data for the time evolution of the MCA channel data within a
        specified range and time range. If nothing is specified, the spectrum will 
        be done using the SELECTED intervals.

        Args:
            MCA_range (list, optional): The range of MCA channels to include.
                                        Defaults to class SELECTED_MCA_RANGE.
            time_range (list, optional): The range of time to include, in minutes.
                                        Defaults to class SELECTED_TIME_RANGE.
            n_timestamp (int, optional): The number of timestamps for the histogram.
                                        If None, the number is determined automatically.

        Returns:
            tuple: A tuple containing the times, rate, and rate error in seconds and hertz
        """
        MCA_range = MCA_range if MCA_range is not None else self.SELECTED_MCA_RANGE
        time_range = time_range if time_range is not None else self.SELECTED_TIME_RANGE

        # print('get_time_evolution')
        # print(MCA_range)
        # print(time_range)
        
        # Apply filters based on MCA range and time range
        mask = (
            (self.mca >= MCA_range[0])
            & (self.mca <= MCA_range[1])
            & (self.runtime > time_range[0])
            & (self.runtime < time_range[1])
        )
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

    def get_base_plot(
        self, MCA_range=None, time_range=None, n_mca=None, n_timestamp=None
    ):
        """Generates and displays base plots for the MCA channel data, including a
        histogram, scatter plot, and error bar plot.

        Args:
            MCA_range (list, optional): The range of MCA channels for the time evolution plot.
                                        Defaults to class SELECTED_MCA_RANGE.
            time_range (list, optional): The range of time for the plots, in minutes.
                                        Defaults to class SELECTED_TIME_RANGE.
            n_mca (int, optional): The number of channels for the histogram.
                                    If None, the number is determined automatically.
            n_timestamp (int, optional): The number of timestamps for the time evolution plot.
                                    If None, the number is determined automatically.
        """
        MCA_range = MCA_range if MCA_range is not None else self.SELECTED_MCA_RANGE
        time_range = time_range if time_range is not None else self.SELECTED_TIME_RANGE

        # print('get_base_plot')
        # print(MCA_range)
        # print(time_range)
        
        # Plotting
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=150)
        axs = axs.flatten()

        # Plot property that depends on the energy_calibration
        if self.energy_calibration is not None:
            label = f"Time evolution in {MCA_range} keV"
            axs[1].set_ylabel("Energy [keV]")
            axs[0].set_xlabel("Energy [keV]")
            _MCA_range = [
                (self.DEFAULT_MCA_RANGE[0] - self.energy_calibration[1]) / self.energy_calibration[0],
                (self.DEFAULT_MCA_RANGE[1] - self.energy_calibration[1]) / self.energy_calibration[0],
            ]
        else:
            label = f"Time evolution in {MCA_range} MCA ch"
            axs[1].set_ylabel("MCA channel")
            axs[0].set_xlabel("MCA channel")
            _MCA_range = self.DEFAULT_MCA_RANGE

        data, mcas = self.get_mca_histogram(
            MCA_range=_MCA_range, time_range=time_range, n_mca=n_mca
        )
        times, rate, rate_err = self.get_time_evolution(
            MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp
        )

        # MCA histogram
        axs[0].plot(mcas, data, ds="steps", color="black")
        axs[0].fill_between(mcas, data, step="mid", color="black", alpha=0.3)
        axs[0].axvspan(*MCA_range, color="pink", lw=0, alpha=0.5)
        axs[0].set_yscale("log")
        axs[0].set_xlim(_MCA_range)
        axs[0].set_ylabel("Counts")
        axs[0].grid()

        # Scatter plot of runtime vs MCA channel
        axs[1].scatter(self.runtime, self.mca, s=3, color="black", alpha=0.3)
        axs[1].axvspan(0, time_range[0], color="grey", lw=0, alpha=0.5)
        if time_range[1] != float('inf') and time_range[1] is not None:
            axs[1].axvspan(time_range[1], max(self.runtime), color="grey", lw=0, alpha=0.5)
        axs[1].axhspan(*MCA_range, color="pink", lw=0, alpha=0.5)
        axs[1].set_xlabel("Runtime [minutes]")
        axs[1].set_ylim(_MCA_range)
        axs[1].grid()

        # Time evolution error bar plot
        axs[2].errorbar(
            times / 60 / 60 / 24,
            rate,
            yerr=rate_err,
            lw=0,
            color="black",
            marker="o",
            ms=2,
            elinewidth=1,
            label=label,
        )
        axs[2].set_xlabel("Times [day]")
        axs[2].set_ylabel("Rate [Hz]")
        axs[2].grid()
        axs[2].legend()

        plt.tight_layout()
        plt.show()
        
    def get_mca_spectrum_fitting_object(
        self,
        model,
        init,
        limits=None,
        fixed=None,
        MCA_range=None,
        time_range=None,
        MCA_counts_limit=5,
        n_mca=None,
        prefit=True,
    ):
        MCA_range = MCA_range if MCA_range is not None else self.SELECTED_MCA_RANGE
        time_range = time_range if time_range is not None else self.SELECTED_TIME_RANGE

        data, mcas = self.get_mca_histogram(
            MCA_range=MCA_range, time_range=time_range, n_mca=n_mca
        )

        # Truncating data biases the tails. Consider dropping this limit in the future
        # if you switch to a Poisson Maximum Likelihood cost function.
        mask = data > MCA_counts_limit
        _data = data[mask]
        _mcas = mcas[mask]

        _parameter_names, _init = self._prepare_init_for_fit(model, init)

        if prefit:
            print("Prefit with scipy for deriving initial values...")
            _bounds = self._prepare_bounds_for_fit(model, init, fixed, limits)
            
            # Protect against 0 counts causing division-by-zero in weights
            _sigma = np.maximum(np.sqrt(_data), 1)
            
            try:
                # Suppress scipy optimize warnings so they don't flood the terminal
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    
                    _pre_init, _init_cov = curve_fit(
                        model.total_model,
                        _mcas,
                        _data,
                        sigma=_sigma,
                        absolute_sigma=True,
                        p0=_init,
                        maxfev=500000,
                        bounds=_bounds,
                    )
                
                # Only overwrite initial values if covariance is valid
                if not np.isinf(_init_cov).any():
                    _init = _pre_init
                    _init_err = np.sqrt(np.diag(_init_cov))
                    self.print_table(_parameter_names, _init, _init_err)
                else:
                    print("Prefit converged but covariance is infinite. Falling back to manual init.")
                    
            except RuntimeError:
                print("Prefit failed to converge. Falling back to manual init.")

        # Note: For low background spectra, consider cost.ExtendedBinnedNLL in the future
        cost_function = cost.LeastSquares(
            _mcas, _data, np.maximum(np.sqrt(_data), 1), model.total_model
        )
        
        m = Minuit(cost_function, *_init)
        
        # Cleaned up dictionary iteration
        if limits:
            for k, v in limits.items():
                m.limits[k] = v
        if fixed:
            for k, v in fixed.items():
                m.fixed[k] = v
                
        return m
    

    def get_time_evolution_fitting_object(
        self,
        model,
        init,
        limits=None,
        fixed=None,
        MCA_range=None,
        time_range=None,
        n_timestamp=None,
        rate_limit=0,
        prefit=True,
    ):
        MCA_range = MCA_range if MCA_range is not None else self.SELECTED_MCA_RANGE
        time_range = time_range if time_range is not None else self.SELECTED_TIME_RANGE

        # print('get_time_evolution_fitting_object')
        # print(MCA_range)
        # print(time_range)
        
        times, rate, rate_err = self.get_time_evolution(
            MCA_range=MCA_range, time_range=time_range, n_timestamp=n_timestamp
        )

        mask = rate > rate_limit
        _times = times[mask]
        _rate = rate[mask]
        _rate_err = rate_err[mask]

        _parameter_names, _init = self._prepare_init_for_fit(model, init)

        if prefit:
            print("Prefit with scipy for deriving inital values")
            _bounds = self._prepare_bounds_for_fit(model, init, fixed, limits)
            _init, _init_cov = curve_fit(
                model.total_model,
                _times,
                _rate,
                sigma=_rate_err,
                absolute_sigma=True,
                p0=_init,
                maxfev=500000,
                bounds=_bounds,
            )
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
    def _radon_db_is_reachable(url):
        """Checks if the URL is accessible."""
        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
            return True
        except (requests.exceptions.RequestException, requests.exceptions.HTTPError):
            print("\n" + "!" * 60)
            print(f"CONNECTION ERROR: Cannot reach {url}")
            print("Please ensure you are connected to the MPIK network or VPN.")
            print("!" * 60 + "\n")
            return False

    @staticmethod
    def _radon_db_name_format(url):
        """Checks if the URL makes sense."""
        folder_name = url.split()[-1].split("/")[-2]
        is_upper = folder_name[0].isupper()
        if not is_upper:
            warnings.warn(
                f"The '{folder_name}' starts with a lowercase letter.\n"
                f"Radon db typically requires capitalization.\n"
                f"This can give an error while fetching the data.",
                UserWarning,
            )

    @staticmethod
    def _scrape_radon_db(url):
        # Sends a GET request to the specified URL and parses the content 
        # using BeautifulSoup.
        RadonAnalysis._radon_db_name_format(url)
        page = requests.get(url).text
        soup = BeautifulSoup(page, "html.parser")
        # Returns a list of URLs, appending the href attribute of each <a> tag 
        # that ends with '.root' to the base URL.
        return [
            url + node.get("href")
            for node in soup.find_all("a")
            if node.get("href").endswith(".root")
        ]

    @staticmethod
    def _prepare_init_for_fit(Model, init):
        # Retrieves the names of the parameters (excluding 'self') of 
        # the Model's constructor.
        parameter_names = list(inspect.signature(Model).parameters.keys())[1:]
        # Returns a list of initial values for these parameters based 
        # on the provided 'init' dictionary.
        return parameter_names, [init[p] for p in parameter_names]

    @staticmethod
    def _prepare_bounds_for_fit(
        model,
        init,
        fixed,
        limits,
        lower_bound=-np.inf,
        upper_bound=np.inf,
        epsilon=1e-9,
    ):
        # Retrieves the names of the parameters (excluding 'x') of 
        # the model's function.
        parameter_names = list(inspect.signature(model).parameters.keys())[1:]
        bounds = []

        # prepare limits and fixed.
        _limits = limits or {}
        _fixed = fixed or {}

        # Creates bounds for each parameter based on whether it is fixed 
        # (using 'epsilon' for tight bounds)
        # or if it is specified in limits
        # or free (using provided 'lower_bound' and 'upper_bound')
        for p in parameter_names:
            if _fixed.get(p, False):
                bounds.append((init[p] - epsilon, init[p] + epsilon))
            elif p in _limits.keys():
                lb = -np.inf if _limits[p][0] is None else _limits[p][0]
                ub =  np.inf if _limits[p][1] is None else _limits[p][1]
                bounds.append((lb, ub))
                # bounds.append((_limits[p][0], _limits[p][1]))
            else:
                bounds.append((lower_bound, upper_bound))
        # Separates the bounds into two lists: one for lower bounds 
        # and one for upper bounds.
        # Unzip the pairs into two lists
        lower_bounds, upper_bounds = zip(*bounds)  
        return lower_bounds, upper_bounds

    @staticmethod
    def print_table(_parameter_names, _init, _init_err):
        print(f"{'| Parameters':12} | {'Value':11} | {'Error':11} |")
        print("-" * 42)
        for l, i, j in zip(_parameter_names, _init, _init_err):
            print(f"| {l:10} | {i:11.3f} | {j:11.3f} |")
        print("-" * 42)
