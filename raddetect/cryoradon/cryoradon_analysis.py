import numpy as np

from ..base_analysis import RadonAnalysis


class CryoRadonAnalysis(RadonAnalysis):
    """A class for analyzing data from ROOT files, specifically designed for CryoRadon
    analysis.

    This class provides methods to retrieve data from a specified URL, generate histograms of MCA channel data,
    plot time evolution of the data, and fit the data using specified models.

    Attributes:
        mca (numpy.ndarray): MCA data.
        timestamp (numpy.ndarray): Timestamp data.
        runtime (numpy.ndarray): Runtime data.
    """
    # Used in the full spectrum plot
    DEFAULT_MCA_RANGE = [0, 800]
    DEFAULT_TIME_RANGE = [0, np.inf]
    
    # To select a specific line to plot the temporal evolution
    SELECTED_MCA_RANGE = [350, 390]
    SELECTED_TIME_RANGE = [0, np.inf]
    
    def __init__(self, filename, compute_runtime_from_timestamp=False, **kwargs):
        # Pass required arguments to the parent class (RadonAnalysis) if needed
        super().__init__(filename=filename) 
        
        self.compute_runtime_from_timestamp = compute_runtime_from_timestamp
        
        # Override the class defaults with passed arguments, or fall back to defaults
        self.DEFAULT_MCA_RANGE = kwargs.get('DEFAULT_MCA_RANGE', self.DEFAULT_MCA_RANGE)
        self.DEFAULT_TIME_RANGE = kwargs.get('DEFAULT_TIME_RANGE', self.DEFAULT_TIME_RANGE)
        self.SELECTED_MCA_RANGE = kwargs.get('SELECTED_MCA_RANGE', self.SELECTED_MCA_RANGE)
        self.SELECTED_TIME_RANGE = kwargs.get('SELECTED_TIME_RANGE', self.SELECTED_TIME_RANGE)
