import numpy as np

from ..base_analysis import RadonAnalysis


class BlumchenAnalysis(RadonAnalysis):
    """
    A class for analyzing data from ROOT files, specifically designed for Blumchen analysis.

    This class provides methods to retrieve data from a specified URL, generate histograms of MCA channel data,
    plot time evolution of the data, and fit the data using specified models.

    Attributes:
        mca (numpy.ndarray): MCA data.
        timestamp (numpy.ndarray): Timestamp data.
        runtime (numpy.ndarray): Runtime data.
    """

    DEFAULT_MCA_RANGE = [910, 1060]
    DEFAULT_TIME_RANGE = [240, np.inf]
