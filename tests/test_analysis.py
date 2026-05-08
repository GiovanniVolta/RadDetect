import os

from raddetect.base_analysis import RadonAnalysis
from raddetect.blumchen.blumchen_analysis import BlumchenAnalysis
from raddetect.monalpha.monalpha_analysis import MonalphaAnalysis

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
# Using one of the test files for all tests, as requested
TEST_DATA_FILE = os.path.join(TEST_DATA_DIR, "rn13032026.root")


def test_radon_analysis_initialization():
    # Initialize the base class with an actual local test file
    analysis = RadonAnalysis(TEST_DATA_FILE)

    # Assert successful initialization
    assert analysis is not None
    assert hasattr(analysis, "mca")


def test_radon_analysis_get_mca_histogram():
    analysis = RadonAnalysis(TEST_DATA_FILE)
    data, mcas = analysis.get_mca_histogram(
        MCA_range=[0, 500], time_range=[0, 50], n_mca=10
    )

    # Verify the histogram outputs
    # 10 bins, 9 values and 9 center bins
    assert len(data) == 9
    assert len(mcas) == 9


def test_monalpha_initialization():
    # Initialize the class with an actual local test file
    analysis = MonalphaAnalysis(TEST_DATA_FILE)

    # Assert successful initialization
    assert analysis is not None


def test_blumchen_initialization():
    # Initialize the class with an actual local test file
    analysis = BlumchenAnalysis(TEST_DATA_FILE)

    # Assert successful initialization
    assert analysis is not None


# An alternative it owuld be to use @patch for places a real
# object or function with a "Mock" object so you can test your
# code in isolation without actually running the original logic.
