# RadDetect
Advanced Spectroscopy and Emanation Analysis for Radon Content Characterization in Diverse Environments

``raddetect`` is a utility package for Radon content characterization. 
It is designed to handle `.root` files containing `channel`, `timestamp`, and `runtime` data, which are the result of the Radon Emanation Database's processing of ORTEC MCA list files. These files originate from detectors like Blumchen (an electrostatic radon monitor), Mona (an alpha spectroscopy setup), and CryoRadon (a cryogenic radon emanation monitor, which is not yet implemented). Due to a bug in the ORTEC file creation process, the runtime can be computed from the timestamps, which are set 60 seconds apart by the DAQ during list mode data acquisition.
 
> [!NOTE]
> If you are not connected to the MPIK network (or VPN), you will not be able to fetch data from the radon database automatically. In this case, you must provide the ROOT files locally to the analysis modules.


## Installation
Ensure your `pip` is up-to-date, as this project uses a `pyproject.toml` build system.

`git clone` this repo and:
```bash
# at top level of raddetect
pip install .

# OR for an editable installation (for development)
pip install -e .
```

## Running tests, linter, and code formatting

Before making a commit, it is important to run the tests and ensure the code is properly formatted. Tests are handled via `pytest`, general linting and formatting via `ruff`, and docstring formatting via `docformatter`. Make sure you have all of them installed (`pip install pytest ruff docformatter`).

To run the test suite:
```bash
pytest tests/
```

`ruff` is used to catch syntax errors, enforce best practices, and format the Python code logic.

- Check the current directory for errors: `ruff check .`
- Check a specific file: `ruff check path/to/file.py`
- Auto-fix fixable errors: `ruff check --fix`
- Format all code: `ruff format .`
- Check if files would be formatted (without changing them): `ruff format --check .`
- **Speed combo:** `ruff check --fix && ruff format`

While `ruff` formats the code, it intentionally ignores long text inside docstrings (avoiding the `E501 Line too long` error). We use `docformatter` specifically to wrap our summaries and descriptions to the 88-character limit.

- Format all docstrings in the project (uses `pyproject.toml` settings): `docformatter .`
- Format a specific file: `docformatter path/to/file.py`
- Check if docstrings need formatting (dry run, makes no changes): `docformatter --check .`
- **Complete Formatting Combo:** `ruff check --fix && ruff format && docformatter .`

## Development & Contribution Notes

To maintain a clean repository, we follow a strict policy regarding **Jupyter Notebooks** and **ROOT data files**.

### 1. Jupyter Notebooks
By default, changes to `.ipynb` files are **ignored** to prevent Git history bloat caused by cell metadata and execution outputs.
* **To contribute a new tutorial or update an existing one:** You must manually "force-add" the file.
* **Requirement:** Please **Clear All Outputs** (`Cell > All Output > Clear`) before staging to keep the diff readable.

```bash
git add -f tutorials/your_notebook.ipynb
```

### 2.Test Data (tests/data/)

Large binary `.root` files are ignored by default to keep the repository size manageable.

To add or update reference datasets for tests:

```bash
git add -f tests/data/your_data.root
```

### 3. ROOT
The `.root` data files are generated using the following ROOT and compiler configuration:

```bash
[radon@lfs1 ~]$ root-config --version
6.30.04
[radon@lfs1 ~]$ root-config --cxx
g++
[radon@lfs1 ~]$ root-config --cflags
-pthread -std=c++17 -m64 -I/cern/root_v6.30.04/include
```
