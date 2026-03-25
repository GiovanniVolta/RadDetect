# RadDetect
Advanced Spectroscopy and Emanation Analysis for Radon Content Characterization in Diverse Environments

``raddetect`` is a utility package for Radon content characterization. 

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

## Running tests, linter and code formatter
Before making a commit it is important to run the tests and be sure the code is properly formmatted. The testsa re done via `pytests` and the formatting via `ruff`. Be sure you have them installed.

```bash
pytest tests/
```

For the `ruff`:
- to check the current directory: `ruff check .`, or `ruff check path/to/file.py` for a specific file
- to auto-fix the errors: `ruff check --fix`
- to format all files: `ruff format .`
- to check if files would be formatted: `ruff format --check .`
- speed combo `ruff check --fix && ruff format`
