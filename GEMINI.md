# Gemini AI Assistant Guidelines

## Project Overview
* **Description:** A Python package for analyzing radon emanation data and alpha spectroscopy.
* **Data Pipeline:** Data is uniformly pre-processed by a standard script and acquired using the same electronics. Do not write code to re-format or pre-process raw data unless explicitly asked.
* **Data Ingestion:** Data is either passed manually or fetched from a private database. 

## Tech Stack
* **Language:** Python 3.x
* **Core Libraries:** [Add your specific libraries here, e.g., NumPy, Pandas, SciPy, Matplotlib]
* **Database Driver:** [Add your DB tool here, e.g., SQLAlchemy, psycopg2]

## Security & Privacy (CRITICAL)
* **No Hardcoded Secrets:** NEVER generate code that hardcodes database passwords, API keys, or sensitive URLs. 
* **Environment Variables:** Always use `os.environ` or the `python-dotenv` package to load database credentials and configurations.
* **Data Protection:** Assume all data fetched from the database is sensitive. Do not write code that logs raw data outputs or exports data to unsecured public directories.

## Coding Standards & Style
* **Style Guide:** Strictly adhere to PEP 8 standards.
* **Type Hinting:** Use Python type hints (`typing` module) for all function arguments and return values to ensure data pipeline consistency.
* **Docstrings:** Use [Choose one: NumPy or Google] style docstrings for all modules, classes, and functions. Clearly document the expected shape and type of the input data.
* **Error Handling:** Implement robust `try/except` blocks, especially around database connection functions and file reading operations. 

## Architecture & File Structure
* `/raddetect`: Contains the core package modules.
* `/tests`: Contains unit tests for the analysis logic.
* `/utility`: Local directory for yml configureation inputs.
* `/tutorials`: Contains tutorial notebooks.
* `/docs`: Contains documentations in principles.

## Important Quirks & Rules
* **Data Format:** Assume all incoming data is already in the standardized, pre-processed format.
* **Modularity:** Keep the database fetching logic strictly separated from the mathematical analysis logic. Analysis functions should accept data objects (like DataFrames or arrays), not database connection objects.