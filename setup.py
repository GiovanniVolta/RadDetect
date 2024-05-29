from setuptools import setup, find_packages

setup(
    name='raddetect',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        "numpy", 
        "scipy", 
        "matplotlib", 
        "pandas", 
        "pyyaml"
    ],
    # Optional metadata
    author='MPIK',
    author_email='your.email@example.com',
    description='Advanced Spectroscopy and Emanation Analysis for Radon Content Characterization in Diverse Environments',
    keywords='Radon',
    url='https://github.com/GiovanniVolta/RadDetect?tab=readme-ov-file',  # project home page
)