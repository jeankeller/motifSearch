from setuptools import setup, find_packages
import os, sys


#Package meta-data
NAME = "motifSearch"
DESCRIPTION = "A python package to find motifs in DNA/protein sequences"
URL = "https://github.com/jeankeller/motifSearch"
EMAIL = "kellerjeanphd@gmail.com, cyril.libourel@gmail.com"
AUTHOR = "Jean Keller, Cyril Libourel"
REQUIRES_PYTHON = ">=3.6.0"
VERSION = "1.0"

#Required packages
REQUIRED = ["biopython", "pandas"]

here = os.path.abspath(os.path.dirname(__file__))

#Import README.md as long-description
try:
    with open(os.path.join(here, "README.md"), encoding = "UTF-8") as freadme:
        ldesc = "\n" + freadme.read()
except FileNotFoundError():
    ldesc = DESCRIPTION


setup(
    name=NAME,
    version="1.0.0",
    long_description=ldesc,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=["motifSearch"],
    install_requires=REQUIRED,
    package_data={"motifSearch": ["motifSearch/motifSearch2_plot.R"]},
    include_package_data=True,

    entry_points={
        "console_scripts": [
            "motifsearch=motifSearch.motifSearch_main:main",
            ],
        },
    license="LGPLv3+",
    classifiers=[
        "License :: OSI Approved :: GNU Lesser Public License v3 or Later",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: R",
        "Programming Language :: R :: 3",
        "Programming Language :: R :: 3.6",
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independant",
        ],
    )
    
    
    
