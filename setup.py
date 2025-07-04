#!/usr/bin/env python
"""Setup script for consensus_framework package."""

from setuptools import setup, find_packages

# Read the content of README.md for the long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="consensus_framework",
    version="0.1.0",
    author="Connor M. Frankston",
    author_email="franksto@ohsu.edu",
    description="A framework for robust statistical locus selection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/frankston/consensus-framework",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.2.0",
        "statsmodels>=0.12.0",
        "click>=8.0.0",
    ],
    entry_points={
        "console_scripts": [
            "consensus-framework=consensus_framework.bed_consensus_cli:construct_consensus_beds",
        ],
    },
)