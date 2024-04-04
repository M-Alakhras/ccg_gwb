#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.rst") as readme_file:
    readme = readme_file.read()


requirements = [
    "appdirs>=1.4.4",
    "astropy>=6.0.0",
    "pint-pulsar>=0.9.8",
    "psrqpy>=1.2.8",
    "tqdm>=4.66.2",
    "ptmcmcsampler>=2.1.1",
    "enterprise_extensions>=2.4.2",
    "enterprise-pulsar @ git+https://github.com/M-Alakhras/enterprise",
    "networkx>=3.2.1",
]

test_requirements = []

setup(
    name="ccg_gwb",
    description="A package developed by Complex & Cosmology Group (CCG) at Shahid Beheshti University (SBU). \
    It uses data-driven techniques to search the pulsar timing arrays (PTAs) data for \
    gravitational wave background (GWB) signal.",
    long_description=readme,
    author="Mohammad Alakhras",
    author_email="mohammadalakhras1989@gmail.com",
    url="https://github.com/M-Alakhras/ccg_gwb",
    packages=find_packages(include=["ccg_gwb", "ccg_gwb.*"]),
    package_dir={"ccg_gwb": "ccg_gwb"},
    python_requires=">=3.9, <3.11",
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords="ccg_gwb",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data=True,
    test_suite="tests",
    tests_require=test_requirements,
)
