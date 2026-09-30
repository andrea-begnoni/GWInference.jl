# GWInference

[![Build Status](https://github.com/andrea-begnoni/GW.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrea-begnoni/GW.jl/actions/workflows/CI.yml?query=branch%3Amain)

![alt text](logo.png)

### Perform accurate and fast parameter estimation of Gravitational Waves sources using the Fisher Information Matrix (alias GWJulia, as [link](https://arxiv.org/abs/2506.21530)).

Main features:

- Generate a catalog of compact binary sources accordingly to GWTC-4 arXiv.2508.18083
- Perform a Fisher Matrix analysis with automatic differentiation
- Analyze the results and obtain the errors on the parameters of the binary
- Very fast, it takes less than a second to compute the Fisher Matrix, independently of the waveform model and the detector(s) configuration

## Where to start?

After the installation there is an extensive tutorial called `easy_GWJulia.ipynb`, we leave to it the description of all the features available in the code.

## Installation

Before proceding with the installation verify that you have Julia installed on your machine and we suggest to use the most recent stable version. You can find instructions on how to install Julia here: https://julialang.org/downloads/.

To install the package you have two possibilities, depending on your needs. If you do not need to add new waveforms, develop new functionalities or extend the parameter space you can download the package from the registry (equivalent as Python's pip install)

### Download from the registry
GWInference is now a Julia package on the registry, so to install it just run
```julia
] #package mode
add GWInference
```
Then just import it when you need in your Julia script
```julia
using GWInference
```
### Git clone repository
If you plan to modify/develop the package, it would be easier to clone this repository and then activate the package.

To clone the repository
```
git clone https://github.com/andrea-begnoni/GWInference.jl.git
```
Then in your Julia script write
```julia
using Pkg
Pkg.activate("/path/to/the/GWInference/folder") 
# alternatively Pkg.develop(path="/path/to/the/GWInference/folder") 
using GWInference
```

You are now ready to go. We suggest you to start from `easy_GWJulia.ipynb`
If you installed the package from Julia command line we suggest to create a new folder and download there the `easy_GWJulia.ipynb` file. Then when you'll run the code it will create automatically the sub-folder needed. 
Otherwise you can clone the whole repository on your machine, then you'd need to pay attention to which GWInference folder you are activating (the package from the registry or the git cloned repository)

## Structure of the repository

Directories:

- `src`: where all the source code is located
- `test`: contains some test functions
- `useful_files`: where the sensitivities of the interferometers and some files needed for the waveforms are kept
- `catalogs`: automatically generated during the first run, used to store generated catalogs
- `output`: automatically generated during the first run, used to store the outputs (e.g. Fisher matrix, forecast uncertainties)

### The `src` directory

The `src` folder contains the 4 modules used by the code:

- `catalog.jl` contains population models and functions to generate and load the catalogs
- `detector.jl` contains all relevant info (e.g., positions, noises, patter functions) to output the signal measured by each detector (for the waveform uses `waveforms`, see below)
- `utils.jl` contains multi-purpose functions and the functions to analyze the Fisher Matrix
- `waveforms` folder containing all available waveforms models, tested to be compatible with their LAL counterparts


### Citating this work

If you use this code in your work, please cite:

```
@article{Begnoni:2025oyd,
    author = "Begnoni, Andrea and Anselmi, Stefano and Pieroni, Mauro and Renzi, Alessandro and Ricciardone, Angelo",
    title = "{Detectability and Parameter Estimation for Einstein Telescope Configurations with GWJulia}",
    eprint = "2506.21530",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "6",
    year = "2025"
}
```

### Contacts

If you need any help, feel free to contact andrea.begnoni@phd.unipd.it
