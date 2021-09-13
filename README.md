# interplay2020acsnano
This project contains codes required to reproduce the publication:  
A. K. Chew, B. C. Dallin, and R. C. Van Lehn. “Interplay of Ligand Properties and Core Size Dictate the Hydrophobicity 
of Monolayer-Protected Gold Nanoparticles.” *ACS Nano*. **2021**, *15*, 3, 4534–4545.

# Initiation and requirements
To download all the code:
```buildoutcfg
# Clone the repository
git clone git@github.com:alexkchew/interplay2020acsnano.git
# Add any submodules crucial to this repository
git submodule update --init
# Pulling submodules
git pull
```
## Installing Python environment
To access the same Python environment, install a conda environment using the Anaconda suite. 


# Overview

In this work, we quantified the hydrophobicity of monolayer-protected gold nanoparticles (GNPs) using atomistic 
molecular dynamics (MD) simulations. First, classical MD simulations were used to model GNPs in solution. Then, we 
compute the local hydration free energies at the nanoparticle-water interface by analyzing the interfacial fluctuations 
of water. These computational tools allow us to analyze interfacial hydrophobicity on non-planar geometries, which can 
yield insights into how GNPs may interact with other biomolecules. 

![Overview of GNP model development](/images/full_system_setup_main.png)



# Zenodo repository

All simulations and large data is stored within the Zenodo repository linked below:

To install, download the repository and decompress using the command below:
```buildoutcfg
# Command to decompress *.tar.gz files
tar -xvf interplay2020acsnano_zenodo.tar.gz
```

This is in a work of progress - still updating this page! Estimated time of completion: October 2021

