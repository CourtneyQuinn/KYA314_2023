# KYA314_2023
Collection of code used in KYA314 Dynamical Systems and Chaos in semester 2 of 2023.

## Miniconda installation
Conda is an open source package/environment management system that we will use to set up and run our Python-based computations.  Miniconda is a free minimal installer for conda from which we will build upon to set up our computing environment.  If you do not have conda installed, please visit 

    https://docs.conda.io/en/latest/miniconda.html
  
and download the latest version suitable for your operating system. Opening the downloaded file will bring up a prompt which you can follow to install Miniconda. *Note: Mac/Linux users will likely need to select "Add Anaconda to my PATH environment variable"*

Open up a window in either **Anaconda Prompt** (Windows) or **Terminal** (Mac/Linux).  The input line should look something like this:

    (base) C:\>
    
If it does not have the `(base)` then you may need to restart your computer. If after restarting there is still no `(base)` then conda is not installed properly in that terminal window. Please come see me to troubleshoot. 

## Creating a conda environment
To create the environment needed for this unit, download the `dynamics_environment.yml` file from this repo. Navigate from your command line to the directory of the file and run

    conda env create -f dynamics_environment.yml
    
This will take a few minutes.  Please ensure you are connected to the internet. To activate the environment simply run

    conda activate dynamics
    
## Initialising Jupyter
To initialise a Jupyter notebook use the command

    jupyter notebook
    
To initialise a Jupyter lab session use the command

    jupyter lab
