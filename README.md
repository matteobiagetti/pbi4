# PBI4

PBI4 is an FFT-based Python code to measure the power spectrum and bispectrum in real and redshift space of a distribution of particles in a box. It  implements up to 4$^{th}$-order grid interpolation and an interlacing technique to reduce the aliasing contribution arising from the sampling of the density field, as introduced in [Sefusatti et al. (2015)](https://arxiv.org/pdf/1512.07295). 

PBI4 is an upgrade of the Fortran code [PowerI4](https://github.com/sefusatti/PowerI4), which was based on a previous version by Roman Scoccimarro which, in turn, used algorithms developed by Hugh Couchmann.

Authors: Matteo Biagetti, Emiliano Sefusatti, Francesco Verdiani

## Installation and Running

The code is written in Python 3. For computationally heavy routines, it calls Fortran-compiled modules.

#### Requirements:

- basic Python packages (`numpy`, `scipy`, `pandas`,`h5py`). 

- a fortran compiler (e.g. [gcc](https://gcc.gnu.org/) )

- [FFTW](https://www.fftw.org/) 

#### PIP Installation:

First, you need to link FFTW libraries in the setup file. To do that, open the `setup.py` file and modify the paths `/usr/local/Cellar/fftw/3.3.10_1/lib` and `/usr/local/Cellar/fftw/3.3.10_1/include` found in `library_dirs` and `include_dirs`, respectively, with the path to FFTW libraries in your system.
 
From the main directory in the repository, you can now run 
 
 `python -m pip install .` 
 
 The installation was also tested within `conda` and `pyenv` virtual environments.

#### Custom Installation:

If the PIP installation fails, you can resort to a brute force installation as follows:

- Go to pbi4/F90_modules/f2pyscriptI4.sh and edit the lib and include path to the correct path to the FFTW libraries
- Run the script using `sh f2pyscriptI4.sh`
- Add `export PYTHONPATH=“${PYTHONPATH}:/path/to/pbi4/pbi4` to your `~/.bashrc` file.

 ### Running

 The code can be run in notebook mode or script mode. An example for both modes is found on the `Tutorials` folder. 
 
 The notebook `test.ipynb` in `Tutorials/Notebook/` contains an hands-on guide to the different parts of the code. 
