# FRET Embedlab

## Table of Contents

- [About](#about)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [License](#license)
- [Contact](#contact)

## About

Computation of electronic energy transfer (EET) rates for three general cases:

   1. Donor to acceptor chromophores EET
   2. Plasmonic substrate to acceptor chromophore EET
   3. Donor to acceptor chromophores EET mediated by plasmonic substrate


### Prerequisites

   - CMake3.14 or higher
   - Python 2.7.
   - fortran compiler (gfortran 9.3.0 or higher suggested)
   - Lapack/Blas libraries (MKL suggested)
   - runtest module python
     (pip install runtest)
   - OpenMP
   - If you have MKL installed (recommended)
      export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin


## Installation

   ./setup.sh -b buildir -omp

   buildir : the name of the directory where the code is built (optional)
   omp     : option to activate OpenMP (Recommended)

## License 

   This project is licensed under the GNU General Public License v3.0

## Contact

   For any issue contact pgrobasillobre@gmail.com
   

