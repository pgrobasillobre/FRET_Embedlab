# FRET Embedlab

## Table of Contents

- [About](#about)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)

## About


Computation of electronic energy transfer (EET) rates for three general cases:

   1. Donor to acceptor chromophores EET
   2. Plasmonic substrate to acceptor chromophore EET
   3. Donor to acceptor chromophores EET mediated by plasmonic substrate


## Installation

    ./setup.sh -b buildir -omp

    buildir : the name of the directory where the code is built (optional)
    omp     : option to activate OpenMP (Recommended)


### Prerequisites

    - CMake3.14 or higher
    - Python 2.7.
    - fortran compiler (gfortran suggested)
    - Lapack/Blas libraries (MKL suggested)
    - runtest module python
      (pip install runtest)
    - OpenMP
    - If you have MKL installed (recommended)
      export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin

```bash
# Example
Node.js
npm
