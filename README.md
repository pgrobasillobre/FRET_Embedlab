# FretLab – Electronic Energy Transfer Computation

## Table of Contents

- [About](#about)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [License](#license)
- [Contact](#contact)

## About

**FretLab** is a high-performance computational tool designed to compute **electronic energy transfer (EET) rates**. It supports three general use cases:

1. **Donor to acceptor chromophore EET**
2. **Plasmonic substrate to acceptor chromophore EET**
3. **Donor to acceptor chromophore EET mediated by a plasmonic substrate**

FretLab is designed for speed and scalability using parallel processing with OpenMP and efficient linear algebra routines.

## Prerequisites

Before building FretLab, ensure the following dependencies are installed:

- CMake 3.14 or higher
- Python 2.7
- Fortran compiler (gfortran 9.3.0 or higher recommended)
- LAPACK/BLAS libraries (MKL suggested)
- Python `runtest` module (`pip install runtest`)
- OpenMP support

If using MKL (recommended), set the following environment variable:

```
export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin
```

## Installation

To build FretLab, run:

```
./setup.sh -b <build-dir> -fc <fort-path> -omp
```

### Options:
- `-b <build-dir>` : Name of the build directory (optional)
- `-fc <fort-path>` : Path to the Fortran compiler (optional)
- `-omp` : Enables OpenMP (recommended)

## License

This project is licensed under the **GNU General Public License v3.0**.

## Contact

For questions, support, or contributions, contact:

**pgrobasillobre@gmail.com**
