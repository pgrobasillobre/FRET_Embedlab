# FretLab – Electronic Energy Transfer Computation

<p align="center">
  <img src="https://raw.githubusercontent.com/pgrobasillobre/FretLab/main/docs/_static/FretLab.png" width="600">
</p>



## Table of Contents

- [About](#about)
- [Installation](#installation)
- [License](#license)
- [Contact](#contact)

## About

**FretLab** is a high-performance computational tool designed to compute **electronic energy transfer (EET) rates**. It supports three general use cases:

1. **Donor to acceptor chromophore EET**
2. **Plasmonic substrate to acceptor chromophore EET**
3. **Donor to acceptor chromophore EET mediated by a plasmonic substrate**

FretLab is designed for speed and scalability using parallel processing with OpenMP and efficient linear algebra routines.

## Theoretical Framework

All quantities in FretLab are expressed in atomic units. The EET rate is 
calculated using Fermi’s Golden Rule:

<p align="center">
  <code>κ_EET = (2π / ℏ) · |V|² · J</code>
</p>

Where:
- V is the total coupling potential between donor and acceptor
- J is the spectral overlap integral

In the presence of a plasmonic substrate (modeled via induced charges q_k),
the total coupling V becomes:

<p align="center">
  <code>
    V = V_Coulomb + V_overlap + V_environment<br><br>
    = ∫ d𝒓 d𝒓′ [ρ_A*(𝒓) ρ_D(𝒓′)] / |𝒓 − 𝒓′|<br>
    − ω₀ ∫ d𝒓 ρ_A*(𝒓) ρ_D(𝒓)<br>
    + ∑ₖ [ ∫ d𝒓 ρ_A*(𝒓) / |𝒓 − 𝒓ₖ| ] · qₖ^ω(𝒓ₖ; ρ_D)
  </code>
</p>

Where:
- ρ_A and ρ_D are the acceptor and donor charge densities
- ω₀ is the incident frequency
- qₖ^ω(𝒓ₖ; ρ_D) are the frequency-dependent induced charges at positions 𝒓ₖ

## Installation

FretLab requires the following dependencies:

- CMake 3.14 or higher
- Python 3.8+
- Fortran compiler (gfortran 9.3.0 or higher recommended)
- LAPACK/BLAS libraries (MKL suggested)
- Python `runtest` module (`pip install runtest`)
- OpenMP support

If using MKL (recommended), set the following environment variable:

```
export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin
```

### Compilation

To build FretLab, run:

```
./setup.sh -b <build-dir> -fc <fort-path> -omp
```

### Options:
- `-b <build-dir>` : Name of the build directory (optional)
- `-fc <fort-path>` : Path to the Fortran compiler (optional)
- `-omp` : Enables OpenMP (recommended)

## License

FretLab is licensed under the **GNU General Public License v3.0**.

## Contact

For issues or contributions:

- Email: **pgrobasillobre@gmail.com**
- Github issues: https://github.com/pgrobasillobre/FretLab/issues
