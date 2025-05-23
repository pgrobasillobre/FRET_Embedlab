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

### Theoretical Framework

All quantities in FretLab are assumed to be given or computed in **atomic units**. The EET rate, denoted \(\kappa_{\text{EET}}\), is computed using **Fermi's Golden Rule**:

\[
\kappa_{\text{EET}} = \frac{2\pi}{\hbar} \, |V|^2 \, J
\]

Where:
- \(V\) is the total coupling potential between donor and acceptor,
- \(J\) is the spectral overlap integral.

If a **plasmonic substrate** is present (modeled as a set of induced charges \(q_k\)), the total potential \(V\) is given by:

\[
\begin{aligned}
V &= V_{\text{Coulomb}} + V_{\text{overlap}} + V_{\text{environment}} = \\
&= \int d\mathbf{r} \, d\mathbf{r'} \, \frac{\rho_A^*(\mathbf{r}) \, \rho_D(\mathbf{r'})}{|\mathbf{r} - \mathbf{r'}|} 
- \omega_0 \int d\mathbf{r} \, \rho_A^*(\mathbf{r}) \, \rho_D(\mathbf{r}) \\
&\quad + \sum_k \left( \int d\mathbf{r} \, \frac{\rho_A^*(\mathbf{r})}{|\mathbf{r} - \mathbf{r}_k|} \right) \, q_k^{\omega}(\mathbf{r}_k; \rho_D)
\end{aligned}
\]

Where:
- \(\rho_A\) and \(\rho_D\) are the acceptor and donor charge densities, respectively,
- \(\omega_0\) is the incident frequency,
- \(q_k^{\omega}(\mathbf{r}_k; \rho_D)\) are the frequency-dependent induced charges at positions \(\mathbf{r}_k\), due to the donor density.


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
