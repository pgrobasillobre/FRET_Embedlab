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

$$
\kappa_{\text{EET}} = \frac{2\pi}{\hbar} \ |V|^2 \ J
$$

Where:
- V is the total coupling potential between donor and acceptor
- J is the [spectral overlap integral](https://github.com/pgrobasillobre/SpectralOverlap)

In the presence of a plasmonic substrate (modeled via induced charges $$q_k$$), the total coupling V becomes:

![Equation](https://latex.codecogs.com/svg.image?V%20=%20V_{\text{Coulomb}}%20+%20V_{\text{overlap}}%20+%20V_{\text{environment}}%20=%20\int%20d\mathbf{r}%20d\mathbf{r'}%20\frac{\rho_{A}^*(\mathbf{r})%20\rho_{D}(\mathbf{r'})}{|\mathbf{r}-\mathbf{r'}|}%20-%20\omega_0%20\int%20d\mathbf{r}%20\rho_{A}^*(\mathbf{r})%20\rho_{D}%20+%20\sum_{k}%20\left(%20\int%20d\mathbf{r}%20\frac{\rho_{A}^*(\mathbf{r})%20}{|\mathbf{r}-\mathbf{r}_{k}|}%20\right)%20q^{\omega}(\mathbf{r}_{k};%20\rho_{D}))


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
./setup.sh -omp -fc <fort-path>

cd build/
make
```

### Options:
- `-omp` : Enables OpenMP (recommended)
- `-fc <fort-path>` : Path to the Fortran compiler (optional)


### Running Tests:

After building, you can run the test suite with:

```
cd build/
ctest
```

## License

FretLab is licensed under the **GNU General Public License v3.0**.

## Contact

For issues or contributions:

- Email: **pgrobasillobre@gmail.com**
- Github issues: https://github.com/pgrobasillobre/FretLab/issues

