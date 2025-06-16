# FretLab – Electronic Energy Transfer Computation

<p align="center">
  <img src="https://raw.githubusercontent.com/pgrobasillobre/FretLab/main/docs/_static/FretLab.png" width="600">
</p>



## Table of Contents

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
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

In the presence of a plasmonic substrate (modeled via induced charges $$q_k$$ and, if applicable, induced dipoles $$\mu_k$$), the total coupling V becomes:

$$
V = V_{\text{coulomb}} \ + V_{\text{overlap}} \ + V_{\text{environment}}^{\text{   q}} + V_{\text{environment}}^{  \mu} 
$$

$$
V_{\text{coulomb}} =  \int d\mathbf{r} \ d\mathbf{r'} \ \frac{\rho_A^*(\mathbf{r}) \rho_D(\mathbf{r'})}{|\mathbf{r} - \mathbf{r'}|}
$$

$$
V_{\text{overlap}} = - \omega_0 \int d\mathbf{r} \ d\mathbf{r'} \ \rho_A^*(\mathbf{r}) \rho_D (\mathbf{r'}) 
$$

$$
V_{\text{environment}}^{\text{   q}} = \sum_k \left( \int d\mathbf{r} \ \frac{\rho_A^*(\mathbf{r})}{|\mathbf{r} - \mathbf{r}_k|} \right) q^\omega(\mathbf{r}_k; \rho_D)
$$


$$
V_{\text{environment}}^{\mu} = \sum_k \left( \int d\mathbf{r} \ -\frac{\mathbf{r} - \mathbf{r}_k}{|\mathbf{r} - \mathbf{r}_k|^3} \cdot \rho_A^*(\mathbf{r}) \right) \mu^\omega(\mathbf{r}_k; \rho_D)
$$

Where:

- $\rho_A$ and $\rho_D$ are the acceptor and donor charge densities  
- $\omega_0$ is the incident frequency  
- $q_k^\omega(\mathbf{r}_k; \rho_D)$ are the frequency-dependent induced charges at positions $\mathbf{r}_k$
- $\mu_k^\omega(\mathbf{r}_k; \rho_D)$ are the frequency-dependent induced dipoles at positions $\mathbf{r}_k$

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

This will generate an executable named `FretLab` inside the `build/` directory.

### Options:
- `-omp` : Enables OpenMP (recommended)
- `-fc <fort-path>` : Path to the Fortran compiler (optional)


### Running Tests:

After building, you can run the test suite with:

```
cd build/
ctest
```



## Usage

To execute a simulation, run:

```
./FretLab input_file.inp
```

To see example input files and different configuration options, refer to the files located in:

```
tests/
```


## License

FretLab is licensed under the **GNU General Public License v3.0**.

## Acknowledgments & Funding

This project has been supported by the **FARE 2020** program — *"Framework per l’attrazione e il rafforzamento delle eccellenze per la ricerca in Italia."*

FretLab has been developed as part of a collaborative scientific effort involving **Tommaso Giovannini**, **Piero Lafiosca**, **Sveva Sodomaco**, and **Chiara Cappelli**.

## Contact

For issues or contributions:

- Email: **pgrobasillobre@gmail.com**
- Github issues: https://github.com/pgrobasillobre/FretLab/issues

