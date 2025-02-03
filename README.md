# 🚀 FRET Embedlab – Electronic Energy Transfer Computation

## 📌 Table of Contents

- [About](#-about)
- [Prerequisites](#-prerequisites)
- [Installation](#-installation)
- [License](#-license)
- [Contact](#-contact)

## ✨ About

FRET Embedlab computes **electronic energy transfer (EET) rates** for three general cases:

1. **Donor to acceptor chromophore EET**  
2. **Plasmonic substrate to acceptor chromophore EET**  
3. **Donor to acceptor chromophore EET mediated by a plasmonic substrate**  

This tool enables precise calculations for energy transfer in molecular and plasmonic systems.

## ⚙️ Prerequisites

Before installing, ensure you have the following dependencies:

- 🔹 **CMake 3.14 or higher**
- 🔹 **Python 2.7**  
- 🔹 **Fortran compiler** (gfortran **9.3.0** or higher recommended)  
- 🔹 **LAPACK/BLAS libraries** (MKL suggested)  
- 🔹 **Python `runtest` module** (install via `pip install runtest`)  
- 🔹 **OpenMP support**  
- 🔹 If using **MKL** (recommended), set:  

```
export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin
```


## 🛠️ Installation

Run the following command to build the project:

```
./setup.sh -b <build-dir> -fc <fort-path> -omp
```


Options:
- **`-b build-dir`** : Name of the directory where the code is built *(optional)*
- **`-fc fort-path`** : Path to the Fortran compiler *(optional)*
- **`-omp`** : Enables OpenMP *(recommended)*

## 📜 License

This project is licensed under the **GNU General Public License v3.0**.

## 📬 Contact

For inquiries or issues, please contact:

📧 **pgrobasillobre@gmail.com**
