# SU(N) Symmetric Heisenberg 3-Column Spin-1 Resummation-Based SSE

This repository contains the implementation of the **Resummation-Based Stochastic Series Expansion (SSE)** method for solving **SU(N) symmetric Heisenberg spin models**. The focus is on 3-column spin-1 systems, leveraging resummation techniques to enhance convergence and efficiency. This code is ideal for exploring quantum phases and transitions in higher-symmetry spin models.

---

## Features

- **SU(N) Symmetry Support:** Fully supports generalized SU(N) Heisenberg models for spin-1 systems.
- **Resummation Techniques:** Implements advanced series resummation for improved numerical accuracy.
- **Columnar Lattice Geometry:** Specializes in 3-column lattice configurations for spin-1 models.
- **Customizable Parameters:** User-defined symmetry order (N), spin coupling constants, and lattice dimensions.
- **Observables Computation:** Calculates thermodynamic properties, correlation functions, and energy spectra.

---

## Getting Started

### **Prerequisites**

Ensure the following are installed:

- **C++ Compiler** (C++17 or later, e.g., GCC, Clang)
- **CMake** (optional for C++ build automation)
- **Python 3.7+** for generating coefficients or preliminary analyses:
  - **NumPy**
  - **Matplotlib** (for visualizations)
  - **SymPy** (for symbolic manipulations)

Install Python dependencies with:

```bash
pip install numpy matplotlib sympy
