# Simple IBM Solver

[![License](https://img.shields.io/github/license/NaokiHori/SimpleIBMSolver)](https://opensource.org/licenses/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SimpleIBMSolver/main)](https://github.com/NaokiHori/SimpleIBMSolver/commits/main)

[![Simulation Snapshot](https://github.com/NaokiHori/SimpleIBMSolver/blob/main/docs/source/snapshot2d.png)](https://youtu.be/nMAyrIYET10)

## Overview

This library numerically simulates the motion of rigid particles in two- and three-dimensional Cartesian domains using the finite-difference and immersed boundary methods.

## Dependency

This solver is built on top of [`SimpleNSSolver`](https://github.com/NaokiHori/SimpleNSSolver). 
Please check its repository for dependency details.

## Quick Start

1. **Set up your workspace**

   ```console
   mkdir -p /path/to/your/directory
   cd /path/to/your/directory
   ```

2. **Clone the repository**

   ```console
   git clone --recurse-submodules https://github.com/NaokiHori/SimpleIBMSolver
   cd SimpleIBMSolver
   ```

3. **Set the initial condition**

   Python 3 is used to conveniently initialize the flow fields. 
   Alternatively, `NPY` files can be provided in a different manner under `initial_condition/output/`.

   ```console
   cd initial_condition
   make output
   bash exec.sh
   cd ..
   ```

4. **Build the solver**

   ```console
   make output
   make all
   ```

5. **Run the simulation**

   ```console
   bash exec.sh
   ```

## Note

The immersed boundary method and the collision model are based on [this publication](https://www.sciencedirect.com/science/article/pii/S0045793021003716) with some modifications.

