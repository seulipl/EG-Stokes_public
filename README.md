# Stokes FEM Codes (based on iFEM)

This repository contains MATLAB codes for solving Stokes problems using pressure-robust Enriched Galerkin (EG) finite element methods. It is built on top of selected routines from [iFEM](https://www.math.uci.edu/~chenlong/programming.html), an open-source MATLAB finite element package by Long Chen (GNU GPL v3).

---

## 📂 Repository Structure

- **Stokes/**  
  Main scripts for running Stokes problem solvers:  
  - `Stokes2.m`, `Stokes3.m`: sample problem information
  - `main_EG2.m`, `main_EG3.m`: core Stokes solvers in 2D and 3D using [pressure-robust EG methods](https://doi.org/10.1016/j.cam.2023.115449).
  - `main_EGWG3.m`, `main_EGWG3.m`: core Stokes solvers in 2D and 3D using [pressure-robust modified EG methods](https://doi.org/10.1016/j.camwa.2024.04.023).

- **iFEM_files/**  
  Essential iFEM routines included here for convenience:  
  - Mesh generation: `squaremesh.m`, `cubemesh.m`  
  - Basis/gradients: `gradbasis.m`, `gradbasis3.m`  
  - Quadrature rules: `quadpts.m`, `quadpts3.m`, `verifyquadpts3.m`  
  - Utilities: `auxstructure.m`, `mycross.m`, `myunique.m`, etc.  

Only the parts of iFEM strictly needed to run these codes are included.  
For the full package, see the [iFEM homepage](https://www.math.uci.edu/~chenlong/programming.html).

---

## ⚖️ License

This repository is released under the **GNU GPL v3 license**.  
Since iFEM is also licensed under GPL v3, redistribution and modifications are permitted under the same terms.  

- Original iFEM copyright:  
  > Copyright (C) Long Chen, University of California, Irvine.  
  > Licensed under the GNU General Public License v3.  

- Extensions, modifications, and Stokes-specific codes:  
  > Copyright (C) Seulip Lee <seulip.lee@tufts.edu>, Lin Mu <linmu@uga.edu>

See the [LICENSE](LICENSE) file for full details.