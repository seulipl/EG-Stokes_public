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

---

## 📖 Citation

If you use this repository in your research, please cite the following paper:

- X. Hu, S. Lee, L. Mu, and S.-Y. Yi (2024). *Pressure-robust enriched Galerkin methods for the Stokes equations*. Journal of Computational and Applied Mathematics. [doi.org/10.1016/j.cam.2023.115449](doi.org/10.1016/j.cam.2023.115449).
- S. Lee and L. Mu (2024). *A low-cost, parameter-free, and pressure-robust enriched Galerkin method for the Stokes equations*. Computer and Mathematics with Applications. [doi.org/10.1016/j.camwa.2024.04.023](doi.org/10.1016/j.camwa.2024.04.023).

BibTeX entries:

```bibtex
@article{hu2024pressure,
  title={Pressure-robust enriched Galerkin methods for the Stokes equations},
  author={Hu, Xiaozhe and Lee, Seulip and Mu, Lin and Yi, Son-Young},
  journal={Journal of Computational and Applied Mathematics},
  volume={436},
  pages={115449},
  year={2024},
  publisher={Elsevier}
}

@article{lee2024low,
  title={A low-cost, penalty parameter-free, and pressure-robust enriched Galerkin method for the Stokes equations},
  author={Lee, Seulip and Mu, Lin},
  journal={Computers \& Mathematics with Applications},
  volume={166},
  pages={51--64},
  year={2024},
  publisher={Elsevier}
}
