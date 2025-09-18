{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 AppleColorEmoji;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Stokes FEM Codes (based on iFEM)\
\
This repository contains MATLAB codes for solving Stokes problems using pressure-robust Enriched Galerkin (EG) finite element methods. It is built on top of selected routines from [iFEM](https://www.math.uci.edu/~chenlong/programming.html), an open-source MATLAB finite element package by Long Chen (GNU GPL v3).\
\
---\
\
## 
\f1 \uc0\u55357 \u56514 
\f0  Repository Structure\
\
- **Stokes/**  \
  Main scripts for running Stokes problem solvers:  \
  - `Stokes2.m`, `Stokes3.m`: core solvers  \
  - `main_EG2.m`, `main_EGWG2.m`: driver scripts  \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0   - `main_EG3.m`, `main_EGWG3.m`: driver scripts  \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 \
- **iFEM_files/**  \
  Essential iFEM routines included here for convenience:  \
  - Mesh generation: `squaremesh.m`, `cubemesh.m`  \
  - Basis/gradients: `gradbasis.m`, `gradbasis3.m`  \
  - Quadrature rules: `quadpts.m`, `quadpts3.m`, `verifyquadpts3.m`  \
  - Utilities: `auxstructure.m`, `mycross.m`, `myunique.m`, etc.  \
\
Only the parts of iFEM strictly needed to run these codes are included.  \
For the full package, see the [iFEM homepage](https://www.math.uci.edu/~chenlong/programming.html).\
\
---\
\
## 
\f1 \uc0\u9878 \u65039 
\f0  License\
\
This repository is released under the **GNU GPL v3 license**.  \
Since iFEM is also licensed under GPL v3, redistribution and modifications are permitted under the same terms.  \
\
- Original iFEM copyright:  \
  > Copyright (C) Long Chen, University of California, Irvine.  \
  > Licensed under the GNU General Public License v3.  \
\
- Extensions, modifications, and Stokes-specific codes:  \
  > Copyright (C) Seulip Lee <seulip.lee@tufts.edu>, Lin Mu <linmu@uga.edu>\
\
See the [LICENSE](LICENSE) file for full details.}