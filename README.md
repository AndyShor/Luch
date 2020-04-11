
[![License: GPL](https://img.shields.io/badge/License-GPL-yellow.svg)](https://opensource.org/licenses/GPL-3.0)
[![build_and_test](https://github.com/AndyShor/Luch/workflows/build_and_test/badge.svg)]
[![codecov](https://codecov.io/gh/AndyShor/Luch/branch/master/graph/badge.svg)](https://codecov.io/gh/AndyShor/Luch)

# About

Luch is a simulation and visualization toolkit for first order charged particle beam optics.
Luch uses transfer matrix approach to trace and visualize beams  in beam lines
made of standardized optical elements such as drift spaces, dipole magnets,
electrostatic accelerating columns, einzel and quadrupole lenses.


The transfer matrixes of optical elements are based on  work of Gillespie [1-4]
Compared  to old realizations Luch is written with OOP approach in mind.
The benefits of such approach for beam optics simulations described by Christofer K Allen [5]
 
Using OOP and making Luch a toolkit rather than an executable 
 prioritizes ease of expansion, future support and flexibility.

Luch represents a beam line as a list of objects, each having a method
for tracing a particle. Thus adding new type of element requires to create
just a new class, defining how it changes particle properties such as
coordinates, energy or charge. 

Luch is written in Python and is meant to be primarily used in interactive manner within
Jupyter notebooks or as interactive web-apps with libraries such as streamlit.

Luch is in the development stage, so unit tests and documentation are not there yet,
and changes to be expected. 

# References

[1] G. H. Gillespie, T. A. Brown, OPTICS ELEMENTS FOR MODELING ELECTROSTATIC LENSES
AND ACCELERATOR COMPONENTS I. EINZEL LENSES,  Proceedings of the 1997 Particle Accelerator Conference,
https://doi.org/10.1109/PAC.1997.751273

[2] G. H. Gillespie, T. A. Brown, OPTICS ELEMENTS FOR MODELING ELECTROSTATIC LENSES
AND ACCELERATOR COMPONENTS II. ACCELERATOR COLUMNS, Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, Volume 427, Issues 1–2, May 1999, Pages 315-320, https://doi.org/10.1016/S0168-9002(98)01543-5

[3] G. H. Gillespie, T. A. Brown, Optics elements for modeling electrostatic lenses and accelerator components: III. Electrostatic deflectors
Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms
Volume 172, Issues 1–4, October 2000, Pages 338-343, https://doi.org/10.1016/S0168-583X(00)00097-5

[4] G. H. Gillespie, Optics Elements for Modeling Electrostatic Lenses and Accelerator Components IV. Electrostatic Quadrupoles and Space Charge Modeling,
19th International Linear Accelerator Conference, Chicago, IL, USA, 23 - 28 Aug 1998, pp.150
https://cds.cern.ch/record/740868

[5] C. K. Allen, A Software Engineering Approach to Particle Beam Simulation,
 Los Alamos National Laboratory Technical Report LA-UR-00-3895
