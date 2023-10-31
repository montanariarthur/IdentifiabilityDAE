# Identifiability of Differential-Algebraic Equations
Codes for parameter identifiability analysis of models described by differential-algebraic equations (DAE).

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Usage

The following codes are direct implementations of the parameter identifiability analysis conducted for the DAE systems studied in Ref. [1].

- `main_idft_chemicalreactor` : Analysis of the identifiability of parameter `Tc` of a chemical reactor model for different types of measurement signals (sensors).

- `main_idft_pendulum` : Analysis of the identifiability of different combinations of parameters of a pendulum equation.

- `main_idft_linearDAE.m`: Analysis of the identifiability of a linear DAE model.

- `DAEobsvmatrix` : Function that computes the observability matrix of a DAE model.


The above codes are dependent on the following functions:

- `DAE models` : Folder contains implementations of the DAEs for the systems investigated in Ref. [1], which are used for numerical integration.
- `sym2str` : Converts symbolic variables to string variables. 


# References
1. A. N. Montanari, F. Lamoline, R. Bereza, J. Goncalves. Identifiability of differential-algebraic equations. *Under review* (**2023**).
