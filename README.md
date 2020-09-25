# PyPrometheus
PyPrometheus generates source code to compute chemical source terms. These source terms appear in the species conservation equations
for reactive flows, so their implementation is needed for reactive-flow simulations. 

# Dependencies
Chemical source terms have models for reaction rates (e.g., Arrhenius rate coefficients) and thermodynamic properties (e.g., NASA polynomials). 
PyPrometheus manages source-term parameterization through Cantera (https://www.cantera.org). We recommend the latest version. There are no other dependencies. 

# Installation
Once you've installed Cantera, follow these steps to get PyPrometheus up and running.
- cd build
- make

# Running PyPrometheus
## Generating Source-term source code

## Using generated code
