# PyPrometheus
PyPrometheus generates source code to compute chemical source terms. These source terms appear in the species conservation equations
for reactive flows, so their implementation is needed for reactive-flow simulations. 

# Dependencies
Chemical source terms have models for reaction rates (e.g., Arrhenius rate coefficients) and thermodynamic properties (e.g., NASA polynomials). 
PyPrometheus manages source-term parameterization through Cantera (https://www.cantera.org). We recommend the latest version. There are no other dependencies. 

# Installation
Once you've installed Cantera, follow these steps to get PyPrometheus up and running.
- `cd build`
- Open `makefile` and substitute your local Cantera path in line 1. 
- `make`

Note: a CMake-based build system is underway.

# Running PyPrometheus
There are two steps to using PyPrometheus: (i) generating the models and (ii) using them in reactive flow simulations.

## Generating source-term source code
Let's start with some preliminary definitions. We call the generated source code **Prometheus thermochemistry kernels** (PTK). To generate a PTK, we need an input config file. This file only needs to include two options:
- `mech`: The name of the mechanism. Note that this has to match the name of the Cantera input file (in cti, xml, or yaml formats).
- `language`: The language in which the source term will be implemented. PyPrometheus supports Python and C++.
We provide an example (inputs/sample_inputs.config). 

After customizinng the input file, we can generate a PTK by typing:
`build/Prometheus.exe <my_input.config>`
Prometheus will then produce a `mech.py` or `mech.H.`

## Using generated PTKs (source code)
At this point, we have a PTK, but how do we even use it? We provide an in-depth tutorial via the jupyter notebook `examples/tutorial_compute_source_terms.ipynb`.
