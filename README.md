# tinker
TINKER Software Tools for Molecular Design

This "amoeba2" branch of Tinker contains energy and gradient routines to implement the improved functional form of AMOEBA 2.0
The major changes are:
  1. Addition of Hydrogen-like Electrostatic Potential for Electrostatic and Polarization interactions
  2. Damped Dispersion PME 
  3. Overlap and Multipole-based Pauli Repulsion
  
This is a beta version and it currently does not support "old" AMOEBA calculations.  This feature is coming soon.
The code also includes virial calculations for NPT simulations

Finally, the first tests have been added to tests/amoeba2test
