# Tinker: Software Tools for Molecular Design

<H2><B>Introduction</B></H2>

The Tinker molecular modeling software is a complete and general package for molecular mechanics and dynamics, with some special features for biopolymers. Tinker has the ability to use any of several common parameter sets, such as Amber (ff94, ff96, ff98, ff99, ff99SB), CHARMM (19, 22, 22/CMAP), Allinger MM (MM2-1991 and MM3-2000), OPLS (OPLS-UA, OPLS-AA), Merck Molecular Force Field (MMFF), Liam Dang's polarizable model, the AMOEBA (2004, 2009, 2013, 2017, 2018) polarizable atomic multipole force field, and a prototype of the new HIPPO force field calibrated against SAPT QM interaction energy components. Parameter sets for other widely-used force fields are under consideration for future releases.

The Tinker software contains a variety of interesting algorithms such as: flexible implementation of atomic multipole-based electrostatics with explicit dipole polarizability, various continuum solvation treatments including several generalized Born (GB/SA) models, generalized Kirkwood implicit solvation for AMOEBA, an interface to APBS for Poisson-Boltzmann calculations, efficient truncated Newton (TNCG) local optimization, surface areas and volumes with derivatives, free energy calculations via the Bennett Acceptance Ratio (BAR) method, normal mode vibrational analysis, minimization in Cartesian, torsional or rigid body space, symplectic RESPA multiple time step integration for molecular dynamics, velocity Verlet stochastic dynamics, pairwise neighbor lists and splined spherical energy cutoff methods, particle mesh Ewald (PME) summation for partial charges and polarizable multipoles, a novel reaction field treatment of long range electrostatics, fast distance geometry metrization with better sampling than standard methods, Elber's reaction path algorithm, potential smoothing and search (PSS) methods for global optimization, Monte Carlo Minimization (MCM) for efficient potential surface scanning, tools for fitting charge, multipole and polarization models to QM-based electrostatic potentials and more....

<H2><B>Current Release</B></H2>

Tinker 8 is a major new release of the Ponder Lab tool set for molecular mechanics and dynamics calculations. An important change in this new version is the switch from old-style common blocks to Fortran modules. Use of modules and greatly increased use of dynamic memory allocation means Tinker can now support very large molecular systems. Tinker 8 also implements improved OpenMP parallelization throughout many parts of the code. Additional big improvements include parallel neighbor list building and updating, and big reduction in iteration needed to converge AMOEBA polarization via an efficient PCG solver. Other changes from the previous Tinker version include new and updated force field parameter sets and numerous minor additions and bug fixes, many of them suggested by users of the package. Please note that as with prior new releases, version 8 is neither backward nor forward compatible with earlier versions of Tinker. In particular, older versions of parameter files should not be used with Tinker 8 executables and vice versa.

While we strongly suggest users switch to Tinker 8 with its many important new features and bug fixes, we provide download links at https://dasher.wustl.edu/tinker/ for prior stable versions, Tinker 7.1.3, Tinker 6.3.3, Tinker 5.1.9 and Tinker 4.3. Tinker 6 and later is OpenMP parallel and written in Fortran 95, while Tinker 4 and 5 are in serial, extended Fortran 77.

