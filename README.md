ChemReader
==========

Library for importing chemical data; focus is on electronic structure and the following file formats:

* Molden
* Gaussian logs and cube
* GAMESS
* ADF

Support for the following data formats will be available as well:

* PDB
* ShellX

The project is dependent on Parsley a C++ parsing framework supporting context-sensitive processing.

Parsed and computed data include:

* unit cell
* affine transform
* atom name, element type and coordinates
* atomic orbitals as GTO or STO basis sets
* molecular obitals as LCAO
* density and spin matrix(computed by the library if required)
* electrostatic potential
* vibrational modes
* scalar fields

Minimal support for biochemical data will be available as well:

* backbone extraction
* residues
* secondary and tertiary structures
