# GUIDEMOL: a Python graphical user interface for molecular descriptors based on RDKit

## Brief description

GUIDEMOL is a Python computer program based on the RDKit software (RDKit:
Open-source cheminformatics, <https://www.rdkit.org>) to process molecular
structures and calculate molecular descriptors with a graphical user interface
using the tkinter package. It can calculate descriptors already implemented in
RDKit as well as grid representations of 3D molecular structures using the
electrostatic potential or voxels. A CLI is also provided for the calculation
of grid representations.

## Installation

GUIDEMOL requires the installation of RDKit. Updated recipes for the
installation of RDKit (e.g. via Conda or PIP) in various operating systems are
provided at <https://www.rdkit.org/docs/Install.html>

## How to use

Download the guidemol_vx.x.x.py file to some directory, navigate via command
line to this directory.
Depending on your installation you may need to activate an rdkit environment by
typing ‘activate my-rdkit-env’.

Then type

```
python guidemol_vx.x.x.py
```

or

```
python3 guidemol_vx.x.x.py
```

### Generate 3D structures

1. Select the operations to be performed before the generation of the 3D models
   (“Remove metals”, “Remove smaller fragments”, or  “Uncharge”).
1. From v1.1.0 there is an option to enable/disable 3D structure optimization
   with the MMFF94 force field. Disabled: 3D structure generated by RDKit ETKDG
   method. Enabled: 3D structure generated by RDKit ETKDG method and further
   optimization with the MMFF94 force field.
1. Click on “3D from SMILES”.
1. Chose the input file with the SMILES strings, then an output file to save
   the 3D structures in the MDL SDFile format.

Note: GUIDEMOL, in its current version, does not explore the generation of
multiple conformers. However, code is publicly available to generate a set of
conformers and get the lowest energy structure, e.g. from
<https://github.com/pstjohn/bde>, which can be easily inserted into GUIDEMOL.

### Calculate molecular descriptors

1. Under “Molecular descriptors” select the descriptors to be calculated.
   For the “Grid of potential” and “Grid of voxels” descriptors you may
   change the default parameters.
1. If the 3D structure misses hydrogen atoms select “Add hydrogens”.
   Hydrogen atoms will be added without changing the 3D coordinates of the other
   atoms.
1. Click on “Calculate descriptors”.
1. Chose the input file with the 3D structures in the MDL SDFile format, then
   an output file to save the descriptors in the CSV format.

A log file is also created in the same directory to register the successful
processing of each molecule, as well as failures.

### Calculation of EP and voxel grids via the command-line interface

The CLI uses the command

```
guidemol_vx.x.x.py [-h] -i I [-o O] -g {pot,vox} \
    -p {Gasteiger,MMFF,MR,LogP,AtomicNumber} \
    -cx CX -cy CY -cz CZ -r R -w W -sx SX -sy SY -sz SZ
```

with the specification of the input file, the output file (optional), the type
of grid (EP or voxels), the atomic property, the xyz coordinates of the grid
center, the cube size, the cutoff distance (in EP grids) or the sigma value (in
voxel grids) and the size of the grid in each direction.

options:

```
-h, --help                                 show this help message and exit
-i I                                       Input file
-o O                                       Output file (optional)
-g {pot,vox}                               Type of grid
-p {Gasteiger,MMFF,MR,LogP,AtomicNumber}   Atomic property
-cx CX                                     x coordinate of grid center
-cy CY                                     y coordinate of grid center
-cz CZ                                     z coordinate of grid center
-r R                                       Cube size (Angstrom)
-w W                                       Atom reach (cutoff distance or sigma)
-sx SX                                     Size of grid (xx axys)
-sy SY                                     Size of grid (yy axys)
-sz SZ                                     Size of grid (zz axys)
```

For example, to calculate an EP grid with Gasteiger charges, centered at
(0,0,0), with size of each side 5 Å, a resolution of 0.5 Å and a cutoff of
3 Å for molecules in the file input.sdf and save the result in the file
outfile.csv:

```
guidemol_vx.x.x.py -i input.sdf -g pot \
    -p Gasteiger -cx 0 -cy 0 -cz 0 -r 0.5 -w 3 -sx 5 -sy 5 -sz 5 -o outfile.csv
```

## Author contact

Joao Aires-de-Sousa (jas_at_fct.unl.pt)

## Acknowledgements

This work was supported by the Associate Laboratory for Green Chemistry –
LAQV, which is financed by national funds from Fundação para a Ciência e
Tecnologia (FCT/MCTES), Portugal, under grant UIDB/50006/2020. It was also
supported by EU ESIF funds (P2020 SI&IDT 17/SI/2019/47212).

Aires De Sousa, J. GUIDEMOL: A Python Graphical User Interface for Molecular
Descriptors Based on RDKit. *Mol. Inf.* **2023**, *42*, e202300190.
<https://doi.org/10.1002/minf.202300190>.
