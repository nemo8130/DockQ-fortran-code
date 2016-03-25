# DockQ-fortran-code
DockQ: A quality measure for PPI complexes based on CAPRI evaluation protocol 

Requires a fortran90 compiler (prefered: ifort)

Installation

git clone https://github.com/nemo8130/DockQ-fortran-code

cd DockQ-fortran-code/DockQ

./install.bash (fortran90 compiler; Default: ifort)

Run with ./DockQ.exe <model> <native>

Example

bash$ ./DockQ.exe ./EXAMPLE_PDBS/model.pdb ./EXAMPLE_PDBS/native.pdb 
 *********************************************************
 *                       DockQ                           *
 *   Scoring function for protein-protein docking models *
 *   Statistics on CAPRI data:                           *
 *    0.00 <= DockQ <  0.23 - Incorrect                  *
 *    0.23 <= DockQ <  0.49 - Acceptable quality         *
 *    0.49 <= DockQ <  0.80 - Medium quality             *
 *            DockQ >= 0.80 - High quality               *
 *   Reference: Sankar Basu and Bjorn Wallner, DockQ:... *
 *   For comments, please email: bjornw@ifm.liu.se       *
 *********************************************************

   Number of equivalent residues at the interface:    103   (receptor:  A)     88   (ligand:   B)
      <=== defined at an all-atom atomic cutoff of  10.00   Å (iRMS) ===>
    fnat:    0.533:   32      correct out of   60  native contacts at the receptor-ligand interface
  fnonat:    0.238:   10   non-native out of   42   model contacts at the receptor-ligand interface
Fnat=   0.533  LRMS_bb=       1.516  IRMS=       1.197  DockQ=     0.70442
CAPRI_class (fnat, iRMS, LRMS): Medium              
CAPRI_class (DockQ)           : Medium              
 outfile: fort.245
 rename the outfile as you wish

