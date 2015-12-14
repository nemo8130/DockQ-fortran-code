# =====================================================================================================================================================
# README file for DockQ
# =====================================================================================================================================================
# THERE ARE TWO VERSIONS OF THE PROGRAM
#
# 1. calDockQNaln.bash : The default 'simpler' version. It assumes that the two input PDB files are already aligned. i.e., the corresponding residues in the two files have identical residue sequence number and chain identifier.
#
# 2. calDockQ.bash : This version can perform a sequence alignment prior the calculation of the DockQ score if specified.
#
# =====================================================================================================================================================
#
# Requirement: A fortran90 compiler (preferred: ifort, gfortran)
# PERL v.5 or higher 
# OS: Linux
# Optional requirement:
# needle as part of EMBOS package (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html)
# If chosen to implement the sequence alignment option prior to the calculation of DockQ.
#
#
#
# =========================================================================================================
# =========================================================================================================
# INSTALLATION:
# =========================================================================================================
# =========================================================================================================
#
# Run ./install.bash with the chosen fortran compiler name (ifort, gfortran) as the only argument ($fcomp)
# The default is: ifort (if nothing provided as a cmd)
# The script will check whether $fcomp is installed and will exit if not found.
# It will prompt 'Installation Successful' for a successful run and prompt the usage.
# =====================================================================================================================================================
#
# =========================================================================================================
# =========================================================================================================
# ./calDockQNaln.bash: 
# =========================================================================================================
# =========================================================================================================
#
#
# USAGE: ./calDockQNoaln.bash model.pdb native.pdb -CA (1/0) 
#
#
# OPTIONAL (flag) ARGUMENTS:
#
# -CA : (0/1) : Choice between CA-only (~6.5X FASTER) or all-atom version of the program
#  0: all-atoms
#  1: CA-only
#  Default: CA-only
#
# DEFAULT-UASGE: ./calDockQNoaln.bash model.pdb native.pdb          (will run the CA-only version)
#
# The program also outputs a RASMOL script to display the interface contact network of the native structure for visualization. 
# The same could also be generated for the model by reverting the order of the two input files.
# To view the network, use: rasmol -script native_filename.spt
# ======================================================================================================================================================
# ./calDockQ.bash: 
#
# =====================================================================================================================================================
#
#
# USAGE: ./calDockQ.bash model.pdb native.pdb -CA (1/0) -aln (1/0) 
#
# OPTIONAL (flag) ARGUMENTS:
#
# -CA : (0/1) : Choice between CA-only (~6.5X FASTER) or all-atom version of the program
#  0: all-atoms
#  1: CA-only
#  Default: CA-only
#
# -aln : (0/1) : Determines whether to perform a SEQUENCE ALIGNMENT (between the MODEL and the NATIVE) prior to the calculation
# 0: The DEFAULT alignment of the equivalent residues is considered based on the combined matching of RESIDUE NUMBER and CHAIN Ids from the two PDB files
# 1: A sequence alignment (using on the Needleman–Wunsch algorithm) is performed from the sequences extracted from the two PDB files 
#    and equivalent residues are defined based on this alignment (needle NEEDS TO BE INSTALLED to avail this option)
#
# DEFAULT_VALUE (for -aln): 0 (Sequence Alignment not to be performed) 
#
#
# DEFAULT-UASGE: ./calDockQ.bash model.pdb native.pdb          (will run the CA-only version, with NO alignment)
# The program also outputs a RASMOL script to display the interface contact network of the native structure for visualization. 
# The same could also be generated for the model by reverting the order of the two input files.
# To view the network, use: rasmol -script native_filename.spt
#
#
# ======================================================================================================================================================
#
