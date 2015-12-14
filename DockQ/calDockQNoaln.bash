#!/bin/bash
#
# The Script to calculate the combined Quality measure 'DockQ' for a given pair of model and native structures. 
# Assuming either an alignemnet based on matching of residue sequence numbers and chain identifiers from the two files (Default)
# Or performing a Needleman–Wunsch global alignment on the two sequences extracted from the two PDB files.
# Simpler version without alignment; rasmol script is generated by default 

path=`echo ${0/\/calDockQNoaln.bash/}`
#echo "==============================================================================================================================="
#echo "CURRENT PATH:" $path
#echo "==============================================================================================================================="

#============================== DEFAULT SETTINGS FOR THE FLAGS =========================
CA=1
#=======================================================================================


model=$1;		# Atomic ccordinates of the model (extension: .pdb / .PDB 
native=$2;		# Atomic coordinates of the native (extension: .pdb / .PDB)
flag1=$3;		# An optional flag (0/1) to input the choice of version (CA-only / all-atoms) 
CA=$4;			# 0: all-atoms; 1: CA-only (Default: 1 ~ CA-only)
rasscript=1
aln=0
if [ "$#" -lt "2" ]; then 
echo
echo
echo
echo
echo "==============================================================================================================================";
cat usageNaln.prompt
echo "==============================================================================================================================";
echo "Also read the README.txt file"
echo "===============================================================================================================================";
echo
echo
echo
echo
echo
echo
exit;
elif [ "$#" -eq "2" ]; then          # Default settings
CA=1
echo "================================================================================="
echo "ONLY THE TWO MANDETORY CMD-ARGUMENTS PROVIDED: MODEL.pdb and NATIVE.pdb"
echo "Default value (0) considered for the -aln flag: SEQUENCE ALIGNMENT NOT TO BE PERFORMED"
echo "================================================================================="
elif [ "$#" -ge "3" ] && [ "$#" -lt "5" ]; then 
#===============================================================================
#===============================================================================
	if [ "$4" -eq "$4" ] 2>/dev/null
	then
		run=1
#	    echo "$4 is an integer !!"
	else
	    echo "ERROR: Option to the flag -CA must be an integer (0/1)."
	    cat usageNaln.prompt
	    exit 1
	fi
#===============================================================================
#===============================================================================
	if [ "$flag1" == "-CA" ]; then
		if [ "$CA" -ne "0" ] && [ "$CA" -ne "1" ]; then
		echo "Incorrect value provided for the -aln flag " $aln
		echo "Should be either 0 or 1"
		echo "The program will exit"
		exit;
		fi
	else 
	echo "Incorrect Input: See the optional flags:"
	cat usageNaln.prompt
	exit;
	fi
elif [ "$#" -ge "5" ]; then
	echo "Toomany Input Arguments: See the optional flags:"
	cat usageNaln.prompt
	exit;
fi

#echo "========================== INPUT TO THE OPTIONAL FLAGS ============================="
#echo 'CA=' $CA
#echo 'rasscript=' $rasscript
#echo "===================================================================================="

#exit;

for i in `cat $path/progs.list`
do
	if [ -e "$path/SRC/$i" ]; then
#	echo "$path/SRC/$i FOUND IN THE PATH";
	run=1
	else 
	echo "$path/SRC/$i NOT FOUND IN THE PATH";
	echo "THE PROGRAM WILL EXCEED";
	exit;
	fi
done

#echo "ALL Files found in the path: The program is RUNNING"

fst1=`$path/SRC/DformcheckNaln.pl $native | grep "FLAG:"`
fst2=`./SRC/DformcheckNaln.pl $model | grep "FLAG:"`

f1st=`echo ${fst1/FLAG: /}`
f2st=`echo ${fst2/FLAG: /}`

#echo $f1st
#echo $f2st

if [ "$f1st" == "0" ] && [ "$f2st" == "0" ] ; then	# A flag (0/1) to check whether both input PDB files contain 'two and exactly two' polypeptide chains.
run=1
#echo "==========================================================================================="
#echo "==========================================================================================="
#echo "Both PDB files contain two and exactly two chains. The Program will proceed"
#echo "==========================================================================================="
#echo "==========================================================================================="
else
echo "==========================================================================================="
echo "==========================================================================================="
echo "At least One of the Two PDB files does contain two and exactly two polypeptide chains. 
                              Program will exit. 
               Both PDB files should contain two and exactly two chains.
                              
                                        OR

The input pdb file contains atleast one instance of unsorted residues [residues sorted in a descending order]
                ===============================================================
                Residue sequence numbers should always follow an ascending order
                ===============================================================
                CHECK THE LOG FILE FOR DETAILS: $model.log / $native.log" 
echo "==========================================================================================="
echo "==========================================================================================="
exit;
fi


#=================================== SPECIFY EXECUTABLE VERSION ==================================

	if [ "$CA" == "0" ]; then 
	echo "==========================================================================================="
	echo "You preferred the 'all-atoms' vesrion"
#	echo "==========================================================================================="
	exe=$path/SRC/DockQ.exe
	elif [ "$CA" == "1" ]; then
	echo "==========================================================================================="
	echo "You preferred the 'CA' vesrion"
#	echo "==========================================================================================="
	exe=$path/SRC/DockQCA.exe
	fi

#=================================================================================================

if [ "$aln" == "1" ]; then 
echo "==========================================================================================="
echo "==========================================================================================="
echo "YOU REQUESTED TO PERFORM A SEQUENCE ALIGNMENT PRIOR TO COMPUTING THE MODEL QUALITY; GREAT!"
echo "YOU NEED TO HAVE 'needle' INSTALLED AS PART OF THE EMBOS PACKAGE IN YOUR WORKSTATION TO AVAIL THIS OPTION"
echo "It could be installed from: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html"
echo "==========================================================================================="
echo "==========================================================================================="
chneedle=`which needle`
	if [ -z "$chneedle" ]; then
	echo "====================================="
	echo "====================================="
	echo "====================================="
	echo "needle not found in the path"
	echo "====================================="
	echo "====================================="
	echo "====================================="
	exit;
	else
	echo "======================================================"
	echo "======================================================"
	echo "======================================================"
	echo "needle FOUND TO BE INSTALLED AS" $chneedle
	echo "======================================================"
	echo "======================================================"
	echo "======================================================"
	fi
echo "==========================================================================================="
echo "==========================================================================================="
#	$path/SRC/renumberF1.pl $model;
#	$path/SRC/renumberF1.pl $native;

	$path/SRC/renumber_pdb.pl $model;
	$path/SRC/renumber_pdb.pl $native;

	modelf1=`echo $model`.renum;
	nativef1=`echo $native`.renum;
	model_renumbered=`echo $modelf1`.fixed;

	echo $model_renumbered;
	echo "=========================================================="
	echo "=========================================================="
	echo "=========================================================="
	echo "PERFORMING THE ALIGNMENT NOW: "
	echo "=========================================================="
	echo "=========================================================="
	$path/SRC/fix_numbering.pl $modelf1 $nativef1;
	echo "=========================================================="
	echo "=========================================================="
	echo "=========================================================="
	echo "=========================================================="		        
	out=`$exe $model_renumbered $nativef1`;
elif [ "$aln" == "0" ]; then
#echo "==========================================================================================================================================="
#echo "==========================================================================================================================================="
#echo "YOU PREFERED THE DEFAULT ALIGNMENT PROVIDED IN THE PDB FILES: "
#echo "THE PROGRAM WILL ASSUME EQUIVALENT RESIDUES IN THE TWO (MODEL AND NATIVE) FILES TO HAVE IDENTICAL RESIDUE SEQUENCE NUMBERS (AND) CHAIN ID's"
#echo "==========================================================================================================================================="
#echo "==========================================================================================================================================="
out=`$exe $model $native`;
fi


echo $out > log.txt;
mv fort.245 $model.DockQ

#echo "==========================================="
#echo 
#echo 
#echo 
#echo 
#echo "==========================================="
echo "CHECK OUTFILE: " $model.DockQ
#echo "==========================================="
#echo 
#echo 
#echo 
#echo 
#echo 
echo "YOUR OUTPUT:"
#echo 
echo "==========================================================================================="
cat $model.DockQ
echo "==========================================================================================="
#echo 
#echo 
#echo 
#echo 
#echo 
#echo 
#echo 
#echo 
#echo 
#echo "==========================================="



	if [ "$rasscript" == 1 ]; then

#                echo "GENERATING RASMOL SCRIPT FOR THE FIRST MODEL";
                scon=$native.scon;
                snet=$native.snet;
                spt=$native.spt;

#		echo $scon
#		echo $snet 

                mv fort.23 $scon;
                mv fort.26 $snet;

			if [ "$aln" == "0" ]; then	
			pdbras=$native
	                $path/SRC/createrasscript.pl $scon $pdbras;
			elif [ "$aln" == "1" ]; then
			pdbras=$nativef1
	                $path/SRC/createrasscript.pl $scon $pdbras;
			fi

		echo "RASMOL SCRIPT TO SHOW THE INTERFACE CONTACT NETWORK FOR THE NATIVE:" $spt

	fi

# CLEANING INTERMEDIATE FILES

rm fort.*



