#!/bin/bash

path=`pwd`

fcomp=$1

if [ "$#" == "0" ]; then
fcomp=ifort
echo "DEFAULT FORTRAN COMPILER CHOSEN:" $fcomp
else
echo "YOUR CHOICE OF FORTRAN COMPILER:" $fcomp
fi

chfort=`which $fcomp`

        if [ -z "$chfort" ]; then
	echo "=============================================="
	echo "=============================================="
	echo "=============================================="
        echo $fcomp "NOT found in the path"
	echo "=============================================="
	echo "=============================================="
	echo "=============================================="
	exit;
	else
	echo $fcomp "FOUND TO BE INSTALLED AT " $chfort
        fi


$fcomp $path/src/DockQfast.f -o $path/DockQ.exe -O3
#$fcomp $path/SRC/DockQCA.f -o $path/SRC/DockQCA.exe -O3

echo "Installation Successful" 
cat README

