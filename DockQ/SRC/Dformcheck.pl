#!/usr/bin/perl

#====================================================================================
#  Check compatibility of the input file
#------------------------------------------------------------------------------
#  It should contain two and exactly two polypeptide chains
#  Residue numbers should always follow an ascending order 
#====================================================================================

$flagUlt=0;

`rm -f formch.out`;

open (OUTF,">formch.out");

$inppdb = $ARGV[0] || die "Enter PDB file\n";
chomp $inppdb;
open (PDBFILE,"<$inppdb");
@fdat = <PDBFILE>;

print "\n\nINPUT PDB : $inppdb\n\n";
$log = $inppdb.'.log';
`rm -f $log`;
open (LOG,">$log");

@pcoords = grep(/^ATOM\s/,@fdat);
%h = ();

foreach $k (@pcoords)
{
chomp $k;
$h{substr($k,21,1)}++;
}

@chain = sort keys %h;
$nch = @chain;
#print $nch,"\n";

if ($nch == 2)
{
print "\n\nThe input pdb file contains \'TWO\' polypeptide chains ($chain[0] & $chain[1]) as required. \n\n";
}
else
{
$flagUlt=1;
print LOG "=======================================================================================================\n";
print LOG "\n\nThe input pdb file does not contain exactly \'TWO\' polypeptide chains and thus not acceptable for calculation \n [Rather it contains $nch chains]\n";
	foreach $c (@chain)
	{
	print LOG $c,"\n";
	}
print "\nProgram will exit: $nch chains\n";
print LOG "It should contain \'TWO and exactly TWO\' polypeptide chains\n\n";
print LOG "=======================================================================================================\n";
}

print OUTF "FLAG: $flagUlt\n";
print "FLAG: $flagUlt\n";

	if ($flagUlt==1)
	{
	print "LOG FILE: $log\n";
	}
	elsif ($flagUlt==0)
	{
	`rm -f $log`;
	}







