#!/usr/bin/perl

# Create pdb file with connections defining the edges in a contact core network

$scon = $ARGV[0] || die "Enter .scon file\n";	# subgraph	/   intsg	/ scon
$pdb = $ARGV[1] || die "Enter .pdb file\n";     # mapped according to .scon

chomp $scon;
chomp $pdb;

open (INP,"<$scon");
@dat = <INP>;
open (PDB,"<$pdb");
$snet = $scon;

if ($scon =~ m/scon/)
{
$snet =~ s/\.scon/-unet\.pdb/g;
}

open (OUT,">$snet");
@dpdb = <PDB>;

$atno = 0;
@coll = ();

foreach (@dpdb)
{
chomp $_;
$atno++;
$atno = sprintf("%5d",$atno);
print OUT substr($_,0,6),$atno,substr($_,11, ),"\n";
$str = substr($_,0,6).$atno.substr($_,11, );
@coll = (@coll,$str);
}

@dpdb = ();
@dpdb = @coll;

@CA = grep(/CA/,@dpdb);

if ($scon =~ m/\.scon/)
{
$nfile = $scon;
$nfile =~ s/\.scon/\.snet/g;
open (ND,"<$nfile") || die "$nfile not found\n";
}

@nun = <ND>;

$sfile = $scon;

if ($scon =~ m/\.scon/)
{
$sfile =~ s/\.scon/\.spt/g;
}

open (SPT,">$sfile");

@anumber = ();
@resn = ();
$c = 0;
foreach $n (@nun)
{
chomp $n;
#print $n,"\n";
@resn = (@resn,int(substr($n,0,4)));
$c++;
}

#print $c,"\n";

%h1 = ();

foreach $k (@CA)
{
chomp $k;
$rest = substr($k,22,4).'-'.substr($k,17,3).'-'.substr($k,21,1);
$anum = int(substr($k,7,4));
#print $rest,"\n";
	foreach $n (@nun)
	{
		if ($rest eq substr($n,0,9))
		{
#		print "$rest -> $anum\n";
		@anumber = (@anumber,$anum);
		$h1{$anum}++;
#		print SPT "select atomno=$anum\n";
		}
	}
}

@keys = sort {$a<=>$b} keys %h1;
@anumber = @keys;
$ll = @anumber;

#print $ll,"\n";

$f10 = int($ll/10);
$f1st = 10*$f10;

$rm = $ll-$f1st;

@ar1 = ();

for $i (0..$f1st-1)
{
@ar1 = (@ar1,$anumber[$i]);
}

for $i ($f1st..$ll-1)
{
@ar2 = (@ar2,$anumber[$i]);
}

foreach (@ar1)
{
#print $_,"\n";
}
#print "\n\n";
foreach (@ar2)
{
#print $_,"\n";
}

$la1 = @ar1;

#print SPT "zap\n";
print SPT "load pdb $snet\n";
print SPT "wireframe off\n";
print SPT "ribbon\n";
print SPT "color chain\n";


for ($i=0;$i<$la1;$i+=10)
{
$j = $i+9;
print SPT "select ";
#print "select ";
	for $k ($i..$j-1)
	{
	print SPT "atomno=$anumber[$k],";
#	print "atomno=$anumber[$k],";
	}
print SPT "atomno=$anumber[$j]";
#print "atomno=$anumber[$j]";
print SPT "\nspacefill 120\n";
#print "\nspacefill 120\n";
#print SPT "label \%r\n";
print SPT "color atom white\n";
#print "color atom white\n";
}

if (scalar(@ar2)>0)
{
print SPT "select ";
#print "select ";
	for $i (0..scalar(@ar2)-2)
	{
	print SPT "atomno=$ar2[$i],";
#	print "atomno=$ar2[$i],";
	}
print SPT "atomno=$ar2[scalar(@ar2)-1]";
#print  "atomno=$ar2[scalar(@ar2)-1]";
print SPT "\n";
#print  "\n";
print SPT "spacefill 120\n";
#print  "spacefill 120\n";
#print SPT "label \%r\n";
print SPT "color atom white\n";
#print "color atom white\n";
}

@col = ();

foreach $k (@dat)
{
chomp $k;
$n1 = int(substr($k,0,4));	# 12
$c1 = substr($k,8,1);	# 12
$n2 = int(substr($k,11,3));	# 21
$c2 = substr($k,19,1);	# 21
#print "$n1  $c1  <=>  $n2  $c2\n";
	foreach $p (@dpdb)
	{
	chomp $p;
	$atomno = int(substr($p,6,5));
	$atype = substr($p,13,3);
	$chain = substr($p,21,1);
	$resno = int(substr($p,22,4));	
#	print "$atype  $resno  $chain\n";
#	print $atomno,"\n";
		if ($n1 == $resno && $atype eq 'CA ' && $c1 eq $chain)
		{
		$atn1 = $atomno;
		}
		if ($n2 == $resno && $atype eq 'CA ' && $c2 eq $chain)
		{
		$atn2 = $atomno;
		}
	}
printf OUT "%6s %4d %4d\n",'CONECT',$atn1,$atn2;
#print  "CONECT 1:$atn1 2:$atn2\n";
$iatn1 = int($atn1);
$iatn2 = int($atn2);
print SPT "select (atomno=$iatn1) or (atomno=$iatn2)\n";
print SPT "wireframe 40\n";
print SPT "color bond yellow\n";
@col = (@col,int($atn1),int($atn2));
}

close OUT;
