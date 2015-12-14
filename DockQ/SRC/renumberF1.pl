#!/usr/bin/perl

$pdbfile = $ARGV[0] || die "Enter a PDB file containing two chains\n";
chomp $pdbfile;

$path=`pwd`;
chomp $path;
print $path,"\n";

`$path/SRC/pdb2resWMchainfM1.pl $pdbfile`;

open (INP,"<$pdbfile");
@dat = <INP>;

@atoms = grep(/^ATOM\s+/,@dat);
@atoms1 = ();
@atoms2 = ();

$chain1 = substr($atoms[$i],21,1);
$chain2 = substr($atoms[scalar(@atoms)-1],21,1);

#print "$chain1  $chain2\n";

foreach $a (@atoms)
{
chomp $a;
	if (substr($a,21,1) eq $chain1)
	{
	@atoms1 = (@atoms1,$a);
	}
	elsif (substr($a,21,1) eq $chain2)
	{
	@atoms2 = (@atoms2,$a);
	}
} 

@atomsc1 = @atoms1;
@atomsc2 = @atoms2;

$l1 = @atomsc1;
$f1 = substr($atomsc1[0],22,4);
$rn1 = $f1;

#print "$l1   $rn1\n";

#print $rn1,"\n";

$outpdb = $pdbfile.'.renum';
$outres = $outpdb.'.res';

open (OUT,">$outpdb");

for $i (0..$l1-1)
{
chomp $atomsc1[$i];
$fp = substr($atomsc1[$i],0,22);
$rn = int(substr($atomsc1[$i],22,4));
$sp = substr($atomsc1[$i],26, );
$rn = $rn - ($rn1-1);
printf OUT "%22s%4d%28s\n",$fp,$rn,$sp;
}

$rnf1 = substr($atomsc1[scalar(@atomsc1)-1],22,4);

#print $rnf1,"\n";

$l2 = @atomsc2;
$rn1 = int(substr($atomsc1[0],22,4));
$rn2 = int(substr($atomsc2[0],22,4));
$rnadd2 = $rnf1-($rn1-1) - ($rn2-1);			

#print "$l1   $rn1\n";

#print $rn1,"\n";
$rnf2 = 0;

for $i (0..$l2-1)
{
chomp $atomsc2[$i];
$fp = substr($atomsc2[$i],0,22);
$rn = int(substr($atomsc2[$i],22,4));
$sp = substr($atomsc2[$i],26, );
$rn = $rn + $rnadd2;
$rnf2 = $rn;
printf OUT "%22s%4d%28s\n",$fp,$rn,$sp;
}

#print $rnf2,"\n";

`$path/SRC/pdb2resWMchainfM1.pl $outpdb`;

$resf = $pdbfile.'.res';
chomp $resf;

open (RES1,"<$resf");
@res1 = <RES1>;
close RES1;

open (RES2,"<$outres");
@res2 = <RES2>;
close RSE2;

$map = $pdbfile.'.map';
open (MAP,">$map");

$l1 = @res1;
$l2 = @res2;

#print "$l1   $l2\n";

if ($l1 != $l2)
{
print "SOMETHING WRONG : $pdbfile\n";
}
else
{
	for $i (0..$l1-1)
	{
	chomp $res1[$i];
	chomp $res2[$i];
	print MAP "$res1[$i] ->  $res2[$i]\n";
	}
}

