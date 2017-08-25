#!/usr/bin/perl -w
use POSIX qw(ceil floor);
use File::Find;

$PeptideAlignment=shift;

open (FILE, $PeptideAlignment)|| die "ERROR: CEBA couldn't open the file $PeptideAlignment!";
while ($line =<FILE>)
{ 
    chomp $line;
    if ($line =~ /^>/)
	{
		#print "SEQ".$index." ".$line."\n";
		@lineContent=split(">",$line);
		$PepName=$lineContent[1];
		#print "HLANamePart1 ".$HLANamePart1."\n";
		#print "HLANameCorrect ".$HLANameCorrect."\n";
    #exit;
    }else{
	if (exists $hashPepSeq{$PepName}){
		#$hashPepSeq{$PepName}=$hashPepSeq{$PepName}.$line;
	}
	else{
		$hashPepSeq{$PepName}=$line;
	}
	}
}
close(FILE);


open (FILEOUT, ">MHC_ligand_selected.txt")|| die "ERROR: CEBA couldn't open the file MHC_ligand_selected.txt!";

@arrayPep=keys(%hashPepSeq);
foreach $pep(@arrayPep)
{
#print $pep."\n";
print FILEOUT ">".$pep."\n".$hashPepSeq{$pep}."\n";
}
close(FILEOUT);


exit;
