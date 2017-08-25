#!/usr/bin/perl -w
use POSIX qw(ceil floor);
use File::Find;

$HLAAlignment=shift;
$PeptideAlignment=shift;
open (FILE, $HLAAlignment)|| die "ERROR: CEBA couldn't open the file $HLAAlignment!";
while ($line =<FILE>)
{ 
    chomp $line;
    if ($line =~ /^>/)
	{
		@lineContent=split(">",$line);
		$HLAName=$lineContent[1];

    }else{
	if (exists $hashHLASeq{$HLAName}){
		 print "ERROR: HLA seq is on one line!! what is going on here? $HLAName\n";
	}else{
		$hashHLASeq{$HLAName}=$line;
	}}
}
close(FILE);

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
		print "ERROR: PEP seq is on one line!! what is going on here? $PepName\n";
	}else{
		$hashPepSeq{$PepName}=$line;
	}}
}
close(FILE);
print "CHECKS......\n";
@arrayPep=keys(%hashPepSeq);
print "scalar PEP=".(scalar @arrayPep)."\n";
foreach $pep(@arrayPep)
{
if(exists $hashHLASeq{$pep}){
}else{
print "HLA sequence not found:".$pep."\n";
}
}



@arrayHLA=keys(%hashHLASeq);
print "scalar HLA=".(scalar @arrayHLA)."\n";

foreach $HLAName(@arrayHLA)
{
if(exists $hashPepSeq{$HLAName}){
}else{
print "HLA sequence has no pepetides :".$HLAName."\n";
}

}

exit;
