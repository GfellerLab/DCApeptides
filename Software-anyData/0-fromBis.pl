#!/usr/bin/perl -w
 use POSIX qw(ceil floor);
$filenameIn=shift;
$filenameOut=shift;
$format=shift;# (FASTA or PHYLIP)

$nbreSeq=0;
$name="";
$seq="";
open (FILE, $filenameIn)|| die "ERROR: CEBA couldn't open the file $filenameIn!";
while ($line =<FILE>)
{ 
 if($line=~ /^(\s)*$/){
   }else{
    chomp $line;
    ($name,$seq)= split(" ",$line);
    $hash{$name}{$seq}=$seq;
    $nbreSeq++;
    }
}
close(FILE);
$seqLength= length $hash{$name};
print "Number of sequences ".$nbreSeq."\n";
	
$totalSeq=0;

if ($format eq "FASTA"){
	open (FILOUT, ">".$filenameOut)|| die "ERROR: CEBA couldn't write into the file $filenameOut!";
	@keysSeq=sort keys  %hash;
	
	foreach $seqName (@keysSeq){
	@keysSeq2=sort keys  %{$hash{$seqName}};
	foreach $peptide (@keysSeq2){
	$totalSeq=$totalSeq+(scalar @keysSeq2);
		$sequence=$hash{$seqName}{$peptide};
		$sequence =~ tr/a-z/A-Z/;
		$sequence =~ tr/\./\-/;
		print FILOUT ">".$seqName."\n".$sequence."\n"
	}
	}
	close(FILOUT);
}
print "TotalSeq ".$totalSeq."\n";

if ($format eq "PHYLIP"){
	open (FILOUT, ">".$filenameOut)|| die "ERROR: CEBA couldn't write into the file $filenameOut!";
	print FILOUT $nbreSeq." ".$seqLength."\n";
	@keysSeq=sort keys  %hash;
	foreach $seqName (@keysSeq){
		$sequence=$hash{$seqName};
		$sequence =~ tr/a-z/A-Z/;
		$sequence =~ tr/\./\-/;
		print FILOUT $seqName." ".$sequence."\n"
	}
	close(FILOUT);
}


exit;
