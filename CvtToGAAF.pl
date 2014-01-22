#!c:/perl64/bin/perl.exe
##
## Version 2003.12.16
## fasta_pro.pl
## Copyright (C) 2013 John Wilkins
##
## Use of this software governed by the Artistic license,
## as reproduced at http://www.opensource.org/licenses/artistic-license.php
##
#
# Convert a text file of modifications in  ????????? format to "gaaf" format.

my $FId = $ARGV[0]; #File name without the extension.
my $MOd = $ARGV[1]; #MOdification name (e.g. Phospho, Acetyl or GlyGly )
my $pathIn = "$FId.txt";
my $pathO  = "$FId.gaaf"; 
open(IFILE,"<$pathIn") or die "Input file path \"$pathIn\" does not exist\n";
open(OFILE,">$pathO"); 

$/="";
$_ = <IFILE>;
#print $_;
my $ID=0;
my @startrow = (); # create table to place data for sorting. this is junk before the sort key
my @posA = ();     # major key is the position in dna string of the beginning of the nucleotide triple
my @posB = ();     #  minor key is the position in dna string of the end of the nucleotide triple
my @endrow = ();   # more junk after the sort key.
my @mid1row = ();  # more junk after the sort key.
my @mid2row = ();  # more junk after the sort key.
my @endrow = ();   # more junk after the sort key.
while ($_) {
	my $protein=""; my $chr=""; my $strand=""; 
	if (m/protein\:(.+?)\s/isg) {$protein=$1;}  #Pluck off protein (peptide) id,
	if (m/chromosome\:(\w+)\s/isg) {$chr=$1;}   #... chromosome (should be same as input arg0),
	if (m/strand\:(\-\d+)\s/isg) {$strand="-";} #  ...  and sign of strand.
	else {$strand="+";}

	# Each row beginning "pos" corresponds to one or more output rows;
	while (/pos\:(\d+).*?res\:([A-Z]).*?gc\:(\d+)\,(\d+)\,(\d+)/gsi) {  # for each pluck off the pos of residue in peptide and residue letter.
                                                                            # as well as positions of nucleotide triplet in dna string.

		# if nucleotide triplet is consecutive, create one row. 
		if (($3+1==$4) && ($4+1==$5)) {
			$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
			$posA[$ID] = $3;
			$posB[$ID]   = $5;
			$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
			$mid2row[$ID] = "\;Reference_aa=$2\;";
			$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
			$ID++;
			#print TFILE "chr$chr\tgpmdb\tPTM\t$3\t$5\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
		}
		# else create a row each for nonconsecutive nucleotide positions.
		else {
			if ($3+1!=$4) {
				$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
				$posA[$ID] = $3;
				$posB[$ID]   = $3;
				$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
				$mid2row[$ID] = "\;Reference_aa=$2\;";
				$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
				$ID++;
				#print TFILE "chr$chr\tgpmdb\tPTM\t$3\t$3\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
				if ($4+1!=$5) {
					$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
					$posA[$ID] = $4;
					$posB[$ID]   = $4;
					$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
					$mid2row[$ID] = "\;Reference_aa=$2\;";
					$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
					$ID++;
					#print TFILE "chr$chr\tgpmdb\tPTM\t$4\t$4\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
					$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
					$posA[$ID] = $5;
					$posB[$ID]   = $5;
					$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
					$mid2row[$ID] = "\;Reference_aa=$2\;";
					$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
					$ID++;
					#print TFILE "chr$chr\tgpmdb\tPTM\t$5\t$5\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
				}
				else {
					$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
					$posA[$ID] = $4;
					$posB[$ID]   = $5;
					$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
					$mid2row[$ID] = "\;Reference_aa=$2\;";
					$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
					$ID++;
					#print TFILE "chr$chr\tgpmdb\tPTM\t$4\t$5\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
				}
			}
			else {
				$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
				$posA[$ID] = $3;
				$posB[$ID]   = $4;
				$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
				$mid2row[$ID] = "\;Reference_aa=$2\;";
				$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
				$ID++;
				#print TFILE "chr$chr\tgpmdb\tPTM\t$3\t$4\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
				$startrow[$ID] = "chr$chr\tgpmdb\tPTM\t";
				$posA[$ID] = $5;
				$posB[$ID]   = $5;
				$mid1row[$ID] = "\t.\t$strand\t.\tID="; 
				$mid2row[$ID] = "\;Reference_aa=$2\;";
				$endrow[$ID] = "Reference_ptm=$protein.pm:$2$1+$MOd\;";
				$ID++;
				#print TFILE "chr$chr\tgpmdb\tPTM\t$5\t$5\t.\t$strand\t.\tID=0000\;Reference_aa=$2\;Reference_ptm=$protein.pm:$2$1+$MOd\;\n";
			}
		}
	}
	$_ = <IFILE>;
}
close IFILE;
#sort the indexes by posB within posA
my @ix; #sorted index;
my $length = scalar(@posA);
@ix = sort { $posA[$a] <=> $posA[$b] || $posB[$a] <=> $posB[$b]}  0 .. $#posA;

#print gff file header.

print OFILE <<End_of_header;
##gaaf-version 1.00
##genome-build GRCh37.p13
##proteome-build ENSEMBL 70

End_of_header

#process in sorted order merging rows to write file records with unique nucleotide tripet positions/ranges.
my $ID=0;
my $i=0;
while ($i <= $#ix) {
	$ID++;
	my $prtline = $startrow[$ix[$i]].$posA[$ix[$i]]."\t".$posB[$ix[$i]].$mid1row[$ix[$i]]."$ID".$mid2row[$ix[$i]].$endrow[$ix[$i]];
	my $prtA = $posA[$ix[$i]];
	my $prtB = $posB[$ix[$i]];
	$i++;
	while (($i <= $#posA) &&  ($posA[$ix[$i]] == $prtA ) && ($posB[$ix[$i]] == $prtB )) {
		$prtline .= $endrow[$ix[$i]];
		$i++
	}
	print OFILE "$prtline\n";
}
close OFILE;
