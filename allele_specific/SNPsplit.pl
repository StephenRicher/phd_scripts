#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Cwd;


## This program is Copyright (C) 2014-17, Felix Krueger (felix.krueger@babraham.ac.uk)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.



## Reading in a BAM or SAM file
my $pipeline_version = '0.3.4';
my $parent_dir = getcwd;
my ($snp_file,$no_sort,$verbose,$samtools_path,$bam,$paired,$hic,$conflict,$output_dir,$singletons,$bisulfite) = process_commandline ();
my %snps; # storing SNP position and sequence information for all SNPs

my $snp_found = 0;
my $no_snp_found = 0;
my $count_SNPs_stored = 0;

my ($no_snp,$genome1,$genome2,$bizarre,$conflicting,$unassigned,$ct_snp,$unassigned_but_ct);

while (@ARGV){
    my $file = shift @ARGV;

    warn "Now starting to process file <<< '$file' >>> \n\n";
    sleep (1);

    my $outfile = open_output_filehandles($file);

    warn "\nSummary of parameters for SNPsplit-tag:\n";
    warn '='x40,"\n";
    warn "SNPsplit infile:\t\t$file\n";
    warn "SNP annotation file:\t\t$snp_file\n";
    warn "Output directory:\t\t>$output_dir<\n";
    warn "Parent directory:\t\t>$parent_dir<\n";
    warn "Samtools path:\t\t\t$samtools_path\n";
    if ($bam){
	warn "Output format:\t\t\tBAM (default)\n";
    }
    else{
	warn "Output format:\t\t\tSAM\n";
    }
    if ($hic){
	warn "Input format:\t\t\tHi-C (by definition paired-end)\n";
    }
    else{
	if ($paired){
	    if ($bisulfite){
		warn "Input format:\t\t\tBismark paired-end (not relevant for tagging process)\n";
	    }
	    else{
		warn "Input format:\t\t\tPaired-end (not relevant for tagging process)\n";
	    }
	}
	else{
	    if ($bisulfite){
		warn "Input format:\t\t\tBismark Single-End\n";
	    }
	    else{
		warn "Input format:\t\t\tSingle-End\n";
	    }
	}
    }
    warn "\n\n";

    sleep (1);


    ### SAM TO BAM CONVERSION - if necessary
    if ($file =~/\.sam$/){
	$file = sam2bam_convert($file);
    }
    elsif ($file =~ /\.bam$/){
	# fine
    }
    else{
	die "Please supply a file in either SAM or BAM format (ending with .sam or .bam). File name given was '$file'\n\n";
    }

    ### Paired-end
    if ($paired){
	unless ($no_sort){
	    warn "File specified as unsorted paired-end BAM file\n";
	    $file = sort_by_name_paired_end($file); # setting this as the new file
	}
	else{
	    warn "File was specified to be sorted (i.e. Read 1 and Read 2 are following each other in the paired-end BAM file), thus skipping the sorting step\n\n";
	}
    }
    else{ ### single-end
	warn "File specified as single-end BAM file\n";
    }

    if (%snps){ # not necessary to read SNPs a second time
	warn "Using already stored SNP information\n";
	sleep (1);
    }
    else{
	### READ IN SNP FILE
	$count_SNPs_stored = read_snps();

	unless (%snps){
	    die "SNP file doesn't appear to contain readable SNP positions. Please check the filename/file format and try again\n\n";
	}
    }

    process_masked_BAM_file($file);

    warn "Finished allele-tagging for file <<< '$file' >>> \n\n";

    ######################
    #####            #####
    #####  tag2sort  #####
    #####            #####
    ######################

    # Calling the allele-specific sorting now
    my @args;
    if ($hic){
	push @args, '--hic';
    }
    if ($samtools_path){
	push @args, "--samtools_path $samtools_path";
    }
    if ($paired){
	push @args, '--paired';
    }
    unless ($bam){
	push @args, '--sam';
    }
    if ($verbose){
	push @args, '--verbose';
    }
    if ($output_dir){ # not passing on empty strings
	push @args, "--output_dir $output_dir";
    }
    if ($conflict){
	push @args, '--conflicting';
    }
    if ($singletons){
	push @args, '--singletons';
    }

    system ("$RealBin/tag2sort @args $outfile");

    ####

    my $sortfile = $outfile;

    if ($sortfile =~ s/allele_flagged\.bam$/bam/){
	if ($sortfile =~ s/bam$/SNPsplit_sort.txt/){
	    if (-e $sortfile){

		open (SORT,"${output_dir}$sortfile") or die "Failed to read sorting report from '${output_dir}$sortfile': $!\n";
		my $doc;
		while(<SORT>){
		    chomp;
		    $_ =~ s/\r//g;
		    $doc .= $_."\n";
		}

		close SORT or die $!;

		### appending the sorting report to the tagging report
		warn "Appending sorting report to tagging report...\n";

		print REPORT $doc;
	    }
	}
    }
    warn "SNPsplit processing finished...\n\n";
}

###

sub process_masked_BAM_file{

    my ($file) = @_;

    open (IN,"$samtools_path view -h $file |") or die "Failed to read from BAM file $file: $!\n";
    warn "Reading from sorted mapping file '$file'\n";

    ### READING FROM MAPPING BAM FILE

    ($no_snp,$genome1,$genome2,$bizarre,$unassigned,$ct_snp,$unassigned_but_ct) = (0,0,0,0,0,0,0);
    my $count = 0;

    my $count_N_containing = 0;
    my $non_N_containing = 0;
    my $contained_N_deletion = 0;
    my $multi_N_deletion = 0;
    my $unmapped = 0;

    while (<IN>){

	if ($_ =~ /^\@/){ # header lines
	    print OUT;
	    next;
	}

	++$count;
	# warn "N-masked line $count\n"; sleep(1);

	if ($count%1000000 == 0){
	    warn "Processed $count lines so far\n";
	}

	my ($id,$flag,$chr,$start,$cigar,$sequence)  = (split /\t/)[0,1,2,3,5,9];

	my ($genome1_count,$genome2_count,$ct_unassigned) = (0 , 0, 0); # will get filled below

	### Checking for unmapped read
	if ($flag & 0x4){
	    # warn "read unmapped:\n$_\n"; sleep(1);
	    $unmapped++;
	    next; # skipping
	}

	my $md_tag = $1 if ($_ =~ /MD:Z:(.+?)\s+/);
	unless ($md_tag){
	    die "Failed to extract MD:Z tag from line:\n$_\n";
	}

	#  if ($verbose){
	#  print $_;
	#  print join ("\t",$id,$flag,$chr,$start,$cigar,$sequence,$md_tag),"\n";
	# }

	my ($XX_tag); # allele-specificity flag

	if ($md_tag =~ /N/){  # skipping these reads for now, may or may not contain additional deletions before the N

	    my $n_dels = 0;
	    while ($md_tag =~ /\^\D*N\D*\d+/g){ # looking for the number of deletions of the N-SNP position in the read
		++$n_dels;
	    }

	    if ($n_dels == 1){
		++$contained_N_deletion;
		++$bizarre;
		$XX_tag = 'XX:Z:CF';     # Reads where the N was deleted are 'conflicting'
		# print $_;
	    }
	    elsif($n_dels > 1){
		++$contained_N_deletion;
		++$multi_N_deletion;
		++$bizarre;
		# warn "Read contained more than 1 deleted N position. MD-tag: $md_tag\n";
		# print join ("\t",$id,$flag,$chr,$start,$cigar,$sequence,$md_tag),"\n";
		$XX_tag = 'XX:Z:CF';     # calling this situation 'conflicting' for the time being
		# print $_;
	    }
	    else{ # no deletion of the N position itself
		++$count_N_containing;
		### Determine whether the N was a known SNP and whether we can assign an allele
		($genome1_count,$genome2_count,$ct_unassigned) = determine_N_position ($chr,$start,$cigar,$sequence,$md_tag,$_);

		if ($genome1_count > 0 and $genome2_count > 0){
		    ++$bizarre;
		    $XX_tag = 'XX:Z:CF'; #  read contained one or more SNPs for both genomes -> dodgy
		}
		elsif ($genome1_count > $genome2_count){
		    $genome1++;
		    $XX_tag = 'XX:Z:G1'; # genome 1 specific
		}
		elsif($genome2_count > $genome1_count){
		    $genome2++;
		    $XX_tag = 'XX:Z:G2'; # genome 2 specific
		}
		elsif($genome1_count == 0 and $genome2_count == 0){
		    # read came back without any assigned number, e.g. the SNP position contained a base different to either genome. The read is thus unassignable
		    if ($bisulfite){
			# for bisulfite reads this may have happened because of C>T SNPs
			if ($ct_unassigned > 0){
			    # we will not score this read as having an unexpected base at the SNP position, but as unassignable because of C>T SNps
			    ++$unassigned_but_ct;
			}
			else{
			    ++$no_snp; # here the base must have been different to ref and SNP, and not involve a C/T of either of them
			}
		    }
		    else{
			++$no_snp;
		    }
		    $XX_tag = 'XX:Z:UA'; # unassignable read
		    ++$unassigned;
		}
		else{
		    warn "g1: $genome1_count\ng2: $genome2_count\n";
		    die "The read could not be assigned an XX-tag\n\n";
		}
	    }
	}
	else{
	    ++$non_N_containing;
	    ### Read did not cover a SNP position and is by definition unassignable
	    $XX_tag = 'XX:Z:UA'; # unassignable read
	    ++$unassigned;
	}

	unless ($XX_tag){
	    die "The XX tag should have been set by now ...\n";
	}

	chomp; # appending the XX tag to the current BAM line
	$_ .= "\t$XX_tag\n";
	print OUT;

    }

    ### Printing Summary reports of the stats on a per-read basis

    my ($perc_no_snp,$perc_genome1,$perc_genome2,$perc_bizarre,$perc_unassigned,$perc_unmapped);
    if ($count == 0){
	$perc_unassigned = $perc_no_snp = $perc_genome1 = $perc_genome2 = $perc_bizarre = $perc_unmapped = 'N/A';
    }
    else{
	$perc_unassigned = sprintf ("%.2f",$unassigned*100/$count);
	$perc_no_snp     = sprintf ("%.2f",$no_snp*100/$count);
	$perc_genome1    = sprintf ("%.2f",$genome1*100/$count);
	$perc_genome2    = sprintf ("%.2f",$genome2*100/$count);
	$perc_bizarre    = sprintf ("%.2f",$bizarre*100/$count);
	$perc_unmapped   = sprintf ("%.2f",$unmapped*100/$count);
    }

    warn "\nAllele-tagging report\n",'='x21,"\n";
    warn "Processed $count read alignments in total\n";
    warn "Reads were unaligned and hence skipped: $unmapped ($perc_unmapped%)\n";
    warn "$unassigned reads were unassignable ($perc_unassigned%)\n";
    warn "$genome1 reads were specific for genome 1 ($perc_genome1%)\n";
    warn "$genome2 reads were specific for genome 2 ($perc_genome2%)\n";
    if ($bisulfite){
	warn "$unassigned_but_ct reads that were unassignable contained C>T SNPs preventing the assignment\n"
    }
    warn "$no_snp reads did not contain one of the expected bases at known SNP positions ($perc_no_snp%)\n";
    warn "$bizarre contained conflicting allele-specific SNPs ($perc_bizarre%)\n\n\n";

    print REPORT "\nAllele-tagging report\n",'='x21,"\n";
    print REPORT "Processed $count read alignments in total\n";
    print REPORT "Reads were unaligned and hence skipped: $unmapped ($perc_unmapped%)\n";
    print REPORT "$unassigned reads were unassignable ($perc_unassigned%)\n";
    print REPORT "$genome1 reads were specific for genome 1 ($perc_genome1%)\n";
    print REPORT "$genome2 reads were specific for genome 2 ($perc_genome2%)\n";
    if ($bisulfite){
	print REPORT "$unassigned_but_ct reads that were unassignable contained C>T SNPs preventing the assignment\n"
    }
    print REPORT "$no_snp reads did not contain one of the expected bases at known SNP positions ($perc_no_snp%)\n";
    print REPORT "$bizarre contained conflicting allele-specific SNPs ($perc_bizarre%)\n\n\n";

    my ($pc_n_found , $pc_no_n_found , $pc_n_deletion , $pc_multi_n_deletion);

    if ($count == 0){
	$pc_n_found = $pc_no_n_found =  $pc_n_deletion = $pc_multi_n_deletion = 'N/A';
    }
    else{
	if ( ($snp_found + $no_snp_found) == 0){
	    $pc_n_found = 'N/A';
	    $pc_no_n_found = 'N/A';
	}
	else{
	    $pc_n_found    = sprintf ("%.2f",$snp_found*100/($snp_found + $no_snp_found));
	    $pc_no_n_found = sprintf ("%.2f",$no_snp_found*100/ ($snp_found + $no_snp_found) );
	}

	$pc_n_deletion   = sprintf ("%.2f",$contained_N_deletion*100/ $count );
	$pc_multi_n_deletion = sprintf ("%.2f",$multi_N_deletion*100/ $count );
    }

    warn "SNP coverage report\n===================\n";
    warn "SNP annotation file:\t$snp_file\n";
    warn "SNPs stored in total:\t$count_SNPs_stored\n";
    warn "N-containing reads:\t$count_N_containing\nnon-N:\t\t\t$non_N_containing\ntotal:\t\t\t$count\n";
    warn "Reads had a deletion of the N-masked position (and were thus called Conflicting):\t$contained_N_deletion ($pc_n_deletion%)\n";
    warn "Of which had multiple deletions of N-masked positions within the same read:\t$multi_N_deletion ($pc_multi_n_deletion%)\n\n";

    warn "Of valid N containing reads,\n";
    warn "N was present in the list of known SNPs:\t$snp_found ($pc_n_found%)\n";
    if ($bisulfite){
	warn "Positions were skipped since they involved C>T SNPs:\t$ct_snp\n";
    }
    warn "N was not present in the list of SNPs:\t\t$no_snp_found ($pc_no_n_found%)\n\n";


    print REPORT "SNP coverage report\n===================\n";
    print REPORT "SNP annotation file:\t$snp_file\n";
    print REPORT "SNPs stored in total:\t$count_SNPs_stored\n";
    print REPORT "N-containing reads:\t$count_N_containing\nnon-N:\t\t\t$non_N_containing\ntotal:\t\t\t$count\n";
    print REPORT "Reads had a deletion of the N-masked position (and were thus dropped):\t$contained_N_deletion ($pc_n_deletion%)\n";
    print REPORT "Of which had multiple deletions of N-masked positions within the same read:\t$multi_N_deletion\n\n";

    print REPORT "Of valid N containing reads,\n";
    print REPORT "N was present in the list of known SNPs:\t$snp_found ($pc_n_found%)\n";
    if ($bisulfite){
	print REPORT "Positions were skipped since they involved C>T SNPs:\t$ct_snp\n";
    }
    print REPORT "N was not present in the list of SNPs:\t\t$no_snp_found ($pc_no_n_found%)\n\n";

    close OUT or warn "Failed to close filehandle OUT: $!. Exit status: $?\n";

    return;

}


sub determine_N_position{

    my ($chr,$start,$cigar,$sequence,$md_tag,$line) = @_;
    my ($genome1,$genome2,$irrelevant,$ct_unassigned) = (0,0,0,0); # determining number of SNPs that can be used to assign a certain allele

    if ($verbose){ warn "Now looking at the following read:\n$chr\t$start\t$cigar\t$sequence\t$md_tag\n";sleep(1);}
    # parsing CIGAR string
    my @len = split (/\D+/,$cigar); # storing the length per operation
    my @ops = split (/\d+/,$cigar); # storing the operation
    shift @ops; # remove the empty first element
    die "CIGAR string contained a non-matching number of lengths and operations\n" unless (scalar @len == scalar @ops);

    my @ns = split /N/,$md_tag;

    if (scalar @ops == 1 and $ops[0] eq 'M'){

	if ($verbose){ warn "Read is a single contiguous match: $cigar\n"; sleep(1);}
	# adjusting the read mapping position the SNP position

	if ($verbose){
	    warn "MD-tag: $md_tag\n";
	    print join ("\t",@ns),"\n";
	}

	die "There were fewer than 2 elements in the MD tag $md_tag! number of elements: ",scalar @ns,"\n" if (scalar @ns < 2);
	my $pos;

	foreach my $index (0..($#ns-1)){ # the last index of @ns is the number of matching strings after the last N, so is not relevant

	    if ($ns[$index] =~/\D/){ # Element is a Non-Digit

		if ($ns[$index] =~/\W/){
		    if ($verbose){ warn "Element contains a Non-Word character: $ns[$index]\n~~~~~~~~!!!\n"};
		}

		if ($verbose){ warn "element contains a Non-DIGIT: $ns[$index]\n"};
		my @array = split /[ATCG]/,$ns[$index];
		if ($verbose){ warn "Elements in the array: ",scalar @array,"\n"};

		foreach my $element(@array){
		    if ($verbose){ warn "Adding $element to pos\n"};
		    $pos += $element;
		    if ($verbose){ warn "pos right now: $pos\n"};
		}
		if ($verbose){ warn "now adjusting for the number of elements in the array, adding ",(scalar @array-1),"\n"};
		$pos += (scalar @array) - 1;
		if ($verbose){ warn "pos right now: $pos\n\n"; sleep(1)};
	    }
	    else{ # element is a digit, so we can just add it to the position
		if ($verbose){ warn "element is $ns[$index]\n"};
		$pos += $ns[$index];
		if ($verbose){ warn "Start: $start\tpos: $pos\n"; sleep(1);};
	    }

	    my $genomic_pos = $pos + $start;

	    ### seeing if the SNP is present in the SNP position list
	    if (exists $snps{$chr}->{$genomic_pos}){
		++$snp_found;

		if ($verbose) { warn "SNP was present: ref: $snps{$chr}->{$genomic_pos}->{ref}\tsnp: $snps{$chr}->{$genomic_pos}->{snp}\n";}
		if ($verbose){ warn "pos: $pos\n$sequence\n";}
		my $read_base = substr($sequence,$pos,1); # The N position has not been adjusted for the N yet so can be used as index
		if ($verbose) { warn "base in read was: $read_base\n";}

		if ($bisulfite){
		    ($genome1,$genome2,$irrelevant,$ct_unassigned) = score_bisulfite_SNPs($line,$chr,$genomic_pos,$read_base,$genome1,$genome2,$irrelevant,$ct_unassigned);
		}
		else{ ## DEFAULT

		    if ($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){
			if ($verbose) { warn "genome 1 specific!\n";}
			++$genome1;
		    }
		    elsif ($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){
			if ($verbose) { warn "genome 2 specific!\n";}
			++$genome2;
		    }
		    else{
			if ($verbose) { warn "The base was different than both genome 1 or genome 2: $read_base\n";}
			++$irrelevant;
		    }
		}
	    }
	    else{
		++$no_snp_found;
		if ($verbose) { warn "SNP was not present in the SNP list\n";}
	    }

	    if ($verbose) { warn "now adjusting the position for the N itself\n";}
	    $pos += 1; # adjusting position for the N itself. irrelevant if this is the last element of the array
	    if ($verbose) { warn "pos right now: $pos\n";}
	}

	if ($verbose) { warn "~~~~~~~~~~~~~~~~~~~~\n\n\n\n";}

    }
    else{ # CIGAR string is not a simple match
	if ($verbose){
	    warn "CIGAR string is not a simple match:\n";
	    print "MD-tag: $md_tag\tIndividual elements are:\n";
	    print join ("\t",@ns),"\n";
	}

	die "There were fewer than 2 elements in the MD tag $md_tag! number of elements: ",scalar @ns,"\n" if (scalar @ns < 2);
	my $pos = 0;
	my %insertions_accounted_for;

	### STEP I: determining the position of the N(s) in the read

	foreach my $index (0..($#ns-1)){ # the last index of @ns is the number of matching strings after the last N, so is not relevant

	    if ($ns[$index] =~/\D/){ # non-digit

		if ($ns[$index] =~/\W/){ # non-word character, e.g. ^ for deletions
		    if ($verbose) { warn "element contains a Non-Word character: $ns[$index]\n";}
		    if ($ns[$index] =~ /\^/){ # contains deletion
			if ($verbose){  warn "This was a deletion\n";}

			# split into single characters
			my @string = split //,$ns[$index];
			if ($verbose){  warn "before starting, \$pos was $pos\n";}
			my $operation; # will store the operation, either a number for matches or a variable number of letters for deletions
			foreach my $char (@string){

			    # setting first element
			    unless (defined $operation){ # can be 0
				if ($verbose){	warn "Setting \$operation to $char\n";}
				$operation = $char;
				next;
			    }
			    if ($verbose){   warn "$char\n";}

			    if ($operation =~ /\d+/){
				if ($verbose){ warn "$operation is a digit\n";}
				if ($char =~ /\d+/){
				    if ($verbose){ warn "appending $char to $operation\n";}
				    $operation .= $char;
				    if ($verbose){ warn "\$operation is now $operation\n";}
				}
				else{
				    if ($verbose){  warn "\$char is now a non-digit character, so need to process \$operation\n";}
				    $pos += $operation;
				    if ($verbose){  warn "\$pos is now $pos\n";}

				    # setting $operation as the new character
				    $operation = $char;
				    if ($verbose){  warn "Setting \$operation to: $operation\n";}
				}
			    }
			    else{ # $operation is a non-digit character
				if ($verbose){ 	warn "$operation is a non-digit character\n";}

				if ($char =~ /\d+/){
				    if ($verbose){ 	warn "\$char is now a digit, so need to process \$operation: $operation\n";}

				    if ($operation =~ /^\^/){ #  if the string starts with a ^ symbol is was a deletion
					if ($verbose){ 	 warn "Deletions do not count towards the position within a read -> skipping...\n";}
				    }
				    else{	                    # else it was a mismatch and the full length needs to be considered
					if ($verbose){ warn "Now adding length of \$operation (just a mismatch)\n";}
					$pos += (length $operation);
				    }
				    if ($verbose){ warn "\$pos is now $pos\n";}

				    # setting $operation as the new character
				    $operation = $char;
				    if ($verbose){ warn "Setting \$operation to: $operation\n";}
				}
				else{
				    if ($verbose){  warn "appending $char to $operation\n";}
				    $operation .= $char;
				    if ($verbose){ warn "\$operation is now $operation\n";}
				}
			    }
			}
			if ($verbose){ 	    warn "Now processing the last string which has to be a number\n";}
			die "\$operation was not a number: '${operation}'!" unless ($operation =~ /\d+/);
			$pos += $operation;
			if ($verbose){   warn "\$pos is now $pos\n\n";}
		    }
		    else{
			# this might give us a clue of which other characters might turn up in the MD-tag
			if ($verbose){ warn "The character was no deletion: $ns[$index]\n"; sleep(1);}
		    }
		}
		else{
		    if ($verbose){ warn "Element contains a Non-DIGIT word character(s), but is no deletion: $ns[$index]\n";}
		    my @array = split /[ATCG]/,$ns[$index];

		    if ($verbose){
			warn "Elements in the array: ",scalar @array,"\n";
			print join ("\t",@array),"\n";
		    }

		    foreach my $element(@array){
			if ($verbose){  warn "Adding $element to pos\n";}
			$pos += $element;
			if ($verbose){ warn "pos right now: $pos\n";}
		    }
		    if ($verbose){ warn "now adjusting for the number of elements in the array, adding ",(scalar @array-1),"\n";}
		    $pos += (scalar @array) - 1;
		    if ($verbose){ warn "pos right now: $pos\n";}
		}
	    }
	    else{ # element is a digit, so we can just add it to the position
		if ($verbose){ warn "element is $ns[$index]\n"; }
		$pos += $ns[$index];
		if ($verbose){ warn "Start: $start\tpos: $pos\n"; }
	    }


	    ### STEP II: Adjusting $pos for insertions
	    # Insertions in the read are omitted from the MD:Z: tag, and can only be seen when processing the CIGAR string. In order not to extract the wrong base from the read we quickly process the CIGAR string and adjust $pos if necessary

	    # Even though insertions should affect the position within a read, they are in fact absent from the MD:Z: field. Here is an example:
	    # HISEQ2500-01:182:H89GMADXX:1:1101:15083:2899    83      17      31124195        42      9M2I40M =       31124087        -157    GAGGCTCATGAGGCCTAGCTGACAGTGTCCCCTGCCCCTTGGTCAGTTTAC     JJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJHHGHHFFFFFCCC AS:i:-12 XN:i:1  XM:i:1  XO:i:1  XG:i:2  NM:i:3  MD:Z:44N4       YS:i:0  YT:Z:CP
	    # The sequence is 51bp long, so is the CIGAR string, but the MD:Z: field is only 49bp long!!

	    if ($verbose){  warn "\nNow adjusting \$pos for potential insertions in the read. \$pos now is $pos\n";}
	    my $total_length_of_ops = 0;
	    if ($verbose){  warn "\$total_length_of_ops is $total_length_of_ops at the start\n";}

	    foreach my $index(0..$#len){
		if ($ops[$index] eq 'M'){  # standard matching bases
		    $total_length_of_ops += $len[$index];
		    if ($verbose){ warn "\$total_length_of_ops is $total_length_of_ops\n";}
		}
		elsif ($ops[$index] eq 'D'){  # deletions do not count towards the position within the read
		    if ($verbose){ warn "Deletion. skipping...\n";}
		}
		elsif ($ops[$index] eq 'N'){  # splice junctions are similar to deletions and do not count towards the position within the read
		    if ($verbose){ warn "Splice junction. skipping...\n";}
		}
		elsif($ops[$index] eq 'I'){
		    if ($total_length_of_ops >= $pos){
			if ($verbose){ warn "\$pos is $pos, but there were already enough matches to be at position $total_length_of_ops! Exiting adjustment\n\n";}
			last; # exiting now
		    }
		    else{
			if (exists $insertions_accounted_for{$total_length_of_ops} ){ # if we have already adjusted the position for this insertion we skip it this time
			    if ($verbose){ warn "This insertion has been accounted for already. Skipping this time...\n";}
			    next;
			}
			else{
			    ++$insertions_accounted_for{$total_length_of_ops};
			    if ($verbose){ warn "Operation is 'I', adding length of $len[$index]\n";}
			    $pos += $len[$index];
			    if ($verbose){ warn "New \$pos is now: $pos\n";}
			}
		    }
		}
		else{
		    die "Found CIGAR operations other than M, I, D or N: '$ops[$index]'. Not allowed at the moment\n";
		}
	    }

	    ### STEP III; Determining the genomic position
	    my $genomic_pos = $start; # this is the start of the read
	    my $temp_pos = $pos;
	    my $endpos = $start; # of the current CIGAR operation
	    my $operation_length = 0;
	    if ($verbose){
		warn "Chromosome is $chr\n";
		warn "Now determining the genomic position of the SNP\n";
		warn "Position within the read is: $temp_pos. MD-tag is: $md_tag\n"; sleep(1);
		warn "To start with, genomic position is: $genomic_pos\n";
		warn "To start with, operation length is $operation_length\n\n";
	    }

	    ### determining end position of a read
	    foreach my $index(0..$#len){

		$operation_length = $len[$index];

		if ($ops[$index] eq 'M'){  # standard matching bases
		    $endpos += $len[$index];
		    if ($verbose){
			warn "Operation is 'M', length $len[$index] bp\n";
			warn "Endpos is now: $endpos\n";
			warn "Operation length is now: $operation_length\n";
		    }

		    # testing to see if the position within the read is contained within this stretch of matching bases
		    if ($temp_pos <= $operation_length){
			if ($verbose){
			    warn "Read overlaps the N position\n";
			    warn "\$temp_pos is: $temp_pos\n";
			}
			$genomic_pos += $temp_pos;
			if ($verbose){ warn "Setting \$genomic_pos to: $genomic_pos\n";}
			last; # once we have got the genomic positions we can continue
		    }
		    else{
			# adding the length of the M operation to genomic position
			if ($verbose){   warn "Read does not yet overlap the N position\n";}
			# adjusting the genomic position
			$genomic_pos += $operation_length;
			$temp_pos    -= $operation_length; # subtracting the operation length from the length within the read
			if ($verbose){
			    warn "Adding $operation_length to \$genomic_pos. Now at: $genomic_pos\n";
			    warn "Subtracting $operation_length from \$temp_pos. Now at: $temp_pos\n";
			}
		    }
		}
		elsif($ops[$index] eq 'I'){ # insertions in the read do not affect the genomic position
		    if ($verbose){
			warn "Operation is 'I', not adjusting genomic position but \$temp_pos\n";
			warn "Operation length is now: $operation_length\n\n";
		    }
		    $temp_pos    -= $operation_length; # subtracting the operation length from the length within the read
		    if ($verbose){  warn "\$temp_pos is now: $temp_pos\n";}
		    # Since we have adjusted $pos in STEP II we can now subtract the length of the insertion from the read (which I find logical)
		}
		elsif($ops[$index] eq 'D'){ # deletions do affect the genomic position but not the $temp_pos within the read
		    if ($verbose){   warn "Operation is 'D',adding $len[$index] bp to the genomic position\n";}
		    ### eaxample read with a deletion, length 51 bp
		    # HISEQ2500-01:182:H89GMADXX:1:1101:15249:2789    99      8       87082089        23      33M6D18M        =       87082185        147     ACAGCCAACTGCTTTCCTCCAGGGGTGCCCATGACCTGCACACCTCAGCTG     CCCFFFFFHHHHHJJJJJJJJJJJJHHIJJJJJJJJJJJJJJJJJJJJJJJ  AS:i:-24        XN:i:1  XM:i:1  XO:i:1  XG:i:6  NM:i:7  MD:Z:33^ACTCGA11N6      YS:i:-1 YT:Z:CP
		    $genomic_pos += $operation_length;
		    if ($verbose){
			warn "Operation length is now: $operation_length\n\n";
			warn "Adding $operation_length to \$genomic_pos. Now at: $genomic_pos\n";
		    }
		}
		elsif($ops[$index] eq 'N'){ # skipped regions (= splice junctions) do affect the genomic position but not the $temp_pos within the read (similar to deletions)
		    if ($verbose){	  warn "Operation is 'N', adding $len[$index] bp to the genomic position\n";}
		    $genomic_pos += $len[$index];
		    if ($verbose){  warn "Endpos is now: $endpos\n";}
		}
		else{
		    die "Found CIGAR operations other than M, I, D or N: '$ops[$index]'. Not allowed at the moment\n";
		}
	    }

	    ### STEP IV: determining if the SNP is present in the SNP position list, and whether the position in the read was the reference or SNP base
	    if ($verbose){   warn "Found a match, genomic position should be: $genomic_pos\n~~~~~~~~~~~~~~\n\n"; sleep(1);}
	    if (exists $snps{$chr}->{$genomic_pos}){
		++$snp_found;

		if ($verbose){	warn "SNP was present: ref: $snps{$chr}->{$genomic_pos}->{ref}\tsnp: $snps{$chr}->{$genomic_pos}->{snp}\n";}
		if ($verbose){	warn "pos: $pos\n$sequence\n";}
		my $read_base = substr($sequence,$pos,1); # The N position has not been adjusted for the N yet so can be used as index

		if ($verbose){	warn "base in read was: $read_base\n";}
		if ($bisulfite){
		    ($genome1,$genome2,$irrelevant,$ct_unassigned) = score_bisulfite_SNPs($line,$chr,$genomic_pos,$read_base,$genome1,$genome2,$irrelevant,$ct_unassigned);
		}
		### Non-bisulfite setting ( = DEFAULT)
		else{
		    if ($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){
			if ($verbose){ warn "genome 1 specific!\n";}
			++$genome1;
		    }
		    elsif ($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){
			if ($verbose){ warn "genome 2 specific!\n";}
			++$genome2;
		    }
		    else{
			if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
			++$irrelevant;
		    }
		}
	    }
	    else{
		++$no_snp_found;
		if ($verbose){
		    warn "SNP was not present in the SNP list\n";
		}

	    }

	    if ($verbose){  warn "now adjusting the position for the N itself\n";}
	    $pos += 1; # adjusting position for the N itself. irrelevant if this is the last element of the array
	    if ($verbose){ warn "pos right now: $pos\n";}
	}

	if ($verbose){  warn "***************************\n\n\n\n";}

    }

    return ($genome1,$genome2,$ct_unassigned);

}



###################################################################################################################
###
### SUBROUTINES
###
###################################################################################################################



sub score_bisulfite_SNPs{
    my ($line,$chr,$genomic_pos,$read_base,$genome1,$genome2,$irrelevant,$ct_unassigned) = @_;
    my $first_read_conversion;
    my $genome_conversion;
    # warn "$line\n"; sleep(1);

    my $strand;
    while ( $line =~ /(XR|XG):Z:([^\t]+)/g ) {
	my $tag = $1;
	my $value = $2;
	chomp $value;

	if ($tag eq "XR") {
	    $first_read_conversion = $value;
	    $first_read_conversion =~ s/\r//;
	}
	elsif ($tag eq "XG") {
	    $genome_conversion = $value;
	    $genome_conversion =~ s/\r//;
	}
    }

    if ($first_read_conversion eq 'CT' and $genome_conversion eq 'CT') {
	$strand = 'OT';		## this is OT
    }
    elsif ($first_read_conversion eq 'GA' and $genome_conversion eq 'CT') {
	$strand = 'CTOT';		## this is CTOT
    }
    elsif ($first_read_conversion eq 'GA' and $genome_conversion eq 'GA') {
	$strand = 'CTOB';		## this is CTOB
    }
    elsif ($first_read_conversion eq 'CT' and $genome_conversion eq 'GA') {
	$strand = 'OB';		## this is OB
    }
    else {
	die "Unexpected combination of read and genome conversion: $first_read_conversion / $genome_conversion\n";
    }
    if ($verbose){  warn "strand is $strand\n"};

    ### if the SNP in question involves a genomic C position we potentially need to allow both C or T as C (G/A for the reverse strand)

    if ($verbose){ warn "Reference base is $snps{$chr}->{$genomic_pos}->{ref}\tSNP base:\t$snps{$chr}->{$genomic_pos}->{snp}\tRead base:\t$read_base\n"; sleep(1);}

    ### REFERENCE STRAIN = C
    if ( $snps{$chr}->{$genomic_pos}->{ref} eq 'C'){    # the position is a C in the reference genome

	### C>T SNP. Can't use the position for top strand alignments
	if ($snps{$chr}->{$genomic_pos}->{snp} eq 'T'){

	    if ($strand eq 'OT' or $strand eq 'CTOT'){
		if ($verbose){ warn "This is a C>T SNP and can't be considered for OT or CTOT alignments!\n";}
		++$ct_snp;
		++$ct_unassigned;
		return ($genome1,$genome2,$irrelevant,$ct_unassigned);
	    }
	}

	# Top strand alignments with C>T SNPs should have been skipped already. For other SNPs we accept both C or T as C.
	if ($strand eq 'OT' or $strand eq 'CTOT'){

	    if ($read_base eq 'C' or $read_base eq 'T'){ # a bisulfite converted position may be a C or T (here we know that the ref base is C)
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){
		if ($verbose){ warn "genome 2 specific!\n";}
		++$genome2;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }

	}
	# Alignments to the bottom strand can still be assigned normally as they do not involve the C directly
	elsif ($strand eq 'OB' or $strand eq 'CTOB'){

	    if ($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){ # == C
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    # added this section 13 June 2017
	    elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'G'){ # a bisulfite converted position on the bottom strand may be a G or A
		if ($read_base eq 'G' or $read_base eq 'A'){
		    if ($verbose){ warn "genome 2 specific!\n";}
		    ++$genome2;
		}
		else{
		    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		    ++$irrelevant;
		}
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){ # only A left...
		if ($verbose){ warn "genome 2 specific!\n";}
		++$genome2;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }

	}
	else{
	    die "Unable to identify read alignment strand\n";
	}
    }

    ### AlTERNATIVE STRAIN = C
    elsif ( $snps{$chr}->{$genomic_pos}->{snp} eq 'C'){ # the position is a C in the alternative/SNP genome

	### C>T SNP. Can't use the position for top strand alignments
	if ($snps{$chr}->{$genomic_pos}->{ref} eq 'T'){

	    if ($strand eq 'OT' or $strand eq 'CTOT'){
		if ($verbose){ warn "This is a C>T SNP between alternative and reference and can't be considered for OT or CTOT alignments!\n";}
		++$ct_snp;
		++$ct_unassigned;
		return ($genome1,$genome2,$irrelevant,$ct_unassigned);
	    }

	}

	# Top strand alignments with C>T SNPs should have been skipped already. For other SNPs we accept both C or T as C.
	if ($strand eq 'OT' or $strand eq 'CTOT'){

	    if ($read_base eq 'C' or $read_base eq 'T'){ # a bisulfite converted position may be a C or T
		if ($verbose){ warn "genome 2 specific!\n";}
		# warn "Alternative base is C\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome2;
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }
	}
	elsif ($strand eq 'OB' or $strand eq 'CTOB'){ # Reads to the bottom strand can still be assigned normally

	    if ($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){ # has to be C
		if ($verbose){ warn "genome 2 specific!\n";}
		# warn "Alternative base is C\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome2;
	    }

	    # added this section 13 June 2017
	    elsif ($snps{$chr}->{$genomic_pos}->{ref} eq 'G'){ # a bisulfite converted position on the bottom strand may be a G or A
		if ($read_base eq 'G' or $read_base eq 'A'){
		    if ($verbose){ warn "genome 1 specific!\n";}
		    ++$genome1;
		}
		else{
		    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		    ++$irrelevant;
		}
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){ # only A left...
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }

	}
	else{
	    die "Unable to identify read alignment strand\n";
	}
    }
    ### REFERENCE STRAIN = G, i.e. bottom strand C might be affected
    elsif ( $snps{$chr}->{$genomic_pos}->{ref} eq 'G'){    # the position is a G in the reference genome

	if ($snps{$chr}->{$genomic_pos}->{snp} eq 'A'){

	    if ($strand eq 'OB' or $strand eq 'CTOB'){
		if ($verbose){ warn "This is a G>A SNP and can't be considered for OB or CTOB alignments!\n";}
		++$ct_snp;
		++$ct_unassigned;
		return ($genome1,$genome2,$irrelevant,$ct_unassigned);
	    }
	}

	# Bottom strand alignments with G>A SNPs should have been skipped already
	if ($strand eq 'OT' or $strand eq 'CTOT'){  # Reads to the top strand can be assigned normally

	    if ($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){ # eq G
		if ($verbose){ warn "genome 1 specific!\n";}
		# warn "Reference base is G\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome1;
	    }
	    # added this section 13 June 2017
	    elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'C'){ # a bisulfite converted position on the bottom strand may be a G or A
		if ($read_base eq 'C' or $read_base eq 'T'){
		    if ($verbose){ warn "genome 2 specific!\n";}
		    ++$genome2;
		}
		else{
		    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		    ++$irrelevant;
		}
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){ # has to be T
		if ($verbose){ warn "genome 2 specific!\n";}
		++$genome2;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }

	}
	elsif ($strand eq 'OB' or $strand eq 'CTOB'){ # Reads to the bottom strand might carry a bisulfite conversion

	    if ($read_base eq 'G' or $read_base eq 'A'){ # a bisulfite converted position may be a G or A
		if ($verbose){ warn "genome 1 specific!\n";}
		# warn "Reference base is G\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome1;
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){
		if ($verbose){ warn "genome 2 specific!\n";}
		++$genome2;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }

	}
	else{
	    die "Unable to identify read alignment strand\n";
	}

    }
    ### AlTERNATIVE STRAIN = G, i.e. bottom strand C might be affected
    elsif ( $snps{$chr}->{$genomic_pos}->{snp} eq 'G'){ # the position is a G in the alternative/SNP genome

	if ($snps{$chr}->{$genomic_pos}->{ref} eq 'A'){

	    if ($strand eq 'OB' or $strand eq 'CTOB'){
		if ($verbose){ warn "This is a G>A SNP between alternative and reference genome and can't be considered for OB or CTOB alignments!\n"; }
		++$ct_snp;
		++$ct_unassigned;
		return ($genome1,$genome2,$irrelevant,$ct_unassigned); # we are also returning 1 so we can discriminate between un-assignable and rejected C>T SNPs
	    }

	}

	# Bottom strand alignments with G>A SNPs should have been skipped already
	if ($strand eq 'OT' or $strand eq 'CTOT'){  # Reads to the top strand can be assigned normally

	    if ($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){ # has to be G
		if ($verbose){ warn "genome 2 specific!\n";}
		# warn "Alternative base is G\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome2;
	    }
	    # added this section 13 June 2017
	    elsif ($snps{$chr}->{$genomic_pos}->{ref} eq 'C'){ # a bisulfite converted position on the bottom strand may be a G or A
		if ($read_base eq 'C' or $read_base eq 'T'){
		    if ($verbose){ warn "genome 1 specific!\n";}
		    ++$genome1;
		}
		else{
		    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		    ++$irrelevant;
		}
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){ # has to be T
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }
	}
	elsif ($strand eq 'OB' or $strand eq 'CTOB'){ # Reads to the bottom strand might carry a bisulfite conversion

	    if ($read_base eq 'G' or $read_base eq 'A'){ # a bisulfite converted position may be a G or A
		if ($verbose){ warn "genome 2 specific!\n";}
		# warn "Alternative base is G\nGenomic base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
		++$genome2;
	    }
	    elsif($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){
		if ($verbose){ warn "genome 1 specific!\n";}
		++$genome1;
	    }
	    else{
		if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
		++$irrelevant;
	    }
	}
	else{
	    die "Unable to identify read alignment strand\n";
	}
    }

    ### The SNP did not affect a C or G position at all
    else{
	# warn "they do exist! Ref base:\t$snps{$chr}->{$genomic_pos}->{ref}\tRead base:\t$read_base\n"; sleep(1);
	if ($read_base eq $snps{$chr}->{$genomic_pos}->{ref}){
	    if ($verbose){ warn "genome 1 specific!\n";}
	    ++$genome1;
	}
	elsif ($read_base eq $snps{$chr}->{$genomic_pos}->{snp}){
	    if ($verbose){ warn "genome 2 specific!\n";}
	    ++$genome2;
	}
	else{
	    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	    ++$irrelevant;
	}
    }

    return ($genome1,$genome2,$irrelevant,$ct_unassigned);

}


sub sam2bam_convert{

    my $file = shift;

    warn "File appears to be in SAM format and needs conversion into BAM format first\n";
    sleep (3);
    my $bamfile = $file;
    $bamfile =~ s/sam$/bam/;
    warn "Now converting '$file' to '$bamfile' ...\n";
    sleep (3);

    system ("$samtools_path view -bS $file > $bamfile");

    warn "BAM conversion finished\n\n";
    sleep (3);
    return $bamfile;  # setting bamfile as new input file

}

sub sort_by_name_paired_end{

    my $file = shift;

    my $sorted = $file;
    $sorted =~ s/bam$/sortedByName.bam/;

    warn "Sorting paired-end BAM file '$file' by read IDs ...\n";
    #   system ("$samtools_path sort -O BAM -n -o $sorted $file"); # added -O BAM on 13 Dec 2016
    system ("$samtools_path sort -n -o $sorted $file");  # removed -O again because it doesn't seem to be working... 09 Jan 2017
    $file = $sorted; # setting sorted BAM file as new $file
    warn "Finished sorting BAM file into new file '$file'\n\n";
    sleep (1);
    return $file;

}

sub read_snps{

    if ($snp_file =~ /\.gz$/){
	open (SNP,"gunzip -c $snp_file |") or die "Failed to read from gzipped $snp_file: $!\n\n";
    }
    else{
	open (SNP,$snp_file) or die "Failed to read from $snp_file: $!\n\n";
    }

    ### READING FROM SNP FILE

    warn "Storing SNP positions provided in '$snp_file'\n";
    sleep(1);

    my $snp_count = 0;

    while (<SNP>){
	++$snp_count;
	$_ =~ s/\r|\n//g; # removing CR and LF line endings
	my ($chr,$pos,$diff)  = (split /\t/)[1,2,4];
	my ($ref,$snp) = (split /\//,$diff);
	# print "$chr\t$pos\t$ref\t$snp\n";

	### for N-masked genomes the SNPs are stored in an. hashes
	$snps{$chr}->{$pos}->{ref} = $ref;
	$snps{$chr}->{$pos}->{snp} = $snp;
    }
    warn "Stored $snp_count positions in total\n\n";
    return ($snp_count);

}

sub read_snps_bisulfite{

    warn "Processing SNPs provided in '$snp_file' in a bisulfite-aware fashion...\n";

    if ($snp_file =~ /\.gz$/){
	open (SNP,"gunzip -c $snp_file |") or die "Failed to read from zgipped $snp_file: $!\n\n";
    }
    else{
	open (SNP,$snp_file) or die "Failed to read from $snp_file: $!\n\n";
    }

    ### READING FROM SNP FILE
    my $snp_count = 0;
    my $total     = 0;
    my $ct_snp    = 0;
    my $ga_snp    = 0;

    while (<SNP>){
	chomp;
	my ($chr,$pos,$diff)  = (split /\t/)[1,2,4];
	my ($ref,$snp) = (split /\//,$diff);
	# print "$chr\t$pos\t$ref\t$snp\n";
	++$total;

	if ($ref eq 'C' and $snp eq 'T'){
	    $ct_snp++;
	    next;
	}
	if ($ref eq 'T' and $snp eq 'C'){
	    $ct_snp++;
	    next;
	}

	if ($ref eq 'G' and $snp eq 'A'){
	    $ga_snp++;
	    next;
	}
	if ($ref eq 'A' and $snp eq 'G'){
	    $ga_snp++;
	    next;
	}

	#sleep(1);

	++$snp_count;

	### for N-masked genomes the SNPs are stored in an. hashes
	$snps{$chr}->{$pos}->{ref} = $ref;
	$snps{$chr}->{$pos}->{snp} = $snp;
    }

    my $perc_CT = sprintf("%.2f",$ct_snp/$total *100);
    my $perc_GA = sprintf("%.2f",$ga_snp/$total *100);
    my $perc_stored = sprintf("%.2f",$snp_count/$total *100);

    warn "Stored positions in total:\t\t$snp_count ($perc_stored%)\n";
    warn "Excluded C>T positions in total:\t$ct_snp ($perc_CT%)\n";
    warn "Excluded G>A positions in total:\t$ga_snp ($perc_GA%)\n\n";

    print REPORT "Stored positions in total:\t\t$snp_count ($perc_stored%)\n";
    print REPORT "Excluded C>T positions in total:\t$ct_snp ($perc_CT%)\n";
    print REPORT "Excluded G>A positions in total:\t$ga_snp ($perc_GA%)\n\n";

    return ($snp_count);

}


sub sort_snps{

    warn "Sorting SNP positions\n";

    foreach my $chr (keys %snps){
	warn "Sorting SNP positions on chromosome $chr...\n";
	@{$snps{$chr}} =  sort {$a->{pos}<=>$b->{pos}} @{$snps{$chr}};
    }
    warn "Finished sorting SNP positions\n\n\n";

}

sub open_output_filehandles{

    my $file = shift;

    # replacing folder information before creating output file names 24 07 2017
    $file =~ s/.*\///;

    my $outfile = my $genome1_file = my $genome2_file = my $unassigned_file = my $impossible_to_distinguish_file = my $report_file = $file;

    # if files were sorted they look like this
    $genome1_file =~ s/sorted\.//;
    $genome2_file =~ s/sorted\.//;
    $unassigned_file =~ s/sorted\.//;
    $impossible_to_distinguish_file =~ s/sorted\.//;
    $report_file =~ s/sorted\.//;

    # if files were supplied with --no_sort then the files might end in .bam only
    $genome1_file =~ s/bam$/genome1.sam/;
    $genome2_file =~ s/bam$/genome2.sam/;
    $unassigned_file =~ s/bam$/unassigned.sam/;
    $impossible_to_distinguish_file =~ s/bam$/undecided.sam/;

    $outfile =~ s/sam$/bam/;
    $outfile =~ s/bam$/allele_flagged.sam/;

    $report_file =~ s/sam$/bam/;
    $report_file =~ s/bam$/SNPsplit_report.txt/;

    if ($bam){
	$outfile =~ s/sam$/bam/;
	open (OUT,"| $samtools_path view -bS 2> /dev/null - > ${output_dir}$outfile") or die "Unable to write to BAM file '$outfile': $!\n";
    }
    else{ # writing out to a SAM file
	open (OUT,'>',"${output_dir}$outfile") or die "Unable to write to file '${output_dir}$outfile': $!\n";
    }

    open (REPORT,'>',"${output_dir}$report_file") or die "Unable to write to file '${output_dir}$report_file': $!\n";

    warn "Input file:\t\t\t\t\t'$file'\n";
    print REPORT "Input file:\t\t\t\t\t'$file'\n";

    warn "Writing SNPplit-tag report to:\t\t\t'$report_file'\n";

    warn "Writing allele-flagged output file to:\t\t'$outfile'\n\n";
    print REPORT "Writing allele-flagged output file to:\t\t'$outfile'\n\n";

    return $outfile;

}



#######################################
###   Bisfulfite Read Detection   #####
#######################################

sub check_for_bs{

    my ($paired_end,$samtools_path) = @_;
    my $isBismark;

    # SAM/BAM format, checking the first file supplied

    my $file = $ARGV[0];

    ### if the user did not specify whether the alignment file was single-end or paired-end we are trying to get this information from the @PG header line in the SAM/BAM file
    if ($file =~ /\.gz$/){
	open (DETERMINE,"gunzip -c $file |") or die "Unable to read from gzipped file $file: $!\n";
    }
    elsif ($file =~ /\.bam$/ ){
	open (DETERMINE,"$samtools_path view -h $file |") or die "Unable to read from BAM file $file: $!\n";
    }
    else{
	open (DETERMINE,$file) or die "Unable to read from $file: $!\n";
    }

    while (<DETERMINE>){
	last unless (/^\@/);
	if ($_ =~ /^\@PG/){
	    #  warn "found the \@PG line:\n";
	    # warn "$_";

	    if ($_ =~ /Bismark/){

		$isBismark = 1;

		unless ($paired_end){ # unless the user specified --paired we try to extract this information from the @PG line

		    if ($_ =~ /\s+-1\s+/ and $_ =~ /\s+-2\s+/){
			warn "Treating file(s) as paired-end data (as extracted from \@PG line)\n\n";
			$paired_end = 1;
		    }
		    else{
			warn "Treating file(s) as single-end data (as extracted from \@PG line)\n\n";
			$paired_end = 0;
		    }
		}
	    }
	    else{
		$isBismark = 0;
	    }
	}
    }
    return ($isBismark,$paired_end);
}


sub test_positional_sorting{

    my $filename = shift;
    my $samtools_path = shift;
    warn "\nNow testing Bismark result file $filename for positional sorting\n";
    sleep(1);

    if ($filename =~ /\.gz$/) {
	open (TEST,"gunzip -c $filename |") or die "Can't open gzipped file $filename: $!\n";
    }
    elsif ($filename =~ /bam$/){ ### this would allow to read BAM files that do not end in *.bam
	if ($samtools_path){
	    open (TEST,"$samtools_path view -h $filename |") or die "Can't open BAM file $filename: $!\n";
	}
	else{
	    die "Sorry couldn't find an installation of Samtools. Either specifiy an alternative path using the option '--samtools_path /your/path/', or use a SAM file instead\n\n";
	}
    }
    else{
	open (TEST,$filename) or die "Can't open file $filename: $!\n";
    }

    my $count = 0;

    while (<TEST>) {
	if (/^\@/) {	     # testing header lines if they contain the @SO flag (for being sorted)
	    if (/^\@SO/) {
		warn "SAM header line '$_' indicates that the Bismark aligment file has been sorted by chromosomal positions. Paired-end files will be sorted by name before proceeding with the SNPsplit tagging process\n\n";
		return 0;
	    }
	    next;
	}
	$count++;

	last if ($count > 100000); # else we test the first 100000 sequences if they start with the same read ID

	my ($id_1) = (split (/\t/));

	### reading the next line which should be read 2
	$_ = <TEST>;
	my ($id_2) = (split (/\t/));
	last unless ($id_2);
	++$count;

	if ($id_1 eq $id_2){
	    ### ids are the same
	    next;
	}
	else{ ### in previous versions of Bismark we appended /1 and /2 to the read IDs for easier eyeballing which read is which. These tags need to be removed first
	    my $id_1_trunc = $id_1;
	    $id_1_trunc =~ s/\/1$//;
	    my $id_2_trunc = $id_2;
	    $id_2_trunc =~ s/\/2$//;

	    unless ($id_1_trunc eq $id_2_trunc){
		warn "The IDs of Read 1 ($id_1) and Read 2 ($id_2) are not the same. This might be a result of sorting the paired-end SAM/BAM files by chromosomal position. Paired-end files will be sorted by name before proceeding with the SNPsplit tagging process\n\n";
		return 0;
	    }
	}
    }

    ### If it hasn't returned so far then it seems the file is in the correct Bismark format (read 1 and read 2 of a pair directly following each other)
    return 1;
}


#######################################
###   Bisfulfite Read Detection   #####
#######################################

sub process_commandline{
    my $help;
    my $version;
    my $snp_file;
    my $no_sort;
    my $verbose;
    my $samtools_path;
    my $sam;
    my $paired;
    my $hic_data;
    my $conflict;
    my $output_dir;
    my $singletons;
    my $bisulfite;

    my $command_line = GetOptions ('help|man'                => \$help,
				   'versions'                => \$version,
				   'output_dir=s'            => \$output_dir,
				   'SNP_file=s'              => \$snp_file,
				   'no_sorting'              => \$no_sort,
				   'verbose'                 => \$verbose,
				   'samtools_path=s'         => \$samtools_path,
				   'sam'                     => \$sam,
				   'paired'                  => \$paired,
				   'hic'                     => \$hic_data,
				   'conflicting|weird'       => \$conflict,
				   'singletons'              => \$singletons,
				   'bisulfite'               => \$bisulfite,
	);

    ### EXIT ON ERROR if there were errors with any of the supplied options
    unless ($command_line){
	die "Please respecify command line options\n";
    }

    ### HELPFILE
    if ($help){
	print_helpfile();
	exit;
    }

    if ($version){
	print << "VERSION";
    SNPsplit - Allele-specific alignment sorting
	        Version: $pipeline_version
	  Copyright 2014-17 Felix Krueger
             Babraham Bioinformatics
       https://github.com/FelixKrueger/SNPsplit
VERSION
	  exit;
    }

    ### no files provided
    unless (@ARGV){
	warn "You need to provide one or more SAM/BAM files to start the allele-specific pipeline. Please respecify!\n";
	print_helpfile();
	exit;
    }


    unless ($snp_file){
	die "You need to provide a text file detailing SNPs with '--SNP_file your.file'. Please respecify!\n";
    }

    if ($no_sort){
	foreach my $file (@ARGV){
	    die "The option --no_sorting requires you to supply a BAM file (ending in .bam), but at least one of the file names did not look like a BAM file: '$file'. Please respecify!\n\n" unless ($file =~ /\.bam$/);
	}
    }

    ### checking to see if all files are indeed present in the folder
    foreach my $file (@ARGV){
	unless(-e $file){
	    die "Input file '$file' doesn't exist in the input folder. Please check filenames and try again!\n\n";
	}
    }

    ### PATH TO SAMTOOLS
    if (defined $samtools_path){
	# if Samtools was specified as full command
	if ($samtools_path =~ /samtools$/){
	    if (-e $samtools_path){
		# Samtools executable found
	    }
	    else{
		die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
	    }
	}
	else{
	    unless ($samtools_path =~ /\/$/){
		$samtools_path =~ s/$/\//;
	    }
	    $samtools_path .= 'samtools';
	    if (-e $samtools_path){
		# Samtools executable found
	    }
	    else{
		die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
	    }
	}

    }
    # Check whether Samtools is in the PATH if no path was supplied by the user
    else{
	if (!system "which samtools >/dev/null 2>&1"){ # STDOUT is binned, STDERR is redirected to STDOUT. Returns 0 if samtools is in the PATH
	    $samtools_path = `which samtools`;
	    chomp $samtools_path;
	}
    }

    unless (defined $samtools_path){
	die "Could not find an installation of Samtools on your system. Please either install Samtools first or provide the path to an existing installation\n\n";
    }

    ### --singletons only make sense for paired-end data
    if ($singletons){
	if ($hic){
	    die "The option --singletons can't be used in Hi-C mode since this expects paired-end data. Please respecify!\n\n";
	}

	unless ($paired){
	    warn "The option --singletons only makes sense for paired-end files. Simply ignoring it in single-end mode...\n"; sleep(1);
	    $singletons = undef;
	}
    }

    if ($hic_data){
	warn "Hi-C data specified. This assumes that the data has been processed with HiCUP, is paired-end and does not require positional sorting.\nSetting '--paired' and '--no_sort'...\n\n";
	$paired = 1;
	$no_sort = 1;
	sleep(1);
    }

    ### OUTPUT DIRECTORY PATH
    chdir $parent_dir or die "Failed to move back to current working directory\n";
    if (defined $output_dir){
	unless ($output_dir eq ''){
	    unless ($output_dir =~ /\/$/){
		$output_dir =~ s/$/\//;
	    }
	}

	if (chdir $output_dir){
	    $output_dir = getcwd(); #  making the path absolute
	    unless ($output_dir =~ /\/$/){
		$output_dir =~ s/$/\//;
	    }
	}
	else{
	    mkdir $output_dir or die "Unable to create directory $output_dir $!\n";
	    warn "Created output directory $output_dir!\n\n";
	    chdir $output_dir or die "Failed to move to $output_dir\n";
	    $output_dir = getcwd; #  making the path absolute
	    unless ($output_dir =~ /\/$/){
		$output_dir =~ s/$/\//;
	    }
	}
	warn "Output will be written into the directory: $output_dir\n";
    }
    else{
	$output_dir = '';
    }
    chdir $parent_dir or die "Failed to move back to current working directory\n";

    if ($sam){
	$bam = 0;
    }
    else{
	$bam = 1;
    }

    ### Checking to see if the input file is a Bismark processed Bisulfite-Seq file
    unless ($hic_data){

	warn "Testing if input file '$ARGV[0]' looks like a Bisulfite-Seq file\n";
	(my $isBismark,$paired) = check_for_bs($paired,$samtools_path);
	if ($isBismark){
	    $bisulfite++;
	    if ($paired){
		warn "File looks like a Bismark paired-end file. Setting '--bisulfite' and '--paired'...\n";

		my ($skip_sorting_by_name) = test_positional_sorting($ARGV[0],$samtools_path);
		if ($skip_sorting_by_name){
		    warn "File appears to be in the right format (Read1 and Read2 following each other), setting '--no_sort'...\n";
		    $no_sort = 1;
		}
		else{
		    warn "File needs sorting by name first (Read1 and Read2 are required to follow each other for SNPsplit processing)...\n";
		}
	    }
	    else{
		warn "File looks like a Bismark single-end file. Setting '--bisulfite'...\n";
	    }
	}
    }

    return ($snp_file,$no_sort,$verbose,$samtools_path,$bam,$paired,$hic_data,$conflict,$output_dir,$singletons,$bisulfite);
}



sub print_helpfile{
    print <<EOF
      SYNOPSIS:
	SNPsplit is designed to read in alignment files in SAM/BAM format and determine the allelic origin of reads
	that cover known SNP positions. For this to work the alignment step must have been carried out against a
	genome that had all SNP positions replaced by Ns, and a list of SNP positions between the two different genomes
	has to be provided using the option --snp_file.
	SNPsplit operates in two stages:
	I) SNPsplit-tag: SNPsplit analyses reads (single-end mode) or read pairs (paired-end mode) for overlaps
        with known SNP positions, and writes out a tagged BAM file in the same order as the original file.
        II) SNPsplit-sort: the tagged BAM file is read in and is being sorted into allele-specific files. This process
        may also be run as a stand-alone module (tag2sort).
        The SNPsplit-tag module determines whether a read can be assigned to a certain allele and appends an additional
        optional field 'XX:Z:' to each read. The tag can be one of the following:
        XX:Z:UA - Unassigned
        XX:Z:G1 - Genome 1-specific
        XX:Z:G2 - Genome 2-specific
        XX:Z:CF - Conflicting
        The SNPsplit-sort module reads in the tagged BAM file and sorts the reads (or read pairs) according to their XX:Z:
        tag (or the combination of tags for paired-end or Hi-C reads) into subfiles.
      USAGE:      SNPsplit  [options] --snp_file <SNP.file.gz> [input file(s)]
    input file(s)          Mapping output file in SAM or BAM format. SAM files (ending in .sam) will first be
    converted to BAM using Samtools.
    --snp_file <file>      Mandatory file specifying SNP positions to be considered; may be a plain text file or
                           gzip compressed. Currently, the SNP file is expected to be in the following format:
                                       SNP-ID     Chromosome  Position    Strand   Ref/SNP
                           example:   33941939        9       68878541       1       T/G
                           Only the information contained in fields 'Chromosome', 'Position' and 'Reference/SNP base'
                           are being used for analysis. The genome referred to as 'Ref' will be used as genome 1,
                           the genome containing the 'SNP' base as genome 2.
    --paired               Paired-end mode. (Default: OFF).
    -o/--outdir <dir>      Write all output files into this directory. By default the output files will be written into
                           the same folder as the input file(s). If the specified folder does not exist, SNPsplit will attempt
                           to create it first. The path to the output folder can be either relative or absolute.
    --singletons           If the allele-tagged paired-end file also contains singleton alignments (which is the
                           default for e.g. TopHat), these will be written out to extra files (ending in _st.bam)
                           instead of writing everything to combined paired-end and singleton files. Default: OFF.
    --no_sort              This option skips the sorting step if BAM files are already sorted by read name (e.g.
                           Hi-C files generated by HiCUP). Please note that setting --no_sort for unsorted paired-end
                           files will break the tagging process!
    --hic                  Assumes Hi-C data processed with HiCUP (www.bioinformatics.babraham.ac.uk/projects/hicup/)
                           as input, i.e. the input BAM file is paired-end and Reads 1 and 2 follow each other. Thus,
                           this option also sets the flags --paired and --no_sort. Default: OFF.
    --bisulfite            Assumes Bisulfite-Seq data processed with Bismark (www.bioinformatics.babraham.ac.uk/projects/bismark/)
                           as input. In paired-end mode (--paired), Read 1 and Read 2 of a pair are expected to follow
                           each other in consecutive lines. SNPsplit will run a quick check at the start of a run to see if the
                           file provided appears to be a Bismark file, and set the flags --bisulfite and/or --paired
                           automatically. In addition it will perform a quick check to see if a paired-end file appears
                           to have been positionally sorted, and if not will set the --no_sort flag.
    --samtools_path        The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to
                           be specified explicitly if Samtools is in the PATH already.
    --verbose              Verbose output (for debugging).
    SNPsplit-sort specific options (tag2sort):
    ==========================================
    --conflicting/--weird  Reads or read pairs that were classified as 'Conflicting' (XX:Z:CF) will be written to
                           an extra file (ending in .conflicting.bam) instead of being simply skipped. Reads may be
                           classified as 'Conflicting' if a single read contains SNP information for both genomes at
                           the same time, or if the SNP position was deleted from the read. Read-pairs are considered
                           'Conflicting' if either read is was tagged with the XX:Z:CF flag. Default: OFF.
    --sam                  The output will be written out in SAM format instead of BAM (default). SNPsplit will attempt to use
                           the path to Samtools that was specified with --samtools_path, or, if it hasn't been
                           specified, attempt to find Samtools in the PATH environment. If no installation of
                           Samtools can be found, the SAM output will be compressed with GZIP instead (yielding a
                           .sam.gz output file).
    --help                 Displays this help information and exits.
    --version              Displays version information and exits.
                                   Script last modified: 20 July 2017
EOF
    ;
  exit 1;
}