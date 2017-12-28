#!/usr/bin/perl
use strict ;
use warnings ;

##################################################################################################################################
### CALCULATES HARD ALLELE FREQUENCY FROM VCF FILE OF TETRAPLOIDS
### REQUIRED INPUTS:
###		1.) Coverage cutoff, or number of reads supporting a base call for an individual
###		2.)	Proportion of missing data allowed per population
##################################################################################################################################

# Input arguments
my $coverage_cutoff = $ARGV[0] ;	# missing data threshold 
my $Prop_missdata = $ARGV[1] ;	# number of individuals allowed to be missing in smallest population
my $VCF_filepath = $ARGV[2] ;
my $PopIDFile = $ARGV[3] ;

# Declare data structures
my %PopIDs ;
my %PopTotals ;
my %PopDownsampled ;
my %MissingDataPerPop ;
my %pop_indices ; #keys = populations, values = array of indicies in VCF
my %AF ; #AF{scaff}{site} = AF ;

open QC, ">./QCFile" ;

#############################################################
# make sure arguments supplied to program
unless($coverage_cutoff){
	print "The coverage cutoff must be greater than 0\nexiting...\n" ;
	exit ;
}
unless( $Prop_missdata > 0 && $Prop_missdata <= 1 ){
	print "The missing data threshold must be greater than 0 and less than or equal to 1\nexiting...\n" ;
	exit ;
}
unless($VCF_filepath){
	print "Please provide a file path to the VCF file\nexiting...\n" ;
	exit ;
}
#############################################################
# open PopIDFile, which assigns individuals to populations
open(my $fh, "<", $PopIDFile) or die "Could not open file ${PopIDFile}. $!";
while (<$fh>){
  chomp $_ ;
	if($_ !~ m/^Individual/){
		my @line = split(/\t/, $_) ;
		$PopIDs{$line[0]} = $line[1] ;
		$PopTotals{$line[1]}++ ;
	}
}
close $fh;
#############################################################
# Calculate downsampled population sizes accounting for missing data
foreach my $pop (sort {$a cmp $b} keys %PopTotals){
	my $x = $PopTotals{$pop}*$Prop_missdata ;
	if(int($x) == $x){
		$PopDownsampled{$pop} = $x ;
	}else{
		$PopDownsampled{$pop} = int($x) + 1 ;
	}
}
print QC "Number of individuals per population:\n" ;
print QC "Population\tNumberIndividuals\tNumberIndividualsDownsampled\n" ;
foreach my $pop (sort{$a cmp $b} keys %PopTotals){
	print QC $pop, "\t", $PopTotals{$pop}, "\t", $PopDownsampled{$pop}, "\n" ;
}
print QC "\n" ;
#############################################################
# Read through VCF file
open(VCFFILE, "<", $VCF_filepath) or die "Could not open VCF file ${VCF_filepath}. When in doubt, just specify the full file path $!";
while(<VCFFILE>){
	chomp $_ ;
	if ( $_ =~ m/^\#\#/ ) { 
			next ; 
	}elsif ( $_ =~ m/^\#CHROM/ ) { 
		my @header = split ( /\t/, $_ ) ; 
		foreach my $index (9..$#header) {
			chomp $header[$index] ;
			my $ind = $header[$index] ;		 
			my $pop = $PopIDs{$ind} ;
			unless($pop){
				print "Warning: Individual \"${ind}\" was identified in the VCF file but not in population identifier file. " ;
				print "Please fix this unless you're intentionally excluding this individual\n" ;
			}
			if($pop){
				push @{$pop_indices{$pop}}, $index ; 
			}				
		}
		next ;
	}else{
		my @vcf_line = split ( /\t/, $_ ) ; 
		my %missing_individuals  ;
		###INITIALIZE MISSING IND HASH
		foreach my $pop ( keys %pop_indices ){
			$missing_individuals{$pop} = 0 ;
		}
		## Get field that corresponds to DP; it sometimes varies
		my @info = split(/:/,$vcf_line[8]) ;
		my $DP_index = "NA" ;
		foreach (0..$#info){
			if($info[$_] eq "DP"){
				$DP_index = $_ ;
			}
		}
		if($DP_index ne "NA"){
			## see how many individuals have low depth or no data
			foreach my $pop ( keys %pop_indices ){
				foreach my $index ( @{$pop_indices{$pop}} ){
					my @ind = split (/:/, $vcf_line[$index]) ;
					#USE 2ND FIELD HERE TO CHECK DP
					if( ($#ind == 0) ){
						$missing_individuals{$pop} ++ ;
					}elsif($#ind > 0){
						if($ind[$DP_index] eq "."){
							$missing_individuals{$pop} ++ ;
						}elsif( $ind[$DP_index] < $coverage_cutoff ){
							$missing_individuals{$pop} ++ ;
						} 	
					}
				}
			}
			## see if any population has more missing data than specified threshold	
			my $pops_with_not_enough_ind = 0 ;
			foreach my $pop ( keys %pop_indices ){
				my $pop_size = scalar @{$pop_indices{$pop}} ;
				if( (($pop_size - $missing_individuals{$pop})/$pop_size) < $Prop_missdata ){
					$MissingDataPerPop{$pop}++ ;
					$pops_with_not_enough_ind ++ ;
				}
			}
			unless( $pops_with_not_enough_ind ){ 
				my $scaffold = $vcf_line[0] ;
				my $position = $vcf_line[1] ;
				if($vcf_line[4] eq "."){ # Site is invariant
					foreach my $pop (keys %pop_indices){
						$AF{$scaffold}{$position}{$pop} = 0 ;	
					}	
				}
				elsif( length($vcf_line[4])==1 && $vcf_line[4] =~ m/[AGTC]/ ){
					foreach my $pop ( keys %pop_indices ){
						my $pop_sample_size = ( $PopTotals{$pop} - $missing_individuals{$pop} ) ;
						my $pop_hard_allele_freq = extract_genos(\$coverage_cutoff, \@vcf_line, \%pop_indices, \$pop) ;
						if( $pop_sample_size > $PopDownsampled{$pop} ){
							# Randomly subsample alleles (without replacement) if there is more data than specified by the downsampled size
							my @allele_counts ;
							foreach(1 .. $pop_hard_allele_freq){
								push @allele_counts, 1 ;
							}
							foreach(1 .. ($pop_sample_size*(4)-$pop_hard_allele_freq)){
								push @allele_counts, 0 ;
							}
							fisher_yates_shuffle(\@allele_counts) ;
							my $sum ;
							foreach ( 0 .. ($pop_sample_size*(4)-1) ){
								$sum += $allele_counts[$_] ;
							}
								$AF{$scaffold}{$position}{$pop} = $sum ;
						}else{
							$AF{$scaffold}{$position}{$pop} = $pop_hard_allele_freq ;
						}
					}			
				}	
			}
		}	
	}
}
close VCFFILE ;

#############################################################
# Printing missing data per population
print QC "Amount of missing Data per population:\n" ;
print QC "Population\tMissingData\n" ;
foreach my $pop (sort{$a cmp $b} keys %MissingDataPerPop){
	print QC $pop, "\t", $MissingDataPerPop{$pop}, "\n" ;
}
print QC "\n" ;
#############################################################
# Printing population allele frequencies
open OUT, ">./AFS_DepthCutoff_${coverage_cutoff}_RequiredPropGenotypes_${Prop_missdata}.txt" ;
	print OUT "scaffold", "\t", "position", "\t" ;
	#print population labels
	foreach my $pop (sort{$a cmp $b}(keys %pop_indices) ){
		print OUT $pop, "\t" ;
	}
	print OUT "\n" ;
	foreach my $scaffold ( sort{$a cmp $b}(keys %AF) ){
		foreach my $pos (sort{$a <=> $b}(keys %{$AF{$scaffold}}) ){
			print OUT $scaffold, "\t", $pos, "\t" ;
			foreach my $pop ( sort{$a cmp $b}( keys %{$AF{$scaffold}{$pos}}) ){
				print OUT $AF{$scaffold}{$pos}{$pop}/($PopDownsampled{$pop}*4), "\t" ;
			}
			print OUT "\n" ; #print newline after each position
		}	
	}
close OUT ;

exit ;
#########################################

sub extract_genos{
	my $coverage_cutoff = $_[0] ;
	my $vcf_line = $_[1] ;
	my $pop_indices = $_[2] ;
	my $pop = $_[3] ;
	
	my $pop_allele_freq = 0 ;
	
	foreach ( @$vcf_line[@{$$pop_indices{$$pop}}] ){
		my @ind = split (/:/, $_) ;
		if( ($#ind < 1) || ($ind[2] < $$coverage_cutoff) ){
			next ;
		}else{	
			if( $ind[0] =~ m/0\/0\/0\/0/){
				$pop_allele_freq += 0 ;
			}elsif( $ind[0] =~ m/0\/0\/0\/1/ ){
				$pop_allele_freq += 1 ;
			}elsif( $ind[0] =~ m/0\/0\/1\/1/ ){
				$pop_allele_freq += 2 ;
			}elsif( $ind[0] =~ m/0\/1\/1\/1/ ){
				$pop_allele_freq += 3 ;
			}elsif( $ind[0] =~ m/1\/1\/1\/1/ ){
				$pop_allele_freq += 4 ;
			}
		}	
	}
	return ($pop_allele_freq) ;
}


sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
