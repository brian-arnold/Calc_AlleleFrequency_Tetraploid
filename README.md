# Calc_AlleleFrequency_Tetraploid.pl

This program calculates allele frequencies from a VCF file containing tetraploid genotypes.
At the moment this script only accommodates tetraploid data and not mixed ploidy

Usage: perl Calc_AlleleFrequency_Tetraploid.pl Argument1 Argument2 Argument3 Argument4

########################################

Input Files:

1.) VCF file with tetraploid genotypes from one or more populations

2.) Population Identification File. This file contains two columns and should be in tab-
delimited format. The first column is a list of names, one for each individual. These 
names must correspond to those present in the VCF file header (specifically the line that
begins with "#CHROM". The second column lists the population to which each individual 
belongs. These population names can be letters, numbers, or words, but should not contain
white spaces (especially tabs). This file does not need a header line, but if you choose 
to supply one, please label the first column "Individual", the second column "Population".
NOTE: if this file does not contain individual sample names present in the VCF file, then
these individuals will be ignored. This could be useful if you want to look at a subset of
the data.

########################################

Input Arguments (specified after script name, separated by spaces; see "Usage" above):

1.) Sequencing depth cutoff. This program only considers genotypes that have a minimum
number of reads. Genotypes with less sequencing depth are treating as missing data. This
argument must be an integer.

2.) Required proportion of genotypes with data. This number represents the fraction of  
individuals required to have data for each population. Loci that have fewer genotypes (more 
missing data) than this threshold are ignored. For example, specifying a proportion of 0.8
indicates that an allele frequency will be calculated for loci in which each population
has at least 80% of its individuals represented with sequence data. This proportion must
be greater than 0 and less than or equal to 1 (no missing data allowed).

3.) VCF file name. You may specify just the name if the file is in the same directory as
the script. Otherwise, provide the full file path in the filename as well.

4.) Population Identification File name. You may specify just the name if the file is in
the same directory as the script. Otherwise, provide the full file path in the filename
as well.

########################################

How script works:

With these files, this program goes through the VCF file and calculates allele frequencies
for loci with sufficient data for each population. For loci in which a population has
*more* genotypes than the required proportion (see input arguments above), alleles are
randomly downsampled to this proportion. For example, if you specify the program to
calculate allele frequencies for loci in which at least 80% of individuals in a population
have high depth genotypes, but the population has 100% of its individuals represented, 
than 20% of its individuals are randomly discarded. This downsampling ensures that all 
loci have an equivalent amount of data and thus power to discover diversity.

If more than one population is used, than this script only calculates allele frequencies
if *all* populations have sufficient data. You may modify the "Population Identification
File" to consider fewer populations, if you wish.
 
########################################

Output Files:

Allele Frequency Spectrum (AFS) file: this file is tab-delimited, with each row
representing a different locus (specific position on a scaffold/contig). The first and
second column represent the scaffold and position, respectively, and the following columns
are the derived allele proportions for each population and range from 0 to 1.

QC file: This file shows some information about the populations, including the number of
individuals per population but also the number of individuals *after* downsampling. It
also shows the number of loci that had too much missing data (as specified by input
argument 2) for each population. This may be useful for identifying populations that may
be ignored from having significantly more missing data than other populations, since
allele frequencies are only reported in the AFS file for loci in which *all* populations
have at least the minimum amount of data.



