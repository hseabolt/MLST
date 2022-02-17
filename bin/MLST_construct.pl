#!/usr/bin/perl

# MLST_construct.pl
# A wrapper script to compute an ad hoc MLST scheme using multiple sequence alignments of reference loci.
# This component is part of a larger pipeline to construct an MLST scheme and write it to a file in pieces or its entirety.
# It also contains the option to write a SQLite database to disk to store the data for piecemeal querying later on.

# Each MLST is instantiated here as an MLST object, from HS custom class MLST.pm
# This class stores a list of Locus objects (HS class Locus.pm) which represent individual loci within the model,
# as well as the list of genomes that we have inserted into the model for analysis.  

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(first shuffle uniq);
use File::Basename; 
use Cwd qw(cwd);

# Essential CPAN modules
use Statistics::R;
use Set::Scalar;
use DBI;
use Parallel::ForkManager;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use MLST;

# Required input parameters
my $input_list;
my $output;
my $input_gff;
my $analyze_list;
my $dir;
my $v;
my $sqlite;
my $threads;

sub usage {
	my $usage = "MLST_construct.pl\n
	PURPOSE:	A Perl pipeline for constructing MLST schemes as Perl objects.
			Requires CPAN modules Statistics::R, DBI, and Set::Scalar.
                        \n
	USAGE: MLST_construct.pl [options]
          
	=== INPUT PARAMETERS ===
	-il    		input list of reference alignment files, one per line.  Files should be in FASTA format.
			(HS script core_genome_reads.pl is a good way to generate these alignments)
			These files are used to construct the model, they are NOT analyzed!
	-ig 		input GFF file containing the genomic coordinates of the reference loci (OPTIONAL)	
	-al 		list of multi-FASTA files to ANALYZE with the model, one per line.  Each file should contain sequences from only one genome.
			The names of the FASTA header for each sequence should match those of loci in the MLST model.
              
	=== OUTPUT PARAMETERS ===
	-out        	output file base name, not including extensions
	-dir		name of the input directory for storing files 
	
	=== OPTIONAL EXTRAS ===
	-v 		verbose
	-t 		INT; number of threads to use for parallel processing? [ Default: 1 ]
	-sqlite		INT flag; write the MLST scheme to a SQLite database? Requires CPAN module DBI.  [ Default: OFF ]
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'input_list|il=s' 	=> \$input_list,
            'input_gff|ig=s'	=> \$input_gff,
            'analyze_list|al=s' => \$analyze_list,
            'verbose|v' 		=> \$v,
            'dir|d=s' 			=> \$dir,
		  'sqlite=i' 			=> \$sqlite,
		  'output|o=s' 			=> \$output,
		  'threads|t=i' 		=> \$threads,
) or die usage(); 

# Parameter setups
$v = defined($v)? 1 : 0;
$sqlite = ( $sqlite && $sqlite == 1 )? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;
$dir = ( $dir )? $dir : cwd();
$output = ( $output )? $output : "OUT.mlst";

###################################################################################################
# Create directories for storing ALOT of files
system("mkdir $dir\\models $dir\\output\\fasta ");
open LOG, ">", "$dir\\output\\$output.log";

# Instantiate the MLST model and initialize a master configuration file that manages the MLST scheme.
my $MLST = MLST->new( 'name' => $output );

# Instantiate the Fork Manager object
my $manager = Parallel::ForkManager->new($threads);

open MASTER, ">", "$dir\\output\\$output.config";
my @colnames = qw( locusID gene_name n_haplotypes length model_file   );
print MASTER join("\t", @colnames), "\n";

# Grab the reference GFF if it's given
open GFF, "$input_gff" or warn "INPUT ERROR -- Your input GFF file doesn't seem to exist!\n" if ( $input_gff && -e $input_gff );
	my @gff = <GFF>;
close GFF;
###################################################################################################

###################################################################################################
# Process all of the reference genes given into the model.
# Open the reference genes file.  Process each alignment and construct a Locus object from the alignment and add it to the MLST master model.
open REF, "$input_list" or die "INPUT ERROR -- Your input list of alignment files doesn't seem to exist!\n";
print STDERR " === Constructing the MLST model from the input sequences === \n" if ( $v == 1 );

my @alis;
my %Haplotypes = ();
my $locusID = 0;
my $position = "";
while ( <REF> )	{

	# Sanity checks
	my $line = $_;
	chomp $line;
	if ( not -e "$line" )	{
		warn " --- ERROR: $line file apparently doesnt exist!  Skipping this file!\n";
		print LOG "Alignment file error: $line file doesnt exist.\n";
	}
	
	# Pre-process the alignment file -- it should have already been dereplicated by distance, all we want to do here is relabel the headers
	my $ali = fileparse($line, qr/\.[^.]*/);
			
	# Store the FASTA sequences in a hash
	$/ = ">";
	open ALI, "$line";
		my @fastas = <ALI>;
		my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
	close ALI;
	$/ = "\n";

	my $alignment = [ ];
	my $haplotype = 1;
	my @headers; my $seq;
	foreach my $record ( @fastas )	{
		my ($header, @seq) = split "\n", $record;
		$seq = join '', @seq;
		@headers = split " ", $header;
		
		$headers[0] =~ s/ /_/g;				# Convert any spaces to underscores
		$headers[0] =~ s/-//g;				# Convert any hyphens in SEQUENCE NAMES to underscores (Paup really doesnt like these)
		$headers[0] =~ s/\s+//g;			# Remove any other odd whitespace
		$headers[0] =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
		$seq =~ s/>//g;						# Remove any lingering ">" symbols
		$seq = uc $seq;						# Convert everything to uppercase 

		# Store the sequences as a matrix and record the haplotype IDs
		push @{$alignment}, $seq;
		$Haplotypes{$ali}->{$headers[0]} = $haplotype;
		$haplotype++;
	}
	
	# Get the chromosomal position from the GFF
	$position = "";
	my @lines = grep { /$ali/ } @gff;
	for( my $i=0; $i < scalar @lines; $i++ )	{
		my @coords = split("\t", $lines[$i]);
		$position .= "$coords[0]:$coords[3]\-$coords[4]" if ( $coords[6] eq "+" );
		$position .= "$coords[0]:complement($coords[3]\-$coords[4])" if ( $coords[6] eq "-" );
		$position .= "," if ( $i < scalar @lines - 1 );
	}
	
	# Create the locus object within the MLST model
	# NOTE: this currently assumes that the reference alignment given to this script as input has already been rigorously curated.  
	# I.e. we will not check for identical sequences at this stage.
	my $successful_add = $MLST->add_locus( $alignment, $ali, $position, 'DNA' );
		
	# Print the model and update the master config file if we successfully added the locus
	if ( $successful_add )	{
		print MASTER "$locusID\t$ali\t", $haplotype - 1, "\t", length $seq, "\t$dir\/models\/$ali.locus$locusID.mod\n";
		print LOG "Added new locus $ali (LocusID=$locusID) to MLST model.\n";
		$locusID++;
	}
	else	{
		print LOG "ERROR: Could not add locus $ali to MLST model.\n";
	}
}
$MLST->print_mlst_model( "$dir\/models\/original.MLST.mod", 1 );
close REF;



# We have alot of loci here, so cache the model with block sizes of 50 genes each (should have 4 blocks)
my $blocks = $MLST->cache_model( $dir, 50 );

exit;

#################################################################
# For each novel genome, traverse the model for each locus.  
# Update models if/when we encounter novel haplotypes at various loci.
# Assign a nomenclature to each genome and save this information off.

print STDOUT " === Analyzing the provided genomes === \n" if ( $v == 1 );

if ( $analyze_list && -e "$analyze_list" )	{
	open ANALYZE, "$analyze_list" or warn " --- ERROR: $analyze_list file unable to be opened.\n $!\n";
	while ( <ANALYZE> )		{
		chomp $_;
		
		# Get the genome's basename
		my $genome = fileparse($_, qr/\.[^.]*/);
		print STDOUT " --- $genome --- \n" if ( $v == 1 );
		
		# Run it through the MLST model to generate our sequence of loci IDs
		my $mlst = $MLST->assign_mlst( $_, $genome );
	}
	close ANALYZE;
}

#################################################################
# Finally, write some results information to a database for storage and export any output files that we want.
# The final collection of MLST loci should be:
# 	i)	 a list of unique FASTA files, one per locus, containing all the unique haplotype sequences within it.
# 	ii)	 a textual HMM model file and/or positional weight matrix describing the model used.
# 	iii)

print STDERR " === Printing the MLST model === \n" if ( $v == 1 ); 

$MLST->print_mlst_model( "$dir\/models\/$output.MLST.mod", 1 );
print LOG "Printed the entire MLST model to file: $dir\/models\/$output.MLST.mod \n";

close MASTER;
close LOG;
exit;
