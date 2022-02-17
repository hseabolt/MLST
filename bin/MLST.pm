#!/usr/bin/perl

# package MLST.pm
# A Perl package for MLST objects

# Author: MH Seabolt
# Last Updated: 8-11-2020 

package MLST; 

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(max min sum lg log10 pow round);			 #Import from other packages

use strict;
use warnings;
use List::Util qw( max min sum first shuffle );
use Scalar::Util;
use File::Basename;
use Storable qw(dclone);
use Data::Dumper;
use Carp;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use Locus; 

# ESSENTIAL CPAN module for Sets
use Set::Scalar;

################################################################################## 

# Class data and methods 
# Attributes 
{
	my %_attribute_properties = (
		_name 					=> '',					# the name of the model, e.g. the species/OTU that the model represents
		_loci					=> [ ],					# a list of loci in the MLST model, these should all be Locus objects
		_genomes					=> [ ],					# a list of the genomes in the model, in the same order that they appear in _alignment
		_alignment				=> [ ],					# an anon hash where the keys are the IDs of an input genome and the values are the haplotype vectors
		_profile					=> [ ],					# 2D matrix representing a positional weight matrix profile of the MLST alignment of haplotype vectors
		_alphabet					=> '',					# the set of possible haplotypes currently represented in the model
	);
	
	# Global variable counter
	my $_count = 0;
	
	# Return a list of all attributes
	sub _all_attributes	{
		keys %_attribute_properties;
	}
	
	# Return the default value for a given attribute
    	sub _attribute_default 	{
     	my( $self, $attribute ) = @_;
        	$_attribute_properties{$attribute};
    	}
    
	# Manage the count of existing objects
	sub get_count	{
		$_count;
	}
	sub _incr_count	{
		++$_count;
	}
	sub _decr_count	{
		--$_count;
	}	
}


############################################################
#                       CONSTRUCTORS                       #
############################################################

# The contructor method
# Construct a new graph (my $node = Markov->new() );
# Returns a scalar reference to a
sub new				{
	my ( $class, %arg ) = @_;
	
	# Create the new object
	my $self = bless {}, $class;

	foreach my $attribute ( $self->_all_attributes() ) {
		# E.g. attribute = "_name",  argument = "name"
		my ($argument) = ( $attribute =~ /^_(.*)/ );
		# If explicitly given
		if (exists $arg{$argument}) 	{
			$self->{$attribute} = $arg{$argument};
		}
		else	{
			$self->{$attribute} = $self->_attribute_default($attribute);
		}
   	}
	
	# Construction code specific to this class
	######################################################################################
   	 
   	my $alphabet = Set::Scalar->new(); 
   	$alphabet->insert( "-" );
   	$self->{_alphabet} = $alphabet;
   	
   	
	# We arent adding any extra constructor code for this class beyond this for the alphabet 
	#   -- these are meant to be instantiated and then built iteratively by the driver code automatically.   
	   
	######################################################################################
	
	
     $class->_incr_count();
	return $self;
}

# The clone method
# All attributes are copied from the calling object, unless specifically overriden
# Called from an existing object ( Syntax: $cloned_obj = $obj->clone(); )
sub clone	{
	my ( $caller, %arg ) = @_;
	# Extract the class name from the calling object
	my $class = ref($caller);
		
	# Create a new object
	my $self = bless {}, $class;
		
	foreach my $attribute ( $self->_all_attributes() )	{
		my ($argument) = ( $attribute =~ /^_(.*)/ );
			
		# If explicitly given
		if ( exists $arg{$argument} )	{
			$self->{$attribute} = $arg{$argument};
		}
			
		# Otherwise, copy attribute of new object from the calling object
		else	{
			$self->{$attribute} = $caller->{$attribute};
		}
	}
	$self->_incr_count();
	return $self;
}

# When an object is no longer being used, garbage collect it and adjust count of existing objects
sub DESTROY	{
	my ( $self ) = @_;
	$self->_decr_count();
}

############################################################
#         ORGANIZATION HASHES AND SUBROUTINES              #
############################################################

# This hash is associated with the _stateToIndex() subroutine, but it must be outside the subroutine scope to be maintainable and accessible.
my $h = 0;			# The initial index in the prebuilt hash of states below, we will increment it if needed
my %Index = ();
my %ReverseIndex = ();		# Associated with the _indexToState() subroutine and automatically updated by _stateToIndex()
my %Names = ();			
my %NamesIndex = ();
my %NamesState = ();

my $g = 0;
my %GenomesIndex = ();
my %ReverseGenomesIndex = ();

# Converts key current character into index
sub _locusToIndex	{
	my ( $self, $locus_ref )	= @_;
	return if ( not $locus_ref );		# Sanity check
	my $index;
	
	if ( exists $Index{$locus_ref}  ) 	{
		$index = $Index{$locus_ref};
	}
	# If the index doesnt exist in the Index hash, then add it at the end and increment the value
	# This should be fine for multiple letter k-mer style indices and some symbols
	else		{
		
		$Index{$locus_ref} = $h;				# Be careful here, as odd symbols may cause errors
		$ReverseIndex{$h} = $locus_ref;
		my $name = ${$locus_ref}->{_name};
		$NamesIndex{$h} = $name;
		$Names{$name} = $h;
		$NamesState{$name} = $locus_ref;
		$index = $Index{$locus_ref};
		$h++;
	}

	return $index;
}

# Converts key current index into character
sub _indexToLocus	{
	my ( $self, $index )	= @_;
	my $locus;
	if ( exists( $ReverseIndex{$index} ) ) 	{
		$locus = $ReverseIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	
	return $locus;
}

# Converts from an index to a Locus NAME
sub _indexToName	{
	my ( $self, $index ) = @_;
	my $name;
	if ( exists $NamesIndex{$index} )	{
		$name = $NamesIndex{$index};
	}
	return $name;
}

# Converts from the name of a Locus to the index
sub _nameToIndex	{
	my ( $self, $name ) = @_;
	my $index;
	if ( exists $Names{$name} )	{
		$index = $Names{$name};
	}
	return $index;
}

sub _nameToLocus	{
	my ( $self, $name ) = @_;
	my $locus;
	if ( exists $Names{$name} )	{
		$locus = $ReverseIndex{ $Names{$name} };
	}
	return $locus;
}

# This hash is associated with the _charToIndex() subroutine, but it must be outside the subroutine scope to be maintainable and accessible.
my $p = 95;			# The last index in the prebuilt hash of chars below, we will increment it if needed
my %charIndex = (
	# Lowercase
	'a' => 0, 'b' => 1, 'c' => 2, 'd' => 3, 'e' => 4,
	'f' => 5, 'g' => 6, 'h' => 7, 'i' => 8, 'j' => 9,
	'k' => 10, 'l' => 11, 'm' => 12, 'n' => 13, 'o' => 14,
	'p' => 15, 'q' => 16, 'r' => 17, 's' => 18, 't' => 19,
	'u' => 20, 'v' => 21, 'w' => 22, 'x' => 23, 'y' => 24,
	'z' => 25,	
	# Uppercase
	'A' => 26, 'B' => 27, 'C' => 28, 'D' => 29, 'E' => 30,
	'F' => 31, 'G' => 32, 'H' => 33, 'I' => 34, 'J' => 35,
	'K' => 36, 'L' => 37, 'M' => 38, 'N' => 39, 'O' => 40,
	'P' => 41, 'Q' => 42, 'R' => 43, 'S' => 44, 'T' => 45,
	'U' => 46, 'V' => 47, 'W' => 48, 'X' => 49, 'Y' => 50,
	'Z' => 51,
	# Numbers
	'0' => 52, '1' => 53, '2' => 54, '3' => 55, '4' => 56,
	'5' => 57, '6' => 58, '7' => 59, '8' => 60, '9' => 61,
	# Other ASCII characters (this is not an exhaustive list, only includes a few symbols)
	           '_' => 62, '!' => 63, '.' => 64, '"' => 65,			# '-' => 62, Removed this value and adjusted subsequent ones due to this being used as the gap character
	'#' => 66, '$' => 67, '%' => 68, '&' => 69, "\'" => 70,
	'(' => 71, ')' => 72, '!' => 73, '.' => 74, '"' => 75,
	'*' => 76, '+' => 77, ',' => 78, '/' => 79, ':' => 80,
	';' => 81, '<' => 82, '>' => 83, '=' => 84, '?' => 85,
	'@' => 86, '[' => 87, ']' => 88, "\\" => 89, '^' => 90,
	'`' => 91, '{' => 92, '}' => 93, '|' => 94, '~' => 95,
);

# This hash is associated with the _indexToChar() subroutine, but it must be outside the subroutine scope to be maintainable and accessible.
# It is automatically updated by _charToIndex()
my %reverseCharIndex = (
	# Lowercase
	'0' => 'a', '1' => 'b', '2' => 'c', '3' => 'd', '4' => 'e',
	'5' => 'f', '6' => 'g', '7' => 'h', '8' => 'i', '9' => 'j',
	'10' => 'k', '11' => 'l', '12' => 'm', '13' => 'n', '14' => 'o',
	'15' => 'p', '16' => 'q', '17' => 'r', '18' => 's', '19' => 't',
	'20' => 'u', '21' => 'v', '22' => 'w', '23' => 'x', '24' => 'y',
	'25' => 'z', 	
	# Uppercase
	'26' => 'A', '27' => 'B', '28' => 'C', '29' => 'D', '30' => 'E',
	'31' => 'F', '32' => 'G', '33' => 'H', '34' => 'I', '35' => 'J',
	'36' => 'K', '37' => 'L', '38' => 'M', '39' => 'N', '40' => 'O',
	'41' => 'P', '42' => 'Q', '43' => 'R', '44' => 'S', '45' => 'T',
	'46' => 'U', '47' => 'V', '48' => 'W', '49' => 'X', '50' => 'Y',
	'51' => 'Z',
	# Numbers
	'52' => '0', '53' => '1', '54' => '2', '55' => '3', '56' => '4',
	'57' => '5', '58' => '6', '59' => '7', '60' => '8', '61' => '9',
	# Other ASCII characters (this is not an exhaustive list, only includes a few symbols)
			   '62' => '_', '63' => '!', '64' => '.', '65' => '"',				# '62' => '-', 
	'66' => '#', '67' => '$', '68' => '%', '69' => '&', '70' => "\'",
	'71' => '(', '72' => ')', '73' => '!', '74' => '.', '75' => '"',
	'76' => '*', '77' => '+', '78' => ',', '79' => '/', '80' => ':',
	'81' => ';', '82' => '<', '83' => '>', '84' => '=', '85' => '?',
	'86' => '@', '87' => '[', '88' => ']', '89' => "\\", '90' => '^',
	'91' => '`', '92' => '{', '93' => '}', '94' => '|', '95' => '~',
);

# Converts key current character into index
sub _charToIndex	{
	my ( $self, $char )	= @_;
	my $index;
	if ( exists( $charIndex{$char} ) ) 	{
		$index = $charIndex{$char};
	}
	# If the index doesnt exist in the Index hash, then add it at the end and increment the value
	# This should be fine for multiple letter k-mer style indices and some symbols
	else		{
		warn " === CHAR $char HAS NO EXISTING INDEX, ADDING IT !\n\n";	
		$p++;
		$charIndex{$char} = $p;				# Be careful here, as odd symbols may cause errors
		$reverseCharIndex{$p} = $char;
		$index = $charIndex{$char};
	}
	return $index;
}

# Converts key current index into character
sub _indexToChar	{
	my ( $self, $index )	= @_;
	my $char;
	if ( exists( $reverseCharIndex{$index} ) ) 	{
		$char = $reverseCharIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	else		{
		warn " === INDEX $index HAS NO EXISTING CHAR!\n\n";	
		$char = ' ';	# Set char to a space just as an empty placeholder.
	}
	return $char;
}

# Converts a genome's name to an index position in the _genomes and _alignment lists
sub _genomeToIndex	{
	my ( $self, $name ) = @_;
	my $index;
	if ( exists $GenomesIndex{$name} )	{
		$index = $GenomesIndex{$name};
	}
	else	{
		$GenomesIndex{$name} = $g;
		$ReverseGenomesIndex{$g} = $name;
		$g++;
		$index = $GenomesIndex{$name};
	}
	return $index;
}

# Returns a genome's name, given it's index in the _genomes list OR _alignment list
sub _indexToGenome	{
	my ( $self, $index ) = @_;
	my $name;
	if ( exists( $ReverseGenomesIndex{$index} ) ) 	{
		$name = $ReverseGenomesIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	else		{
		warn " === INDEX $index HAS NO EXISTING GENOME!\n\n";	
	}
	return $name;
}

############################################################
#            LOCUS OBJECT HANDLING SUBROUTINES             #
############################################################

# Add a new locus to the MLST model.
# The constructor for the Locus object will check that the provided data is up to snuff
sub add_locus	{
	my ( $self, $data, $new_locus, $position, $type ) = @_;
	
	# Sanity check
	return if ( scalar @{$data} == 0 );
	
	# Initialize a new Locus object
	$new_locus = ( $new_locus )? $new_locus : $h;		# ALL states must have a name! Use $h as a placeholder if no name is supplied.
	my $node = Locus->new( "name" => "$new_locus", "DNA" => $data, "position" => $position, "type" => $type );
	
	# Update the requisite hashes and lists with the new node
	# If this node already exists in the model, then here we will just update it with new information
	# WARNING: be careful updating nodes this way to prevent overwriting data unintentionally!
	my $index = $self->_locusToIndex( \$node );
	return if ( not defined $index );
	$self->{_loci}->[$index] = $node;
	
	# Update the MSLT alphabet if needed
	my $nhaplotypes = $node->{_t};
	if ( $nhaplotypes > $self->{_alphabet}->size )	{
		for ( my $i=0; $i < $nhaplotypes; $i++ )	{
			$self->{_alphabet}->insert( $reverseCharIndex{$i} );
		}
	}
	
	# Lastly, update the profile
	$self->{_profile} = $self->update_profile;
	
	return 1; 
}

# Delete a locus from the MLST model
sub delete_locus	{
	my ( $self, $name ) = @_;
	my $kill_state_ref = $self->_nameToLocus( $name );
	
	# Update _loci
	my $index = $self->_locusToIndex( $kill_state_ref );
	splice( @{ $self->{_loci} }, $index, 1 );
	${$kill_state_ref}->DESTROY;
	
	# Update the alignment, if it exists
	foreach ( @{$self->{_alignment}} )	{
		my @sequence = split('', $self->{_alignment}->[$_]);
		splice( @sequence, $index, 1);
		my $sequence = join('', @sequence);
		$self->{_alignment}->[$_] = $sequence;
	}
	
	# Lastly, update the profile
	$self->{_profile} = $self->update_profile;
	
	# Update the organizational hashes --> requires resetting and rehashing them all
	$h = 0;		
	%Index = ();
	%ReverseIndex = ();		
	%Names = ();			
	%NamesIndex = ();
	%NamesState = ();
	for ( my $i=0; $i < scalar @{$self->{_loci}}; $i++ )	{
		$self->_locusToIndex( \$self->{_loci}->[$i] );
	}
}

# A generic function to build the model in its entirety, given a list of input alignment FASTA files.
# You MUST instantiate an empty MLST object first to call this function!
sub build_mod		{
	my ( $self, $input_list, $input_gff ) = @_;
	
	open MASTER, ">", "MLST.build_mod.config";
	my @colnames = qw( locusID gene_name n_haplotypes length model_file   );
	print MASTER join("\t", @colnames), "\n";

	# Open and read the reference genes file.  Process each alignment and construct a Locus object from the alignment and add it to the MLST master model.
	open REF, "$input_list" or die "MLST::build_mod INPUT ERROR -- Your input list of alignment files doesn't seem to exist!\n";
		my @input = <REF>;
	close REF;

	# Grab the reference GFF if it's given
	my @gff;
	open GFF, "$input_gff" or die "INPUT ERROR -- Your input GFF file doesn't seem to exist!\n" if ( -e "$input_gff" );
		@gff = <GFF>;
	close GFF;

	my @alis;
	my $locusID = 0;
	my %Haplotypes = ();
	foreach my $line ( @input )	{
		chomp $line;
	
		# Sanity checks
		if ( not -e "$line" )	{
			warn "MLST::build_mod ERROR: $line file apparently doesnt exist!  Skipping this file!\n";
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
		my $position = "";
		my @lines = grep { /$ali/ } @gff;
		for( my $i=0; $i < scalar @lines; $i++ )	{
			my @coords = split("\t", $lines[$i]);
			$position .= "$coords[0]:$coords[3]\-$coords[4]" if ( $coords[6] eq "+" );
			$position .= "$coords[0]:complement($coords[3]\-$coords[4])" if ( $coords[6] eq "-" );
			$position .= "," if ( $i < scalar @lines - 1 );
		}
			
		# Create the locus object within the MLST model
		$self->add_locus( $alignment, $ali, $position );
	
		# Print the model and update the master config file
		print MASTER "$locusID\t$ali\t", $haplotype - 1, "\t", length $seq, "\tlocus$locusID.MLST.mod\n";
		$locusID++;
	}
	return $self;
}

############################################################
#               SEQUENCE HANDLING SUBROUTINES              #
############################################################

# Add a new sequence to this locus, including all requisite updates 
sub add_sequence	{
	my ( $self, $new_sequence, $name ) = @_;
	
	# Make sure that we have a name
	$name = ( $name )? $name : ("Genome" . scalar @{$self->{_genomes}} );
	
	# Sanity check -- is the new sequence the same length as the other sequences?
	if ( length $new_sequence == length $self->{_alignment}->[0] )		{
		;  # If they are the same, all good.  Carry on.
	}
	else	{
		return "MLST::add_sequence ERROR --- cannot add the requested sequence as it is not the same length as the other MLST sequences!\n";
	}
	
	# Check if we need to update the alphabet with new haplotypes found in this new sequence
	# Be very careful with this block, since unexpected symbols will totally corrupt our alphabet
	my @haplotypes = split('', $new_sequence );
	foreach my $haplotype ( @haplotypes )	{
		$self->{_alphabet}->insert( $haplotype );
	}
	
	# Update the organizational hash and push the new sequence to the _alignment and _genomes attributes
	my $gindex = $self->_genomeToIndex( $name );
	$self->{_genomes}->[$gindex] = $name;
	$self->{_alignment}->[$gindex] = $new_sequence;
	
	# Lastly, update the profile
	$self->{_profile} = $self->update_profile;
}

# Delete an MLST sequence from the MLST model, given the genome name
sub delete_sequence		{
	my ( $self, $kill_genome ) = @_;
	my $kill_index = $self->_genomeToIndex( $kill_genome );
	
	# Splice out the sequence from the _DNA attribute and update _t
	splice( @{$self->{_genomes}}, $kill_index, 1 );
	splice( @{$self->{_alignment}}, $kill_index, 1 );
	
	# Update the profile
	$self->{_profile} = $self->update_profile;
	
	# Update the organizational hashes --> requires resetting and rehashing them all
	$g = 0;			
	%GenomesIndex = ();
	%ReverseGenomesIndex = ();		
	for ( my $i=0; $i < scalar @{$self->{_genomes}}; $i++ )	{
		my $genome = $self->{_genomes}->[$i];
		$self->_genomeToIndex( $genome );
	}
}

# Updates/replaces a sequence in place within the _alignment attribute 
sub update_sequence	{
	my ( $self, $genome, $replacement_sequence ) = @_;
	my $index = $self->_genomeToIndex( $genome );
	$self->{_alignment}->[$index] = $replacement_sequence;
	
	# Update the alphabet if needed for the new sequence
	my @haplotypes = split('', $replacement_sequence );
	foreach my $haplotype ( @haplotypes )	{
		$self->{_alphabet}->insert( $haplotype );
	}
}

# Checks if a given sequence is present in the model already or if it is a novel sequence
sub check_sequence	{
	my ( $self, $sequence ) = @_;
	
	print "MLST::check_sequence --- Length of given sequence: ", length $sequence, " | Length of model: ", length $self->{_alignment}->[0], "\n";
	
	if ( length $sequence != length $self->{_alignment}->[0] )	{	
		print STDERR "MLST::check_sequence ERROR -- the given sequence is not the same length as the others in the model.\n";
		return;
	}
	
	if ( grep { $sequence eq $_ } @{$self->{_alignment}} ) 	{
		return "Sequence is in the model ";
	}	
	else	{
		return "Novel sequence";
	}
}

### Below are aliases for the above sequence handling subroutines, if the user prefers to call them by genome instead
# Alias for add_sequence subroutine
sub add_genome		{
	my ( $self, $name, $new_sequence ) = @_;
	
	# Make sure that we have a name
	$name = ( $name )? $name : ("Genome" . scalar @{$self->{_genomes}} );
	
	# Sanity check -- is the new sequence the same length as the other sequences?
	if ( scalar @{$self->{_alignment}} > 0 && (length $new_sequence == length $self->{_alignment}->[0]) )		{
		;  # If they are the same, all good.  Carry on.
	}
	elsif ( scalar @{$self->{_alignment}} == 0 )	{
		;  # If this is the first genome to be added, carry on.
	}
	else	{
		print STDERR "MLST::add_genome ERROR --- cannot add the requested genome sequence as it is not the same length as the other MLST sequences!\n";
		return;
	}
	
	# Update the alphabet if needed for the new sequence
	my @haplotypes = split('', $new_sequence );
	foreach my $haplotype ( @haplotypes )	{
		$self->{_alphabet}->insert( $haplotype );
	}
	
	# Update the organizational hash and push the new sequence to the _alignment and _genomes attributes
	my $gindex = $self->_genomeToIndex( $name );
	$self->{_genomes}->[$gindex] = $name;
	$self->{_alignment}->[$gindex] = $new_sequence;
	
	# Lastly, update the profile
	$self->{_profile} = $self->update_profile;
}

sub delete_genome		{
	my ( $self, $kill_genome ) = @_;
	my $kill_index = $self->_genomeToIndex( $kill_genome );
	
	# Splice out the sequence from the _DNA attribute and update _t
	splice( @{$self->{_genomes}}, $kill_index, 1 );
	splice( @{$self->{_alignment}}, $kill_index, 1 );
	
	# Update the profile
	$self->{_profile} = $self->update_profile;
	
	# Update the organizational hashes --> requires resetting and rehashing them all
	$g = 0;			
	%GenomesIndex = ();
	%ReverseGenomesIndex = ();		
	for ( my $i=0; $i < scalar @{$self->{_genomes}}; $i++ )	{
		my $genome = $self->{_genomes}->[$i];
		$self->_genomeToIndex( $genome );
	}
}

sub update_genome	{
	my ( $self, $genome, $replacement_sequence ) = @_;
	my $index = $self->_genomeToIndex( $genome );
	$self->{_alignment}->[$index] = $replacement_sequence;
	
	# Update the alphabet if needed for the new sequence
	my @haplotypes = split('', $replacement_sequence );
	foreach my $haplotype ( @haplotypes )	{
		$self->{_alphabet}->insert( $haplotype );
	}
	
	# Lastly, update the profile
	$self->{_profile} = $self->update_profile;
}

sub check_genome	{
	my ( $self, $sequence ) = @_;
	if ( length $sequence != length $self->{_alignment}->[0] )		{
		print STDERR "MLST::check_genome ERROR -- the given sequence is not the same length as the others in the model.\n";
		return;
	}
	if ( grep { $sequence eq $_ } @{$self->{_genomes}} ) 	{
		return "Genome is in the model ";
	}	
	else	{
		return "Novel genome name";
	}
}

############################################################
#                  TRAVERSE THE MODEL!                     #
############################################################

# Accepts a Hash reference to a FASTA file of new sequences to search in the model, plus the name of the genome
# (keys are locus names which MATCH THOSE IN THE MODEL, and values are the sequences themselves).
# The most difficult technical part here will be to make sure that the names in the hash match the ones in the model
# Generates a vector of haplotype IDs encoded by ASCII characters (not the ASCII values, but the same concept) and returns the vector.
sub assign_mlst	{
	my ( $self, $new_sequences_file, $genome ) = @_;
	
	# Initialize the MLST alignment as all gaps
	my $mlst = [ ];
	for ( my $i=0; $i < scalar @{$self->{_loci}}; $i++ ) 	{ 
		$mlst->[$i] = "-";
	}

	# Store the FASTA sequences in a hash
	$/ = ">";
	open FASTA, "$new_sequences_file" or warn;	
	while ( <FASTA> )	{
		chomp $_;
		my ($header, @seq) = split "\n", $_;
		next if ( not $header );
		my $seq = join '', @seq;
		my @headers = split " ", $header;
		
		$headers[0] =~ s/ /_/g;				# Convert any spaces to underscores
		$headers[0] =~ s/-//g;				# Convert any hyphens in SEQUENCE NAMES to underscores (Paup really doesnt like these)
		$headers[0] =~ s/\s+//g;			# Remove any other odd whitespace
		$headers[0] =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
		$seq =~ s/>//g;						# Remove any lingering ">" symbols
		$seq = uc $seq;						# Convert everything to uppercase 
		
		# Analyze the sequences in the file:
		# Get the index of the Locus model in our MLST object
		my $model_index = $self->_nameToIndex( $headers[0] );
		if ( not defined $model_index )		{
			next;		# Skip this if we don't have any Locus model for this sequence
		}
		# If we do locate the correct Locus, assign the haplotype for the new sequence
		else	{
			print STDERR "\t~~~>> MLST::assign_mlst::$self->{_loci}->[$model_index]->{_name}";
			my $haplotype = $self->{_loci}->[$model_index]->add_sequence( $seq );
			$mlst->[$model_index] = $self->{_loci}->[$model_index]->_indexToChar( $haplotype );
		}
	} 	
	close FASTA;
	$/ = "\n";
	
	# Update the _genome and _alignment lists
	$mlst = join('', @{$mlst});
	$genome = ( $genome )? $genome : "Genome" . scalar @{$self->{_genomes}};
	$self->add_genome( $genome, $mlst );
	
	return $mlst;
}

# Alias for the above code
sub analyze	{
	my ( $self, $new_sequences, $genome ) = @_;
	my $mlst = $self->assign_mlst( $new_sequences, $genome );
	return $mlst;
}

############################################################
#                  PROFILING SUBROUTINES                   #
############################################################

# Generates a normalized profile from an aligned sequence matrix
sub create_profile 		{
	my ( $self, $matrix ) = @_;
	my %Profile = ();
	my @alphabet = $self->{_alphabet}->members;
	
	# Initialize a hash where keys are each possible symbol in our alphabet,
	# and values are an anon array, where each index represents a position in the sampled motif
	foreach my $haplotype ( @alphabet )	{
		for ( my $i=0; $i < @{$self->{_loci}}; $i++ )	{
			push @{$Profile{$haplotype}}, 0;
		}
	}
	
	# Count the occurences of each symbol
	for ( my $i=0; $i < @{$self->{_alignment}}; $i++ )	{
		my @sequence = split('', $matrix->[$i]);		
		for ( my $j=0; $j < $self->{_loci}; $j++ )	{
			$Profile{ $sequence[$j] }->[$j] += 1;
		}
	}
	
	# Normalize the counts
	foreach my $haplotype ( @alphabet )	{
		for ( my $j=0; $j < $self->{_loci}; $j++ )	{
			$Profile{$haplotype}->[$j] = sprintf("%3.3f", $Profile{$haplotype}->[$j] / scalar @{$self->{_alignment}});
		}
	}
	return \%Profile;	
}

# Updates the profile from the _DNA attribute.  
# The code here is the same as create_profile, just doesn't require passing an argument.
sub update_profile 		{
	my ( $self ) = @_;
	my %Profile = ();
	my @alphabet = $self->{_alphabet}->members;	
	my $matrix = $self->{_alignment};
	
	# Initialize a hash where keys are each possible symbol in our alphabet,
	# and values are an anon array, where each index represents a position in the sampled motif
	foreach my $haplotype ( @alphabet )	{
		for ( my $i=0; $i < scalar @{$self->{_loci}}; $i++ )	{
			push @{$Profile{$haplotype}}, 0;
		}
	}
	
	# Count the occurences of each symbol
	for ( my $i=0; $i < scalar @{$self->{_genomes}}; $i++ )	{
		my @sequence = split('', $matrix->[$i]);
		for ( my $j=0; $j < scalar @{$self->{_loci}}; $j++ )	{
			next if ( not defined $sequence[$j] );		# be careful here with this line...
			$Profile{ $sequence[$j] }->[$j] += 1;
		}
	}
	
	# Normalize the counts
	foreach my $haplotype ( @alphabet )	{
		for ( my $j=0; $j < scalar @{$self->{_loci}}; $j++ )	{
			$Profile{$haplotype}->[$j] = ( scalar @{$self->{_genomes}} > 0  )? ( sprintf("%3.3f", $Profile{$haplotype}->[$j] / scalar @{$self->{_genomes}}) ) : "0.000";
		}
	}
	$self->{_profile} = \%Profile;	
}

# Reads a profile and returns the consensus sequence
sub consensus_sequence 	{
	my ( $self, $profile ) = @_;
	my @alphabet = $self->{_alphabet}->members;
	my @consensus;
	
	# Get the most probable letter at each position in the profile
	for ( my $i=0; $i < scalar @{$self->{_loci}}; $i++ )	{
		my @candidates;
		foreach my $haplotype ( @alphabet )	{
			my @tmp = ( $haplotype, $profile->{$haplotype}->[$i] );
			push @candidates, \@tmp;
		}
		push @consensus, $self->max_score( \@candidates )->[0];
	}
	# Join the consensus array into a string and return
	return join("", @consensus);
}

# Calculate the score of a given profile -- the higher the score, the more homogenous the sequences are that make up the profile
sub get_profile_score		{
	my ( $self, $profile ) = @_;
	my @list;
	
	# The score of the profile is the sum of the most frequent letter in every position in the profile's indices
	# For each position in the motif, get the column from the profile.
	# Then take the maximum value of the column and sum it with all the other colMaxes.
	for (my $i=0; $i < scalar @{$self->{_loci}}; $i++ )	{
		my @col;
		foreach my $freq ( values %{$profile} )	{
			push @col, $freq->[$i];
		}
		push @list, max @col;
	}
	return sum @list;
}

# Returns the maximum valued tuple from a list of tuples, where the value to maximize is the second element in the tuple.		
sub max_score	{
	my ( $self, $tuples ) = @_;
	my $max = -10**10;			# Approximate a very low number as the initial max
	my $tup;
	
	for ( my $i=0; $i < scalar @{$tuples}; $i++ )	{
		if ( $tuples->[$i]->[1] >= $max )	{
			$max = $tuples->[$i]->[1];
			$tup = $tuples->[$i];
		}
		else	{
			next;
		}
	}
	return $tup;
}

############################################################
#                      I/O SUBROUTINES                     #
############################################################

# Parse an MLST model, including the loci, and regenerate the entire model from a text file
# Returns a populated MLST object from a given empty model instance.
sub parse_mlst_model	{
	my ( $self, $model ) = @_;
	my $flag = 0;
	
	# Read the model with the record separator since we can have multiple models in the same file
	$/ = "\/\/";
	open MODEL, "$model";
		my @input = <MODEL>;
	close MODEL;
	$/ = "\n";
	
	# Here we only need to parse the basics of the model (name, sequence data, and position) and the instantiation will regenerate the rest automatically.
	foreach my $model ( @input )	{
		my @model = split("\n", $model);
		my $name; my $data; my $position; my $type;
		
		foreach my $line ( @model )	{
			chomp $line;
			
			########## PARSE THE MLST RECORD #############
			if ( $line =~ /^MLST NAME/ )		{
				my @tmp = split(": ", $line);
				$self->{_name} = $tmp[1];
			}
			
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			elsif ( $line =~ /^GENOME/ )	{
				my @tmp = split(": ", $line);
				push @{ $self->{_genomes}}, $tmp[1];
			}
			
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			elsif ( $line =~ /^MLST SEQUENCE/ )	{
				my @tmp = split(":\t", $line);
				push @{$self->{_alignment}}, $tmp[1];
			}
			
			########## PARSE THE LOCI RECORDS #############
			# Get the model name
			elsif ( $line =~ /^MODEL NAME/ )	{
				my @tmp = split(": ", $line);
				$name = $tmp[1];
				$flag = 1;
			}
			# Get the position, if it's documented
			if ( $line =~ /^POSITION/ )		{
				my @tmp = split(": ", $line);
				$position = $tmp[1];
			}
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			if ( $line =~ /^SEQUENCE/ )	{
				my @tmp = split(":\t", $line);
				push @{$data}, $tmp[1];
			}
			# Get the type of model
			if ( $line =~ /^TYPE/ )	{
				my @tmp = split(": ", $line);
				$type = $tmp[1];
			}
		}
		
		# Instantiate the Locus object and add it to the list, if we have the flag is active
		$self->add_locus( $data, $name, $position, $type ) if ( $flag == 1 );
		$flag = 0;
	}
	
	# Generate the MLST profile and return 
	$self->update_profile();
	return $self;
}


# Parse an existing Locus model file a regenerate Locus objects from the file --> SKIPS ANY MLST MODELS PRINTED IN THE SAME FILE
sub parse_locus_model	{
	my ( $self, $model_file ) = @_;

	# Read the model with the record separator since we can have multiple models in the same file
	$/ = "\/\/";
	open MODEL, "$model_file";
		my @input = <MODEL>;
	close MODEL;
	$/ = "\n";

	# Here we only need to parse the basics of the model (name, sequence data, and position) and the instantiation will regenerate the rest automatically.
	foreach my $model ( @input )	{
		my @model = split("\n", $model);
		next if ( scalar @model == 0 );
		my $name; my $data; my $position; my $type;
		
		foreach my $line ( @model )	{
			chomp $line;
			
			# We might have the complete MLST model in any given record, so if we hit that, then break this inner loop and go to the next record.
			last if ( $line =~ /^MLST NAME/ );
			
			# Get the model name
			if ( $line =~ /^MODEL NAME/ )	{
				my @tmp = split(": ", $line);
				$name = $tmp[1];
			}
			# Get the position, if it's documented
			if ( $line =~ /^POSITION/ )		{
				my @tmp = split(": ", $line);
				$position = $tmp[1];
			}
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			if ( $line =~ /^SEQUENCE/ )	{
				my @tmp = split(":\t", $line);
				push @{$data}, $tmp[1];
			}
			# Get the type of model
			if ( $line =~ /^TYPE/ )	{
				my @tmp = split(": ", $line);
				$type = $tmp[1];
			}
		}
		
		# Instantiate the Locus object and add it to the list, if we have the flag is active
		$self->add_locus( $data, $name, $position, $type );
	}
	return $self;
}

# Read a FASTA file and return a hash reference with the sequence headers as keys and the sequences themselves as values.
# WARNING: this subroutine does include some HS-typical processing of the headers and sequences for formatting, be aware of this.
sub hash_ref_from_fasta_file	{
	my ( $self, $file ) = @_;
	my %Fasta;
	
	# Store the FASTA sequences in a hash
	$/ = ">";
	open FASTA, "$file";
		my @fastas = <FASTA>;
		my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
	close FASTA;
	$/ = "\n";

	foreach my $record ( @fastas )	{
		my ($header, @seq) = split "\n", $record;
		my $seq = join '', @seq;
		my @headers = split " ", $header;
		
		$headers[0] =~ s/ /_/g;				# Convert any spaces to underscores
		$headers[0] =~ s/-//g;				# Convert any hyphens in SEQUENCE NAMES to underscores (Paup really doesnt like these)
		$headers[0] =~ s/\s+//g;			# Remove any other odd whitespace
		$headers[0] =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
		$seq =~ s/>//g;						# Remove any lingering ">" symbols
		$seq = uc $seq;						# Convert everything to uppercase 
		
		# Store the sequences as a matrix and record the haplotype IDs
		$Fasta{$headers[0]} = $seq;
	}
	return \%Fasta;
}


# Print the MLST alignment of haplotype vectors as a FASTA format
sub print_mlst_fasta		{
	my ( $self, $filename, $mode ) = @_;
	my $seqs = $self->{_alignment};
	$mode = ( $mode && $mode eq ">>" )? ">>" : ">";
	open OUT, "$mode", "$filename";
	for ( my $i=0; $i < scalar @{$seqs}; $i++ )	{
		my $header = $self->{_genomes}->[$i];
		print OUT ">$header\n$seqs->[$i]\n";
	}	
	close OUT;
}

# Print out the profile in a nice format
sub print_mlst_profile	{
	my ( $self, $filename, $mode ) = @_;
	my $profile = $self->{_profile};
	$mode = ( $mode && $mode eq ">>" )? ">>" : ">";
	my @alphabet = $self->{_alphabet}->members;
	
	open OUT, "$mode", "$filename";
	foreach my $haplotype ( @alphabet )	{
		print OUT "$haplotype";
		for ( my $j=0; $j < scalar @{$self->{_loci}}; $j++ )	{
			print OUT " $profile->{$haplotype}->[$j]";
		}
		print OUT "\n";
	}
	close OUT;
}

# Print the entire MLST model into a text file --> this is the same as printing each locus in the model into the same text file
sub print_mlst_model	{
	my ( $self, $output, $split ) = @_;
	$split = ( $split && $split != 0 )? 1 : 0;
	
	# Get our alphabet, which is the highest number of sequences in any given Locus
	my @alphabet = $self->{_alphabet}->members;
	
	# Print the basic model metadata
	open MODEL, ">", "$output";
		my $date = 	`date`; chomp $date;
	print MODEL "MLST NAME: $self->{_name}\nDATE WRITTEN: $date\nN_SEQUENCES: ", scalar @{$self->{_genomes}}, "\nMLST_LENGTH: ", scalar @{$self->{_loci}}, "\n";
	
	# Print the genomes present in the model in the format "GENOME 1: Giardia_duodenalis_46504", one per line
	for ( my $i=0; $i < scalar @{$self->{_genomes}}; $i++ )	{
		print MODEL "GENOME ", $i + 1, ": $self->{_genomes}->[$i]\n";
	}
	
	# Print the sequences in the format "MLST 1: ATGGCTATCGGCGCTAA", one per line (note that the alignments will be encoded in ASCII chars)
	for ( my $i=0; $i < scalar @{$self->{_genomes}}; $i++ )	{
		print MODEL "MLST SEQUENCE ", $i + 1, ":\t$self->{_alignment}->[$i]\n";
	}
	
	# Print the profile in the format "A:	0.24	0.5		0.8	...", one per line
	foreach my $haplotype ( @alphabet )	{
		print MODEL "$haplotype:";
		for ( my $j=0; $j < scalar @{$self->{_loci}}; $j++ )	{
			print MODEL " $self->{_profile}->{$haplotype}->[$j]";
		}
		print MODEL "\n";
	}
	
	# The record separator // 
	print MODEL "\/\/\n";	
	close MODEL;
	
	# Finally, print each Locus model out with our "//" record separator
	for ( my $i=0; $i < scalar @{$self->{_loci}}; $i++ )	{
		my $name = $self->{_loci}->[$i]->{_name};
		if ( not $name )	{
			next;
		}
		my $dirname = dirname($output);
		$self->{_loci}->[$i]->print_model( "$dirname\/$name.locus$i.MLST.mod" ) if ( $split == 1 );
		$self->{_loci}->[$i]->print_model( "$output", ">>" );
	}
}

# Wipe out all data and loci in the current model and return the empty shell.
sub empty_model		{
	my ( $self ) = @_;
	
	# Reset all the attributes except the name
	$self->{_loci} = [ ]; 		
	$self->{_genomes} = [ ];					
	$self->{_alignment} = [ ];
	$self->{_profile} = [ ];
	$self->{_alphabet} = Set::Scalar->new();
	
	# Reset the organizational hash tables
	$h = 0;			
	%Index = ();
	%ReverseIndex = ();		
	%Names = ();			
	%NamesIndex = ();
	%NamesState = ();
	$g = 0;
	%GenomesIndex = ();
	%ReverseGenomesIndex = ();
	
	return $self;
}

# Split an MLST model into chunks of loci of a user-defined size, printing each chunk out
# You might want to do this if you have a large number of loci or are running into RAM issues
# The subroutine expects you to tell it where to write the block files, and how many loci to include in each block (200 loci by default)
sub cache_model		{
	my ( $self, $dir, $block_size ) = @_;
	$block_size = ( $block_size && $block_size > 0 )? $block_size : 200;	
	$dir = ( $dir && -d $dir )? $dir : ".";
	my $block_list;
	my $j=1;
	for ( my $i=0; $i < scalar @{$self->{_loci}}; $i+=$block_size )	{
		for ( my $k=$i; $k < $i+$block_size; $k++ )		{
			my $locus = $self->{_loci}->[$k];
			next if ( not $locus );
			$locus->print_model( "$dir\/mlst_chunk.$j.mod", ">>" );
		}
		push @{$block_list}, "$dir\/mlst_chunk.$j.mod";
		$j++;
	}
	return $block_list;
}	

# Given a block file, load in the correct cache file and return this object
sub load_cache	{
	my ( $self, $cache_file ) = @_;
	return " --> MLST::load_cache ERROR :: your cache file doesnt exist!!\n" if ( not -e $cache_file );
	$self->parse_locus_model( $cache_file );
	return $self;
}

# Given a list of schema files, read them and merge them into a single large schema.
# This function operates the same as if we were concatenating a FASTA file.
sub merge_schemas	{
	my ( $self, $list_of_schemas ) = @_;
	
	# Sanity check: do we have a valid list of schemas?
	return if ( not $list_of_schemas );

	# Loop through the list of schemas, read each one into a tmp object, and add it to the overall MLST
	for ( my $i=0; $i < scalar @{$list_of_schemas}; $i++ )	{
		$self = $self->_merge( $list_of_schemas->[$i] );
	}
	return $self;
}	

# The actual workhorse for merging two schemas.  Merges a 2nd MLST schema with the input object instance.
# If you want to merge more than two, this is done in an iterative loop.
sub _merge	{	
	my ( $self, $second_schema ) = @_;
	return if ( not $second_schema && -e $second_schema );
	
	# Read the model with the record separator since we can have multiple models in the same file
	$/ = "\/\/";
	open MODEL, "$second_schema";
		my @input = <MODEL>;
	close MODEL;
	$/ = "\n";

	# Parse the model chunks in $list_of_schemas->[$i]
	# Note: we are doing this in reverse order, from the last element in the model to the first.
	# This is done because we need to add the loci to the overall MLST model FIRST before we do the MLST genomes/sequences
	my $NewGenomes = [ ]; my @list_of_loci;
	my $flag = 0;
	foreach my $model ( reverse @input )	{
		my @model = split("\n", $model);
		my $name; my $data; my $position; my $type;
		
		foreach my $line ( @model )	{
			chomp $line;
			
			########## PARSE THE MLST RECORD #############
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			if ( $line =~ /^GENOME/ )	{
				my @tmp = split(": ", $line);
				push @{$NewGenomes->[0]}, $tmp[1];
			}
			
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			elsif ( $line =~ /^MLST SEQUENCE/ )	{
				my @tmp = split(":\t", $line);
				push @{$NewGenomes->[1]}, $tmp[1];
			}
			
			########## PARSE THE LOCI RECORDS #############
			# Get the model name
			elsif ( $line =~ /^MODEL NAME/ )	{
				my @tmp = split(": ", $line);
				$name = $tmp[1];
				push @list_of_loci, $name;
				$flag = 1;
			}
			# Get the position, if it's documented
			if ( $line =~ /^POSITION/ )		{
				my @tmp = split(": ", $line);
				$position = $tmp[1];
			}
			# Get the sequence data --> this maintains the order of the sequences as they are written in the model.
			if ( $line =~ /^SEQUENCE/ )	{
				my @tmp = split(":\t", $line);
				push @{$data}, $tmp[1];
			}
			# Get the type of model
			if ( $line =~ /^TYPE/ )	{
				my @tmp = split(": ", $line);
				$type = $tmp[1];
			}
		}
		
		# Instantiate the Locus object and add it to the list, if we have the flag is active
		if ( $flag == 1 )	{
			$self->add_locus( $data, $name, $position, $type ) if ( not defined $self->_nameToIndex($name) );
		}

		# Reset the flag
		$flag = 0;
	}

	# Add the new genomes and their sequences to the model
	for ( my $a=0; $a < scalar @{$NewGenomes->[0]}; $a++ )	{
		my $genome_idx = $self->_genomeToIndex( $NewGenomes->[0]->[$a] );
		
		# Derive the correct sequence, including any redundant loci
		# Start by initializing an empty mlst sequence of all gaps equal to the number of loci in the model currently.
		my @sequence;
		for ( my $s=0; $s < scalar @{$self->{_loci}}; $s++ )	{	
			$sequence[$s] = "-";
		}
		# If we already have some existing sequence data, slot it into place
		if ( $self->{_alignment}->[$genome_idx] )	{
			for ( my $b=0; $b < length $self->{_alignment}->[$genome_idx]; $b++ )	{
				$sequence[$b] = substr($self->{_alignment}->[$genome_idx], $b, 1);
			}
		}
		# Loop through the loci in the merging model to generate the correct MLST sequence for this genome
		for ( my $l=0; $l < scalar @list_of_loci; $l++ )		{
			$sequence[ $self->_nameToIndex( $list_of_loci[$l] ) ] = substr( $NewGenomes->[1]->[$a], $l, 1 );
		}
	
		# Finally, add the new genome to the model
		if ( $self->{_alignment}->[$genome_idx] )	{
			$self->update_genome( $NewGenomes->[0]->[$a], join('', @sequence) );
		}
		else	{	
			$self->add_genome( $NewGenomes->[0]->[$a], join('', @sequence) );		# Genome name, derived MLST sequence
		}
	}

	# Update the profile and reset the tmpMLST object for the next loop
	return $self;
}

############################################################
#          ALIGNMENT DISTANCE SUBROUTINES                  #
############################################################

# Compute the Levenshtein distance between two strings.
# This is a dynamic programming algorithm much like sequence alignment,
# except that here, we expect the two strings to already be aligned in the MLST object container.
# This computes their respective distance.
# CREDIT: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Perl
# I shameslessly jacked this algorithm from the above source, blame them if it ain't right.
sub levenshtein_distance	{
	my ( $self, $string1, $string2 ) = @_;
    my @ar1 = split('', $string1);
    my @ar2 = split('', $string2);

    my @dist;
    $dist[$_][0] = $_ foreach (0 .. @ar1);
    $dist[0][$_] = $_ foreach (0 .. @ar2);

    foreach my $i (1 .. @ar1){
        foreach my $j (1 .. @ar2){
            my $cost = $ar1[$i - 1] eq $ar2[$j - 1] ? 0 : 1;
            $dist[$i][$j] = min(
                        $dist[$i - 1][$j] + 1, 
                        $dist[$i][$j - 1] + 1, 
                        $dist[$i - 1][$j - 1] + $cost );
        }
    }

    return $dist[@ar1][@ar2];
}

# Generates a pairwise distance matrix of Levenshtein distances
# for all possible pairs in the MLST object.
sub pairwise_distance_matrix	{
	my ( $self ) = @_;
	my $matrix = [ ];
	for ( my $i=0; $i < scalar @{$self->{_alignment}}; $i++ )	{
		for ( my $j=$i; $j < scalar @{$self->{_alignment}}; $j++ )	{
			my $distance = ( $i == $j )? 0 : $self->levenshtein_distance( $self->{_alignment}->[$i], $self->{_alignment}->[$j] );
			$matrix->[$i]->[$j] = $distance;
			$matrix->[$j]->[$i] = $distance;
		}
	}
	return $matrix;
}






# Le fin
1;
