#     GUAVA — an rapid analysis of diversity in gene-by-gene datasets.
#     Copyright (C) 2017 Sofia Hauck
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/perl

#--------------------------------------------------------
# PROGRAM: GUAVA.pl
# AUTHOR:  Sofia Hauck
# CREATED: 05.02.2016
#--------------------------------------------------------

use strict;
use warnings;
use Benchmark; 
use LWP::Simple; # for downloading API pages 
use List::Util qw( min max sum ); # for measuring lengths

#$| = 1; # for dynamic output
my $version = "1.0";

# Declares subroutines
sub Usage( ; $ );


# Get current time in a nice format
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();


# Get Command line options, exits if conditions don't look right
my $fTable;			# table with isolates and loci
my $transpose;		# whether the table needs to be transposed or no
my $dFAS;			# directory where all the FASTA files for each locus are, if they don't need to be exported
my $dbname;			# database name if want to export the files
my $FASTAoption = 0;# 1 = have complete dir, 2 = make complete dir, 3 = take straight to Observed
my $align;			# whether to complete alignments
my $dOut;			# directory where the results will go
my $dup = 1;		# the smallest number of duplicates necessary for locus to be considered "seen", default is 1
my $locuscat;		# file with categories of loci for graphs
my @mafftarg = ("--clustalout","--quiet"); # start with at least these, can ask for more, add "-mafft --auto" to keep silent

my $i = 0; 
my $arg_cnt = 0; 
for ($i=0; $i<=$#ARGV; $i++)
{
 	if($ARGV[$i] eq "-help")		{ Usage("You asked for help"); exit; }
 	if($ARGV[$i] eq "-table")		{ $fTable		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-transpose")	{ $transpose 	= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-FASTA")		{ $dFAS 		= $ARGV[$i+1] || ''; $FASTAoption = "A"; $arg_cnt++; }
 	if($ARGV[$i] eq "-dbname")		{ $dbname 		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-dboption")	{ $FASTAoption 	= $ARGV[$i+1] || ''; $arg_cnt++; }
	if($ARGV[$i] eq "-align")		{ $align	 	= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-out")			{ $dOut  		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-dup")			{ $dup  		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-locuscat")	{ $locuscat		= $ARGV[$i+1] || ''; $arg_cnt++; }
 	if($ARGV[$i] eq "-mafft")		{ push (@mafftarg, $ARGV[$i+1]) ; $arg_cnt++; }
}


# Hello!
	print << 'EOU';

GUAVA  Copyright (C) 2017 Sofia Hauck
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions.

============================================================
Welcome to the 'GUAVA' processing script!
============================================================

EOU


# =================================
# Receiving input options
# =================================

# a table of loci and isolates 
if(! defined $fTable) 
{ 
	print "1. Where is the tab-separated table of isolates and loci?\n"; 
	$fTable = <STDIN>; 
	chomp $fTable; $fTable =~ s/\s+$//; # removes white spaces and line breaks
}
if(! -e $fTable)  { Usage("Input table file does not exist: $fTable"); exit; }

# are the loci split into categories?
if(! defined $locuscat) 
{ 
	print "\n2. Do you have a tab-separated file with loci categories? If not, leave blank.\n"; 
	$locuscat = <STDIN>; 
	chomp $locuscat; $locuscat =~ s/\s+$//; # removes white spaces and line breaks
	
	if( $locuscat !~ /\w/ ) # if left blank
	{ 	$locuscat = "notdef"	}
	else
	{ if(! -e $locuscat)  { Usage("Input category file does not exist: $locuscat"); exit; } }
}

# and whether to transpose that table or no
if (! defined $transpose )
{
	print "\n\n3. Does the table have loci as columns? \nIf yes, please enter 'Yes' as the table needs transposing. " .
	"\nIf no, that is, if the each row of your table is a locus, transposing is not needed, so enter 'No' below.\n"; 
	$transpose = <STDIN>; 
	if ( ! $transpose =~ /^[YyNn]/ ) # as long as it starts with a "y" or "n"
	{ 
		print "That didn't look like a yes or a no... Try again? \n";
		$transpose = <STDIN>; 
	}
	elsif ( $transpose =~ /^[Yy]/ )
	{ $transpose = "Yes" }
	elsif ( $transpose =~ /^[Nn]/ )
	{ $transpose = "No" }
	else 
	{ Usage("Something went wrong with the transposing options"); exit; }
}

# FASTA files or where to grab them
if( ! defined $dFAS && ! defined $dbname)  # don't have anything defined 
{ 
	if ( $FASTAoption !~ /^[ABC]/ ) # if doesn't start with 1, 2 or 3
	{ 
		print 	"\n\n4. Where is your deposit of FASTA sequences? Choose a letter (A, B or C).\n".
		"   A. You have a folder with all the sequences that you can point me to (fastest)\n".
		"   B. You want to get all the sequences from PubMLST to have your own deposit (slow now, fast later)\n".
		"   C. You want to grab only the relevant sequences for this (best if you're doing just this run)\n";
		$FASTAoption = <STDIN>;

		if ( $FASTAoption !~ /^[ABC]/ )
		{
			print "That didn't look like an A, B or C... Try again? \n";
			$FASTAoption = <STDIN>; 
		}
	}

	if ( $FASTAoption =~ /^A/ ) # You have a folder with all the sequences that you can point me to
	{
		print "\n4A. Where is your folder with your complete FASTA sequence deposit then?\n";
		
		$dFAS = <STDIN>; 
		chomp $dFAS; $dFAS =~ s/\s+$//; # removes white spaces and line breaks
		
		if(! -e $dFAS)
		{ Usage("Input FASTA directory doesn't exist: $dFAS"); exit; }
	}
	elsif	( $FASTAoption =~ /^[BC]/ ) # You want to get all the sequences from PubMLST to have your own deposit
	{
		print "\n4BC. What's the name of the database in PubMLST? For example: 'pubmlst_mycobacteria_seqdef'.\n";
		$dbname = <STDIN>; 
		chomp $dbname; $dbname =~ s/\s+$//; # removes white spaces and line breaks
	}
	
	else
	{ Usage("Something went wrong with the FASTA deposit options"); exit; }
	
}

elsif ( defined $dFAS ) # had defined a folder in the command line options, so going for option 2
{	$FASTAoption = "A";	}

elsif ( defined $dbname ) # had defined a PubMLST database name in the command line, so options B or C
{	
	if ( $FASTAoption !~ /^[BC]/ )
	{
		print "4bc. Do you want to grab the complete FASTA sequences (enter 'B') or just the necessary ones for this run ('C')?\n";
		$FASTAoption = <STDIN>;
	
		if ( $FASTAoption !~ /^[BC]/ ) 
		{	Usage("Something went wrong with the FASTA deposit options, maybe you didn't choose B or C?"); exit; }
	}	
}

my @tableaddress = split ("/", $fTable); #so that @tableaddress is usable later too
	$tableaddress[-1] =~ s/\.[^.]+$//; # remove extension from file name of input table


# and whether to align or not
if (! defined $align )
{
	print "\n\n5. Do you want to complete alignments? This is a very time-consuming step, but it is needed for measuring" .
	"\nproportions of variable sites. Enter 'yes' or 'no' below.\n"; 
	$align = <STDIN>; 
	if ( $align =~ /^[Yy]/ )
	{ $align = "Yes" }
	elsif ( $align =~ /^[Nn]/ )
	{ $align = "No" }
	else 
	{ print "Defaulting to no alignments."; $align = "No"; }
}


# adding any MAFFT arguments?
if( $align eq "Yes" && scalar(@mafftarg) == 2) 
{ 
	print "\n\n5B. Do you have any MAFFT arguments to include?\n" . 
		"You can just leave this blank if not, and MAFFT will run using its automatic option.\n" .
		"If you are using a large data set, try '--retree 1 --maxiterate 0' for the fast FFT-NS-1 method.\n"; 
	my $addtoMAFFT = <STDIN>; 
	chomp $addtoMAFFT; # removes white spaces and line breaks
	
	push ( @mafftarg, $addtoMAFFT );
}


# a directory where results can go
if(! defined $dOut)  
{ 	
	print 	"\n\n6. Name the directory where you'd like the results to go.\n".
			"If you can't decide, leave this blank and a new folder, in the same place as your ".
			"table is held will be created, named '$tableaddress[-1]'\n";
	$dOut = <STDIN>; 
	chomp $dOut; $dOut =~ s/\s+$//;

	if ( $dOut =~ /\w/ )
	{	if(  -e $dOut)  { Usage("Output directory already exists: $dOut"); exit; } }
	else 
	{	$dOut = join ("/", @tableaddress); }
}


# =================================
# Check before starting the work
# =================================

print "\n\n\nLet's check if everything is correct before we start...\n\n";

# Table location
if ( $transpose =~ /^[Nn]/ )
{ print	"1. Your table of isolates as columns and loci as rows is held at $fTable.\n"  }
elsif ( $transpose =~ /^[Yy]/ )
{ print	"1. Your table of loci as columns and isolates as rows is held at $fTable.\n"  }

# Categories included
if ( $locuscat =~ /\w/ )
{ print	"2. Your table of locus names and categories is held at $locuscat.\n"  }
else
{ print	"2. You have not chosen to include locus categories.\n"  }

# Transposing
if ( $transpose =~ /^[Nn]/ )
{ print	"3. Your table will not be transposed.\n"  }
elsif ( $transpose =~ /^[Yy]/ )
{ print	"3. Your table will be transposed.\n"  }


# Where the sequences are 
if ( $FASTAoption =~ /^A/ )
{	print "4. Your directory of FASTA sequences is at $dFAS.\n";	}
elsif ( $FASTAoption =~ /^B/ )
{	print "4. The complete FASTA sequences will be downloaded from the PubMLST API from database $dbname, and put in /PubMLST-FASTA/ in your results folder.\n"; 	}
elsif ( $FASTAoption =~ /^C/ )
{	print "4. There won't be a deposit of FASTA sequences; we'll grab them straight from $dbname. \n";	}

# Alignments
if 		( $align eq "Yes" )	{ print "5. Alignments will be complete, with MAFFT arguments: @mafftarg\n"; } 
elsif 	( $align eq "No" )	{ print "5. You chose not to perform alignments.\n"; } 

# Output directory
if( -e $dOut) { Usage("Output directory already exists: $dOut"); exit; }

if ( mkdir $dOut )
{ print "6. Your results will be in the new directory $dOut\n"; }
else { Usage("Could not create output directory: $dOut"); exit; }

# Others
print "7. Alleles must be seen at least $dup time(s) in order to be included.\n\n\n"; 


# Give the option to get out if it goes in fact look dodgy
print "\nLeave blank to continue, enter any letters to exit.\n\n";
my $confirmation = <STDIN>;
if ( $confirmation =~ /\w/ )
{ Usage("You chose to exit: $confirmation"); exit;  }


# First benchmark point
my $t0 = Benchmark->new;



if ( $FASTAoption =~ /^B/ )
{ mkdir $dOut."/PubMLST-FASTA/" or Usage("Could not create FASTA exports folder: $dOut/PubMLST-FASTA/"); }



# =================================
# Actually starting the work now
# =================================


# If we have a category table
# =================================

my %categories = (); 
my $catrowcount = 0;

if ( $locuscat ne "notdef" )
{
	open(CATTABLE, $locuscat) or die "Cannot open $locuscat\n";

		while ( my $line = <CATTABLE> )
		{
			chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		
			if ( my @row = split ('\t', $line) ) # split row by tabs into array of cells
			{
				# first column is isolate unique ID, second column is whatever group, ignore the rest
				$categories{$row[0]} = $row[1]; 
				my $catrowcount ++; 
			}
			else { next; } # if line is blank or doesn't have at least two elements, skip 
				
		} 

		#if ( $catrowcount <= 1 ) { die "You only have one column so something has probably gone wrong with line breaks.\n"; }
	
	close(CATTABLE);

	my @obscats = do { my %seen; grep { !$seen{$_}++ } values %categories };

	print "The locus categories imported are (" . scalar(@obscats). " total): " . join ("; ", @obscats) . "\n\n";
}


# Reading the table in
# =================================

my $tableref = ReadTableIn ( $fTable );
my @aoaTable = @$tableref; #puts table, in array of arrays format, into aoaTable

if ( $transpose =~ /^[Yy]/ ) # transposes the table and puts it back in aoaTable
{
	print "Transposing your table...";
	$tableref = TransposeTable ( \@aoaTable );
	@aoaTable = @$tableref;
	print " and done.\n";
}

my $headerrow = $aoaTable[0];
my $isolatecount = scalar (@$headerrow) - 1; # since the first column is the header

# Printing the header of the results file
open ( RESULTS, '>', $dOut."/ResultsTable.txt" ) or die "$dOut /ResultsTable.txt"; # then also open where the reduced one will go

	my @headerlist = ("Locus","Missing","MultipleCopies",
					  "CountNuc","CountAA",
					  "MinLength","MaxLength","AvgLength");

	if ( $locuscat ne "notdef" )
	{	splice @headerlist, 1, 0, "Category";	}

	if ( $align eq "Yes" )
	{	splice @headerlist, 9, 0, "NonVarNuc\tNonVarAA\tAliLenNuc\tAliLenAA";	}


	print RESULTS	"GUAVA results table from '". $tableaddress[-1] ."' records table.\n" . 
    				"Created at ".sprintf("%02d", $hour).":".sprintf("%02d", $min).", on $days[$wday] $mday $months[$mon], ". (1900 + $year) . "\n" . 
    				"Software version: " . $version . "\n" . 
    				"Isolate records read from table: ". $isolatecount . "\n\n".
					join ("\t", @headerlist ) . "\n" ;
	
	my @resultsline = (); 


# Going through table and keeping the observed alleles, plus counting missing and unique nucleotide sequences 
# ============================================================================================================

print "Now reading your table in and grabbing the relevant FASTA sequences...\n\n";

mkdir $dOut."/Observed-FASTA/" or die "Could not create filtered FASTA folder: $dOut/Observed-FASTA/";

foreach my $locusrow (@aoaTable) # loop per locus
{
	# empty hash to find unique values, values to be found
	my %unique_alleles = ();
	my $paralogous = 0; 
	my $missing = 0;
	my $countnuc = 0;

	# set locus name
	my $locusname = shift(@{$locusrow});
	$locusname =~ /^(\S+)/; $locusname = $1; # if there are aliases or anything else after the locus ID
			
	# populates that hash with unique allele numbers as keys and their number of appearances as values
	foreach my $allele (@{$locusrow}) #go through each allele in locus
    {
    		if ( length ( $allele ) == 0 || $allele eq "X" || $allele eq "I" )
    		{	$unique_alleles{"0"}++;	} # if cell is empty or has "X", file as "missing locus"
    		
    		else
    		{
    			my @alleles = split (/;/, $allele); # in case of paralogous loci
				if ( scalar(@alleles) >= 2 )
				{	$paralogous++; }

				foreach (@alleles) # go through each allele in locus
				{	$unique_alleles{$_}++;	} #increase the frequency count	
			}
	}

	if ( exists($unique_alleles{"0"}) ) # move count away from empty if any missing alleles were seen
	{	$missing = $unique_alleles{"0"}	} # gives the value in the frequency hash when the key is allele "0", the missing allele		

	if ( $missing == $isolatecount )
	{ print "Locus '$locusname' is being excluded since there were no allele designations for it in any genomes.\n"; next; }

	# find the file where the original FASTA sequences are	
	
	if ( $FASTAoption =~ /^[AB]/ ) # if there is or will be a directory
	{
		my $fullFAS; # where the directory with all possible sequence is / will be
		
		if ( $FASTAoption =~ /^B/ ) # making directory by exporting from BIGS
		{
			GetFASTASeqs ( $dbname, $locusname ) ; # the bit where the file is copied from PubMLST
			$fullFAS = $dOut."/PubMLST-FASTA/".$locusname.".FAS";
		}
		elsif ( $FASTAoption =~ /^A/ ) # if directory already exists and we were pointed to it
		{	$fullFAS = $dFAS."/".$locusname.".FAS"; }
		
		if ( open(FULLFASTA, $fullFAS) ) # if can open file in directory 
		{
			my $save = 1; # whether to copy the sequence following a >identifier line
		
			my $reducedFAS = $dOut."/Observed-FASTA/".$locusname.".FAS"; # where the "seen" alleles get copied to  
			open(REDFASTA, '>', $reducedFAS) or die "Cannot write to $reducedFAS\n"; # then also open where the reduced one will go

			while ( my $line = <FULLFASTA> ) # reading through original FASTA
			{
				if ( $line =~ /^>/ ) # if line is >identifier
				{
					chomp $line; 
				
					$line =~ /([^>_]+)$/; $line = $1; #removes > and anything before the last _
				
					if ( exists($unique_alleles{$line}) && $unique_alleles{$line} >= $dup ) # exists in hash of wanted loci
					{
						print REDFASTA "\n>", $line, "\n";
						$save = 1;
						$countnuc ++;
						$unique_alleles{$line} = 0; # set value to zero (normally would be undef or 1+) as check that has been read
					} 
					else # if not in the hash created from the table, indicating whether this allele is observed
					{ $save = 0; } 
				}
			
				elsif ( $line =~ /^([A-Za-z*-])/ ) # if is a sequence line (starts with a letter, - or *) 
				{
					chomp $line; 
					if ($save == 1) # copy only if "save" is turned on because last ID was wanted
					{ print REDFASTA $line; }
				}	
			}
			
			# check if all alleles seen in table were copied 
			foreach my $key ( sort keys %unique_alleles ) 
			{ 
				if ( $unique_alleles{$key} >= $dup && $key ne 0 ) # if frequency is still above cut-off (copied ones are reset to 0) and isn't the missing allele, 0
				{ print "Did not find sequence for locus $locusname, allele $key.\n"; }
			}
			
			
 			@resultsline = ($locusname, $missing, $paralogous, $countnuc);
 			
 			if ( $locuscat ne "notdef" )
			{	splice @resultsline, 1, 0, $categories{$locusname}; }
 			
 			print RESULTS join ("\t", @resultsline ), "\n";	
			
			close(FULLFASTA);
		}
		
		else # if couldn't open the file
		{ 
			my @seenalleles = sort keys %unique_alleles;

			print "'$locusname' did not exist as a FASTA file. Alleles in table were..."; 
		
			my $sampleallele = 0; 
			foreach my $key ( @seenalleles ) 
			{ 
				if ( $sampleallele < 5 )
				{ print " $key"; $sampleallele++; }
				else # if already printed 5 allele numbers, just leave it
				{ print " etc."; last; }
			 } print "\n"; 
			
			next; # since don't need to include results for a non-locus
		} 
		
		
	}
	
	elsif ( $FASTAoption =~ /^C/ )
	{ 		 
 		my $locusurl = "http://rest.pubmlst.org/db/".$dbname."/loci/".$locusname ;
 		 		
 		my $locusJSON = get($locusurl);
 		
 		if ( ! defined $locusJSON ) # $locusJSON =~ /"status":404/ && $locusJSON =~ /"message":"Locus/) ?
 		{
			print "Locus '$locusname' was not found in PubMLST database '$dbname' and will not be included in results table.\n";
 		}
 		else # not a locus based 404, so it exists 
 		{
 		 	my $reducedFAS = $dOut."/Observed-FASTA/".$locusname.".FAS"; # where the "seen" alleles get copied to  
			open(REDFASTA, '>', $reducedFAS) or die "Cannot write to $reducedFAS\n"; # then also open where the reduced one will go
 		 	
 		 	foreach my $allele ( sort { $a <=> $b } keys %unique_alleles )
			{
	 			if ( $allele eq "0" || $allele eq "" )
	 			{ next; }
	 			
	 			my $url = $locusurl . "/alleles/".$allele;
	 			my $alleleJSON = get($url);
				
				if ( ! defined $alleleJSON )
				{ print "Could not find locus $locusname, allele $allele"; next; }
				
				
				if ( $alleleJSON =~ /"status":"Allele/ ) # then just couldn't find this allele
				{
					my $message =~ /"message":"(.+?)"/ ;
					print $1 . "\n";
				}
				else
				{
					$alleleJSON =~ /"sequence":"([a-zA-Z]+?)"/; # assuming sequence will only have letters
					my $sequence = $1;
					$alleleJSON =~ /"allele_id":"([0-9]+?)"/; # assuming allele ID will only be integers (in BIGS anyway)
					my $allele_id = $1;
	
					print REDFASTA ">" . $allele_id . "\n" . $sequence . "\n";
					$countnuc ++; 
				}

 			} # closes foreach allele loop
 			
 			@resultsline = ($locusname, $missing, $paralogous, $countnuc);
 			
 			if ( $locuscat ne "notdef" )
			{	splice @resultsline, 1, 0, $categories{$locusname}; }
 			
 			print RESULTS join ("\t", @resultsline ), "\n";	

 		} # closes found JSON else
	} # closes opt 3 elsif
	
	close(REDFASTA);
		
	
} # close per locus loop

close RESULTS; # close results file to adding one line per locus as it goes through table-reading loop
print "\nTable reading complete.\n";

my %results = (); # creates a hash where the rest of results go under the locusname as the key, to be added to the results table at the end



# =================================
# Translations
# =================================
mkdir $dOut."/Translated-FASTA/" or die "Cannot create /Translated-FASTA/ folder";

# opening the directory where the observed alleles are and getting all the file names
opendir (OBSDIR, $dOut."/Observed-FASTA/" ) or die "Cannot open directory: $!";
	my @files = readdir OBSDIR;
	@files = grep(/^([A-Za-z0-9])/,@files); # remove any hidden files
	if ($#files < 1)
	{ die "Couldn't find any FASTA files. Maybe their names don't begin with a letter or number?\n";}
	@files = sort @files; # so they're in order and the "Currently up to" is a bit more informative
closedir OBSDIR;


# Loop that translates all the files and puts them in new directory
print "\n\nTranslated up to...\n";

foreach my $file (@files)
{
 	my $infile  = $dOut . "/Observed-FASTA/"   . $file; # file with FASTA sequences in nucleotide form
 	my $outfile = $dOut . "/Translated-FASTA/" . $file; # file with FASTA sequences in amino acid form
 	
	my %unique = (); # empty hash to check for synonymous mutations 

	# for counting GC content in locus
	my $countGC  = 0; 	my $countAT  = 0;


	open(NUCLEOTIDE, $infile) 		or die "Cannot open $infile\n";
	open(AMINOACID, '>', $outfile) 	or die "Cannot open $outfile\n";

 	while ( my $line = <NUCLEOTIDE> )
	{
		if ( $line =~ /^>/ ) # if line is title of sequence
		{ 
			print AMINOACID $line;
		}
		elsif ( $line =~ /^[A-Z]/ ) # if line is sequence
		{
			chomp $line;
			my $peptide = Translate($line, "1");
			
			# counting GC content 
			$countGC = $countGC + ($line =~ tr/GC//); 
			$countAT = $countAT + ($line =~ tr/AT//); 

			# for later checking the count of unique amino acids
			if ( !exists($unique{$peptide})) #if allele doesn't exist in unique-hash
			{ $unique{$peptide} = 1; } #then add it with allele number as key, value as 1
			
			print AMINOACID $peptide, "\n";
		}	
	}
	
	# counting number of unique amino acids and measuring the length of the locus
 	my ($locusname, $extension) = split (/\./, $file);
 	my $countaa = 0; 
 	my @lengths = ();
 	
 	foreach my $key ( sort keys %unique ) # goes through hash where keys are the amino acid sequences for each allele
 	{ 	
 		$countaa++;
 		push (@lengths, (3 * length($key) ) );
 	}

	my $min = "0";
	my $max = "0";
	my $avg = "0";

	if ( $countaa != 0 ) # if no alleles for that locus, just keep it all at zero
	{
		$min =  &min ( @lengths ); 
		$max =  &max ( @lengths ) ;
		$avg =  sprintf ('%.4g', (&sum ( @lengths )) / $countaa) ;
 	}
 	
 	my $GCpercent = 0; 
 	
 	if ( $countGC != "0" ) # if no alleles for that locus, avoid dividing by zero
 	{ $GCpercent = $countGC / ( $countGC + $countAT ); }
 	 	
 	close NUCLEOTIDE; close AMINOACID;
 	 	
 	$results{$locusname} = join ("\t", ($countaa, $min, $max, $avg, $GCpercent) ) ;

 	#So you have something to watch while it runs... 
 	print "\r$locusname";
} #close translating loop

print "\nTranslation complete.\n";


# =================================
# Alignments
# =================================

if ( $align eq "Yes" ) # only if requested
{
	mkdir $dOut."/AlignedNuc-Clustal/" or die "Cannot create /AlignedNuc-Clustal/ folder";
	mkdir $dOut."/AlignedAA-Clustal/"  or die "Cannot create /AlignedAA-Clustal/ folder";

	print "\nCurrently aligning...\n";

	foreach my $file (@files)
	{
		my ($locusname, $extension) = split (/\./, $file);
		
		# where all my files at
		my $inNuc	=	$dOut."/Observed-FASTA/"		. $file; 
		my $inAA	=	$dOut."/Translated-FASTA/"		. $file;
		my $outNuc	=	$dOut."/AlignedNuc-Clustal/"	. $locusname.".aln"; 
		my $outAA	=	$dOut."/AlignedAA-Clustal/"		. $locusname.".aln"; 

		# runs MAFFT commands 
		system ( "mafft " . join (" ", @mafftarg) . " " . $inNuc . " > " . $outNuc  ); 
		system ( "mafft " . join (" ", @mafftarg) . " " . $inAA . " > " . $outAA  ); 

		# now counting the variable sites
		my $nonvarNuc = 0;
		my $alilenNuc   = 0;
		my $nonvarAA  = 0;
		my $alilenAA    = 0;
	
		open (ALIGNEDNUC , $outNuc ) or die "Cannot open $outNuc";
			while (my $line = <ALIGNEDNUC>)
			{	
				next if $. < 2; # skips the first three lines of CLUSTAL file
				$alilenNuc  = $alilenNuc + ($line =~ tr/-atcgATCG//);
				$nonvarNuc  = $nonvarNuc + ($line =~ tr/\*//);
			}
		close ALIGNEDNUC;

		open (ALIGNEDAA , $outAA ) or die "Cannot open $outAA";
			while (my $line = <ALIGNEDAA>)
			{	
				next if $. < 2; # skips the first three lines of CLUSTAL file
				$alilenAA   = $alilenAA + ($line =~ tr/-ACDEFGHIKLMNPQRSTVWY//);
				$nonvarAA   = $nonvarAA + ($line =~ tr/\*//); 
			}
		close ALIGNEDAA;

		$results{$locusname} = $results{$locusname} . "\t" . join ("\t", ($nonvarNuc, $nonvarAA, $alilenNuc, $alilenAA) ) ;
	
		# So you have something to watch while it runs... 
		print "\r$file";
	} # closes per-locus loop

	print "\nAlignments complete.\n";
}


# =================================
# Printing results to file, plus showing time taken
# =================================

open ( RESULTSIN,  '<', $dOut."/ResultsTable.txt" ) or die "$dOut /ResultsTable.txt";
open ( RESULTSOUT, '>', $dOut."/ResultsTable-tmp.txt" ) or die "$dOut /ResultsTable.txt";	
	
	my $linecount = 1; 
	while ( my $line = <RESULTSIN> )
	{
		if ( $linecount <= 5 )
		{
				print RESULTSOUT $line ;
				$linecount ++ ;
		}
		
		elsif ( $linecount == 6 )
		{
			my @headerline = (	"Locus","Missing","MultipleCopies",
								"CountNuc","CountAA",
								"MinLength","MaxLength","AvgLength",
								"GCContent","AllelicDiv", "RatioCount");
			
			if ( $align eq "Yes"  )
			{	push (@headerline, ("AliLenNuc", "AliLenAA", 
									"VSitesNuc", "VSitesAA", "RatioVS") ); }			
			
			if ( $locuscat ne "notdef" ) 
			{ 	splice @headerline, 1, 0, "Category";	}	
			
			print RESULTSOUT join ("\t", @headerline) . "\n"; 
			
			$linecount ++ ;
		}
		
		else 
		{
			$linecount ++ ;
			
			chomp $line; 
			
			# getting the relevant stats for that locus from existing printout and hash
			my @printedresults = split (/\t/, $line);
			
			# getting the locus name from the line we're at
			my $locusname = shift ( @printedresults ); 

			# finding the category name (if used)
			my $category; 
			if ( $locuscat ne "notdef" )
			{ $category = shift ( @printedresults ) ; }
			
			my $missing = shift ( @printedresults );  
			my $paralogous = shift ( @printedresults ); 
			my $countnuc = shift ( @printedresults ); 
									
			my @hashedresults  = split (/\t/, $results{$locusname} ); 
				my $countaa = $hashedresults[0];  
				my $lenmin  = $hashedresults[1];  
				my $lenmax  = $hashedresults[2];  
				my $lenavg  = sprintf ('%.4g', ($hashedresults[3]));  
				my $GCcont  = sprintf ('%.4g', ($hashedresults[4]));  
				
			# in case of alignments
			my $nvNuc; my $nvAA; 
			my $alNuc; my $alAA;
			my $VSitesAA; my $VSitesNuc; my $RatioVS; 
			
			if ( $align eq "Yes" )
			{	
				$nvNuc   = $hashedresults[5];
				$nvAA    = $hashedresults[6]; 
				$alNuc   = $hashedresults[7];
				$alAA    = $hashedresults[8];
			}


			# making calculations of new parameters
			my $AllelicDiv = 0 ;
			my $RatioCount = 0 ;

			if ( $lenavg != 0 )
			{
				$AllelicDiv = sprintf ('%.4g', ($countnuc / $lenavg  )) ;
				$RatioCount = sprintf ('%.4g', ($countaa  / $countnuc)) ;
			
			}
			
			if ( $align eq "Yes"  )
			{
				if ( $alNuc == 0 )
				{	$alNuc = "NA"; 
					$VSitesNuc = "NA" }
				else
				{	$alNuc     = sprintf ('%.4g', ($alNuc / $countnuc));
					$VSitesNuc = sprintf ('%.4g', (( $alNuc - $nvNuc ) / $alNuc ));	}

				if ( $alAA == 0 )
				{	$alAA     = "NA"; 
					$VSitesAA = "NA"; 
					$RatioVS  = "NA" }
				else
				{	$alAA     = sprintf ('%.4g', (  $alAA / $countnuc )); # alignments include redundant sequences
					$VSitesAA = sprintf ('%.4g', (( $alAA - $nvAA ) / $alAA)) ;
					$RatioVS  = sprintf ('%.4g', ($VSitesAA / $VSitesNuc))	}				
			}


			my @finalresults = ($locusname, $missing, $paralogous, 
								$countnuc, $countaa, 
								$lenmin, $lenmax, $lenavg, $GCcont, 
								$AllelicDiv, $RatioCount );
			
			# adding in special variables if alignments or categories were included
			if ( $align eq "Yes"  )
			{	push (@finalresults, ($alNuc, $alAA, $VSitesNuc, $VSitesAA, $RatioVS) ); }
			
			if ( $locuscat ne "notdef" ) 
			{ 	splice @finalresults, 1, 0, $category;	}

			print RESULTSOUT join ("\t", @finalresults) . "\n";
		}
	}
	
close RESULTSIN;
close RESULTSOUT; 

# after making changes by adding in more info, replace the original with the improved temporary file
rename ( $dOut."/ResultsTable-tmp.txt" , $dOut."/ResultsTable.txt" ) or die "Cannot rename temporary results over older.";

# Second Benchmark point
my $t1 = Benchmark->new;

my $td1 = timediff($t1, $t0);


print "\n\nGUAVA all complete! Processing took ",timestr($td1),".\n";

#---------------------------------------------------------------
# Subroutines
#---------------------------------------------------------------


sub ReadTableIn
{
	my $infile = $_[0];
	my @aoaTable = ();
	my $rowcount = 0; 
	
	open(INTABLE, $infile) or die "Cannot open $infile\n";

		while ( my $line = <INTABLE> )
		{
			chomp($line); $line =~ s/\r\n?|\n//g; # clean up line breaks
		  
			my @row = split ('\t', $line); # split row by tabs into array of cells
			$rowcount ++;
		
			push ( @aoaTable, [ @row ] );
		} 

		if ( $rowcount <= 1 ) { die "You only have one column so something has probably gone wrong with line breaks.\n"; }
	close(INTABLE);
	
	return \@aoaTable;
}

sub TransposeTable
{		
	my $aoaRef = $_[0]; 
	my @original = @$aoaRef; # @original is now array of references to arrays
		
	my @transposed = ();
	my $columncount = 0; 
		
	for my $row (@original) 
	{
  		for my $column (0 .. $#{$row}) 
  		{
    		push(@{$transposed[$column]}, $row->[$column]);
    		$columncount = $#{$row};
  		}
	}
		
	@transposed = splice (@transposed, 0, scalar(@{$original[0]})); #removes empty rows if more isolates than loci
	
	return \@transposed; 
}

sub GetFASTASeqs
{
	my $dbname	= $_[0];
	my $locus	= $_[1];

	my $url = "http://rest.pubmlst.org/db/".$dbname."/loci/".$locus."/alleles_fasta";
	my $file = $dOut."/PubMLST-FASTA/".$locus.".FAS";

	if ( ! getstore($url, $file) )
	{ print "Could not get complete FASTA sequences for '$locus' in database '$dbname'.\n"; }
	}

sub Translate 
{

# Author: John Nash
#
# Copyright (c) Government of Canada, 2000-2012,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged 
#   for use, and as long as the author/copyright attributions 
#   are not removed.


# Does:  Translates incoming string from nucleic acid to peptide
# Expects: nucleic acid string ($sequence_str); 
#          reading frame value: +1, +2, +3, -1, -2, -3
# Returns: string of translated peptide
# Uses:
    
# vars for translation:
    my ($codon, $peptide, $seq_str, $frame);
    
# reading frame:
    my $startpoint;
    
# incoming:
    $seq_str = uc $_[0];
    
# assess the reading frame:
    $frame = $_[1];
    $seq_str =~ tr/natgcbdkrvhmyxwsNATGCBDKRVHMYXWS
    /ntacgvhmybdkrxswNTACGVHMYBDKRXSW/ if ($frame < 0);
    $startpoint = (abs $frame) - 1;
    
# Perform the translation:
    for (my $j=$startpoint; $j < length $seq_str; $j+=3) {
	$codon=(substr($seq_str, $j, 3));
	if (($codon eq "GCA") || 
	    ($codon eq "GCC") ||
	    ($codon eq "GCG") ||
	    ($codon eq "GCT") ||
	    ($codon eq "GCR") ||
	    ($codon eq "GCY") ||
	    ($codon eq "GCN")) {
	    $peptide.="A";
	}   # Ala A
	
	elsif (($codon eq "CGA") || 
	       ($codon eq "CGC") ||
	       ($codon eq "CGG") ||
	       ($codon eq "CGT") ||
	       ($codon eq "CGR") ||
	       ($codon eq "CGY") ||
	       ($codon eq "CGN") ||
	       ($codon eq "AGA") ||
	       ($codon eq "AGG") ||
	       ($codon eq "AGR")) {
	    $peptide.="R";
	}   # Arg R
	
	elsif (($codon eq "AAC") || 
	       ($codon eq "AAT") ||
	       ($codon eq "AAY")) {
	    $peptide.="N";
	}   # Asn N
	
	elsif (($codon eq "GAC") || 
	       ($codon eq "GAT") ||
	       ($codon eq "GAY")) {
	    $peptide.="D";
	}   # Asp D
	
	elsif (($codon eq "TGC") || 
	       ($codon eq "TGT") ||
	       ($codon eq "TGY")) {
	    $peptide.="C";
	}   # Cys C
	
	elsif (($codon eq "CAA") || 
	       ($codon eq "CAG") ||
	       ($codon eq "CAR")) {
	    $peptide.="Q";
	}   # Gln Q
	
	elsif (($codon eq "GAA") || 
	       ($codon eq "GAG") ||
	       ($codon eq "GAR")) {
	    $peptide.="E";
	}   # Glu E
	
	elsif (($codon eq "GGA") || 
	       ($codon eq "GGC") ||
	       ($codon eq "GGG") ||
	       ($codon eq "GGT") ||
	       ($codon eq "GGR") ||
	       ($codon eq "GGY") ||
	       ($codon eq "GGN")) {
	    $peptide.="G";
	}   # Gly G
	
    elsif (($codon eq "CAC") || 
	   ($codon eq "CAT") ||
	   ($codon eq "CAT")) {
	$peptide.="H";
    }   # His H
	
	elsif (($codon eq "ATA") || 
	       ($codon eq "ATC") ||
	       ($codon eq "ATT") ||
	       ($codon eq "ATH")) {
	    $peptide.="I";
	}   # Ile I
	
	elsif (($codon eq "CTA") || 
	       ($codon eq "CTC") ||
	       ($codon eq "CTG") ||
	       ($codon eq "CTT") ||
	       ($codon eq "CTR") ||
	       ($codon eq "CTY") ||
	       ($codon eq "CTN") ||
	       ($codon eq "TTA") ||
	       ($codon eq "TTG") ||
	       ($codon eq "TTR")) {
	    $peptide.="L";
	}   # Leu L
	
	elsif (($codon eq "AAA") || 
	       ($codon eq "AAG") ||
	       ($codon eq "AAR")) {
	    $peptide.="K";
	}   # Lys K
	
	elsif (($codon eq "ATG")) {
	    $peptide.="M";
	}   # Met M
	
	elsif (($codon eq "TTC") || 
	       ($codon eq "TTT") ||
	       ($codon eq "TTY")) {
	    $peptide.="F";
	}   # Phe F
	
	elsif (($codon eq "CCA") || 
	       ($codon eq "CCC") ||
	       ($codon eq "CCG") ||
	       ($codon eq "CCT") ||
	       ($codon eq "CCR") ||
	       ($codon eq "CCY") ||
	       ($codon eq "CCN")) {
	    $peptide.="P";
	}   # Pro P
	
	elsif (($codon eq "TCA") || 
	       ($codon eq "TCC") ||
	       ($codon eq "TCG") ||
	       ($codon eq "TCT") ||
	       ($codon eq "TCR") ||
	       ($codon eq "TCY") ||
	       ($codon eq "TCN") ||
	       ($codon eq "AGC") ||
	       ($codon eq "AGT") ||
	       ($codon eq "AGY")) {
	    $peptide.="S";
	}   # Ser S
	
	elsif (($codon eq "ACA") || 
	       ($codon eq "ACC") ||
	       ($codon eq "ACG") ||
	       ($codon eq "ACT") ||
	       ($codon eq "ACR") ||
	       ($codon eq "ACY") ||
	       ($codon eq "ACN")) {
	    $peptide.="T";
	}   # Thr T
	
    elsif (($codon eq "TGG")) {
	$peptide.="W";
    }   # Trp W
	
	elsif (($codon eq "TAC") || 
	       ($codon eq "TAT") ||
	       ($codon eq "TAY")) {
	    $peptide.="Y";
	}   # Tyr Y
	
	elsif (($codon eq "GTA") || 
	       ($codon eq "GTC") ||
	       ($codon eq "GTG") ||
	       ($codon eq "GTT") ||
	       ($codon eq "GTR") ||
	       ($codon eq "GTY") ||
	       ($codon eq "GTN")) {
	    $peptide.="V";
	}   # Val V
	  
	elsif (($codon eq "TAA") || 
	       ($codon eq "TAG") ||
	       ($codon eq "TAR") ||
	       ($codon eq "TGA")) {
	    $peptide.="*";
	}   # Stop *
	
	else {
	    ;
	} # do nothing for now...
    } # end of  for ($j=0; etc...
    
# return the value:
    return $peptide;
}


sub Usage( ; $ )
{
	my $message = $_[0] || '';

	print << 'EOU';

Gene-by-gene Diversity

Description: Examines the diversity of a set of genomes over any number of tagged loci. 
	  
Usage:
GUAVA.pl

Command line options:
-table: Table as tab-separated, with loci and isolates as rows and columns, or reverse if "transpose" option used.
-transpose: "Yes" if the table has loci as columns and isolates as rows, "No" otherwise. 
-FASTA: Directory with FASTA files where "locusname.FAS" is the file name, like BACT000001.FAS
-dbname: If grabbing FASTA sequences from PubMLST, the name of the database. 
-align: "Yes" if alignments are wanted for porportion of variable sites (slow process). "No" otherwise.
-out: Directory where all results will be saved.
-dboption: 1 if complete FASTA directory exists, 2 if creating the same, 3 if pulling just what is needed
-dup: Value for the mininum frequency for allele appearance to be included. Default is "1".
-locuscat: Can point to a file with categories for loci.
-mafft: Can be used to add more parameters to MAFFT. 
EOU

print "\n\nThis script has quit because of a fatal error, namely: ' $message '\n\n";
print "Arguments entered by command line were: " . join (", ", @ARGV) . "\n";
}
