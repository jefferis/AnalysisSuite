#!/usr/bin/perl  
#use strict;

# Script to read in csv file containing registration scores
# and use it to parse a directory containing amira model files
# producing as output file lists called <glomerulus>filelist.txt
# for each glomerulus
# these lists are of the form
#   relative_path_to_model_file1 score1
#   relative_path_to_model_file2 score2
#   ...

require 5.004;
use File::Basename;   # to extract a filename from a full path
use File::Find;       # to build a list of files
use File::Spec;       # Platform independent way to build paths

use vars qw/ %opt /;  # for command line options - see init()

##################
# ENTRY POINT #
##################
# which sets things up and figures out what kind of input we have
###################

# PROCESS THE COMMAND LINE
# ########################

init(); # process command line options
 
# FIGURE OUT WHAT KIND OF INPUT WE HAVE
# #####################################
# ie see if we're working on STDOUT, a directory tree or a file

if (@ARGV == 2){
	$csvfile=$ARGV[0];
	$indir=$ARGV[1];
	
	print "Info file = $csvfile and Directory to parse = $indir\n" if $opt{v};
	
	# read csv file
	# 1 read first line
	die "Can't open $csvfile \: $!\n"  unless open(INFILE, "<$csvfile");
	my @csvlines = <INFILE>;
	close(INFILE);
	
	@colnames=split(/\"*\t\"*/,$csvlines[0]);
	print @colnames if $opt{v};	
	
	$braincol=$scorecol=$glomeruluscol=0;
	for ($i=0; $i<=$#colnames;$i++){
		$col=$colnames[$i];
		$braincol=$i if($col=~/^brain$/i);
		$scorecol=$i if($col=~/score$/i);
		$glomeruluscol=$i if($col=~/^glomerulus$/i);
	}
	
	# 2 
	print "braincol is $braincol\n";
	print "scorecol is $scorecol\n";
	print "glomeruluscol is $glomeruluscol\n";
	my @brains,@scores,@glomeruli;
	foreach $line (@csvlines[1..$#csvlines]){
		@items=split(/\"*\t\"*/,$line);
		#print join(",",@items),"\n";
		if(@items>2){
			push (@brains, $items[$braincol]);
			push (@scores, $items[$scorecol]);
			push (@glomeruli, $items[$glomeruluscol]);
		}
	}
	#print join(",",@brains),"\n";
	#print join(",",@glomeruli),"\n";

	# read in matching model file paths
	@modelfiles={};
	find({ wanted => \&handleFind, follow => 1 },$indir);
	#print join(", ",@modelfiles);
	print "There are ",scalar(@modelfiles)," modelfiles\n";
	
	# which glomeruli do we have?
	%seen = ();
	foreach $glom (@glomeruli) {
    	push(@uniqueglomeruli, $glom) unless $seen{$glom}++;
	}
	@uniqueglomeruli = sort grep { /\w/ } @uniqueglomeruli;
	
	print "There are ",scalar(@uniqueglomeruli)," unique glomeruli in $csvfile :\n",join(",",@uniqueglomeruli),"\n";
	
	# iterate over these glomeruli
	foreach $glom (@uniqueglomeruli){
		$filename = $glom."_listfile.txt";
		open(OUTFILE, ">$filename") unless $opt{n};
		
		#%glomscores={};		
		my $numModelsInListFile=0;		
		for ($i=0;$i<@brains & $numModelsInListFile<$MaxNeuronsPerFile;$i++){
			next if($glom ne $glomeruli[$i]);
			#print "Matching glom $glom to ",$glomeruli[$i];
			$regex=$brains[$i];
			# trim the image specifier off in case the brain col doesn't
			# include it.
			$regex=~s/([1-4]|LH)$//;			
			#print ": regex is $regex\n";
			@matchingmodelfiles=grep {/[\\\/:]$regex/i} @modelfiles;
			if(@matchingmodelfiles>0){
				$model=$matchingmodelfiles[0];
				print STDERR $model,"\n" if $opt{v};
				if (@matchingmodelfiles>0){
					print STDERR "regex $regex matches ",scalar(@matchingmodelfiles),
						" modelfiles - choosing modelfile $matchingmodelfiles[0]\n";						
				}
				print OUTFILE $model," ",$scores[$i],"\n" unless $opt{n};
				$numModelsInListFile++;
			}
		}
		close(OUTFILE) unless $opt{n};
		unlink($filename) if($numModelsInListFile==0);		
	}
} else {
	usage();
}


# this handles the case where we are traversing a directory tree
sub handleFind{
	# get the file name
	my $FullFoundFile = $File::Find::name;
	# first of all see if it ends in .asc case insensitive)
	my $suffix=".t*am";
	if ($opt{s}){
		# a suffix has been specified
		$suffix =$opt{s};		
	}
	$suffix=~s/\./\\./; # replace . with \\. (to end up with \. in the regex)

	#print STDERR $FullFoundFile ,"\n";	
	
	if ($FullFoundFile =~ /$suffix$/i){
		# check if matches exclude pattern
		if ($opt{x}){
			return if /$opt{x}/;			
		}
		print STDERR "handleFind: potential match",$FullFoundFile,"\n" if $opt{v};		
		push(@modelfiles,$FullFoundFile);
	}
	#print "THere are ",scalar(@modelfiles)," entries in modelfiles inside handleFind";
}

	

sub init()
# copied from: http://www.cs.mcgill.ca/~abatko/computers/programming/perl/howto/getopts
{
	use Getopt::Std;      # to handle command line options
	my $opt_string = 'hnvs:x:m:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
}

# Print out usage information
sub usage()
{
	print STDERR << "EOF"; 
Usage: $0 [OPTIONS] <CSVFILE> <MODELFILEDIR>

	-h print this help
	-n dry run
	-v verbose
	-x a pattern specifying files to exclude
	-s suffix of model files - default is .t*am
	-m max number of files to put in each individual list file
	
	This utility parses a directory (or direcrtory tree) containing a set of
	neuron tracings and compares the names of those tracings with a CSV file 
	containing 3 columns called brain, *score, glomerulus.  The csv file can
	be saved from Excel or generated using xlsfilter.pl.  For each glomerulus
	for which it finds some neuron tracings it will write out a list file.
	This list file contains one line for each tracing consitsing of the relative
	path to the file and the registration score.
	
	The list files can be used as input to the Amira masterObject.scro and 
	subsidiary loadAllFiles.scro scripts which allow convenient visualisation
	of sets of neurons.
EOF
	exit();  
}