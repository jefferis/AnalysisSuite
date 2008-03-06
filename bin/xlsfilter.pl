#!/usr/bin/perl 

# xlsfilter.pl
# 2005-02-02
# Script to convert an excel spreadsheet to a tab delimited txt file
# designed to be a filter for use with filemerge
use strict;
use Spreadsheet::ParseExcel;
require 5.004;

use vars qw/ %opt /;  # for command line options - see init()
init(); # process command line options

my $oExcel = new Spreadsheet::ParseExcel;

my $quote=""; $quote = qw/"/ if $opt{q};
my $delim = "\t"; $delim =$opt{d} if $opt{d};
my $suffix = ".txt"; $suffix=$opt{x} if $opt{x};

die "You must provide a filename to $0 to be parsed as an Excel file" unless @ARGV;
my $fileStem=$ARGV[0];
# trim terminal . suffix
$fileStem=~s/\.[^\.]*$//;

my $oBook = $oExcel->Parse($ARGV[0]);
die "Cannot parse Excel file $ARGV[0]\n" if not ( $oBook->{File});

my $outFile=(defined $ARGV[1])?$ARGV[1]:$fileStem.$suffix;
if($outFile eq "-") {
	#OUTFILE=STDERR;	
} else {
	die "could not create output file $outFile" unless open(STDOUT, ">$outFile");	
}

my($iR, $iC, $oWkS, $oWkC);
print STDERR "FILE  :", $oBook->{File} , "\n";
print STDERR "SHEET COUNT :", $oBook->{SheetCount} , "\n";
print STDERR "AUTHOR :", $oBook->{Author} , "\n"
 if defined $oBook->{Author}; 

for(my $iSheet=0; $iSheet < $oBook->{SheetCount} ; $iSheet++)
{
	$oWkS = $oBook->{Worksheet}[$iSheet];
	# skip this sheet if it's not what we're after
	if($opt{s}){
		next if $opt{s} ne $iSheet+1;		
		
		#next if($oWkS->{Name} ne "AllWarpedImageNotes") ; 
	}
	
	# note the sheet name to stderr
	print STDERR "--------- SHEET:", $oWkS->{Name}, "\n";
	# Print to Stdout as well
	#print "--------- SHEET:", $oWkS->{Name}, "\n" unless $opt{s};
	
	for(my $iR = $oWkS->{MinRow} ; defined $oWkS->{MaxRow} && $iR <= $oWkS->{MaxRow} ; $iR++)
	{
		for(my $iC = $oWkS->{MinCol} ;
			defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
			$iC++)
		{
			$oWkC = $oWkS->{Cells}[$iR][$iC];
			#print "( $iR , $iC ) =>", $oWkC->Value, "\n" if($oWkC);
			if($oWkC) {
				# this cell contains something
				print STDOUT $quote,$oWkC->Value,$quote;
			} else {
				# empty cell
				print STDOUT $quote x 2;				
			}
			# if this isn't the last cell in the row print a tab
			if ($iC < $oWkS->{MaxCol}) {
				print STDOUT $delim;				
			}
		}
		print STDOUT "\n";		
	}
}

sub init()
# copied from: http://www.cs.mcgill.ca/~abatko/computers/programming/perl/howto/getopts
{
	use Getopt::Std;      # to handle command line options
	my $opt_string = 'hqd:s:x:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
}

# Print out usage information
sub usage()
{
	print STDERR << "EOF"; 
	Usage: $0 [OPTIONS] <XLSFile> <OutFile>

	-h print this help
	-s number of sheet to parse
	-q quote cell contents
	-d delimiter (default = tab)
	-x suffix (default =.txt)
	
	NB to use as a filter give - (ie stdout) as OutFile
EOF
	exit();  
}
