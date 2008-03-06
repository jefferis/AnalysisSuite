#!/usr/bin/perl 

# File to send to parse the header of a biorad PIC file and print out the
# notes to STDOUT

# GJ 031030
# Made changes so that:
# 1 - it can either print a verbose or short (default) output
# 2 - calculates the voxel size, correcting for lens types
#     based on my previous calculations for Biorad_Microscope_Scale
# GJ 031031
# Now uses getopts
# 1 - can process directories
# 2 - recursively
# 3 - 
# ########
# GJ 040428 
# Successfully added gzip opening to handle .pic.gz files

require 5.004;
use File::Basename;   # to extract a filename from a full path
use File::Find;       # to build a list of files
use File::Spec;       # Platform independent way to build paths



use vars qw/ %opt /;  # for command line options - see init()
                      
init(); # process command line options

if(-d $ARGV[0]){
	# nb it is necesary to convert the directory specification
	# to an absolute path to ensure that the open in &readheader
	# works properly during multi directory traversal
	my $InDir=File::Spec->rel2abs($ARGV[0]);	
	find(\&handleFind,$InDir);	
} elsif (-f $ARGV[0]) {
	&readheader(File::Spec->rel2abs($ARGV[0])) ;
} else {
	die usage();
}

sub init()
# copied from: http://www.cs.mcgill.ca/~abatko/computers/programming/perl/howto/getopts
{
	use Getopt::Std;      # to handle command line options
	my $opt_string = 'hvlxzd:n';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
}

# Print out usage information
sub usage()
{
	print STDERR << "EOF"; 
Usage: $0 [OPTIONS] <PICFILE/DIR>

	-h print this help
	-v verbose ouput
	-l print file names with full path
	-z NO Z axis correction (for refractive index mismatch)
	-x eXclude .pic.gz files from directory search
	-n make a NRRD detached header file (.nhdr) for the image(s) 

Parse the header of a BioRad PIC file or (recursively) parse a directory of PIC files and print a summary of the notes to STDOUT
EOF
  
	exit();  
}

# This extracts the file info from the footer
sub handleFind{
	# get the file name
	my $FullFoundFile = $File::Find::name;
	#if it ends in .pic (case insensitive) read the header

	if( $opt{x}) {
		# only look for .pic files
		if ($FullFoundFile =~ /.*\.pic$/i){
			readheader($FullFoundFile);		
		}
	} else {
		# look for .pic and .pic.gz files
		if ($FullFoundFile =~ /.*\.pic(\.gz){0,1}$/i){
			readheader($FullFoundFile);		
		}		
	}
}


# This does the actual reading in of the header (and footer info)
sub readheader{
	my ($InFile)=@_; 
	# The guts
	
	# Handle gzipped files
	$gzFileOpen=0;	# true if gz file is open
	my $gz;	
	
	if (substr($InFile, -3) eq ".gz") {
		# gzipped pic file
		use Compress::Zlib;   # Allow opening of gzipped PIC files
		$gz=gzopen($InFile,"rb");		
		$gzFileOpen=1;		
	} else {
		# plain pic file
		die "Can't open $InFile \: $!\n"  unless open(PICFILE, "<$InFile");
	}
	
	binmode(PICFILE);
	my $readSuccess=0;
	$readSuccess = $gz->gzread( $in, 76) if $gzFileOpen;		
	$readSuccess = read(PICFILE, $in, 76) if not $gzFileOpen;		
	if(!$readSuccess){
		print STDERR "can't read header of $InFile\n" if $opt{v} ; # the header length
		return;			
	}
	
	if ($opt{l}){
		print "FILE = ",$InFile,"\n";	
	} else {
		print "FILE = ",basename($InFile),"\n";		
	}
	
	# See http://www.perldoc.com/perl5.8.0/pod/func/pack.html for details
	# of unpack format strings
	my $HeaderFormat='vvv@8v@54vh4';
	# looks like mag isn't really used any more
	# See http://www.cs.ubc.ca/spider/ladic/text/biorad.txt for details
	# of biorad header
	($nx, $ny, $npics,$byte_format,$lens,$mag) = unpack($HeaderFormat, $in);
	# removed LENS = since the information appears more accurately in notes
	print "WIDTH = $nx\nHEIGHT = $ny\nNPICS = $npics\n";
	#print "WIDTH = $nx\nHEIGHT = $ny\nNPICS = $npics\nLENS = $lens\n";
	
	# Move to end of raw data
	my $datalength=$nx*$ny*$npics;
	if($gzFileOpen) {
		# This is a pain - there is no seek method for gz - just have to decompress
		# the whole thing
		$gz->gzread( $trashData, $datalength);
	} else {
		seek(PICFILE,$datalength+76,SEEK_SET);
	}

	# Now parse the note information
	# 	Bytes   Description     Details
	# 	------------------------------------------------------------------------------
	# 	0-1     Display level of this note
	# 
	# 	2-5     =0 if this is the last note, else there is another note (32 bit integer)
	# 
	# 	10-11   Note type := 1 for live collection note,
	# 					  := 2 for note including file name,
	# 					  := 3 if note for multiplier file,
	# 					  := 4, 5, etc.,; additional descriptive notes
	# 
	# 	16-95   Text of note (80 bytes)

	my $notes='';  # var to store all of the notes in the footer
	my $morenotes=1;	
	my $level;
	my $endOfFile=0;	
	until ($endOfFile || $morenotes==0){
		# Read in the 16 chars at the head of each line
		if($gzFileOpen) {
			$gz->gzread($lineheader,16);
		} else {
			read(PICFILE,$lineheader,16);
		}
		# get a short and a long (2 and 4 bytes respectively)
		($level,$morenotes) = unpack('vN', $lineheader);

		# seem to be 80 chars of real stuff
		if($gzFileOpen) {
			$gz->gzread($thisline,80);
			$endOfFile=$gz->gzeof();			
		} else {
			read(PICFILE,$thisline,80);		
			$endOfFile=eof(PICFILE);			
		}
		
		
		
		# remove all the extra null characters since they seem to confuse
		#  replacing with an arbitrary invisble character
		$thisline=~ s/[\0]/\11/gs;
		# at the first non word or punctuation character, remove
		# all subsequent characters
		# Problem which took me an hour to solve! 
		# Need the s switch at the end so that all 80 chars are
		# treated as a single line since \0 occurs in the data and is otherwise
		# considered an end of line character!
		#$thisline=~ s/[\0]/\01/g;
		$thisline=~ s/[^\40-\176^].*$//s;
		# (apparently those \ indicate octal ASCII chars)
		$notes.="$thisline\n";
	}
	if($gzFileOpen) {
		$gz->gzclose;		
		
	} else {
		close(PICFILE);
	}

	# print out all or a selection of the notes
	if ($opt{v}){
		print $notes;
		print "---\n";		
	}
	@ToMatch=("PIXEL_BIT_DEPTH","LENS_MAGNIFICATION","INFO_OBJECTIVE_NAME","INFO_OBJECTIVE_ZOOM");		
	#$ToMatch="INFO_OBJECTIVE_NAME|INFO_OBJECTIVE_MAGNIFICATION|INFO_OBJECTIVE_ZOOM";		
	foreach(@ToMatch){
		if($notes=~/($_.*)/){
			print "$1\n";
		}
	}
	%LensFactorHash=("20x Dry",0.978406268,"40x Dry",0.509765625/0.5,
	  "40x Oil", 0.494472656/0.5,"100x Oil",0.2001779/0.2);			

	$notes=~ /INFO_OBJECTIVE_NAME = ([\w ]+).*/;		
	$Lens=$1;

	# Some lenses do not have exactly the expected magnification
	# This is a correction factor for that
	$LensScale=1; # in case we can't find the lens info			
	if ($LensFactorHash{$Lens}){
		$LensScale=$LensFactorHash{$Lens};	
	}
	
	#$PixelWidth= $LensFactorHash{$Lens}*$BoxFactorHash{$nx."x".$ny}/$Zoom;		

	my $ZCorr=1;		
	# Get the Z axis correction - should divide the voxel depth by this
	if(!$opt{z} && $notes=~/Z_CORRECT_FACTOR = ([\d.]+).*/){
		$ZCorr=$1 if $Zcorr!=0;
	}
			
	$notes=~/AXIS_4\s+[0-9.]+\s+[\-+0-9e.]+\s+([\-+0-9e.]+)\W+(\w+)/;		
	$PixelDepth=$1*1.0/$ZCorr;
	$DepthUnits=$2;	

	$notes=~/AXIS_2\s+[0-9.]+\s+[\-+0-9e.]+\s+([\-+0-9e.]+)\W+(\w+)/;		
	$PixelWidth=int($1*$LensScale*1e6)/1e6;
	$WidthUnits=$2;
	$notes=~/AXIS_3\s+[0-9.]+\s+[\-+0-9e.]+\s+([\-+0-9e.]+)\W+(\w+)/;		
	$PixelHeight=int($1*$LensScale*1e6)/1e6;
	$HeightUnits=$2;
	print("Voxel Size (WxHxD) = $PixelWidth $WidthUnits x $PixelHeight $HeightUnits x $PixelDepth $DepthUnits\n","-"x25,"\n");		
	if($opt{n}){
		# make Nrdd header
		my $nhdrFile=$InFile;		
		$nhdrFile.=".nhdr";
		die "Can't open $nhdrFile \: $!\n"  unless open(NHDRFILE, ">$nhdrFile");
		print "Writing $nhdrFile\n";		
		
		print NHDRFILE "NRRD0004\nencoding: ".($gzFileOpen?"gzip":"raw")."\nendian: little\n";		
		print NHDRFILE "dimension: 3\nsizes: $nx $ny $npics\ntype: ".($byte_format?"uint8":"uint16")."\n";		
		print NHDRFILE "spacings: $PixelWidth $PixelHeight $PixelDepth\n";		
		print NHDRFILE "units: $WidthUnits $HeightUnits $DepthUnits\n";		
		print NHDRFILE "byte skip: 76\ndata file: ".($opt{l}?$InFile:basename($InFile))."\n\n";		
	}
}

