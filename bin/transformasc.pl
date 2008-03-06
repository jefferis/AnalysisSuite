#!/usr/bin/perl  
use strict;

# Script to take a Neurolucida .asc file, pick out all the cases
# of lines with 3 or 4 co-ordinates, strip off surrounding
# text, send the bare co-ord list to a function en bloc
# and then put transformed co-ords back into their usual location

require 5.005;
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
my $commandLine=join(" ",@ARGV); # keep a copy of the command line for later
# and tidy it up so that there are no spaces between options and their values
# (if they take one)
$commandLine =~ s/(-\w)\s*/$1/g;
# pass on any commands that begin with icfwan
my $trimmedCommandLine=$commandLine; # make a copy
# ... then remove things starting in -^icwfan
# OR without a dash (ie the filename/dir)
$trimmedCommandLine=~s/(-[^icfawn]{1}|\s[^\-]{1})\S*//g;
$trimmedCommandLine=~s/-/ -/g;  # add spaces back between multiple options

init(); # process command line options
print STDERR "Original command line:",$commandLine,"\n" if $opt{v};

my @infoFile;
# GJ: Hard code this for now!
$opt{i}="/Users/jefferis/projects/PN2/Notes/info.txt";
if ($opt{t}){
	# we are going to apply some kind of transformation
	# so we need a file containing dimensions of corresponding images
	if($opt{i}){
		@infoFile=readInfoFile($opt{i})
	} else {
		if (-f "info.txt"){
			# use info.txt in the current directory
			@infoFile=readInfoFile("info.txt");			
		} else {
			die "If transforming (-t), you must either provide -i option or\n info.txt in current dir\n";						
		}
	}
}

# this global var will store commands to carry out on multiple files
my @commandList;
 
# FIGURE OUT WHAT KIND OF INPUT WE HAVE
# #####################################
# ie see if we're working on STDOUT, a directory tree or a file
if (!defined($ARGV[0])){
	die "no input specified, use $0 -h for help\n";
}

my $inFile=$ARGV[0];
if($inFile eq "-"){
	# ie used STDIN for input
	die ("If using STDIN for input, must specify output file\n") unless $opt{o};	
	&mungeFile("-", $opt{o});	
} elsif(-d $inFile){
	#  if we are processing a directory, normally need a file list
	#  to decide which files to choose and what to do with them
	print STDERR "Trimmed command line:",$trimmedCommandLine,"\n" if $opt{v};
	if ($opt{l}) {
		print STDERR "Reading in command list file $opt{l}\n" if $opt{v};
		@commandList=&readCSV($opt{l});
		print STDERR $#commandList+1," commands in command list\n" if $opt{v};
	}
	# nb it is necessary to convert the directory specification
	# to an absolute path to ensure that the open in &mungeFile
	# works properly during multi directory traversal
	my $InDir=File::Spec->rel2abs($inFile);	
	find(\&handleFind,$InDir);	
} elsif (-f $inFile) {
	print STDERR "Processing $inFile\n" if $opt{v};	
	#my @outFileParts=split(/\.(asc|am)/,$inFile);
	#my $outFile=$outFileParts[0].".t$1";
	my $outFile=$inFile;	
	if ($opt{o}){
		$outFile=$opt{o};		
	} else {
		$outFile=~s/^(.*\.)(asc|am)/$1t$2/;	
	}
	&mungeFile(File::Spec->rel2abs($inFile),$outFile) ;
} else {
	die "no input specified, use $0 -h for help\n";
}


##############################
# MAIN SUBROUTINE, mungeFile #
##############################
# which is run for each file to be processed
##########################
sub mungeFile{
	# open the file
	my ($inFile,$outFile)=@_;

	die "Can't open $inFile \: $!\n"  unless open(INFILE, "<$inFile");
	my @lines = <INFILE>;
	close(INFILE);
	
	# are we being asked to trim the file down to the axon?
	if($opt{n}){
		my @headerLines=grep{/^\s*\;.*$/} @lines;		
		
		my $lines=join("",@lines);		
		#$lines=s/\n\r//;
		$lines=~s/^.*\n(.*?\n\s*\(Axon\).*$)/$1/s;
		@lines=split(/\n/,$lines);
		@lines=(@headerLines,@lines);		
		
		#print STDERR join("\n",@lines);
	}
	
	# extract the co-ords
	my @coords;
	if(&extractCoords(\@lines,\@coords)<1){
		# no co-ords found
		return;
	}
	
	# subtract any deltas and apply scaling
	if ($opt{x}||$opt{X} || $opt{y}||$opt{Y} || $opt{z}||$opt{Z}){
		&scaleAxes(\@coords);
	}

	# apply transformations 
	if ($opt{t}){
		my ($x,$y,)=&getXYZCentre($inFile);
		# split the transforms to carry out into separate elements
		my @transforms=split(/[^lrhv]*/,$opt{t});
		foreach (@transforms){
			&transformCoords((\@coords,$_,$x,$y));
		}
	}
	
	# call command line tool to do further transformations
	if($opt{c}){
		@coords=&callCommand(\@coords);
	}
	
	# write out coords
	# find the name of the output file
	&writeCoords(\@lines,\@coords,$outFile);	
}

########################
# THE GUTS SUBROUTINES #
########################
# which are called by mungeFile to do the real work
# ########################

sub extractCoords{
	# get the co-ords out of the file we are processing
	# expects references to two arrays
	my ($lines,$coords)=@_;
	my ($i,$l)=(0,0);	
	
	# First off determine if this an amira file or a 
	# Neurolucida file
	my $firstLine=@{$lines}[0];	
	
	if($firstLine=~/AmiraMesh/i){
		# we're dealing with an amira file
		my @VertexDefinitions=grep {/^\s*Vertices\s*{.*}\s*={0,1}\s*@\d+\s*$/} @{$lines};		
		#print STDERR "There are ",$#VertexDefinitions+1," VertexDefinitions\n";
		#print STDERR join("",@VertexDefinitions);		
		
		my @CoordDef=grep {/Coordinates/} @VertexDefinitions;		
		if(@CoordDef!=1){
			# ambiguous coordinate definitions
			print STDERR "ambiguous coordinate definitions (",scalar(@CoordDef),") in amiramesh file";			
			return 0;
		} 
		
		my $atNum=0; # amira mesh defines each var by @1, @2 etc.
		if($CoordDef[0]=~/.*\@(\d+)\s*$/) {
			$atNum=$1; #print STDERR "atNum = $atNum\n";	
		}
		my $foundCorrectAt=0;
		foreach (@{$lines}){
			# remove any remaining carriage returns or newlines
			s/\r//g; s/\n//g;
			
			if (/^\s*\@(\d+)/){
				if($1 eq $atNum){
					$foundCorrectAt=1;
					#print STDERR "Found start of @",$atNum," block\n";
					next;					
				} else {
					# if we have already found the correct @ block
					# then we must now be entering a new one, so we're done
					last if ($foundCorrectAt==1);					
					#	print STDERR "Found end of @",$atNum," block\n";					
					
					# ... otherwise just keep looking for correct @block
				}
			}
			# now process candidate lines to get numbers
			if($foundCorrectAt && /([\-+.0-9]+)\s*([\-+.0-9]+)\s*([\-+.0-9]+)/ ) {
				# note the requirement for square brackets around the row x,y,z
				# and the dereferencing of the array ref $coords
				${$coords}[$i++]=[$1,$2,$3];
			}
		}
	} else {
		# a neurolucida file
		foreach(@{$lines}){
			# remove any remaining carriage returns or newlines
			s/\r//g; s/\n//g;
			# nb match to start of line to ensure not a comment			
			#if(/$simpleLineRegex/){
			if(/^[^;\(]*\(\s*([\-+.0-9]+)\s*([\-+.0-9]+)\s*([\-+.0-9]+)(\s*[\-+.0-9]+)*\)/) {
				my $x=$1; my $y=$2; my $z=$3;
				
				# note the requirement for square brackets around the row x,y,z
				# and the dereferencing of the array ref $coords
				${$coords}[$i++]=[$x,$y,$z];
				#print STDERR "( $x,$y,$z )  ";				
			}
		}
	}
	# to let the caller know how many coords we found
	return $i;	
	
}

sub scaleAxes{
# scale the co-ords of the specified axes
	my ($coords)=@_;
	print STDERR "$opt{x} $opt{X} $opt{y} $opt{Y} $opt{z} $opt{Z}\n" if $opt{v};	
	my $dx=0; if ($opt{x}){$dx=$opt{x}}
	my $sx=1; if ($opt{X}){$sx=$opt{X}}
	my $dy=0; if ($opt{y}){$dy=$opt{y}}
	my $sy=1; if ($opt{Y}){$sy=$opt{Y}}
	my $dz=0; if ($opt{z}){$dz=$opt{z}}
	my $sz=1; if ($opt{Z}){$sz=$opt{Z}}
	print STDERR $dx+$sx," ",$dy+$sy," ",$dz+$sz,"\n" if $opt{v};	
	
	# Note that $coords is dereferenced to get the 
	# original array (and not a copy)
	# Also note how in the map the $_ reference is
	# dereferenced as an array to get each element
	if ($opt{d}){
		# if the delta AFTER scale option is set do things
		# in that order 
		@{$coords} = map { [($dx+$sx*@{$_}[0],
			$dy+$sy*@{$_}[1],
			$dz+$sz*@{$_}[2])] } @{$coords};
	} else {
		# rewrote this so that addition comes BEFORE scaling by default
		@{$coords} = map { [  (   $sx*($dx+@{$_}[0]),
			$sy*($dy+@{$_}[1]),
			$sz*($dz+@{$_}[2])   )  ] } @{$coords};		
	}

}

sub transformCoords{
# apply transformation $instruct to coordinates
	my ($coords,$instruct,$x,$y)=(@_);

	# width and height of image
	my ($w,$h)=($x*2,$y*2);	
	
	if ($instruct eq "h"){
		@{$coords} = map { [($w-@{$_}[0],@{$_}[1],@{$_}[2])] } @{$coords};
	}
	if ($instruct eq "v"){
		@{$coords} = map { [(@{$_}[0],$h-@{$_}[1],@{$_}[2])] } @{$coords};
	}
	if ($instruct eq "l"){
		@{$coords} = map { [($x-$y+@{$_}[1],$x+$y-@{$_}[0],@{$_}[2])] } @{$coords};
	}
	if ($instruct eq "r"){
		@{$coords} = map { [($x+$y-@{$_}[1],$y-$x+@{$_}[0],@{$_}[2])] } @{$coords};
	}
}

sub callCommandOld{
# run an external command on the output
# (from man perlipc or:
# http://www.perldoc.com/perl5.6.1/pod/perlipc.html#Bidirectional-Communication-with-Another-Process
	my ($coords)=@_;	
	use FileHandle;
	use IPC::Open2;
	my $pid = open2(*Reader, *Writer, "cat -u -n" );
	print STDERR "length: ",$#{$coords},"\n";
 	foreach (@{$coords}){
 		#note that each element of coords is actually a reference
 		print Writer join(" ",@{$_}),"\n";		
 	}
	close(Writer);
	my @newCoords = <Reader>;
	close(Reader);	
	
	return @newCoords;
}
# Hmm seems like it is better to write to a temporary file
sub callCommand{
	# run an external command on the output
	my ($coords)=@_;	
	# get a temp file
	# from http://perlmonks.thepen.com/128876.html	use Fcntl;
	use POSIX;
	my $tmpname;
	do { $tmpname = tmpnam();
	} until open(TMP,">$tmpname");
	
	# write out coords space-delimited line by line
	foreach (@{$coords}){
		# note that each element of coords is actually a reference
		# to a 3 element array
		print TMP join(" ",@{$_}),"\n";		
	}
	close(TMP);  # but I do I ever delete it? no!
	#print STDERR "$tmpname\n";	
	
	# NOW RUN THE EXTERNAL COMMAND and collects its output
	my $param="";	
	$param=$opt{a} if $opt{a};	# the command line arguments to pass to the tool
	my $cmd=$opt{c}." ".$param." < $tmpname";	
	my @newLines=`$cmd`;	
	print STDERR "read in ",$#newLines+1, " coordinates\n" if $opt{v};
	# now delete the temporary file
	unlink $tmpname;	
	
	# Now make a new coordinate list from that output
	my @newCoords;	
	my $i=0;	
	foreach (@newLines){
		$newCoords[$i++]=[split];		
	}
	return @newCoords;
}

sub writeCoords{
	my ($lines,$coords,$outFile)=@_;
	# nb $coords is a ref to external array @coords, so need to write
	# {$coords} wherever one would write coords
	# when referring to @coords
	
	# open outfile
	die "Can't open $outFile \: $!\n"  unless open(OUTFILE, ">$outFile");

	# now write out new coordinates
	if($opt{w}){
		# output in raw format
		for my $i ( 0 .. $#{$coords} ) {
			# print line number if requested
			print OUTFILE $i+1," " if $opt{b};			
			print OUTFILE join(" ",@{$coords->[$i]}),"\n";
		}
		print STDERR basename($outFile),": $#{$coords} coordinates written out\n";
		return;		
	} 
	# GJ 2005-02-02 FIXME
	# Would now like to make it possible to replace in an amiramesh file
	# Basically amira files are pretty simple since co-ord, width and
	# connectivity are separated 
	
	# replace in original file
	# replace the old co-ords with the newly transformed co-ord array
	# also make a note of what we have done in comments in the header!
	# First off determine if this an amira file or a 
	# Neurolucida file
	my $firstLine=@{$lines}[0];	
	my $i=0;  # this will be the counter used for either file type
	if($firstLine=~/AmiraMesh/){
		# we're dealing with an amira file
		#print STDERR "Recognised an amira mesh file, but don't yet know how to handle it\n";
		my @VertexDefinitions=grep {/^\s*Vertices\s*{.*}\s*={0,1}\s*@\d+\s*$/} @{$lines};		
		#print STDERR "There are ",$#VertexDefinitions+1," VertexDefinitions\n";
		#print STDERR join("",@VertexDefinitions);		
		
		my @CoordDef=grep {/Coordinates/} @VertexDefinitions;		
		if(@CoordDef!=1){
			# ambiguous coordinate definitions
			print STDERR "ambiguous coordinate definitions in amiramesh file";			
			return 0;			
		}
		my $atNum=0; # amira mesh defines each var by @1, @2 etc.
		if($CoordDef[0]=~/.*\@(\d+)\s*$/) {
			$atNum=$1; #print STDERR "atNum = $atNum\n";	
		}
		my $foundCorrectAt=0;
		print STDERR $#{$lines}," lines in output file\n";

		# Print out first line followed by a comment.
		print OUTFILE $firstLine."\n";
		print OUTFILE "# $0 ",$commandLine,"\n";
		
		foreach (@{$lines}[1..$#{$lines}]){
			# remove any remaining carriage returns or newlines
			s/(\r|\n)+//;
			
			if (/^\s*\@(\d+)/){
				if($1 eq $atNum){
					$foundCorrectAt=1;
					#print STDERR "Found start of @",$atNum," block\n";
					print OUTFILE $_,"\n" ;
					next;					
				} else {
					# if we have alread found the correct @ block
					# then we must now be entering a new one
					# so we're done
					$foundCorrectAt=0;					
					
					# else we want to keep going
					#if ($foundCorrectAt==1){
					#	print STDERR "Found end of @",$atNum," block\n";					
					#	last;						
					#	}
				}
			}
			# now process candidate lines to get numbers
			if($foundCorrectAt && /([\-+.0-9]+)\s*([\-+.0-9]+)\s*([\-+.0-9]+)/ ) {
				# note the requirement for square brackets around the row x,y,z
				# and the dereferencing of the array ref $coords
				my @coord=(@{$coords->[$i++]});
				print OUTFILE $coord[0]."  ".$coord[1]."  ".$coord[2]."\n";
			} else {
				print OUTFILE $_,"\n" ;
			}
		} # end of foreach
		# end of processing amira file
	} else {
		# process Neurolucida file
		
		# Keep header lines until we get to a line that is neither blank
		# nor a comment
		my $lineIdx=0;		
		while(${$lines}[$lineIdx]=~/^\s*(;|$)/){
			print OUTFILE ${$lines}[$lineIdx]."\n";
			$lineIdx++;			
		}
		print OUTFILE "; $0 ",$commandLine,"\n";		
		
		#print OUTFILE $firstLine."\n";
		foreach my $line (@{$lines}[$lineIdx..$#{$lines}]){
			$line=~s/\r//g; $line=~s/\n//g;# remove any newlines			
			# ^ stuff ( 3 numbers stuff ) more stuff $ 
 			# nb match to start of line and ensure no comment [^;]
			if($line =~ /^([^;\(]*\(\s*)([\-+.0-9]+)\s+([\-+.0-9]+)\s+([\-+.0-9]+)([^\)]*\)[^\)]*)$/) {
			#if(/$fullLineRegex/){
				# there was a match, so fetch coords, increment index
				# note how we get ith element of array referenced by $coords
				# thie element is also a ref, so we must then deref that!
				my @coord=(@{$coords->[$i++]});
				# and replace line with new version
				my $newline=$1.$coord[0]."   ".$coord[1]."   ".$coord[2]."".$5;
				#print OUTFILE "old:$line\nnew:$newline\n";
				print OUTFILE "$newline\n";
			} else {
				print OUTFILE "$line\n";
			}
		}
	} # finish of if amira / else Neurolucida
	
	if ($i eq @{$coords}) {
		print STDERR basename($outFile),": $i coordinates replaced\n";			
	} else {
		print STDERR "$outFile: $i coordinates replaced but ",($#{$coords})+1," available\n";
	}
}

####################
# HELPER FUNCTIONS #
####################
#ie functions to do small tasks related to major tasks above
sub readCSV{
	# read in the list of asc files to process
	my ($fileName)=@_;
	die "Can't open csv list file $fileName \: $!\n"  unless open(INFILE, "<$fileName");
	local $/="";# deal with nasty returns by setting para mark to nothing
	my $lines = <INFILE>; # reading in whole file as one line
	close(INFILE);
	my @lines= split(/[\n\r]+/,$lines);
	my @commands; my $i=0; my @command;	
	foreach (@lines){
		if(/BrainKey/){ # skip the line containing field names
		    next;		    
		}
		# split by separating commas
		@command=split(/,/);
		#print STDERR "@command";		
		
		# remove the filename at the beginning of the line
		my $commandLine=join(" ", @command[1..$#command]);
		#print STDERR "command:$commandLine\n";                                        #
		# now split by th dashes that make each option                                       
		my @options=split(/(-[a-zA-Z]{1})/,$commandLine);			
		#print STDERR join(",",@options);
		
		my $cleanCommandLine="";		
		my $option;		
		
		for my $i ( 0 .. $#options ) {
		#foreach my (@options){
			my $option=$options[$i];			

			if($option=~/-[txyzXYZa]{1}/){ # a valid option
#				if(/[txyzXYZa]{1}\s*\w+\s*/){ # a valid option
				# the argument of the option is the next array element
				my $arg=$options[++$i];
				if ($arg =~ /\S+/){  # if there is an argument
					$option.=$arg; 		
					$cleanCommandLine.=" ".$option;
				}
			}
		}
		if ($cleanCommandLine ne "") {
			# there's something to do with this file so do it
			# note that the file name is switched to the end
			$commands[$i++]=$cleanCommandLine." ".$command[0];
			print STDERR $commands[$i-1],"\n" if $opt{v};			
		} else {
			#print STDERR "No command:",join("+",@options),"\n";
		}
	}
	return @commands;
}

sub findImageBaseName{
# given a name lies JD3L201.PIC or JD3L2.clm2.asc
# find the basename - JD2L
	# nb first thing is to remove the path
	my $in=basename(@_);	
	if($in =~ /(\w+\d{1,2}[LRTB]{1}).*/){
		return $1;		
	} else {
		return "ERROR:$in";		
	}
}

sub readInfoFile{
	my ($infoFile,)=@_;
	# read in the file containing information about PIC files
  	die "Can't open $infoFile \: $!\n" unless open(INFILE, "<$infoFile");
	# note the use of the record separator
	local $/="-------------------------";	
	# read in all the records
	my @infoFile=<INFILE>;	
	print STDERR "Read $#infoFile records from $infoFile\n" if $opt{v};
	return @infoFile;
}

sub getXYZCentre{
# find the centre in the X,Y,Z volume of the image corresponding to a trace file
	my ($fileName,)=@_;	
	my $baseFileName=&findImageBaseName($fileName);	
	my $rec;
	# initialise to the most commonly used values
	my ($x,$y,$z)=(168.78/2,168.78/2,87/2);
	my $foundFlag=0;	
	foreach $rec (@infoFile){
		if ($rec =~ /$baseFileName/){
			$foundFlag=1;			
			
			my ($dx,$dy,$dz)=(0,0,0);
			if($rec =~/([\d+\-.]+).*?x.*?([\d+\-.]+).*?x.*?([\d+\-.]+).*?/){
				($dx,$dy,$dz)=($1,$2,$3);
				#print STDERR "$dx $dy $dz";	
			}
			my ($w,$h,$d)=(0,0,0);			
			if($rec =~/WIDTH = (\d+)/){
				$w=$1;
			}
			if($rec =~/HEIGHT = (\d+)/){
				$h=$1;
			}
			if($rec =~/NPICS = (\d+)/){
				$d=($1-1); # take one fewer than the number of slices
			}
			#print STDERR "w: ",$w,"h: ",$h,"\n";	
			
			# find the centre - ie half height and width
			($x,$y,$z)=($w*$dx/2,$h*$dy/2,$d*$dz/2);
		}
	}
	
	print STDERR "Did not find record in info file for $fileName\n";	
	
	print STDERR "($x,$y,$z)";	
	
	return ($x,$y,$z);
}
			

#############
#HOUSEKEEPING SUBROUTINES#
#############
# this handles the case where we are traversing a directory tree
sub handleFind{
	# get the file name
	my $FullFoundFile = $File::Find::name;
	# first of all see if it ends in .asc case insensitive)
	my $suffix=".a(m|sc)";  # changed so that default is to search for
	                        # either amiramesh or neurolucida
	if ($opt{s}){
		# a suffix has been specified
		$suffix =$opt{s};		
	}
	$suffix=~s/\./\\./; # replace . with \\. (to end up with \. in the regex)
	if ($FullFoundFile =~ /$suffix$/i){
		#print STDERR "handleFind: potential match",$FullFoundFile,"\n";		
		
		# Figure out what the output filename should be
		my ($outFile,$originalSuffix);
		if (basename($FullFoundFile)=~/^(.*)\.([^\.]*)$/){
			$outFile=$1;# ie filename up to last period
			$originalSuffix=$2;		
			print STDERR "stem = $outFile suffix = $originalSuffix\n";		
		} else {
			print STDERR "Problem parsing filename $FullFoundFile\n";			
			print STDERR "Skipping this file\n";			
			return;			
		}
		
		my $outSuffix = ".t".$originalSuffix;		
		if ($opt{o}){
			$outFile.=$opt{o};
		} else {
			$outFile.=$outSuffix;
		}
		
		# if there is a list file use that to decide what to do
		if ($opt{l}) {
			my $baseName=findImageBaseName($FullFoundFile);		
			foreach my $thisCommand (@commandList){
				if($thisCommand=~s/$baseName\S*/$FullFoundFile/){
					# note that we have replaced basefilename since
					# it is not specific enough
					my $cmd="perl $0 $trimmedCommandLine -o $outFile $thisCommand";
					if ($opt{m}){  # simulate directory traversal
						print STDERR "Would have run $cmd\n";
					} else { # actually run the command!
						my $tmp=`$cmd` ;				
					}
				} else {
					#print STDERR "no match in command: $thisCommand\n";					
				}
			}
		} else {
			# simple version - just go ahead and munge it
			&mungeFile($FullFoundFile,$outFile);		
		}
	}
}

	

sub init()
# copied from: http://www.cs.mcgill.ca/~abatko/computers/programming/perl/howto/getopts
{
	use Getopt::Std;      # to handle command line options
	my $opt_string = 'hHi:t:c:fa:x:y:z:X:Y:Z:vweo:l:s:mbnd';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
	axesHelp() if $opt{H};
}

# Print out usage information
sub usage()
{
	print STDERR << "EOF"; 
Usage: $0 [OPTIONS] <ASCFILE/DIR/->

	-h print this help
	-H print help about axes
	-x,y,z <delta> add delta to one of these axes
	-X,Y,Z <scale> scale one of these axes
	-t<lrhv> simple transformation instructions in XY plane
			 (left, right, horizont, vert)
	-d apply delta AFTER scale (default is BEFORE)
  
	-n trim file down to axoN only
		 
	-i <TXTFILE> specify info file containing information about pic file(s)
	-l <CSVFILE> specify list of brain names and conversion paramaters
	-s <SUFFIX> specify the suffix of files to check when doing
		directory traversal - default is .a(m|sc) ie both file types
	-m siMulate what dir traversal will do

	-c name of conversion program
	-f pass filename to conversion program
	-a arguments to pass to conversion program (in quotes ''?)

	-o output file name (otherwise -> .tasc/.tam)
	   (use -o - to send output to STDOUT)
	-w produce output in raw format rather than ASC format
	-b add line numBers to output
	-v verbose mode

Script to take a Neurolucida .asc file, pick out all the cases
of lines with 3 or 4 co-ordinates, strip off surrounding
text, send the bare co-ord list to a function en bloc
and then put transformed co-ords back into their proper location
EOF
  
	exit();  
}
sub axesHelp()
{
	print STDERR << "EOF"; 
Axes:
	The preferred co-ordinate system is as follows.  Looking at a stack of
	frontal sections through a fly brain
	(eg confocal stack of a brain mounted posterior face down on a slide):
	0,0,0 is the top left of the most anterior section
	(typically dorsomedial anterior corner)
	then going positive in x means going lateral 
	going positive in y means going ventral
	BUT going negative in z means going posterior
	Hmm actually have redefined Z to be +ve going posterior
EOF
	exit();  
}
