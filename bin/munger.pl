#!/usr/bin/perl 

# munger.pl
# script to process all of the images in the specified directory tree
# 1) Run affine transform
# 2) Run Warp transform
# 3) Reformat

# v1.01 041126 to allow inverse consistent switch and to change 
#       default ref brain for HPCF and change bin directory
# v1.02 modified 050111 to allow
#       1) conserver images dir hierarchy in ouput dirs
#       2) fix problem with cwd preventing from 
# v1.03 050124 added version string
# v1.04 050124 fixed paths for vesalius
# v1.05 050124 added support for deep dirs to reformat
# v1.06 050127 fixed luog5 paths and allowed affine reformat
# v1.07 051927 fixed bug in "rootDir" and now copes with .study dirs for
#       both ref or target
# v1.08 060130 added option to specify Registration subdirectory
# v1.09 060207 fixed ability to run on gjz5 remotely (when it is called notched4)
# v1.10 2006-10-13 - added directory locking for warp since on biox cluster
#       multiple processes are started simultaneously and a directory mod time
#       check was not enough
#       Also a great deal of tidying up including reducing host dependence
# v1.11 2006-10-14 - fixed a locking bug for reformat and significantly
#       improved status function (faster, prints dirs with -v option)

require 5.004;
my $version= 1.11;
use vars qw/ %opt /;  # for command line options - see init()
use File::Find;
use File::Basename;
use File::Spec;
use strict;

# Autoflush stdout - so that it always keeps up with STDERR
$|=1;

# hvtawuc:r:s:b:f:E:X:M:C:G:R:
init(); # process command line options
my $energyweight=$opt{E}?$opt{E}:"1e-1";
my $exploration=$opt{X}?$opt{X}:"16";
my $metric=$opt{M}?$opt{M}:"0"; # 0=NMI, 1= MI
my $coarsest=$opt{C}?$opt{C}:"4";
my $gridspacing=$opt{G}?$opt{G}:"40";
my $refine=$opt{R}?$opt{R}:"3";
my $jacobian=$opt{J}?$opt{J}:"0";

# this will be used to name the warp output files
my $warpSuffix="warp_m".$metric."g".$gridspacing."c".$coarsest."e".$energyweight."x".$exploration."r".$refine;
my $icweight=$opt{I}?$opt{I}:"0";

# RECOGNISE SPECIFIC MACHINE TO SET MACHINE-SPECIFIC OPTIONS
# EDIT TO RECOGNISE YOUR MACHINE IF REQUIRED
my $hostName=`hostname`; chomp($hostName); print STDERR "hostname = $hostName";
$hostName = "gjpb" if ($hostName =~/gjpb/i or $hostName eq "jefferis.joh.cam.ac.uk" or $hostName eq "notched1.zoo.cam.ac.uk" or $hostName eq "rechaffe.joh.cam.ac.uk" );
$hostName = "gjz5" if ($hostName =~/(gjz5|notched4)/i);
$hostName = "biox" if ($hostName =~/(bioxcluster|compute)/i);
print STDERR "; short hostname = $hostName\n";
print `echo User path is \$PATH` if $opt{v};

# If we are on vesalius threads defaults to 12 otherwise auto
my $threads="auto";
$threads = 12 if $hostName =~/vesalius/i;
$threads=$opt{T}?$opt{T}:$threads;

# DEFAULT LOCATION OF REFERENCE IMAGE
# EDIT IF REQUIRED
# default
my $referenceImage=$ENV{"HOME"}."/projects/ReferenceImages/average-goodbrains-warp40-5_e1e-2.study";
# ... alternative if that doesn't exist
$referenceImage=$ENV{"HOME"}."/projects/PN2/ReferenceImages/average-goodbrains-warp40-5_e1e-2.study" unless (-d $referenceImage);
# machine specific settings
$referenceImage="/Users/jefferis/projects/ReferenceImages/average/average-goodbrains-warp40-5_e1e-2.study" if $hostName =~/(gsxej2|gjg5|spiracle)/;
$referenceImage="/Users/LuoG5/Greg/PN2/ReferenceImages/average/average-goodbrains-warp40-5_e1e-2.study" if $hostName eq "Bio-LL13.Stanford.EDU";
$referenceImage="/u/jefferis/torsten/ReferenceImages/average/average-goodbrains-warp40-5_e1e-2.study" if $hostName =~/vesalius/i;
# if supplied on command line, thent that overrides
$referenceImage=$opt{s}?$opt{s}:$referenceImage;

print "Reference brain is $referenceImage\n" if $opt{v};
die "Unable to read reference brain $referenceImage" unless -r($referenceImage);

my $referenceStem="average-goodbrains";
if($opt{s}){
	$referenceStem=$opt{s};	
	# remove any terminal slashes 
	# (which would cause basename to return empty)
	$referenceStem=~s/^(.*)\/$/$1/; 
	$referenceStem=basename($referenceStem);

	# if this is a reformatted brain then change the underscore to a hyphen 
	# between brain id and warp/9d0f
	$referenceStem=~s/[_](warp|9dof)/-$1/; 

	# remove up to first dot or underscore
	$referenceStem =~ s/^([^._]*)[._].*$/$1/;
}
print "Reference brain stem is $referenceStem\n" if $opt{v};

# SET UP LOCATION OF WARPING TOOLS
# EDIT IF REQUIRED
my $binDir=$ENV{"HOME"}."/bin/";
$binDir="/u/rohlfing/projects/FlyBrain/bin/" if $hostName =~/vesalius/i;
$binDir=$opt{b}?$opt{b}:$binDir;
my $warpCommand=File::Spec->catdir($binDir,"warp");
my $affCommand=File::Spec->catdir($binDir,"registration");
my $reformatCommand=File::Spec->catdir($binDir,"reformat");

my $regChannels=$opt{c}?$opt{c}:"01";
my $reformatChannels=$opt{r}?$opt{r}:"01";
my $reformatLevels=($opt{l} ne "")?$opt{l}:"5";
my $referenceImageChannel=$opt{f}?$opt{f}:"01";

my $regRoot=$opt{d}?$opt{d}:"Registration";
my $imageRoot="images";
my $reformatRoot="reformatted";
my $rootDir="";

# Set default lock message
#my $lockmessage=$opt{k}?$opt{k}:"";
my $lockmessage=$opt{k}?$opt{k}:getpgrp(0);
print "JOB ID = $lockmessage\n";

my $affineTotal=0;	
my $warpTotal=0;	
my $reformatTotal=0;	

my %found;	# hash to store filenames for status function

if ($opt{p}){
	print "Generating script file $opt{p}" if $opt{v};	
	open SCRIPT, "> $opt{p}";	
}

# Make a note of ARGV0 - the target file or dir
my $inputFileSpec=$ARGV[0];
# remove any terminal slashes 
$inputFileSpec=~s/^(.*)\/$/$1/; 

# NB There was a daft error in which $rootDir was local to the if
# clauses and therefore not found in munge
chomp(my $rootDir=`pwd`);
# so that we can do munger.pl . instead of munger.pl `pwd`
$inputFileSpec=$rootDir if $inputFileSpec eq ".";

if(-f $inputFileSpec || $inputFileSpec=~m/\.study/ ) {
	# find the root dir even if we specifed images dir or subdir
	# $rootDir=findRootDir($inputFileSpec);	

	# 2005-10-18 Actually would rather just use current as root dir
	# that isn't too much of a hardship and if the image 
	# is off in some other location ...
	
	print "Root directory is $rootDir\n";	
	# Hmm not sure that this will work
	&munge($inputFileSpec) ;
	
} elsif(-d $inputFileSpec){
	# find the root dir even if we specifed images dir or subdir
	$rootDir=findRootDir($inputFileSpec);
	print "Changing to root directory: $rootDir\n";	
	chdir($rootDir);	
	
	if($opt{u}){
		status();
		exit 0;
	}
	
	# GJ 130105 - I think I prefer $imageRoot to $rootDir
	# ie only look for images in specifed images dir.
	# which will be the images subdir of $rootDir
	# nb follow=1 implies that we will follow symlinks
	# GJ 2006-10-22 - I want to be able to specify a subdir of image
	# dir and restrict action to that
	if ($inputFileSpec=~/$imageRoot/){
		$imageRoot=$inputFileSpec;	
		print "Setting image root to: ",$imageRoot,"\n";
	} else {
		print "image root is: ",$imageRoot,"\n";
	}
	
	find({ wanted => \&handleFind, follow => 1 },$imageRoot);	
	print "\nRan $affineTotal affine registrations\n";	
	print "Ran $warpTotal warp registrations\n";	
	print "Reformatted $reformatTotal images\n";
} else {
	die usage();
}

sub findRootDir {
	# returns the root directory for a working tree
	# by looking for the dir which has an images subdirectory
	# or returning the original path if no luck
	my $fullpath=shift;	

	# nb it is necesary to convert the directory specification
	# to an absolute path to ensure that the open in &readheader
	# works properly during multi directory traversal
	# since vesalius doesn't have an up to date File::Spec module
	# have to make my own
	if(!File::Spec->file_name_is_absolute($fullpath)){
		my $curDir=cwd();
		$rootDir=File::Spec->catdir($curDir,File::Spec->canonpath($fullpath));
	} else {
		$rootDir=File::Spec->canonpath($fullpath);
	}

	my $partialpath=$rootDir;
	my $sysroot=File::Spec->rootdir();
	while ($partialpath ne $sysroot){
		# Chop off the last directory in the path we are checking
		$partialpath=dirname($partialpath);
		# check if we have a dir called images
		last if (-d File::Spec->catdir($partialpath,"images"));
	}
	# if we have a valid partial path, return that, else return what we were given
	return ($partialpath eq File::Spec->rootdir()?$fullpath:$partialpath);	
}

sub findRelPathToImgDir {
	# returns the relative path of an image filepath to the images directory
	# this just involves removing the first and last elements
	# if we assume that input paths are relative to root dir
	my $filepath=shift;	
	
	# the quick way, but not all systems have File::Spec->abs2rel
	# my($volume,$directories,$file)= File::Spec->splitpath(File::Spec->abs2rel($filepath,$imageRoot));	
	#my($volume2,$directories2,$file2)= File::Spec->splitpath($filepath)[1];	
	
	# doesn't work on vesalius
	#my $rdirectories2 = (File::Spec->splitpath(File::Spec->canonpath($filepath)))[1];
#	my @rdirs = File::Spec->splitdir($rdirectories2);	
#	my @rdirs = File::Spec->splitdir((File::Spec->splitpath($filepath))[1]);	
#	print STDERR "\n",File::Spec->catdir(@rdirs),"\n";
	my ($volume,$directories,$file) = File::Spec->splitpath( $filepath );
	my @dirs = File::Spec->splitdir($directories);	
	# check this works if @dirs has one element
	# not clever, just removes the first and last element
	@dirs=@dirs[1..($#dirs-1)];
	my $dirs=File::Spec->catdir(@dirs);
#	print STDERR "$dirs\n";	
	return ($dirs);	
}
 

sub handleFind {
	# check if file ends in .pic or .pic.gz case insensitive
	if ($File::Find::name =~ /.*\.pic(\.gz){0,1}$/i){
		munge($File::Find::name);		
	}
}

sub munge {
	my $filepath=shift;
	my $filename=basename($filepath);	
	
	# get the brain name
	my $brain=$filename;	
	$brain=~s/(_raw)?(0\d)?\.pic(\.gz)?//i;	
	# the channel of the image
	my $channel=($filename=~/(0\d)\.pic/i)?$1:"";
	
	print  "Found brain name $brain $channel ($filepath)\n" if $opt{v};
	
	# change the working dir to the root dir 
	# (rather than wherever we are in the image dir hierarchy)
	use Cwd;
	my $old_directory = cwd;
	chdir($rootDir) or die "Can't change directory: $!";
	print "New working directory is: $rootDir, Old one was $old_directory \n" if $opt{v};	
	
	# nb by using local we only affect downstream stuff
		 
	if (!$opt{i} && $filename eq basename($referenceImage)){
		print STDERR "Bailing out because target: ",$filename," and ref: ",$referenceImage," are the same" if $opt{v};	
		return 0;
	}
	# run registrations if this file is of the correct channel
	if($channel=="" || $regChannels=~/$channel/){
		runAffine( $filepath,$brain,$channel) if $opt{a};		
		# run the warp transformation
		runWarp($filepath,$brain,$channel) if $opt{w};
	} 
	if ($channel=="" || $reformatChannels=~/$channel/) {
		foreach (split(//,$reformatLevels)){
			runReformat( $filepath,$brain,$channel,$_) if $opt{r};	
		}
	}
	# unset the dir change
	chdir($old_directory) or die "Can't change directory: $!";
}

sub runWarp {
	my ($filepath,$brain,$channel) = @_;	
	my $inlist=File::Spec->catdir($regRoot,"affine",&findRelPathToImgDir($filepath),$referenceStem."_".$brain.$channel."_9dof.list");
	print "inlist = $inlist\n" if $opt{v};	
	
	# new version has relative filenames in output dir depending on input hierarchy
	my $outlist=File::Spec->catdir($regRoot,"warp",&findRelPathToImgDir($filepath),$referenceStem."_".$brain.$channel."_".$warpSuffix.".list");
	print "W: outlist = $outlist\n" if $opt{v};	
	
	my $args="-v --metric $metric --jacobian-weight $jacobian";	
	$args.=" --threads $threads" unless $hostName =~/gjpb/i;			
	$args.=" --spline --fast -e $exploration --grid-spacing $gridspacing ";	
	$args.=" --energy-weight $energyweight --adaptive-fix --refine $refine --coarsest $coarsest";	
	# vesalius doesn't  yet have this parameter
	$args.=" --ic-weight $icweight" if !($hostName =~/vesalius/i);
	$args.=" --output-intermediate" unless $opt{0};	
	# add any extra arguments
	$args.=" ".$opt{W};	
	
	# bail out if infile doesn't exist (or is empty)
	return 0 unless (-s "$inlist/registration.gz") || (-s "$inlist/registration") ;	
	
	my $outfile;
	if($opt{0}){
		# test for the final registration output file ... 
		$outfile = File::Spec->catfile($outlist,"registration");
	} else {
		# ... BUT where possible prefer to test for the level-04 file which is identical
		# to the one that is then saved in the root directory (if all went well)
		# this will avoid an incomplete terminated registration which does 
		# generate a registration file in the root dir from blocking reregistration
		$outfile = File::Spec->catfile($outlist,"level-04.list","registration");		
	}
	# Continue if outdir doesn't exist OR
	# if age of indir > age of outdir & there is a registration
	if( (! -d "$outlist") ){
		myexec("mkdir -p \"$outlist\"") unless $opt{t};
	} else {
		print "outdir exists\n" if $opt{v};
		# there is an output dir ... is there a registration?
		$outfile="${outfile}.gz" if(-f "${outfile}.gz"); # check for a zipped one
		if ( (-s $outfile) && # there is a non-zero registration file already
			(-M $outfile < -M "$inlist/registration") ) { # and it's newer than the input affine reg
			return 0; # then bail out			
		}
	}

	# bail out if somebody else is working on this
	return 0 if (-f "$outlist/registration.lock");
	# now committed, so imediately lock output dir
	makelock("$outlist/registration.lock");	
	
	my $cmd="\"$warpCommand\" $args -o \"$outlist\" \"$inlist\"";
	if($opt{v}){
	    print  "W: Running warp reg with command: $cmd\n";
	} else {
	    print "Warp $brain,";		
	}
	    
	if (!$opt{t}){
		my $cmdnote="# Command run by process with lock ID".$opt{k};
		myexec("echo $cmdnote > \"$outlist/cmd.sh\"") if !$opt{t};
		myexec("echo $cmd >> \"$outlist/cmd.sh\"") if !$opt{t};
	    myexec("$cmd");
		$warpTotal++;
		# the glob here didn't work in perl (only sh)
		#`gzip -f -9 \"${outlist}/registration\" \"${outlist}/level*/registration\"`;
		# so try this:
		# nb opt z indicates that we don't want to gzip; nb -f to over-write
		myexec( "find \"${outlist}\" -name \"registration\" -exec gzip -f -9 {} \";\"") unless $opt{z};		
		myexec("rm \"$outlist/registration.lock\"");
		return $outlist;
	} else {
		return 0;		
	}
}

sub runReformat {
	my ($inputimgfilepath,$brain,$channel,$level) = @_;
	# nb note that the input registration may not be in the same channel as the image to be reformatted
	# therefore use $referenceImageChannel to specify the channel of the input registration
	# (not $channel of current image)
	
	my ($baseinlist,$inlist);
	if($level eq "a"){
		# reformat from affine
		$baseinlist=File::Spec->catdir($regRoot,"affine",&findRelPathToImgDir($inputimgfilepath),$referenceStem."_".$brain.$referenceImageChannel."_"."9dof.list");
		$inlist=$baseinlist;
	} else {
		$baseinlist=File::Spec->catdir($regRoot,"warp",&findRelPathToImgDir($inputimgfilepath),$referenceStem."_".$brain.$referenceImageChannel."_".$warpSuffix.".list");
		$inlist=$baseinlist;
		if($level>=0 && $level<=4){
			# registration will be in a subdir in this case
			$inlist.="/level-0".$level.".list";	
		}
	}	
	
	print "Reformat:inlist would be: $inlist\n" if $opt{v};		
	# bail out if input registration itself doesn't exist
	return 0 if (! -s "$inlist/registration"  && ! -s "$inlist/registration.gz");
	print "inlist exists\n" if $opt{v};
	
	# Construct the outlist - basename gets the name of the file from full path
	my $outlist=basename($baseinlist);  #eg averagegoodbrain_brainame_warp40-5_e1e-1_c4.list
	# remove up to (and including) second underscore and trim terminal .list
	$outlist=~s/^[^_]+_[^_]+_(.*)\.list$/$1/i;
	# nb registration channel may be different from channel of current image
	# which is what we want here
		
	# "a" seems to be trapped by the numeric comparison somehow
	if($level ne "a" && $level>=0 && $level<=4){
		# registration will be in a subdir in this case
		$outlist=$referenceStem."_".$brain.$channel."_".$outlist."_lev".$level;
	} else {
		$outlist=$referenceStem."_".$brain.$channel."_".$outlist;			
	}

	$outlist=File::Spec->catdir($reformatRoot,&findRelPathToImgDir($inputimgfilepath),$outlist);	
	my $outfile = File::Spec->catfile($outlist,"image.bin");
	my $testoutfile=$outfile;	
	
	print "outlist is: $outlist\n" if $opt{v};

	# nb -M gives time since last modification of file
	# it's in days but it's a float with lots of digits	
	# Bail out if we have already reformatted 
	if ( ! -d $outlist){
		# no output dir, so make one and continue
		myexec("mkdir -p \"$outlist\"") unless $opt{t};
	} else {
		print "outdir exists\n" if $opt{v};
		# there is an output dir ... is there an image file?
		# check for a zipped one
		$testoutfile="${outfile}.gz" if(-f "${outfile}.gz");		
		if ( (-f $testoutfile) && # there is an image file already
			(-M $testoutfile < -M $inlist) && # and it's newer than the input reg
			(-M $testoutfile < -M $inputimgfilepath) ) { # and newer than input img
			return 0; # then bail out			
		}
	}
	
	# bail out if somebody else is working on this
	return 0 if (-f "${outlist}.lock");
	# now committed, so immediately lock this directory
	makelock("${outlist}.lock");
	
	# make command 
	my $args="-v --set-null 0";	# makes null pixels black instead of white
	my $cmd="\'$reformatCommand\' $args -o RAW3D:${outfile} --study0 $referenceImage --study1 ${inputimgfilepath} ${inlist}";

	if($opt{v}){
		print  "Running reformat with command: $cmd\n";
	} else {
		# print full name of the file being reformatted
		print "Reformat ".basename($inputimgfilepath).",";		
	}
	if(!$opt{t}){
		print myexec("$cmd");	
		$reformatTotal++;
		# note -f forces overwrite of existing gz 
		myexec ("gzip -f -9 \'${outlist}/image.bin\'") unless $opt{z};
		myexec ("rm \'${outlist}.lock\'");
		return $outlist;
	} else {
		return 0;
	}
}

sub runAffine {
	my ($filepath,$brain,$channel) = @_;	
	
	my $args="-i -v --dofs 6 --dofs 9";
	$args.=" --threads $threads" unless $hostName =~/gjpb/i;
	# add any extra arguments
	$args.=" ".$opt{A};	
	# new version has relative filenames in output dir depending on input hierarchy
	my $outlist=File::Spec->catdir($regRoot,"affine",&findRelPathToImgDir($filepath),$referenceStem."_".$brain.$channel."_9dof.list");

	# Continue if an output file doesn't exist or
	# -s means file exists and has non zero size
	if( ! -s File::Spec->catfile($outlist,"registration") ){
		# no output file, so just continue
	} elsif ( -M "$filepath" > -M File::Spec->catfile($outlist,"registration") ) {
		# ok age of indir > age of outdir so no need to rerun
		return 0;		
	}
	
	# bail out if somebody else is working on this
	return 0 if (-f "$outlist/registration.lock");
	# now committed, so immediately make output dir and lock it
	myexec ("mkdir -p \"$outlist\"") unless $opt{t};			
	makelock("$outlist/registration.lock");

	my $cmd="\'$affCommand\' $args -o \'$outlist\' \'$referenceImage\' \'$filepath\'";
	if( $opt{v}){
		print  "A: Running affine reg with command: $cmd\n";
	} else {
		print  "Aff:$brain$channel ";
	}
	# if not in test mode
	if(!$opt{t}){
		# keep a copy of the commandline
		myexec ("echo $cmd > \'$outlist/cmd.sh\'");	
		# run the command
		#print "Actually running cmd\n";		
		#print `$cmd`;	
		myexec ($cmd);	
		#print "Actually finished cmd\n";		
		myexec ("rm \"$outlist/registration.lock\"");
		$affineTotal++;
		return $outlist;			
	} else {
		return 0;			
	}
}
	
sub status {
	# Displays number of images
	# affine registrations, warp registations, and reformatted
	# images (separated into the two channels)

	# nb follow=1 implies that we will follow symlinks
	print "Searching directory tree ..." if $opt{v};	
	find({ wanted => \&findAllFiles, follow => 1 },$rootDir);	
	print " Finished!\n" if $opt{v};	

	
	my @paths=keys %found;
	my @filenames=values %found;
	my @images=grep /\.pic(\.gz)$/i, @paths;	
	my @channel1images=grep /01\.pic(\.gz)/i, @images;	
	my @channel2images=grep /02\.pic(\.gz)/i, @images;	
  
	print "Total Images: ".scalar(@images)."\n";	
	print "Channel 1 images: ".scalar(@channel1images)."\n";	
	print "Channel 2 images: ".scalar(@channel2images)."\n";
  
	my @affineRegistrations=grep /affine.*9dof\.list$/i, @paths;	
	my @lockedAffineRegistrations=grep /affine.*registration.lock$/i, @paths;
	my @lockedAffineIDs=map {&getidfromlockfile($_)} @lockedAffineRegistrations;	
	
	my @finishedAffineRegistrations=grep /affine.*9dof\.list\/registration$/i, @paths;	
	# make a hash containing the directory name of all finished affines
	my %finished = map { dirname($_) => $_ } @finishedAffineRegistrations;
	# Now go through the array of all registration dirs tossing those that
	# are in the finished hash
	my @unfinishedAffineRegistrations = grep { !exists $finished{$_}  } @affineRegistrations;	
	
	my @warpRegistrations=grep /\/warp\/.*warp[^\/]*\.list$/, @paths;
	my @finishedWarpRegistrations=grep /\/warp\/.*warp[^\/]*\.list\/registration(.gz)$/, @paths;
	my @lockedWarpRegistrations=grep /\/warp\/.*warp[^\/]*\.list\/registration.lock$/i, @paths;	

	# make a hash containing the directory name of all finished warps
	my %finished = map { dirname($_) => $_ } @finishedWarpRegistrations;
	# Now go through all registration dirs tossing those that
	# are in the finished hash
	my @unfinishedWarpRegistrations = grep { !exists $finished{$_}  } @warpRegistrations;	
	
	my @reformattedImages=`find $rootDir/$reformatRoot/ -type d -name \'*.study\'`;
	@channel1images=grep /^[^_]+01_/i, @reformattedImages;	
	@channel2images=grep /^[^_]+02_/i, @reformattedImages;
	
	print "\nAffine registration directories: ".scalar(@affineRegistrations)."\n";
	print "Locked affine registration directories: ".scalar(@lockedAffineRegistrations)."\n";

	if($opt{v}){
		for (my $i=0;$i<$#lockedAffineRegistrations;$i++) {
			print "\t",$lockedAffineRegistrations[$i],"\t",$lockedAffineIDs[$i],"\n"
		}
	}
	
	#print "\t",join("\n\t",sort @lockedAffineRegistrations),"\n" if $opt{v} && @lockedAffineRegistrations;	

	print "Unfinished affine registration directories: ".scalar(@unfinishedAffineRegistrations)."\n";
	print "\t",join("\n\t",sort @unfinishedAffineRegistrations),"\n" if $opt{v} && @unfinishedAffineRegistrations;	

	print "\nWarp registration directories: ".scalar(@warpRegistrations)."\n";
	print "Unfinished warp registration directories: ".scalar(@unfinishedWarpRegistrations)."\n";
	print "\t",join("\n\t",sort @unfinishedWarpRegistrations),"\n" if $opt{v} && @unfinishedWarpRegistrations;	

	print "Locked warp registration directories: ".scalar(@lockedWarpRegistrations)."\n";
#	print "\t",join("\n\t",sort @lockedWarpRegistrations),"\n" if $opt{v} && @lockedWarpRegistrations;	
	if($opt{v}){
		foreach (@lockedWarpRegistrations) {
			print "\t",$_,"\t",getidfromlockfile($_),"\n"
		}
	}

	print "Reformatted image directories: ".scalar(@reformattedImages)."\n";
}

sub findAllFiles {
	# check if file ends in .pic or .pic.gz case insensitive
	$found{$File::Find::name}=$_;
}

# sub selectImageFiles {
# 	# check if file ends in .pic or .pic.gz case insensitive
# 	if ($File::Find::name =~ /.*\.pic(\.gz){0,1}$/i){
# 		$images{$File::Find::name}=1;
# 	}
# }

sub myexec {
	my ($cmd) = @_;	
	if ($opt{p}){
		print SCRIPT $cmd,"\n";
	} else {
		# should get to see output with system
		my $rval = system $cmd;
		if ($? == -1) {
			  print "MYEXEC: CMD = $cmd failed to execute: $!\n";
		  }
		  elsif ($? & 127) {
			  printf "MYEXEC: CMD = $cmd died with signal %d, %s coredump\n",
				  ($? & 127),  ($? & 128) ? 'with' : 'without';
		  }
		  else {
			  printf "MYEXEC:  CMD = $cmd exited with value %d\n", $? >> 8;
		  }
						  
		print STDERR "MYEXEC-DEBUG: RVAL = $rval: CMD = $cmd\n" if($opt{g});
		return $rval;		
	}
}

sub makelock {
	my ($lockfile)=@_;
	# bail if we are in test mode
	return if $opt{t};	
	if($lockmessage){
		open(FH,"> $lockfile") or die ("Can't make lockfile at $lockfile: $!\n");
		print FH $lockmessage;
		close FH;		
	} else {
		myexec("touch \"$lockfile\"");
	}	
}

sub getidfromlockfile {
	my ($file)=@_;	
	open (FH,"$file") or return "NULL";
	my $line=<FH>;
	close(FH);
	return($line);	
}

sub usage {
	print STDOUT << "EOF"; 
Usage: $0 [OPTIONS] <PICFILE/DIR>
Version: $version

	-h print this help
	-v verbose (provide extra feed back at runtime)
	-t test run (no commands are actually run)
	-g debug: prints every command run by myexec and the return value
	-p make a scriPt from commands that would be run 
	   (nb cannot produce commands that depend on earlier commands)
	-u statUs - display number of images, registrations etc
	-z turn gzip off (on by default)
	-k lock message ie contents of lock file (defaults to proces id)

	-a run affine transform
	-w run warp transform
	-c [01|02|..] channels for registration (default 01 or "")
	-r [01|02|..] run reformat on these channels
	-l [a|0|1|2|3|4|5] run reformat on these levels (default 5, a=affine)
	-f [01|02|..] channel of the images used for registration - default is 01

	[nb use -f to specify the channel of the images that were previously used to
	 generate a registration if you now want to reformat a different channel using
	 that registration information]

	-i register brain to itself if possible (default is to skip)
	-0 Don't output intermediate warp registration levels

	-s [file|fileStem] Reference brain (average e-2 by default)
	-b [path] bin directory
	-d [stem] registration subdirectory (default ./Registration)

	-I inverse consistent warp weight (--ic-weight) default 0, try 1e-5
	-E [energy] energy of warp transform (default e-1)
	-X [exploration] (default 16)
	-M [metric] (default 0 (NMI), options 1=MI)
	-C [coarsest] (default 4)
	-G [grid-spacing] (default 40)
	-R [refine] (default 3)
	-J [0 to 1] jacobian-weight volume constraining param (default 0)
	-T [threads] (default auto or 12 on vesalius)
	-A [option] additional options for affine transformation
	-W [option] additional options for warp transformation
  
Munge a BioRad PIC file or (recursively) parse a directory of PIC files
by running the affine and warp registrations and reformatting images
as required.  Final argument must be the images directory or a single image.
EOF
  
	exit();  
}

sub init {
# copied from: http://www.cs.mcgill.ca/~abatko/computers/programming/perl/howto/getopts
	use Getopt::Std;      # to handle command line options
	my $opt_string = 'hvtawuic:r:l:s:b:f:E:X:M:C:G:R:T:J:I:zp:d:k:g0A:W:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if $opt{h};
	
}
