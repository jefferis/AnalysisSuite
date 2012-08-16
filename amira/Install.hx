# Amira Script
# Script to store the location of this script directory in a global variable 
# $GJAnalysisSuiteAmiraScriptDir that is set on Amira startup by reading
# the .Amira file in the user's home folder

set myScriptDir ${SCRIPTDIR}
set amiraUser [file join [file normalize ~] ".Amira"]
echo "Amira User file is $amiraUser"

if { [info exists GJAnalysisSuiteAmiraScriptDir] == 0 } {
    set amiraUserExists [file exists $amiraUser]
    set fileId [open $amiraUser "a"]
    if { $amiraUserExists != 1 } {
        # Seems like .Amira needs to source the bundled startup script
        # or Amira will crash
        source $AMIRA_ROOT/share/resources/Amira/Amira.init
        puts $fileId "# Source bundled startup script or Amira crashes"
        puts $fileId "# You may need to point this to a different script"
        puts $fileId "source $AMIRA_ROOT/share/resources/Amira/Amira.init"
    }
    puts $fileId "# Set the path to GJAnalysisSuiteAmiraScriptDir"
    puts $fileId "# see https://github.com/jefferis/AnalysisSuite/tree/master/amira for details"
    puts $fileId "set GJAnalysisSuiteAmiraScriptDir $SCRIPTDIR"
    puts $fileId "echo \"GJAnalysisSuiteAmiraScriptDir is set to \$GJAnalysisSuiteAmiraScriptDir\""
    close $fileId
} else {
    # we have already set this
    echo "GJAnalysisSuiteAmiraScriptDir is already set to $GJAnalysisSuiteAmiraScriptDir"
    echo "edit $amiraUser file to change this"
}
