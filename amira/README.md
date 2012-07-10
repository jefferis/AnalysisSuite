Installation
============

  * Start a clean Amira session
  * Run Install.hx script in this directory (easiest to do this by drag and drop)
  * Restart Amira
  * Check that you see a startup message like `GJAnalysisSuiteAmiraScriptDir is set to /GD/projects/AnalysisSuite/amira`
  * On multi-user machines, each user must do this.

It will store the path to this directory in a `~/.Amira` file in your home folder. Code in this file is executed when you start Amira. Henceforth the variable `$GJAnalysisSuiteAmiraScriptDir` will be available to your scripts.

Use
===
If you want to use e.g. ShowClones.scro in a script you can do something like

<code>
    load "$GJAnalysisSuiteAmiraScriptDir/masterObject.scro"
    # or (renaming the script object)
    [load "$GJAnalysisSuiteAmiraScriptDir/masterObject.scro"] setLabel IRPNs
</code>