This package provides numerous functions for reading/writing and analysing neuroanatomical data in R.  It is an update to accompany the paper:

Cachero, Ostrovsky et al., Sexual Dimorphism in the Fly Brain, Current Biology (2010), doi:10.1016/j.cub.2010.07.045

http://dx.doi.org/10.1016/j.cub.2010.07.045

of software first released for the paper:

"Comprehensive Maps of Drosophila Higher Olfactory Centers: 
Spatially Segregated Fruit and Pheromone Representation"
Cell (2007), doi:10.1016/j.cell.2007.01.040
by Gregory S.X.E. Jefferis*, Christopher J. Potter*
Alexander M. Chan, Elizabeth C. Marin
Torsten Rohlfing, Calvin R. Maurer, Jr., and Liqun Luo

http://dx.doi.org/10.1016/j.cell.2007.01.040

R Install
=========
Make sure you have a recent version of R installed (R>=2.15.1 recommended). Choose from one of the 3 install
methods below and then see **Running AnalysisSuite in R**.
Easy Install
------------
This is the recommended installation method if you just want to try this out. 
This code changes frequently, so for regular use, the Developer install is recommended. 
  * Start R
  * `source("http://flybrain.mrc-lmb.cam.ac.uk/R/AnalysisSuite/WebInstall.R")`
    * You will be asked to choose an install/download location for the source code
    * R unfortunately doesn't have a function to choose a directory, so once you have
      selected your directory you will have to type a fake filename and then click "Save"
  * R will download my code and install a lot of required packages
  * Now go to **Running AnalysisSuite in R**

Manual Install
--------------
You should only need do this if the Easy Install fails
  * Download this repository as [zip file](https://github.com/jefferis/AnalysisSuite/zipball/master)
  * Unzip the zip file in /path/to/some/sensible/directory
  * Start R
  * source `/path/to/some/sensible/directory/AnalysisSuite/R/Install.R`
    * either use the File ... Source File menu and choose the file in the GUI
    * or `source('/path/to/some/sensible/directory/AnalysisSuite/R/Install.R')` at the command line
  * Now go to **Running AnalysisSuite in R**
  
Developer install
-----------------
This is strongly recommended if you plan to use this code on a regular basis

  * Make sure you have a recent version of R installed (R>=2.15.1 recommended)
  * install [git](http://git-scm.com/) if you don't already have it
  * In a shell/terminal session
    * `cd /path/to/some/sensible/directory`
    * `git clone https://github.com/jefferis/AnalysisSuite.git`
  * Start R
  * `source('/path/to/some/sensible/directory/AnalysisSuite/R/Install.R')`
  * Now go to **Running AnalysisSuite in R**

Running AnalysisSuite in R
==========================
  * You can now start AnalysisSuite (in this and future R sessions) by doing 
     `source(options()$gjanalysissuite.startup)`
  * See [our wiki](http://flybrain.mrc-lmb.cam.ac.uk/dokuwiki/doku.php?id=warping_manual:start) for further details.

Moving the location of AnalysiSuite
-----------------------------------
If you move the AnalysisSuite directory on your filesystem just source `AnalysisSuite/R/Install.R` again to update
the path to AnalysisSuite.


Amira Installation
==================
There are a few useful Amira scripts located in the amira subdirectory.
  * see instructions in the readme in amira subdirectory
  * Essentially run Install.hx script