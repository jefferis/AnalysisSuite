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

INSTALLATION
============
R
-

  * Download the source code to a suitable location on your hard drive
    * you can either get a [zip snapshot of the code](https://github.com/jefferis/AnalysisSuite/zipball/master)
    * or git clone the repository (strongly recommended if you plan to use this code regularly)
  * Start R
  * source `/path/to/AnalysisSuite/R/Code/Startup.R`
  * See [our wiki](http://flybrain.mrc-lmb.cam.ac.uk/dokuwiki/doku.php?id=warping_manual:start) for further details.

Note that as of commit 0879feb91c (September 11, 2011), it is no longer necessary to edit Startup.R to indicate its location. However the code that identifies the location of the AnalysisSuite source code depends on a recent version of R (>2.12 definitely works).

Amira
-----
  * see instructions in the readme in amira subdirectory
  * Essentially run Install.hx script