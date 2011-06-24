================================================================================
                                 biceps README
================================================================================

1. Summary
----------


2. Installation
===============
Installation from the Windows binary package should be straightforward.

On Linux, building biceps has a few requirements:
*gcc(>=4.4) 
*CMake(>=2.6) 
*Boost(>=1.40)

Additionally, 
*Doxygen
*dot
are recommended though not required to build the Documentation.

Using these packages, building biceps should work with
cd builds
mkdir biceps && cd biceps
cmake "PATH TO BICEPS SOURCE" e.g. cmake ~/code/biceps
make && sudo make install

Installation Notes
=====
Installing on either windows or linux will result in a binary file biceps, installed to {PREFIX_PATH}/bin, e.g. /usr/local/bin
A configuration directory called .biceps will also be created as a subfolder of your home directory, containing a
* directory  "Models", which is the Model for the pepnovo component (at the moment, the name of the model to be used is
hardcoded)
* file in_AAmodifications.param, which describes the internal modifications or rather mutations pepsplice will try, e.g.
the different modifications of Lysine.



