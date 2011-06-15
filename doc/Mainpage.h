/*!

\mainpage A Novel Error Tolerant Search Strategy for Cross-Species Proteomics 


\section sec_intro Introduction

\c TODO 

\section sec_license Licensing
\c TODO 

\section sec_download Download
The source repository,including packaged windows binaries is available from 
http://hci.iwr.uni-heidelberg.de/MIP/Software/

\section sec_install Installation
Using the binary, just place it in a convenient location, take care that the \c Models
folder is in the same place.


\subsection sec_install_src Building from Source
The following illustrates how to build the examples, tests and
distribution packages. Bundling \c libfbi requires a working CMake build system
(available from http://cmake.org/) and CMake >= 2.6.

With cmake in the system path, the build process is
\verbatim
 Either
 tar xvzf biceps-xxxxxxx.tar.gz
 or
 git clone git://github.com/buote/BICEPS.git

 mkdir biceps
 cd biceps
 ccmake ../biceps
 make && sudo make install
\endverbatim
\verbatim
 Example line:
 biceps --mgf pickedpeaks_no_pp_Helatest_with_none_resultRawList.mgf --fasta ipi.CHICK.v3.54.fasta 
 For other parameters, just execute biceps.
\endverbatim
*/

