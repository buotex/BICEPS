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


\section results Output

The output of BICEPS consists of several elements as detailed below.
The documentation is shown on the left, example entries on the right
@brief Output Example
<table border="1">
<tr>
<th>Value</th>
<th>Interpretation</th>
<\tr>
<tr>
<td>spectrum number (in the mgf-file)</td><td>8 </td></tr><tr>
<td>spectrum title (in the mgf-file)</td><td>Title=File: rd_1_070822094221, Sample: Space_For_Rent (sample number 1), Elution: 1.06 min, Period: 1,Cycle(s): 20 (Experiment 1)</td></tr><tr>
<td>sequence as identified by BICEPS (small case letters indicate modifications)</td><td>Sequence: AGAHLQGGAK</td></tr><tr>
<td>original sequence as found in the fasta file</td><td>OrigSequence: AGAHLQGGAK</td></tr><tr>
<td>entry of the fasta-file containing the sequence(multiple hits possible)</td><td>fastaId: 18785</td></tr><tr>
<td>header of that fasta entry </td><td>fastaId: >IPI:IPI00219018.5| SWISS-PROT:P04406| TREMBL:Q5ZEY3; Q5D0F4; Q16768| REFSEQ_NP:NP_002037| ENSEMBL:ENSP00000229239| H-INV:HIT000034167; HIT000036204; HIT000040230; HIT000009799; HIT000030508; HIT000040030; HIT000031710; HIT000040848 Tax_Id=9606 Glyceraldehyde 3-phosphate dehydrogenase</td></tr><tr>
<td>maximum number of possible matching ions</td><td>n: 18</td></tr><tr>
<td>number of matching ions</td><td>k: 15</td></tr><tr>
<td>bic score (combining penalty and score)</td><td>bic: 2.5114</td></tr><tr>
<td>penalty used</td><td>penalty: 0</td></tr><tr>
<td>score for the PSM</td><td>score: 12.3221</td></tr><tr>
<td>Tool describes the tag generation used </td><td>Tool: 1</tr><tr> 
<td>pmax the maximum allowed penalty for this search </td><td>pmax: 2.31 </td></tr><tr>
<td>whether the tag itself was allowed to carry a mutation or not </td><td>mutation: 0</td></tr><tr>
<td>indicating the component of the mixture model as described in the curve FDP tool</td><td>Label: 2</td></tr><tr>
<td>confidence gives the curveFDP estimate </td><td>Confidence: 0.25891
</tr>
</table>
*/
