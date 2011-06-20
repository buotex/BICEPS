/*
* Copyright (c) <2007>, <Franz Roos, ETH Zurich>
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the copyright holders nor the names of 
*       its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef SPECTRA_H_
#define SPECTRA_H_
#include <cfloat>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <tuple>
#include "Spectrum.h"
#include "Services.h"
#include "Distribution.h"
#include "Protein.h"
#include "Match.h"
#include "PeptideConfidenceComparator.h"
#include "ProteinConfidenceComparator.h"
#include "UnreachableBins.h"
#include "Results.h"

using namespace std;

//class Scoring; //include did not work otherwise
namespace Pepsplice{
class Spectra
{
private:

	vector<double> parentmasses;
	
public:
	Spectra(Services *se0, DnaAA *dnaAA1);
	virtual ~Spectra();

	//int outputlevel;
	//int identifiedproteins, identifiedspectra;
	

	///holding the parsed spectra from mgfparser
	vector<Spectrum*> spectra;
	
    ///pointer to main services
	Services *se1;

	UnreachableBins *unreachablebins1;
    ///Handles the results which will be written to files/memory
	Results *results1;
    
    const std::vector<std::tuple<unsigned int, std::string, std::string> >* currentfasta;
    vector<string> * oldSequences;
    ///This function parses the spectrum
	void addSpectrum(Spectrum *spectrum1);
    
    ///sort, rearrange in memory, initialize best matches

	void sortAndPrepareSpectra(); 
    
    void writeSpectrumPositionsInMemory(string filename);
	void initializeResults();
    
    ///tests how good the results are,
    ///the pepsplice-main chooses the sequence with the
    ///best bic.
    bool returnBIC(vector<PepspliceResult> & pep);
	void writeResults();
};
}
#endif /*SPECTRA_H_*/
