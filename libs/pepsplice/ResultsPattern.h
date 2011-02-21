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

#ifndef RESULTSPATTERN_H_
#define RESULTSPATTERN_H_

#include "Spectrum.h"
#include "Services.h"
#include "Distribution.h"
#include "Protein.h"
#include "Match.h"
#include "BestMatches.h"

using namespace std;
namespace Pepsplice{
class ResultsPattern
{
private:
	//series pattern statistics
	float *bions;
	float *yions;
	int *bionsint;
	int *yionsint;
	int *expspec;
	Services *se1;

public:
	ResultsPattern(Services *se0);
	virtual ~ResultsPattern();

	//series pattern statistics
	void calculateBYIons(unsigned char *aasequence, int len);
	void calculateBYIonsInt(unsigned char *aasequence, int len);
	void calcMassDeviationAndSeriesPattern(vector<BestMatches*> b, int i);
	void getPatternAroundMatch(Match *mt, Distribution *d);
	
	//write theoretical spectra
	void writeBestTheoreticalSpectra(vector<BestMatches*> b);
	void writeTheoreticalSpectrum(Match *mt, ofstream outfile);

	
};
}
#endif /*RESULTSPATTERN_H_*/
