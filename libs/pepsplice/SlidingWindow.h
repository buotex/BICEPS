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

#ifndef SLIDINGWINDOW_H_
#define SLIDINGWINDOW_H_

#include <iostream>	//cout
#include <fstream>
#include <deque>
#include <cmath>	//logarithms
#include <cstdlib>	//abs
#include "Timer.h"
#include "DnaAA.h"
#include "Distribution.h"
#include "Services.h"
#include "Scoring.h"
#include "Tuple.h"
#include "Tuples.h"
#include "Chromosome.h"

using namespace std;
namespace Pepsplice{
/* 
 * A sliding window travels along the DNA and enumerates all potential spliced peptides within the window.
 * Corresponding to biological rules, each spliced peptide corresponds to a tuple:
 * 1. prefix start
 * 2. prefix end
 * 3. suffix start
 * 4. suffix end
 * 5. sequence
 * 6. parent mass
 * 
 * - No stop codons must occur within the final peptide
 * - There is a threshold for the maximum gap length
 * 
 * - A DNA parser feeds nucleotides into the sliding window
 * - The sliding window outputs peptides which are then matched against spectra
 * 
 */


class SlidingWindow
{
private:
	//deque<long> deqKR, deqGT, deqAG, deqStop;
	long circlelength, circleextra;
	long plenmax, slenmax;
	long i, imodz, imoda, imodb, imodc, imodd;
	long illegalnucleotides;
	
	double knownprefixmass, knownsuffixmass;
	double monoparentmassMH, peptideresiduemass, aamass_temp;
	float peptideresiduemass_float1, peptideresiduemass_float2;

	long ps, pe, ss, se;
	long ps3, pe3, ss3, se3;
	long ps_old, pe_old, ss_old, se_old;
	long plen, gaplen, slen, totlen;
	long pexcessnt, sexcessnt;


	int outputlevel;
	bool chr_reverse;


	Timer timer1;
	DnaAA dnaAA1;
	Services *se1;
	Distribution *distribution1; // delete afterwards? not if statistics over more than one chromosome are required
	
	Chromosome *activechromosome; // do not delete!

	float *aamass;
	double *massruler;
	long *globalposition;
				
	char *nt;
	char *ntnum;
	char *readingframe;
	char *ntsequence;
	unsigned char *aa;
	unsigned char *aasequence;
	float *aamasssequence;


	long *sincestop;
	long *sinceprefixstart;
	long *sinceprefixend;
	long *sincesuffixstart;
	long *sincesuffixend;
	
	bool *isstop;
	bool *isprefixstart;
	bool *isprefixend;
	bool *issuffixstart;
	bool *issuffixend;
	


	void enumerateSplicedPeptides();
	void enumerateEntirePeptides();
	long cheapMod(long x);
	long cheapUpMod(long x);
	long cheapDownMod(long x);
	long cheap3Mod(long x);
	double getBrokenTripletMass();
	unsigned char getBrokenTripletAA();
	string getAASequence();
	string getAASequenceFast();
	unsigned char* getAASequenceFastChar();
	string getNTSequence();
	void showAAMassSequenceSlow();



public:
	
	bool hot;
	
	Tuples *tuples;
	
	SlidingWindow(Services *se0);
	virtual ~SlidingWindow();
	
	bool addNextNucleotide(char newNucleotide);
	void writeProgress();
	void showParameters();
	void writeCircleContent();
	void setTuples(Tuples *tps);
	void setChromosome(Chromosome *c, bool rev);
	double j, nextout;
};
}
#endif /*SLIDINGWINDOW_H_*/
