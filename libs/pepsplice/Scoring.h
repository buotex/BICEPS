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

#ifndef SCORING_H_
#define SCORING_H_

#include <ctime>
#include <iostream>
#include <vector>
#include <math.h>
#include "DnaAA.h"
#include "Spectrum.h"
#include "Spectra.h"
#include "Tuple.h"
#include "Match.h"
#include "BestMatches.h"
#include "Services.h"
#include "Timer.h"
#include "Distribution.h"
#include "Hypergeometric.h"

using namespace std;

class Spectra;
namespace Pepsplice{
class Scoring
{
private:
	int aamaxlen;
	int maxparentmass;
	//int cumulatedpeaksmax;
	double pmtolneg, pmtolpos, pmmin, pmmax;
	double discretization_factor;
	
	//objects
	DnaAA dnaAA1;
	Hypergeometric * hypergeometric1;
	Services *se1;
	Distribution *d1;
	
	//tuples (for spectra see public:)
	vector<Tuple*> *vtuples;
	
	//cache-aware approach, stripwise approach
	void scoreTuplesAmortizedCacheAware2(bool slicealways);
		
	//recursively score tuples vs spectra	
	/*
	void recursivelyScoreTuples(bool onebyone);
	bool recursiveIntersectsArea(int x1, int x2, int y1, int y2);
	void recurse1by1(int x1, int x2, int y1, int y2);
	void scoreTupleVsSpectrum1by1(int sp1, int tp1);
	void recurse1byx(int x1, int x2, int y1, int y2);
	void scoreTupleVsSpectrum1byx(int tp1, int sp1, int sp2);
	*/
	
    //int *spcumpeaks;
	//int *tpcumpeaks;
	
	//ion calculations per tuple
	//start changes BYR, added b/yions for charge 3
	void calculateBYIons(Tuple *tup, int length);
	void calculateBYIonsInt(Tuple *tup, int length);
   	void calculateBYIons3(Tuple *tup, int length);
        void calculateBYIonsInt3(Tuple *tup, int length);
	        
	float *bions   , *yions; 
	int   *bionsint, *yionsint;
	float *bions3   , *yions3;
	int   *bionsint3, *yionsint3;



	int *thspec;
	int *thspec_prefix;
	int *thspec_suffix;
    float *thspecdisc2; 
	
	//chief functions for scoring
	float calculateThSpec(Tuple *tup, bool erase, int *thspec);
	void calculateScore(Tuple *tup, Spectrum *spectrum, int *thspec, float thspec_avgint);
	
	//score 2: Rikos idea, detailed negative score
	void calculateThSpec2(Tuple *tup);
	void setCurrentIon2(int pos, int weight);
	float calculateScore2(Spectrum *spectrum, int seqlen);
		
	//score 3: Rikos idea, uniform negative score; factor 2 faster, one single += per peak only
	float calculateThSpec3(Tuple *tup, bool erase, int *thspec);
	float calculateScore3(Spectrum *spectrum, int *thspec, float thspec_avgint);
	//void setCurrentIon3(int pos, int weight);
	//float thspec3_negweight;
	//bool  thspec3_cleaned;
	
	//scores 4-6: hypergeometric model, Sadygov and Yates 2003
	void calculateThSpec4(Tuple *tup, bool erase, int *thspec);
	float calculateScore4(Tuple *tup, Spectrum *spectrum, int *thspec);
	void calculateThSpec5(Tuple *tup, bool erase, int *thspec);
	float calculateScore5(Tuple *tup, Spectrum *spectrum, int *thspec);
	void calculateThSpec6(Tuple *tup, bool erase, int *thspec);
	float calculateScore6(Tuple *tup, Spectrum *spectrum, int *thspec);
	//calculateThSpec7 is not necessary, b- and y-ions are sufficient
	float calculateScore7(Tuple *tup, Spectrum *spectrum, int *thspec, int *sharedpeakcount);

		
	//debugging
	void checkParentMass(Tuple *tup);
	
public:
	Scoring(Services *se0);
	virtual ~Scoring();

	Spectra *spectra1; //is initialized in main; access needed for progress reports that are initiated by Chromosomes

	void setSpectra(Spectra *sp1);
	void addTuples(vector<Tuple*> *vtuples);
	//end changes BYR
		
};
}
#endif /*SCORING_H_*/
