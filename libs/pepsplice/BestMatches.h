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

#ifndef BESTMATCHES_H_
#define BESTMATCHES_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include "Spectrum.h"
#include "Match.h"

//there are circular dependencies and compilation problems with BestMatches-Match-Spectrum

using namespace std;
namespace Pepsplice{
class Match;
class Spectrum;

class BestMatches
{
private:
	//int i; //spectrum number after sorting
	
public:

	//long length, length2;
	//float lowerBoundScore;
	double pvalue;
	float discriminationscore;
	float precision;
	float recall;
	float fraction;
	float tp;
	float fp;
	
	//values from local analysis
	double thpm;
	float bestscore;
	float deltascore12;
	//float deltascores[22];
	float *bestscores;
	
	bool random;
	int show_hits;
	int tot_best_hits;
	int nb_prot;
	int nb_chromunspliced;
	int nb_chromspliced;
	
	//from discrimination function
	int confidence_bin;
	double confidence;
	
	//helper variables
	int matches;

	BestMatches(Spectrum *ownerspec);
	virtual ~BestMatches();
	Spectrum *ownerSpectrum;
	vector<Match* > bm;
	
	//methods
	//start changes BYR
	void resize();
	void addMatch(Match * mt1);
	Match* getBestMatch();
	void sortBestMatches();
	void analyzeBestHits(); //extracts them and stores them within class BestMatches
	void reset(); //Because we need additional runs, starting with a clean BestMatches
	//end changes
};
}
#endif /*BESTMATCHES_H_*/

