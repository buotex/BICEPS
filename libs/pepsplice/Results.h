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

#ifndef RESULTS_H_
#define RESULTS_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <set>

#include "Spectrum.h"
#include "Services.h"
#include "Distribution.h"
#include "Protein.h"
#include "Match.h"
#include "BestMatches.h"
#include "PeptideConfidenceComparator.h"
#include "ProteinConfidenceComparator.h"
#include "ResultsPattern.h"

using namespace std;
namespace Pepsplice{
class Results
{
private:
	Services *se1;

	ResultsPattern *resultspattern1;
	double sf; //scaling_factor
public:
	Results(Services *se0);
	virtual ~Results();
	
	void writeFiles(vector<BestMatches*> b);
	
	void writeResultsPeptideWise(vector<BestMatches*> b, int dl);
	void writePepXMLLite(vector<BestMatches*> b, int dl);
	void writeResultsProteinWise(vector<BestMatches*> b, int dl, bool details);
	void calcScoreDistribution(vector<BestMatches*> b, int discriminate_by, bool writedistribution);	
	void calcLocalRecallPrecision(vector<BestMatches*> b, int delta_i, Distribution *deltascore_performances);
	


};

class PepspliceResult
{
  public:

    string Sequence;
    string OrigSequence;
    //unsigned int seqid;
    int n;
    int k;
    float bic;
    float score;
    float penalty;
    bool tool;
    float penalty_max;
    bool mutation;
    std::set<unsigned int> fastaIds;
    //    list<unsigned int> * seqIds;

    PepspliceResult():Sequence(""),OrigSequence(""),n(0),k(0),bic(-99999.0f),score(0.0f),penalty(0.0f){}
    void reset(){
      Sequence ="";
      OrigSequence="";
      //seqid=0;
      n=0;
      k=0;
      bic=-99999.0f;
      score=0.0f;
      penalty=0.0f;
      fastaIds.clear();
    }

};

ostream& operator << (ostream& os, const PepspliceResult& pep);


struct PepspliceResultComparator
{ //TODO the Match version is derived from public binary_function<Match, Match, bool>, use that approach when this breaks.
  bool operator()(const PepspliceResult a, const PepspliceResult b) const
  {return b.bic < a.bic;}
};

} //namespace

#endif /*RESULTS_H_*/
