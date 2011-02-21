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

#ifndef TUPLES_H_
#define TUPLES_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

#include "DnaAA.h"
#include "Protein.h"
#include "Scoring.h"
#include "Services.h"
#include "Tuple.h"
#include "TupleComparator.h"


//class Scoring;

using namespace std;
namespace Pepsplice{
class Tuples
{
private:
	vector<double> tuples1;

public:
	Tuples(Services *se0, DnaAA *d1);
	virtual ~Tuples();
	
	//int outputlevel;
	long buffersize;

	double starttrigger, endtrigger, interval;
	double starttuples;
	double starttime;
	double bufferincrease, lastbufferincrease;
	double tuples_per_second, last_tuples_per_second;
	double performanceincrease, lastperformanceincrease;
	
	//Protein *randomprotein;
	Scoring *scoring1;
	Services *se1;
	DnaAA *dnaAA1;

    vector<Tuple*> *tuples;
	list<Tuple*> *oldtuples11; //enumeration
	list<Tuple*> *oldtuples12; //enumeration2
    list<pair<Tuple*, bool> > *oldtuples13; //enumeration3
    list<pair<Tuple*, bool> > *oldtuples14; //enumeration4
    list<Tuple*> *oldtuples31; //terminal
    list<Tuple*> *oldtuples32; //terminal2

    list<pair<Tuple*,int> > *oldtuples41; //internal
    list<pair<Tuple*,int> > *oldtuples42; //internal2
    //internal2...


	void setScoring(Scoring *sc1);
	void addTuple(Tuple* tup);
	void addTuple2(Tuple* tup);
	
    //void addTupleDeleted(Tuple* tup);

    void enumerateSNPs(Tuple *tup);
    
    void enumerateTerminalModifications(Tuple *tup);
    
	void enumerateInternalModifications(Tuple *tup, int pos);
	

    
    void setPenalty(Tuple *tup);
	void addPM(double pm);
	void sortTuples();
	void forwardTuples();
	void adaptTupleBuffer();
	Tuple *dummyTuple(Tuple *tup);
	

	
};
}
#endif /*TUPLES_H_*/

