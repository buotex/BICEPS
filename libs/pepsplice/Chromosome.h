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

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Services.h"
#include "DnaAA.h"

using namespace std;
namespace Pepsplice{
//class SlidingWindow;

class Chromosome
{
private:

	
public:
	//Chromosome();
	Chromosome(string f, Services *se0);
	virtual ~Chromosome();
	
	string file;
	string *NTsequence;
	string test2;
	string description_line;
	int fullchromsize;
	bool unloaded;
//	DnaAA dnaAA1;
	Services *se1;
	//start changes BYR
	char getNucleotide(long pos);
	char getReverseNucleotide(long pos);
	void loadFullSequence();
	void unloadFullSequence();
	int getFullChromSize();
	char* getNTSeqForTuple(long len, long ps, long pe, long ss, long se, bool rev);
	//end changes
};
}
#endif /*CHROMOSOME_H_*/
