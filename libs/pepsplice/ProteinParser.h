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

#ifndef PROTEINPARSER_H_
#define PROTEINPARSER_H_

#include <iostream>
#include <fstream>
#include "Services.h"
#include "Protein.h"
#include "Tuples.h"
#include "Tuple.h"
#include "DnaAA.h"
#include <map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <tuple>
using namespace std;
namespace Pepsplice{
class ProteinParser
{
public:
	ProteinParser(Services *se0, Tuples *tp0);
	virtual ~ProteinParser();
    Tuples *tuples1;
	Protein *currentprotein;
	Services *se1;
	DnaAA dnaAA1;
	int proteincount;
    int oldSequenceCount;
    vector<string> * oldSequences;
    
//	void parseFASTA(const char* file);
    void parseFASTA(const std::vector<std::tuple<unsigned int, std::string, std::string> > & currentfasta);    
	void finishProtein();
	void generateTuples();
	void doAnotherRun(const std::vector<std::tuple<unsigned int, std::string, std::string> > & currentfasta);
	unsigned char* getAASeqChar(string seq);
	
};
}
#endif /*PROTEINPARSER_H_*/
