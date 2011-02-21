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

#ifndef UNREACHABLEBINS_H_
#define UNREACHABLEBINS_H_

#include <iostream>
#include "Services.h"

using namespace std;
namespace Pepsplice{
class UnreachableBins
{

public:
	
	UnreachableBins(Services *se0);
	virtual ~UnreachableBins();
	
	void initialize();
	int estimateK(bool *binaryB, int parentmassMHdisc, int trypticends);
	void estimateNK(bool *binaryB, int parentmassMHdisc, int trypticends, int *N, int *K);
	bool isBinReachable(int bin, int trypticends);
	int getTotal(int trypticends);	
	
private:
	
	bool *lowmassbins_tryptic;
	bool *lowmassbins_nontryptic;
	bool *highmassbins_tryptic;
	bool *highmassbins_nontryptic;
	bool *Nspectrum;
	int total_tryptic;
	int total_nontryptic;
	bool initialized;
	
	Services *se1;
	int fillBins(bool tryptic, bool* unreachablebins);
	int fillBins2(bool tryptic, bool* unreachablebins);
	int smearBins(bool *sourcebins, bool *targetbins, int length, int negtol, int postol);

};
}
#endif /*UNREACHABLEBINS_H_*/
