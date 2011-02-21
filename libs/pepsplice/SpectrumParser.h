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

#ifndef SPECTRUMPARSER_H_
#define SPECTRUMPARSER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Spectrum.h"
#include "Spectra.h"
#include "DnaAA.h"
#include "Services.h"
#include "Peak.h"
#include "HotSpectra.h"

using namespace std;

namespace Pepsplice{
class SpectrumParser
{
private:

public:
	SpectrumParser(Spectra *s1, Services *se1);
	virtual ~SpectrumParser();

	//int outputlevel;
	//double sf;
	int spectrawritten;
	
	Spectra *spectra1;
	//DnaAA *dnaAA1;
	Services *se1;
	HotSpectra *hotspectra1;

	//void addFile(string *file);
	void parseFiles(vector<string*> files);
	void parseFile(string *file);
	void parseDTA(string *file);
	void parseMGF(string *file);
	double string_to_double(const string & s);


};
}
#endif /*DTAPARSER_H_*/

