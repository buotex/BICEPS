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

#ifndef __PEPSPLICE_H__
#define __PEPSPLICE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

//#include <gsl/gsl_mode.h>
//#include <gsl/randist/gsl_randist.h>

#include "Services.h"
#include "Chromosomes.h"
#include "ProteinParser.h"
#include "SpectrumParser.h"
#include "SlidingWindow.h"
#include "Tuples.h"
#include "Spectra.h"
#include "Scoring.h"
#include "HotSpots.h"
#include "HotSpectra.h"
#include <vector> //Added by BX, needed for max_penalties and sequences_array.
#include <stdexcept>

namespace Pepsplice{
	///the pepsplice main function, it should use the pepsplice_options
	///given by the ets main function and fill the PepspliceResult vector
	
	void pepsplice_func(int argc, vector<string> & pepsplice_options, float penalty_mutation, std::vector<float> & max_penalties, std::vector<PepspliceResult>& results, const map<string,string> & currentfasta);
}

#endif
