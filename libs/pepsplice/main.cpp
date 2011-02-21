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

int main( int argc, char * argv[] )
{	
	//argument containers
	string* *argstrings = new string*[argc]; //dynamic array of string pointers
	vector<string*> specfiles;	
	vector<string> params;

	//parse arguments (parameters and spectrum files mixed)
	cout << "\n\nmain.cpp: Reading " << argc - 1 << " arguments\n";
	string paramstring = "";
	for(int i = 1; i < argc; i++){
		argstrings[i] = new string();
		(*argstrings[i]) = argv[i];

		if((*argstrings[i])[0] == '-'){
			//parse parameter
			params.push_back( (*argstrings[i]) ); //parse parameters all at once afterwards
			paramstring += (*argstrings[i]); //for output file names later on
		}else{
			//parse spectrum
			specfiles.push_back( argstrings[i] );
		}
	}
	cout << "\n";
	
	//define THE parameter and service class, available to all classes
	Services *se1 = new Services();
	cout<<se1->masstol_belowDa;
	cout<<"\nposition nach Services";
	se1->parseParameters(params); //vector of parameters
	cout<<se1->masstol_belowDa;
	cout<<"\nposition nach parseParameters";
	cout<<se1->masstol_below;

	//define unique objects
	DnaAA *dnaAA1 = new DnaAA();
	//dnaAA1->checkInitialization();
	Spectra *spectra1 = new Spectra(se1, dnaAA1);
	SpectrumParser *spectrumparser1 = new SpectrumParser(spectra1, se1);
	Tuples *tuples1 = new Tuples(se1, dnaAA1);
	ProteinParser *proteinparser1 = new ProteinParser(se1, tuples1);
	Chromosomes *chromosomes1 = new Chromosomes(se1);
	SlidingWindow *slidingwindow1 = new SlidingWindow(se1);
	Scoring *scoring1 = new Scoring(se1);
	
	//connect unique objects, define flow
	chromosomes1->setSlidingWindow(slidingwindow1);
	slidingwindow1->setTuples(tuples1);
	tuples1->setScoring(scoring1);
	scoring1->setSpectra(spectra1);
	
	if(se1->hotspectra == true) spectrumparser1->hotspectra1->parseHotSpectrumFiles("in_hotspectra.param");

	//load spectra, learn spectrum intensities, discard spectra
	se1->learning_spectrum_intensities = true;
	spectrumparser1->parseFiles(specfiles);
	se1->invertLearnedPeakDistribution();
	
	//load spectra, preprocess them, keep them
	se1->learning_spectrum_intensities = false;
	se1->spectra = 0;
	spectrumparser1->parseFiles(specfiles);

	//add number of spectra to output filename
	paramstring += "_" + se1->intToString(se1->spectra);	
	se1->outfileparams = paramstring;

	//last preparations
	se1->peakdistribution_rel->writeDistribution(se1->getOutFileName("peakdistribution_rel")); //needs parameter string
	se1->peakdistribution_abs->writeDistribution(se1->getOutFileName("peakdistribution_abs")); //needs parameter string
	se1->parentmassdistribution->writeDistribution(se1->getOutFileName("parentmassdistribution")); //needs parameter string
	spectra1->sortAndPrepareSpectra(); //needs parameter string
	slidingwindow1->showParameters(); //show all sliding window parameters
	se1->globaltime->reset();
			
	//PARSE FASTA FILES AND INSTRUCTIONS
	fstream chromprotfile;
	chromprotfile.open("in_fastafiles.param");
	string parsedline;
	
	//start changes BYR
	cout << "\n parsedline:" << parsedline;
	
	int nbresults = 0;

	//bool tag = false;
	bool chromosomes = false;
	bool hotspots = false;
	bool proteins = false;
	//bool write = false;


	while(getline(chromprotfile, parsedline)){
//		//cout << "\nparsedline: " << parsedline;
		if(parsedline.size() > 0 && parsedline[0] != '#'){
			if(parsedline[0] == '<'){
				
				//write tag
				if(parsedline.substr(0, 15) == "<WRITE RESULTS>"){
					tuples1->forwardTuples();
					cout<<"\n this is what the tuples look like in main" << tuples1;
					se1->doProgressReport();
					nbresults++;
					if(se1->outputlevel > 2) cout << "\n\nmain.cpp Writing result files of round " << nbresults << ".\n\n" << flush;				
					se1->outfilenumber = se1->intToString(nbresults);
					spectra1->writeResults();
                    //spectra1->getValues();
					cout << "\n";	
				
				//chromosome tag
				}else if(parsedline.substr(0, 19) == "<CHROMOSOMES START>"){
					chromosomes = true;
					hotspots = false;
					proteins = false;
		
				//hotspots tag
				}else if(parsedline.substr(0, 16) == "<HOTSPOTS START>"){
					chromosomes = false;
					hotspots = true;
					proteins = false;
					
				//proteins tag
				}else if(parsedline.substr(0, 16) == "<PROTEINS START>"){
					chromosomes = false;
					hotspots = false;
					proteins = true;
				}
				
			}else{
				//IS NO TAG
				if(proteins == true){
			//		cout << "\nProcessing protein file " << parsedline << flush;
					proteinparser1->parseFASTA(parsedline,true);
					tuples1->forwardTuples();
					//Standard up to here
				//	cout << "pass done with penalty_max=" << se1->penalty_max;
					double test[] = {3.0,5.0,15.0};
					vector<double> penalty_max(test,test+3);
					for (unsigned int i = 0; i < penalty_max.size();++i)
					{
						se1->penalty_max = penalty_max[i];
						proteinparser1->doAnotherRun();
						tuples1->forwardTuples();
						cout << "pass done with penalty_max=" << se1->penalty_max;

					}
					se1->doProgressReport();

				}
				
				else if(chromosomes == true){
					cout << "\nProcessing chromosome file " << parsedline << flush;
					chromosomes1->parseAndFeedChromosome(parsedline); //c_str converts a cpp-string to a c-string
					tuples1->forwardTuples();
					se1->doProgressReport();
				}else if(hotspots == true){
					cout << "\nProcessing hot spots file " << parsedline << flush;
					chromosomes1->hotspots1->parseHotSpots(parsedline);
				}				
			}
			
		} //line without #
	} //end while
	
	delete se1;
	delete tuples1;
	delete proteinparser1;
	delete chromosomes1;
	delete slidingwindow1;
	delete spectra1;
	delete spectrumparser1;
	delete scoring1;

}
