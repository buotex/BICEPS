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

#ifndef SYSTEMENVIRONMENT_H_
#define SYSTEMENVIRONMENT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "Timer.h"
#include "DnaAA.h"
#include "Distribution.h"
#include "Peak.h"
#include "Hypergeometric.h"

using namespace std;
namespace Pepsplice{
class Services
{

public:
    bool changedMass;//you get 2 mgf files, one for pepmass-1/charge
	//parameters
	//bool parametersparsed;
	int	eachXspec;
	int tupspec;
	double masstol_above;
	double masstol_below;
	double masstol;
	int safety_margin_disc;
	double min_monoparentmassMH, max_monoparentmassMH;
	int trypticendsrequired;
	double realtuples_per_randomtuple;
	int randomtype;
	bool spliced;
	bool rearrangespectra;
	int processorcacheL2;
	int discriminationscore;
	int scorepval;
	int scoring;
	int extraload;
	bool adapttuplebuffer;
	double masstol_belowDa;
	double masstol_aboveDa;
	bool dismisstuples;
	bool force_tb_report;
	int scorebins;
	bool firstscore;
	double firstpeaksperhundred;
	double firstpdfcutoff;
	bool writespecpos;
	int bestseries;
	bool noprefixend;
	bool nosuffixstart;
	//int modifications;
	int blockedbins;
	bool determineNK;
	bool writethspec;
	int bestmatches;
	bool writebestmatches;
	bool write20bestscores;
	int resizeeveryX;
	bool chargestate2;
	double minpepconfidence;
	int hotspottolerance;
	bool hotspectra;
	bool doWholeGenome;
	bool writespecsubsetasmgf;
	bool writespecsubsetasdta;
	bool doModifications;
	bool saveOld; //to save tuples instead of deleting them
	//int modificationvalue;
	bool doSNPnt;
	bool doSNPaa;
	
	float penalty_mutation; //added by BX
	float penalty_max;
    vector<float> max_penalties;
    float omax_penalty; //added by BX
	float penalty_ptm;
	float penalty_snp;
	float penalty_methox;
	float penalty_tryptic;
	float penalty_genomic;
	float penalty_spliced;
	float penalty_misscleav;

    
	//modifications
	bool mod_m16;
	bool mod_s80;
	
	//output parameters
	int outputlevel;
	string outfileprefix;
	string outfilesuffix;
	string outfilenumber;
	string outfileparams;
		
	//preprocessing parameters
	bool learning_spectrum_intensities;
	int specglobintpeaksper100;
	int specwinintsize;
	int specwinintpeaks;
	int specwinisosize;
	int specwinisopeaks;
	
	//splicing parameters
	int gaplenmin;		//minimum gap length 
	int gaplenmax;		//maximum gap length
	int totlenmax;		//maximum (prefix + suffix)
	int plenmin;		//minimum prefix length (DNA)
	int slenmin;		//minimum suffix length	
	
	//Spectrum.cpp
	vector<Peak> *temppeaks1; //give the same to each spectrum for preprocessing, prevents heap fragmentation, i.e. bad memory allocation
	vector<Peak> *temppeaks2;
    //SpectrumParser.cpp
	double nexttime;
	//bool results_desired;
	

	//counters
	double nucleotides, triggernucleotide;
	double aminoacids;
	double tuples, randomtuples, realtuples, triggertuple;
	unsigned long tuplebuffersize;
	double tuples_per_second;
	double matches, triggermatch;
	double scores;
	int spectra;
	double prescored;
	double postscored;
	ofstream os;
    
	Timer *globaltime;
	DnaAA *dnaAA1;
	Hypergeometric *hypergeometric1;
	//Distribution *distr_specwithin;
	Distribution *parentmassdistribution;
	Distribution *peakdistribution_abs;
	Distribution *peakdistribution_rel;
	Distribution *progressreport;
    int reportnumber;
	
	Services();
	void showSizeOfTypes();
	virtual ~Services();
	void setMinMaxPM(double min_meas_PM, double max_meas_PM);
//	void updateParentMassTol();
	void suggestProgressReport();

	void writeLogFiles();
	void incrementMatches();
	void incrementTuples(bool FP);
	void incrementSpectra();
	void invertLearnedPeakDistribution();
	//double oneOutOfXTuplesIsFP();
	void incrementNucleotides();
	string getOutFileName(string x);	
	string intToString(int i);
	double string_to_double(string s);
	void doProgressReport();
	void parseParameters(vector<string> & args);
	void parseParameter(string arg);
	
private:

};
}
#endif /*SYSTEMENVIRONMENT_H_*/
