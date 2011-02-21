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

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <iostream>	//cout
#include <string>
#include <vector>
#include <cmath> //pow
#include "BestMatches.h"
#include "Distribution.h"
#include "Peak.h"
#include "PeakComparator.h"
#include "Services.h"

//there are circular dependencies and compilation problems with BestMatches-Match-Spectrum

#include "DnaAA.h"

using namespace std;

namespace Pepsplice{
class BestMatches;


        
    class Spectrum
    {
    private:
        //long i;

        
    public:
        //construction
        int chargestate;
        double parentmassMH;
        int parentmassMHB;
        float lowerBoundScore;

        int skippedpeaks;
        long lengthgood;
            
        //peak lists
        vector<Peak> *temppeaks1;
        vector<Peak> *temppeaks2;
        int *binsB;
        //float *intensitiesB;
        bool *binaryB;
        int lengthB;
        int scored;
        int N, K; //variable for hypergeometric score
        bool *arrayN;
        int empiricalN, empiricalK;
        //float ticB;

        string filename, mgfname;
        
        //OBJECTS
        BestMatches *bestmatches1;
        //DnaAA *dnaAA1;
        Services *se1;
        
        //METHODS
        
        //constructing empty spectrum
        Spectrum(double pmMH, int cs, DnaAA *d1, Services *se0);
        Spectrum(Spectrum *s);
        //Spectrum(DnaAA *d1, BestMatches *bestm, double pmMH, int cs, int lenB, int *bnsB, int lengood, int skp, string *filen, string *mgfn);
        virtual ~Spectrum();
        void setFileName(const string & fn);
        void setMGFName(const string & mn);
        string getFileName();
        string getMGFName();
        
        //add peaks
        void addDataPoint(float mz, float intensity);
        void initializeBestMatches();
        
        //preprocessing
        void preprocessSpectrum();
        void normalizeByInvertedDistribution(vector<Peak> *p);
        void sortPeaksByIntensity(vector<Peak> *p);
        void bestGlobalIntensityPeaks(vector<Peak> *p, int bestX);
        void bestLocalIntensityPeaks(vector<Peak> *p, int tolerance, int bestX);
        void bestLocalIsotopePeaks(vector<Peak> *p, int tolerance, int bestX);
        vector<Peak> *shortenPeakList(vector<Peak> *p);
        void peakListToArray(vector<Peak> *p);		
        void normalizeToLocalAvgIntensityPerDa(vector<Peak> *p, int tolerance);
            
        //deisotoping
        void isotopeSimilarity(vector<Peak> *p);
        float computePoisson(int x, float mass);
        float compareDistributions(float idealDistribution[], float realDistribution[], int len);

        //initialize parameters for hypergeometric score, is called by Spectra::sortAndPrepareSpectra()
        void initializeHypergeometricNK();
        int getEmpiricalN();
        int getEmpiricalK();

        
    };
}

#endif /*SPECTRUM_H_*/
