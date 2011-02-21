#include "Spectrum.h"

/*
 * Contains the information of a single spectrum.
 * 
 * There is one BestMatches per spectrum, as a child of spectrum.
 * 
 * Beware: the copy constructor must be complete! In case you forget to
 * include variables, you risk segmentation faults! 
 * 
 */
namespace Pepsplice{
Spectrum::Spectrum(double pmMH, int cs, DnaAA *d1, Services *se0)
{
	//constructor arguments

	parentmassMH = pmMH;
	chargestate = cs;

	//dnaAA1 = d1;
	parentmassMHB = (int)(parentmassMH/se0->dnaAA1->scaling_factor + 0.5) +1;
	se1 = se0;
	se1->incrementSpectra();
	bestmatches1 = NULL;

	//initializations
	filename = "";
	mgfname = "";
	skippedpeaks = 0; //invalid temppeaks1
		
	//peak list
	lengthgood = 0;
	scored = 0;
	N = 0;
	K = 0;
	//temppeaks1 = new vector<Peak>;
	
	//recycle same temporary peak list vector over and over again to prevent heap fragmentation
	temppeaks1 = se1->temppeaks1;
	temppeaks1->clear();
	temppeaks2 = se1->temppeaks2;
	temppeaks2->clear();
    

}


//similar to copy constructor
Spectrum::Spectrum(Spectrum *s)
{
	//helper objects
	//dnaAA1 = s->dnaAA1;
	se1 = s->se1;
	bestmatches1 = s->bestmatches1;

	//parent mass
	parentmassMH = s->parentmassMH;
	parentmassMHB = (int)(parentmassMH/se1->dnaAA1->scaling_factor + 0.5)+1;
	chargestate = s->chargestate;
		
	//peaks
	lengthB = s->lengthB;
	binsB = new int[lengthB];
	for(int i = 0; i < lengthB; i++){
		binsB[i] = s->binsB[i];
	}
	binaryB = new bool[parentmassMHB];
	for(int i = 0; i < parentmassMHB; i++){
		binaryB[i] = s->binaryB[i];
	}
	
	if(se1->determineNK){
		arrayN = new bool[parentmassMHB];
		for(int i = 0; i < parentmassMHB; i++){
			arrayN[i] = s->arrayN[i];
		}
	}
	
	lengthgood = s->lengthgood;
	skippedpeaks = s->skippedpeaks;
	scored = s->scored;
	N = s->N;
	K = s->K;
	
	//file names
	filename = s->filename;
	mgfname = s->mgfname;
}


Spectrum::~Spectrum()
{
	delete[] binsB;
	//delete[] *intensitiesB;
	delete[] binaryB;
    delete bestmatches1;
}

void Spectrum::addDataPoint(float mz, float intensity)
{
	//peak distribution statistics for normalization later on
	if(se1->learning_spectrum_intensities == true){
		se1->peakdistribution_rel->addElement(0, mz/parentmassMH, intensity);
		se1->peakdistribution_abs->addElement(0, (int)(mz/se1->dnaAA1->discretization_factor + 0.5), 1); //frequency only
	}
	
	//reasons not to add a data point
	if(se1->learning_spectrum_intensities || mz < 0 || mz >= parentmassMH - 5 || intensity <= 0){
		skippedpeaks++;
		return;
	}
	
	//add data point, discretize immediately
	int bin = (int)(mz/se1->dnaAA1->discretization_factor + 0.5);
	if(temppeaks1->size() > 0 && temppeaks1->back().bin == bin){
		temppeaks1->back().intensity += intensity; //merge temppeaks1 that fall into same bin immediately
	}else{
		temppeaks1->push_back(Peak(bin, mz, intensity)); //add peak in regular way
	}
}

void Spectrum::setFileName(const string & fn)
{
	filename = fn;
}

void Spectrum::setMGFName(const string & mn)
{
	mgfname = mn;
}

string Spectrum::getFileName()
{
	return filename;
}

string Spectrum::getMGFName()
{
return mgfname;
}

void Spectrum::initializeBestMatches(){
	bestmatches1 = new BestMatches(this);
}

void Spectrum::preprocessSpectrum()
{
	lengthgood = temppeaks1->size();

	isotopeSimilarity(temppeaks1);
	bestLocalIsotopePeaks(temppeaks1, se1->specwinisosize, se1->specwinisopeaks);
	
	normalizeByInvertedDistribution(temppeaks1);
	bestLocalIntensityPeaks(temppeaks1, se1->specwinintsize, se1->specwinintpeaks);
	bestGlobalIntensityPeaks(temppeaks1, (int)(parentmassMHB/100*se1->specglobintpeaksper100));
	
	temppeaks1 = shortenPeakList(temppeaks1);
	sortPeaksByIntensity(temppeaks1);
	//normalizeToLocalAvgIntensityPerDa(temppeaks1, 100);
	peakListToArray(temppeaks1);
}

void Spectrum::sortPeaksByIntensity(vector<Peak> *p){
	//a Peak is NOT on the heap, only the vector of Peaks
	sort(p->begin(), p->end(), PeakComparator());
	//cout << "\n" << (*p)[0].intensity << " " << (*p)[1].intensity << " " << (*p)[2].intensity;
}

void Spectrum::normalizeByInvertedDistribution(vector<Peak> *p){
	float relmz = 0;
	int binpos = 0; //distribution=histogram bin
	float factor = 0;
	for(int i = 0; i < p->size(); i++){
		relmz = (*p)[i].mz / parentmassMH;
		binpos = se1->peakdistribution_rel->getBinPos(3, relmz);
		(*p)[i].intensity = (*p)[i].intensity * se1->peakdistribution_rel->getBinValue(3, binpos);
	}
}

void Spectrum::bestGlobalIntensityPeaks(vector<Peak> *p, int bestX){
	
	if(bestX == 0) return;
	
	//determine cutoff intensity value
	vector<float> tempintensities;
	for(int i = 0; i < p->size(); i++){
		tempintensities.push_back((*p)[i].intensity);
	}
	sort(tempintensities.begin(), tempintensities.end());
	int cutoffpos = tempintensities.size() - 1 - bestX;
	double cutoffvalue = 0;
	if(cutoffpos >= 0) cutoffvalue = tempintensities[cutoffpos];
	
	//label peaks: keep or delete
	for(int i = 0; i < p->size(); i++){
		if((*p)[i].intensity >= cutoffvalue){
			(*p)[i].above_glob_intensity_cutoff = true;
		}else{
			(*p)[i].above_glob_intensity_cutoff = false;
		}
	}
}

void Spectrum::bestLocalIntensityPeaks(vector<Peak> *p, int tolerance, int bestX){
	
	//a peak is accepted if it is among the top X temppeaks1 within the tolerance region

	if(tolerance == 0 || bestX == 0) return;
		
	for(int i = 0; i < p->size(); i++){
		
		int j = 0;
		int othersbetter = 0;
		
		//left half window
		j = i;
		while(j >= 0 && (*p)[j].bin > (*p)[i].bin - tolerance){
			if((*p)[i].intensity < (*p)[j].intensity) othersbetter++;
			if(othersbetter >= bestX) break;
			j--;
		}
		
		//right half window
		j = i;
		while(j < p->size() && (*p)[j].bin < (*p)[i].bin + tolerance){
			if((*p)[i].intensity < (*p)[j].intensity) othersbetter++;
			if(othersbetter >= bestX) break;
			j++;
		}

		//if the peak has survived till here, add to club of top temppeaks1
		if(othersbetter < bestX){
			(*p)[i].above_intensity_cutoff = true;
		}else{
			(*p)[i].above_intensity_cutoff = false;
		}
		
	}
}

void Spectrum::bestLocalIsotopePeaks(vector<Peak> *p, int tolerance, int bestX){
	
	//a peak is only accepted if it is among the top X temppeaks1 within the tolerance region
	
	if(tolerance == 0 || bestX == 0) return;

	for(int i = 0; i < p->size(); i++){
		
		int j = 0;
		int othersbetter = 0;
		
		//left half window
		j = i;
		while(j >= 0 && (*p)[j].bin > (*p)[i].bin - tolerance){
			if((*p)[i].isotopesim < (*p)[j].isotopesim)	othersbetter++;
			if(othersbetter >= bestX) break;
			j--;
		}
		
		//right half window
		j = i;
		while(j < p->size() && (*p)[j].bin < (*p)[i].bin + tolerance){
			if((*p)[i].isotopesim < (*p)[j].isotopesim) othersbetter++;
			if(othersbetter >= bestX) break;
			j++;
		}
		
		//if the peak has survived till here, add to club of top temppeaks1
		if(othersbetter < bestX){
			(*p)[i].above_isotope_cutoff = true;
		}else{
			(*p)[i].above_isotope_cutoff = false;
		}
	}
}


vector<Peak> *Spectrum::shortenPeakList(vector<Peak> *p){

	//all conditions must be true, but default is true anyway
	
	temppeaks2->clear();
	for(int i = 0; i < p->size(); i++){
		//cout << (*p)[i].bin << " " << (*p)[i].intensity;
		if((*p)[i].above_intensity_cutoff == true && (*p)[i].above_isotope_cutoff == true && (*p)[i].above_inv_intensity_cutoff == true && (*p)[i].above_glob_intensity_cutoff == true){
			temppeaks2->push_back((*p)[i]);
		}
	}
	
	p->clear();
	lengthB = temppeaks2->size();
	return temppeaks2;
}

void Spectrum::peakListToArray(vector<Peak> *p){
	
	int psize = p->size();
	//ticB = 0;
	
	//initialize dynamic arrays
	binsB = new int[psize + 1];
	//intensitiesB = new float[psize + 1];
	binaryB = new bool[parentmassMHB];
	for(int i = 0; i < parentmassMHB; i++){binaryB[i] = 0;}
	
	if(se1->determineNK){
		//cout << "-";
		arrayN = new bool[parentmassMHB];
		for(int i = 0; i < parentmassMHB; i++){arrayN[i] = 0;}
	}

	//fill in arrays
	for(int i = 0; i < psize; i++){
		binsB[i] = (*p)[i].bin;
		//intensitiesB[i] = (*p)[i].intensityNormLocAvg;
		if((*p)[i].bin < parentmassMHB) binaryB[(*p)[i].bin] = 1;
		//ticB += intensitiesB[i];
	}

}


 /**
   *  Poisson model for isotopes
   *  This method calculates the theoretical relative height of an isotope peak within an isotope pattern
   *  normalization: the abundance of all the isotope temppeaks1 in the pattern adds up to 1
   *  formula: p(x, lambda) = e^(-lambda)*lambda^x/x!   (exclamation mark means factorial here)
   * 
   **/
   
void Spectrum::isotopeSimilarity(vector<Peak> *p){
	
    float similarityScore = 0;
    float idealDistribution[4]; // Theoretical Poisson distribution for temppeaks1 0..3
    float realDistribution[4]; // Real Poisson distribution for temppeaks1 0..3

    for (int i = 0; i < p->size(); i++) {

    	//fill in real distribution
    	int basepeak = (*p)[i].bin;
    	
    	//initialize distributions
    	for(int j = 0; j < 3; j++){
    		idealDistribution[j] = computePoisson(j, basepeak);
    		realDistribution[j] = 0;
    	}

    	for(int j = 0; j < 3; j++){
    		int k = i + j;
    		if(k < p->size()){
	        	int offset = (*p)[k].bin - basepeak;
	        	if(offset <= 3) realDistribution[offset] = (*p)[k].intensity;
    		}
        }
        
        // compare real distribution with ideal distribution
        (*p)[i].isotopesim = 1./ (compareDistributions(idealDistribution, realDistribution, 4) + 0.1); //changes by BX

    }
}


float Spectrum::computePoisson(int x, float mass) {
      int factorial[] = {1, 1, 2, 6}; // factorials precalculated for fact(0)..fact(3)
      double excess = 0.6760/1000; // this amount of heavy isotopes occurs per dalton; per 1000 Dalton, a peptide will on average contain 1.2582 H2 or C13 or O18 or S34, ...
      double lambda = mass * excess; // lambda as in the Poisson formula
      double density = 0; // result, i.e. relative frequency 0..1 for an isotope peak x (e.g. 0.7 for x=0 which is the main peak)
      float density2 = 0; // to return a float

      density = exp(-lambda) * pow(lambda, x) / (float)factorial[x];
      density2 = (float)density;
      return density2;
}


float Spectrum::compareDistributions(float idealDistribution[], float realDistribution[], int len) {

      double score1 = 0;
      float realSum = 0;
      double t_chisquare = 0; // root(square of differences between distributions)

      // getTotal the content of the bins to normalize afterwards
      for (int xid = 0; xid < len - 1; xid++) {
        realSum += realDistribution[xid]; // sum up real distribution
      }

      // normalize real distribution to 1
      for (int xid = 0; xid < len - 1; xid++) {
        realDistribution[xid] = realDistribution[xid] / realSum;
      }
      realSum = 0;

      // compare ideal && real distributions, calculate squared differences
      for (int xid = 0; xid < len - 1; xid++) {
        t_chisquare += (pow(idealDistribution[xid] - realDistribution[xid], 2) / (realDistribution[xid] + 0.001)); // square differences
      }
      score1 = t_chisquare; // t_chisquare is the smaller the better
      return (float)score1;
}

void Spectrum::normalizeToLocalAvgIntensityPerDa(vector<Peak> *p, int tolerance){
	
	if(tolerance == 0) return;
		
	int lowerbound_i = 0;
	int upperbound_i = -1;
	int lowerbound_bin = 0;
	int upperbound_bin = 0;
	int psize = p->size();
	
	float peakswithin = 0;
	float intensityOfPeaksWithin = 0;
	float avg = 0;
	
	for(int i = 0; i < psize; i++){
	
		//find upperbound; position of upperbound is included; test if next peak should be included
		while(upperbound_i < psize-1 && (*p)[upperbound_i + 1].bin < (*p)[i].bin + tolerance){
			upperbound_i++;
			intensityOfPeaksWithin += (*p)[upperbound_i].intensity;
			peakswithin++;			
		}
		
		//find lowerbound; test if current lower boundary peak should be evicted
		while(lowerbound_i < psize-1 && (*p)[lowerbound_i].bin < (*p)[i].bin - tolerance){
			lowerbound_i++;
			intensityOfPeaksWithin -= (*p)[lowerbound_i].intensity;
			peakswithin--;
		}
		
		//cout << "\nlowerbound_i: " << lowerbound_i << " upperbound_i: " << upperbound_i << " peakswithin: " << peakswithin << " intensityOfPeaksWithin: " << intensityOfPeaksWithin;
		
		//final result for each peak
		float avg = intensityOfPeaksWithin/(tolerance*2);
		(*p)[i].intensityNormLocAvg = (*p)[i].intensity / avg;
		
		//if((*p)[i].intensityNormLocAvg > 10) cout << "\n" << (*p)[i].bin << " " << (*p)[i].intensityNormLocAvg << "  " << (*p)[i].intensity;
	}
}

void Spectrum::initializeHypergeometricNK(){

	N = parentmassMHB - se1->blockedbins;
	K = lengthB;
	if(N < K) N = K;
	if(K < 1) K = 1;

}

int Spectrum::getEmpiricalN(){
	if(se1->determineNK == false) return 0;
	int N = 0;
	for(int i = 0; i < parentmassMHB; i++){
		if(arrayN[i] == 1) N++;
	}
	return N;
}

int Spectrum::getEmpiricalK(){
	if(se1->determineNK == false) return 0;
	int K = 0;
	for(int i = 0; i < parentmassMHB; i++){
		//cout << arrayN[i] << ":" << (int)binaryB[i] << " "; 
		if(arrayN[i] == 1 && binaryB[i] == 1) K++;
	}
	return K;
}

/*

void Spectrum::discretize(int inverse_granularity, float *sourcearray, int *destinationarray, long len){
	int pos = 0;
	for(int i = 0; i < len; i++){
		pos = (int)(sourcearray[i] / dnaAA1->discretization_factor + 0.5) * inverse_granularity;
		destinationarray[i] = pos;
	}
}

void Spectrum::normalizeRegionWise(int nbwin, float *mzs, float *source_intens, float *dest_intens, float *TIC, long len)
{
	//cout << "normalizeRegionWise entered";
	float windowsize = parentmassMH / nbwin;
	float max = 0;
	int j = 0;
	int startj = 0;
	
	j = 0;
	for(int i = 1; i <= nbwin; i++)
	{
		max = 1;
		startj = j;
		while(j < len && mzs[j] < i * windowsize){
			if(source_intens[j] > max) max = source_intens[j];
			j++;
		}
		//cout << "\nstartj: " << startj << " j: " << j;
		for(int k = startj; k < j; k++){
			dest_intens[k] = source_intens[k] * 1 / max;
		}
        
	}
	
	//in case there are temppeaks1 beyond the last window (should not happen)
	if(j < len){
		dest_intens[j] = 0;
		j++;
	}
	
	//calculate TIC for scoring3
	(*TIC) = 0;
	for(int i = 0; i < len; i++){
		(*TIC) += dest_intens[i];
	}
		
}

*/
}
