#include "ResultsPattern.h"
namespace Pepsplice{
ResultsPattern::ResultsPattern(Services *se0)
{
	se1 = se0;

	//theoretical ions
	bions = new float[1000];
	yions = new float[1000];
	bionsint = new int[1000];
	yionsint = new int[1000];
	
	//for series pattern statistics
	expspec = new int[100000];
	for(int i = 0; i < 100000; i++){expspec[i] = 0;} //initialize
}

ResultsPattern::~ResultsPattern()
{
    delete[] bions;
    delete[] yions;
    delete[] bionsint;
    delete[] yionsint;
    delete[] expspec;
}

void ResultsPattern::calcMassDeviationAndSeriesPattern(vector<BestMatches*> b, int dl){
	
	Distribution massdeviations = Distribution(-10, 10, 1, 200);
	Distribution *seriespatterns = new Distribution(-50.5, 50.5, 100, 101);

	//collect statistics for reliable matches	
	for(int i = 0; i < b.size(); i++){	
		if(b[i]->precision > 0.90){
			massdeviations.addElement(0, (b[i]->ownerSpectrum->parentmassMH - b[i]->thpm) / se1->dnaAA1->scaling_factor);
			if(b[i]->getBestMatch() != NULL){
				getPatternAroundMatch(b[i]->getBestMatch(), seriespatterns);
			}
		}
	}

	massdeviations.writeDistribution(se1->getOutFileName("_dl" + se1->intToString(dl) + "_massdeviations"));
		
	//write and delete seriespatterns
	seriespatterns->writeDistribution(se1->getOutFileName("_dl" + se1->intToString(dl) + "_seriespattern"));
	delete seriespatterns;
		
}

void ResultsPattern::calculateBYIonsInt(unsigned char *aasequence, int len)
{
	//ions are shifted by 1 with respect to aasequence
	bionsint[0] = se1->dnaAA1->integermassH;
	yionsint[0] = se1->dnaAA1->integermassHOHH;
	for(int i = 1; i <= len; i++){
		bionsint[i] = bionsint[i - 1] + se1->dnaAA1->aa_integermasses256[aasequence[i - 1]];
		yionsint[i] = yionsint[i - 1] + se1->dnaAA1->aa_integermasses256[aasequence[len - i]];
	}
}

void ResultsPattern::calculateBYIons(unsigned char *aasequence, int len)
{
	//ions are shifted by 1 with respect to aasequence
	//float pm1 = 0;
	bions[0] = se1->dnaAA1->monomassH;
	yions[0] = se1->dnaAA1->monomassHOHH;
	for(int i = 1; i <= len; i++){
		bions[i] = bions[i - 1] + se1->dnaAA1->aa_monomasses256float[aasequence[i - 1]];
		yions[i] = yions[i - 1] + se1->dnaAA1->aa_monomasses256float[aasequence[len - i]];
	}
}


void ResultsPattern::getPatternAroundMatch(Match *mt, Distribution *d){
	Spectrum *spec = mt->spectrum;
	Tuple *tup = mt->tuple;
	
	calculateBYIonsInt(tup->aasequence, tup->length);
	
	//expand experimental spectrum
	for(int i = 0; i < spec->lengthB; i++){expspec[spec->binsB[i]] = 1;}

	//getTotal peaks around b- and y-ions
	for(int i = 1; i < tup->length; i++){
		for(int j = -50; j < 50; j++){
			if(expspec[bionsint[i] + j] > 0) d->addElement(i, j);
			if(expspec[yionsint[i] + j] > 0) d->addElement(i + 40, j);
		}
	}
	
	//erase experimental spectrum
	for(int i = 0; i < spec->lengthB; i++){expspec[spec->binsB[i]] = 0;}	
}

void ResultsPattern::writeBestTheoreticalSpectra(vector<BestMatches*> b){
	ofstream outfile;
	string file = se1->getOutFileName("theoreticalspectra");
	outfile.open(file.c_str());
	outfile.precision(7);
	for(int i = 0; i < b.size(); i++){
		if(b[i]->getBestMatch() != NULL && b[i]->precision > 0.95){
			Match *mt = b[i]->getBestMatch();
			Spectrum *spec = mt->spectrum;
			Tuple *tup = mt->tuple;
			vector<float> v;
			
			calculateBYIons(tup->aasequence, tup->length);
			
			//getTotal peaks around b- && y-ions
			for(int i = 1; i < tup->length; i++){
				v.push_back(bions[i]/se1->dnaAA1->scaling_factor);
				v.push_back(yions[i]/se1->dnaAA1->scaling_factor);
			}
			
			sort(v.begin(), v.end());
			
			outfile << "\nBEGIN IONS";
			outfile << "\nTITLE=" << spec->filename;
			outfile << "\nCHARGE=" << spec->chargestate;
			double parentmassMHX = spec->parentmassMH + se1->dnaAA1->monomassH*(spec->chargestate-1);
			double parentmassMmeas = parentmassMHX / spec->chargestate / se1->dnaAA1->scaling_factor;
			//parentmassMHX = (parentmassMmeas * chargestate * dnaAA1->scaling_factor); // - (dnaAA1->monomassH * chargestate) + dnaAA1->monomassH;
			//parentmassMH = parentmassMHX - dnaAA1->monomassH*(chargestate-1);
			outfile << "\nPEPMASS=" << parentmassMmeas;
			
			for(int i = 0; i < v.size(); i++){
				outfile << "\n" << v[i] << "\t" << 1;
			}
			outfile << "\nEND IONS\n";
		}
	}
	outfile.close();
}


}