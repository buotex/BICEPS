
#include "Spectra.h"
#include "SpectrumComparator.h"
namespace Pepsplice{
///Initializer, takes the main instances of the
///services and dna classes

Spectra::Spectra(Services *se0, DnaAA *d)
{
	se1 = se0;
	unreachablebins1 = new UnreachableBins(se1);
	results1 = new Results(se1);

}

Spectra::~Spectra() //Changes by BX
{
    for (int i  = 0; i < spectra.size(); i ++)
    {
        delete spectra[i];
    }
    delete unreachablebins1;
    delete results1;
    
}



void Spectra::addSpectrum(Spectrum *spectrum1){
	if(se1->learning_spectrum_intensities == true){	
		delete spectrum1;
	}else{
		spectra.push_back(spectrum1);
		se1->parentmassdistribution->addElement(0, spectrum1->parentmassMHB, 1);
	}
	if(se1->outputlevel > 5) cout << "\nsize of spectra1: " << spectra.size() << "  parent mass added: " << (*spectrum1).parentmassMH;
	//parentmasses_size = spectra.size();
}



void Spectra::sortAndPrepareSpectra(){
	
	long totskippedpeaks = 0;
	long totlen = 0;
	long totlenB = 0;

	//sort by parent mass
	sort(spectra.begin(), spectra.end(), SpectrumComparator());  //greater<Spectrum*>());
	
	//rearrange spectra
	if(se1->rearrangespectra == true){
		if(se1->writespecpos == true) writeSpectrumPositionsInMemory(se1->getOutFileName("specposfirst"));
		if(se1->outputlevel > 2) cout << "\nRearranging spectra in memory to increase performance for long-running searches. This may take seconds to minutes." << flush;
		//rearrange spectra in memory
		Spectrum *src, *destinationspectrum;
		vector<Spectrum*> spectraold;
		for(int i = 0; i < spectra.size(); i++){
			src = spectra[i];
			destinationspectrum = new Spectrum(src);
			spectra[i] = destinationspectrum;
			spectraold.push_back(src);
		}
		for(int i = 0; i < spectraold.size(); i++){
			delete spectraold[i];
		}
		spectraold.clear();
		if(se1->writespecpos == true) writeSpectrumPositionsInMemory(se1->getOutFileName("specposafter"));
	}
	
	//initialize BestMatches, index them for writing results later on
	for(int i = 0; i < spectra.size(); i++){
		spectra[i]->initializeBestMatches();
		//spectra[i]->bestmatches1->i = i;
	}
	
	//collect spectrum statistics
	for(int i = 0; i < spectra.size(); i++){
		parentmasses.push_back(spectra[i]->parentmassMH);
		totlen += spectra[i]->lengthgood;
        if (se1->outputlevel >2 ) cout << spectra[i] -> lengthgood << endl;
		totlenB += spectra[i]->lengthB;
		totskippedpeaks += spectra[i]->skippedpeaks;
		if(se1->outputlevel > 5) cout << "\nSpectra::sortSpectraAndExtractPM: pointer " << (long)spectra[i] << "  parent mass " << spectra[i]->parentmassMH;
	}

	//initialize spectrum parameters for hypergeometric score
	unreachablebins1->initialize();
	for(int i = 0; i < spectra.size(); i++){
		Spectrum *spec = spectra[i];
		spec->initializeHypergeometricNK(); //naive estimation for N and K
		
		//reestimation for N
		spec->N = spec->parentmassMHB - unreachablebins1->getTotal(se1->trypticendsrequired) - se1->blockedbins;
		
		//reestimation for K
		spec->K = unreachablebins1->estimateK(spec->binaryB, spec->parentmassMHB, se1->trypticendsrequired);
		
		//avoid exceptions in hypergeometric calculations
		if(spec->N < spec->K) spec->N = spec->K;
	}
	
	//set parent mass limits for tuple generation
	if(parentmasses.size() > 0) se1->setMinMaxPM(parentmasses[0], parentmasses[spectra.size() - 1]);
	
	//progressreport statistics
	double sf = se1->dnaAA1->scaling_factor;
	if(se1->outputlevel > 2){
		if(spectra.size() > 0){
			cout << "\nSpectra::sortSpectraAndExtractPM:";
			cout << "\nThere are " << spectra.size() << " spectra which were sorted according to their parent mass.";
			cout << "\nA spectrum contains on average " << (double)totlen / (double)spectra.size() << " correct peaks.";
			cout << "\nAfter preprocessing, the spectrum contains on average " << (double)totlenB / (double)spectra.size() << " peaks.";
			cout << "\nA spectrum contains on average " << (double)totskippedpeaks / (double) spectra.size() << " peaks of m/z < 0 or m/z > parentmass. Those were skipped.";
			cout << "\nLowest parent mass in spectra : " << parentmasses[0]/sf << "   Lowest theoretical parent mass considered: " << se1->min_monoparentmassMH / sf;
			cout << "\nHighest parent mass in spectra: " << parentmasses[spectra.size() - 1] / sf << "   Highest theoretical parent mass considered: " << se1->max_monoparentmassMH / sf;
		}else{
			cout << "\nSpectra::sortSpectraAndExtractPM: No spectra were parsed!";
		}
	}
}

void Spectra::writeSpectrumPositionsInMemory(string filename){
	
	ofstream outfile;
	outfile.open(filename.c_str());
	outfile << "i \tsp_pointer \tpkl_pointer \tdiff_spec_pkl \tdiff_spec_spec \tparent_mass";
	int diff_spec_spec = 0;
	for(int i = 0; i < spectra.size(); i++){
		if(i == 0){
			diff_spec_spec = 0;
		}else{
			diff_spec_spec = (long)spectra[i] - (long)spectra[i-1];
		}
		outfile << "\n" << i << "\t" << (long)spectra[i] << "\t" << (long)spectra[i]->binsB << "\t" << (long)spectra[i]->binsB - (long)spectra[i] << "\t" << diff_spec_spec << "\t" << spectra[i]->parentmassMH/se1->dnaAA1->scaling_factor;
	}
	outfile.close();
}

//begin changes BXu
void Spectra::initializeResults(){
	//analyze best matches, retrieve values
	for(int i = 0; i < spectra.size(); i++){
		spectra[i]->bestmatches1->sortBestMatches();
		spectra[i]->bestmatches1->analyzeBestHits();
	}
}
    
    bool Spectra::returnBIC(vector<PepspliceResult> & results){
        
        initializeResults();
        for (int j = 0; j < spectra.size(); ++j)
        {
            BestMatches *bm1 = spectra[j]->bestmatches1; //the best match of spectrum 1, without changes to pepmass    
            
            for(int i = 0; i < bm1->bm.size(); ++i)
            {
                PepspliceResult pep;
                
                pep.penalty = bm1->bm[i]->tuple->penalty;
                
                pep.score=bm1->bm[i]->score;
                
                float likelihood1 = log((pep.score > 0)?pep.score:FLT_MIN);
                pep.n = (bm1->bm[i]->tuple->length-1) * 2;
                pep.k = bm1->bm[i]->sharedpeakcount;		
                
                pep.bic = likelihood1-(pep.penalty+se1->penalty_mutation)/2*log((float)pep.n/2);
                
                pep.Sequence = bm1->bm[i]->tuple->getAASeq();
                //pep.OrigSequence = bm1->bm[i].tuple->getAASeqOriginal();
                pep.OrigSequence = (*oldSequences)[bm1->bm[i]->tuple->aasequenceorig];
                pep.penalty_max = se1->penalty_max;
                results.push_back(pep);
                delete bm1->bm[i];
                //bm1->bm[i].tuple->decreasePointerCount();
            }
            bm1->bm.clear();
        }
        return 1;
    }
    




void Spectra::writeResults() 
{
	
	vector<BestMatches*> b = vector<BestMatches*>(spectra.size());
	//analyze best matches, retrieve values
//	for(int i = 0; i < spectra.size(); i++){
//		spectra[i]->bestmatches1->sortBestMatches();
//		spectra[i]->bestmatches1->analyzeBestHits();
//		b[i] = spectra[i]->bestmatches1;
//	}

	initializeResults();
	for(int i = 0; i < spectra.size(); i++){
		b[i] = spectra[i]->bestmatches1;
	}
	results1->writeFiles(b);
	
}

}//namespace pepsplice

