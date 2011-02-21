#include "Chromosomes.h"
#include "Tuples.h"

namespace Pepsplice{
/*
 * The Chromosomes is unique. It creates chromosomes, parses them and feeds them to SlidingWindow.
 * Chromosome should not know SlidingWindow to avoid circular dependencies.
 *  
 */

Chromosomes::Chromosomes(Services *se0)
{
	se1 = se0;
	hotspots1 = new HotSpots(se0);
}

Chromosomes::~Chromosomes()
{
    delete hotspots1;
}

void Chromosomes::setSlidingWindow(SlidingWindow *sw1)
{
	slidingwindow1 = sw1;
}

void Chromosomes::parseAndFeedChromosome(string file)
{
	
	
	Chromosome *chr = new Chromosome(file, se1);
	chromosomes.push_back(chr);
	chr->loadFullSequence();
	
	//hotspots
	vector<int> hs;
	int hssize = 0;
	int h = 0;
	int tolerancepast = se1->hotspottolerance;
	int toleranceahead = se1->hotspottolerance + se1->gaplenmax;
	
	//FORWARD FEED, COUNT UP
	if(se1->outputlevel > 2) cout << "\nChromosome forward parsing starts. " << chr->description_line;
	resetSlidingWindow();
	slidingwindow1->setChromosome(chr, false); //bool chromosome_reverse
	
	hs = hotspots1->getHotSpots(chr->description_line, false);
	hssize = hs.size();

	for(long i = 0; i < chr->getFullChromSize(); i++){

		//check for hotspots, if there are any
		if(hssize > 0){
			//for reporting
			bool hotbefore = slidingwindow1->hot;
			//check if still true
			if(slidingwindow1->hot == true){
				if(i > hs[h] + toleranceahead){ //past upper bound
					slidingwindow1->hot = false;
				}
			}
			//check if still false
			if(slidingwindow1->hot == false){
				while (i > hs[h] + toleranceahead && h < hssize - 1) h++; //
				if(i >= hs[h] - tolerancepast && i <= hs[h] + toleranceahead){ //entered region
					slidingwindow1->hot = true;
				}
			}
			if(slidingwindow1->hot == false && hotbefore == true) cout << "\n Hot spot ending, prefix end at " << i; 
			if(slidingwindow1->hot == true && hotbefore == false) cout << "\n Hot spot starting, prefix end at " << i;
		}
		
		
		slidingwindow1->addNextNucleotide(chr->getNucleotide(i));
		se1->incrementNucleotides();
	}

	//REVERSE FEED, COUNT DOWN
	if(se1->outputlevel > 2) cout << "\nChromosome reverse parsing starts. " << chr->description_line;
	resetSlidingWindow();
	slidingwindow1->setChromosome(chr, true); //bool chromosome_reverse

	hs = hotspots1->getHotSpots(chr->description_line, true); //bool chromosome_reverse
	hssize = hs.size();

	for(long i = chr->getFullChromSize() - 1; i >= 0; i--){
		
		//check for hotspots, if there are any
		if(hssize > 0){
			//for reporting
			bool hotbefore = slidingwindow1->hot;
			//check if still true
			if(slidingwindow1->hot == true){
				if(i < hs[h] - toleranceahead){ //past lower bound
					slidingwindow1->hot = false;
				}
			}
			//check if still false
			if(slidingwindow1->hot == false){
				while (i < hs[h] - toleranceahead && h > 0) h--; // if i is past bound, adjust h
				if(i >= hs[h] - toleranceahead && i <= hs[h] + tolerancepast){ //entered region
					slidingwindow1->hot = true;
				}
			}
			if(slidingwindow1->hot == true && hotbefore == false) cout << "\n Hot spot starting at " << i;
			if(slidingwindow1->hot == false && hotbefore == true) cout << "\n Hot spot ending at " << i; 
		}
		
		slidingwindow1->addNextNucleotide(chr->getReverseNucleotide(i));
		se1->incrementNucleotides();
	}

}

void Chromosomes::resetSlidingWindow()
{
	Tuples *tuplestemp = slidingwindow1->tuples;
	slidingwindow1 = new SlidingWindow(se1);
	slidingwindow1->setTuples(tuplestemp);
}
}