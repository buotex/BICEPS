#include "BestMatches.h"
#include "MatchComparator.h" //integration was cumbersome; omit first in BestMatches.h, compile, then add it
namespace Pepsplice{
BestMatches::BestMatches(Spectrum *ownerspec)
{
	ownerSpectrum = ownerspec;
	ownerSpectrum->lowerBoundScore = 0;
	bm.reserve((ownerSpectrum->se1->bestmatches + 4) * 4);
	pvalue = 0;
	matches = 0;
	bestscores = new float[ownerSpectrum->se1->bestmatches + 1];
	//length = len;
	//lowerBoundScore = 0;
	//length2 = 0;
}

BestMatches::~BestMatches(){
	delete[] bestscores;
}

void BestMatches::reset()
{
    for (int i = 0; i < bm.size(); i++)
        delete bm[i];
	bm.clear();
    
	//bm.reserve((ownerSpectrum->se1->bestmatches + 4) * 4);
	pvalue = 0;
	matches = 0;
	for (unsigned int i = 0; i < ownerSpectrum->se1->bestmatches+1;++i)
	{
		bestscores[i] = 0.0f;
	
	}
	
	
	
}


void BestMatches::addMatch(Match * mt1){
	bm.push_back(mt1);
	matches++;
	//if(matches % ownerSpectrum->se1->resizeeveryX == 0) resize();
	if(matches % (ownerSpectrum->se1->bestmatches + 4) == 0) resize();
	//if(bm.size() > ownerSpectrum->se1->bestmatches * 2) resize();
}

void BestMatches::resize(){
	
	sortBestMatches();
	
	float bestscore = bm[0]->score;
	float currentscore = bestscore;
	string bestsequence = bm[0]->tuple->getAASeq();
	string currentsequence = bestsequence;
	int i = 0;
	int delta_i = 0;
	
	while(i < bm.size()){
		if(bm[i]->tuple->getAASeq() != currentsequence) delta_i++; //delta_i = 1 for second best match
		if(delta_i > ownerSpectrum->se1->bestmatches) break;
		//if(bm[i]->score < bestscore) break;
		currentsequence = bm[i]->tuple->getAASeq();
		i++;
	}
    
	while(bm.size() > i + 1 && bm.size() >=2){
		delete bm.back(); //delete match; match will decrease pointercount of tuple
		//bm.back().tuple->decreasePointerCount();
		bm.pop_back(); //shorten vector by one
	}
	
	ownerSpectrum->lowerBoundScore = bm.back()->score;
	//tot_best_hits = bm.size();
}

Match* BestMatches::getBestMatch(){
	Match *rt = NULL;
	sortBestMatches();
	if(bm.size() > 0){
		rt = bm[0];
	}
	return rt;
}

void BestMatches::sortBestMatches(){
	stable_sort(bm.begin(), bm.end(), PMatchComparator());
}

void BestMatches::analyzeBestHits(){

	/*
	 * NORMAL CASE
	 * 1 best match, 1 second best match
	 * several best matches, 1 second best match
	 * 
	 * SPECIAL CASES
	 * 0 matches
	 * 1 match
	 * several best matches
	 * different AA sequences correspond to same best score
	 * 
	 * SOLUTION
	 * Number of best matches means number of different sequences
	 *  
	 ****/
	
	sortBestMatches();

	//initalize variables
	thpm = 0;
	bestscore = 0;
	deltascore12 = 0;
	//for(int i=0; i < 20; i++){deltascores[i]=0;}
	for(int i=0; i <= ownerSpectrum->se1->bestmatches; i++){bestscores[i]=0;}
	random = false;
	show_hits = 0;
	tot_best_hits = 0;
	nb_prot = 0;
	nb_chromunspliced = 0;
	nb_chromspliced = 0;

	//helper variables	
	int delta_i = 1;
	string bestsequence = "";
	string currentsequence = "";
	
	//analyze matches
	for(int i = 0; i < bm.size(); i++){
		
		//set reference values, requesting here in for loop is safe regarding null pointer exception
		if(i == 0){
			bestscore = bm[0]->score;
			bestsequence = bm[0]->tuple->getAASeq();
			currentsequence = bestsequence;
			thpm = bm[0]->tuple->thparentmassMH;
			if(bm[0]->tuple->random == true) random = true; //causes a (desired) representation of 50% of random hits if the score does not work at all	
		}
		
		if(bm[i]->tuple->getAASeq() == bestsequence){
			/* 
			 * COUNT MULTIPLE BEST MATCHES
			 * one new rank per sequence, even if it is the same score
			 * 
			 */
			tot_best_hits++;
			show_hits = tot_best_hits + 1;
			if(bm[i]->tuple->protein != NULL) nb_prot++;
			if(bm[i]->tuple->chromosome != NULL){
				if(bm[i]->tuple->isSpliced()==true){
					nb_chromspliced++;
				}else{
					nb_chromunspliced++;
				}
			}
		}else if(bm[i]->tuple->getAASeq() != currentsequence){
			/* 
			 * COUNT INFERIOR MATCHES
			 * one new rank for every new sequence
			 */
			delta_i++; //different sequence matters, may be same score as previous
			if(delta_i == 2) deltascore12 = (bestscore - bm[i]->score);
			if(delta_i <= ownerSpectrum->se1->bestmatches){
				//deltascores[delta_i] = bestscore - bm[i]->score;
				bestscores[delta_i] = bm[i]->score;
			}
			currentsequence = bm[i]->tuple->getAASeq();
		}				
	}
}
}
