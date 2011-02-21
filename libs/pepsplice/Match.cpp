#include "Match.h"
namespace Pepsplice{
Match::Match(Spectrum *sp1, Tuple *tp1, float sc1, int spc1)
{
	spectrum = sp1;
	tuple = tp1;
	score = sc1;
	sharedpeakcount = spc1;
	tuple->increasePointerCount();
	if(tuple->pointercount < 2) cout << "\nMatch::Match() tuple was added to match and has now pointercount: " << tuple->pointercount;
	confidence = 0;
	competitors = 0;
}

    
    
Match::~Match()
{
	//if(tuple->pointercount == 2) cout << " 2to1 " << flush;
	tuple->decreasePointerCount(); 
}

/*
void Match::setNKnk(int N, int K, int n, int k){
	this->N = N;
	this->K = K;
	this->n = n;
	this->k = k;
}
*/
}
