#include "Protein.h"
namespace Pepsplice{
Protein::Protein(string id)
{
    
	protid = id;
	protseq = "";
	initializeConfidence();
    
}

//
//Protein::Protein(string id, unsigned int fastaid_){
//
//
//  //  fastaid = fastaid_;
//	protid = id;
//	protseq = "";
//	initializeConfidence();
//
//}
// 
    
    Protein::Protein(string id_, string sequence_): protid(id_), protseq(sequence_)
{
    initializeConfidence();
}


Protein::~Protein()
{
}

void Protein::appendSubSeq(string & parsedline)
{
    //the checks should already be done during the main function and the loading fasta stage, so we can improve the speed here.
    protseq = parsedline;
}

int Protein::getLength()
{
    return protseq.size();	
}

void Protein::initializeConfidence(){
    risk = 1;
    matches.clear();
}

void Protein::adjustProteinConfidence(double peptideconfidence, int matchindex){
    risk = risk * (1 - peptideconfidence);
    matches.push_back(matchindex);
}

double Protein::getProteinConfidence(){
    return 1-risk;
}
}