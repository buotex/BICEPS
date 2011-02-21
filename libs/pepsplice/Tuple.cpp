#include "Tuple.h"
#include <algorithm>
#include <cstring>

namespace Pepsplice{
    
    
    extern Pool<Pepsplice::Tuple, POOLSIZE> __POOL__;    
    
    
    void * Tuple::operator new(size_t)
    {
        return __POOL__.getNewObject();
    }
    void Tuple::operator delete(void * oldpointer)
    {
        __POOL__.deleteObject(oldpointer);
    }
    
    
    Tuple::Tuple(double thpm, long len, long ps, long pe, long ss, long se, string & aaseq, int aaseqorig,Chromosome *c, Protein *p, bool rev, bool rnd, bool trypstart, bool trypend)
    {
        //IMPORTANT: CONSTRUCTOR MUST BE SYNCHRONIZED WITH COPY CONSTRUCTOR Tuple::Tuple(Tuple *tup)
        
        
        pointercount = 0; //0 means no one needs the tuple
        thparentmassMH = thpm;
        length = len;
        //	aasequence = aaseq;
        aaseq.copy((char* )aasequence, len, ps);
        //changes BYR added aasequence original
        aasequenceorig=aaseqorig;
        ntsequence = NULL;
        prefixstart = ps;
        prefixend = pe;
        suffixstart = ss;
        suffixend = se;
        chromosome = c;
        protein = p;
        //start changes BYR
        chromosome_reverse = false; //rev;
        //chromosome_reverse=rev;
        //end changes BYR
        random = rnd;
        modifications = 0;
        oxidizedmethionines = 0;
        penalty = 0;
        trypticstart = trypstart;
        trypticend = trypend;
        protNterm = false;
        protCterm = false;
        charNterm = 0;
        charCterm = 0;
        if(protein != NULL && prefixstart == 0) protNterm = true;
        if(protein != NULL && protein->getLength() - 1 == prefixend){
            //cout << "\nProtein: " << protein->getLength() << " " << getAASeq() << " " << prefixend << "\n" << protein->protseq;
            protCterm = true;
        }
        //if (protein != NULL) fastaid = protein->fastaid;
        isSNP = false;
        isSNP2 = false;
        SNPoldAA = '-';
        SNPnewAA = '-';
        SNPpos = -1;
        //IMPORTANT: CONSTRUCTOR MUST BE SYNCHRONIZED WITH COPY CONSTRUCTOR Tuple::Tuple(Tuple *tup)
    }
    
    
    Tuple::Tuple(double thpm, long len, long ps, long pe, long ss, long se, unsigned char *aaseq, int aaseqorig,Chromosome *c, Protein *p, bool rev, bool rnd, bool trypstart, bool trypend)
    {
        //IMPORTANT: CONSTRUCTOR MUST BE SYNCHRONIZED WITH COPY CONSTRUCTOR Tuple::Tuple(Tuple *tup)
        
        
        pointercount = 0; //0 means no one needs the tuple
        thparentmassMH = thpm;
        length = len;
        //	aasequence = aaseq;
        memcpy(aasequence, aaseq, length * sizeof(unsigned char));
        
        //changes BYR added aasequence original
        aasequenceorig=aaseqorig;
        ntsequence = NULL;
        prefixstart = ps;
        prefixend = pe;
        suffixstart = ss;
        suffixend = se;
        chromosome = c;
        protein = p;
        //start changes BYR
        chromosome_reverse = false; //rev;
        //chromosome_reverse=rev;
        //end changes BYR
        random = rnd;
        modifications = 0;
        oxidizedmethionines = 0;
        penalty = 0;
        trypticstart = trypstart;
        trypticend = trypend;
        protNterm = false;
        protCterm = false;
        charNterm = 0;
        charCterm = 0;
        if(protein != NULL && prefixstart == 0) protNterm = true;
        if(protein != NULL && protein->getLength() - 1 == prefixend){
            //cout << "\nProtein: " << protein->getLength() << " " << getAASeq() << " " << prefixend << "\n" << protein->protseq;
            protCterm = true;
        }
        //if (protein != NULL) fastaid = protein->fastaid;
        isSNP = false;
        isSNP2 = false;
        SNPoldAA = '-';
        SNPnewAA = '-';
        SNPpos = -1;
        //IMPORTANT: CONSTRUCTOR MUST BE SYNCHRONIZED WITH COPY CONSTRUCTOR Tuple::Tuple(Tuple *tup)
    }
    
    Tuple::Tuple(const Tuple *tup)
    {
        // fastaid = tup->fastaid;
        pointercount = 0; //0 means no one needs the tuple
        thparentmassMH = tup->thparentmassMH;
        length = tup->length;
        //	aasequence = new unsigned char[length];
        //memcpy(aasequence, tup->aasequence, length * sizeof(unsigned char)); //to copy the array values
        memcpy(aasequence, tup->aasequence, length * sizeof(unsigned char)); //has to be same length as the complete array

        //copy(tup->aasequence, tup->aasequence+length, aasequence);
        //for(int i = 0; i < length; i++){
        //	aasequence[i] = tup->aasequence[i];
        //}
        //changes BYR added original sequence
        //	aasequenceorig = new unsigned char[length];
        //    copy(tup->aasequenceorig, tup->aasequenceorig+length, aasequenceorig);
        aasequenceorig = tup->aasequenceorig;
        //for (long i = 0; i < length; ++i)
        //{
        //   aasequenceorig[i] = tup->aasequenceorig[i];
        //}
        ntsequence = tup->ntsequence;
        prefixstart = tup->prefixstart;
        prefixend = tup->prefixend;
        suffixstart = tup->suffixstart;
        suffixend = tup->suffixend;
        chromosome = tup->chromosome;
        protein = tup->protein;
        chromosome_reverse = tup->chromosome_reverse;//tup->chromosome_reverse;
        random = tup->random;
        modifications = tup->modifications;
        oxidizedmethionines = tup->oxidizedmethionines;
        penalty = tup->penalty;
        trypticstart = tup->trypticstart;
        trypticend = tup->trypticend;
        protNterm = tup->protNterm;
        protCterm = tup->protCterm;
        charNterm = tup->charNterm;
        charCterm = tup->charCterm;
        isSNP = tup->isSNP;
        isSNP2 = tup->isSNP2;
        SNPoldAA = tup->SNPoldAA;
        SNPnewAA = tup->SNPnewAA;
        SNPpos = tup->SNPpos;
    }
    
    Tuple::Tuple(const Tuple& tup)
    {
        // fastaid = tup.fastaid;
        pointercount = 0; //0 means no one needs the tuple
        thparentmassMH = tup.thparentmassMH;
        length = tup.length;
        //aasequence = new unsigned char[length];
        //copy(tup.aasequence, tup.aasequence+length, aasequence);
        memcpy(aasequence, tup.aasequence, length * sizeof(unsigned char));
        
        //	for(long i = 0; i < length; i++){
        //		aasequence[i] = tup.aasequence[i];
        //	}
        //changes BYR added original sequence
        //aasequenceorig=tup.aasequenceorig;
        aasequenceorig = tup.aasequenceorig;
        //for (long i = 0; i < length; ++i)
        //{
        //   aasequenceorig[i] = tup->aasequenceorig[i];
        //}
        ntsequence = NULL;
        prefixstart = tup.prefixstart;
        prefixend = tup.prefixend;
        suffixstart = tup.suffixstart;
        suffixend = tup.suffixend;
        chromosome = tup.chromosome;
        protein = tup.protein;
        chromosome_reverse = tup.chromosome_reverse;//tup->chromosome_reverse;
        random = tup.random;
        modifications = tup.modifications;
        oxidizedmethionines = tup.oxidizedmethionines;
        penalty = tup.penalty;
        trypticstart = tup.trypticstart;
        trypticend = tup.trypticend;
        protNterm = tup.protNterm;
        protCterm = tup.protCterm;
        charNterm = tup.charNterm;
        charCterm = tup.charCterm;
        isSNP = tup.isSNP;
        isSNP2 = tup.isSNP2;
        SNPoldAA = tup.SNPoldAA;
        SNPnewAA = tup.SNPnewAA;
        SNPpos = tup.SNPpos;
    }
    
    
    //Tuple::~Tuple()
    //{
	//if(pointercount > 0) cout << "\nTuple is getting deleted although pointer getTotal is : " << pointercount << flush;
	//if(aasequence != NULL){
	//	delete[] aasequence;
	//	aasequence = NULL;
	//}else{
	//	cout << "\nTuple destructor: aaseq is NULL" << flush;
	//}
	
    //Changes by BX: For us, ntsequence is always NULL
	//if(ntsequence != NULL){
	//	delete[] ntsequence;
	//	ntsequence = NULL;
	//}
	
	/*
     if(thspec != NULL){
     delete[] thspec;
     thspec = NULL;
     }
     */
    //}
    
    
    bool Tuple::decreasePointerCount2(){
        pointercount--;
        if (pointercount <=0)
        {
            delete this;
            return true;
        }
        else return false;
    }
    
    
    
    
    void Tuple::setChromosome(Chromosome *c, bool rev){
        chromosome = c;	
        chromosome_reverse = rev;
    }
    
    void Tuple::setProtein(Protein *p){
        protein = p;	
    }
    
    void Tuple::getNTSeqFromChromosome()
    {
        if(ntsequence == NULL && chromosome != NULL){
            ntsequence = chromosome->getNTSeqForTuple(length, prefixstart, prefixend, suffixstart, suffixend, chromosome_reverse);
            
            /*
             int cs = chromosome->size;
             int ps = prefixstart;
             int se = suffixend;
             ps = prefixend - 47;
             se = suffixstart + 47;
             if(ps < 1) ps = 1;
             if(se > cs - 1) se = cs - 1;
             int plen = prefixend - ps + 1;
             int slen = se - suffixstart + 1;
             excessntseq = chromosome->getNTSeqForTuple((plen + slen)/3, ps, prefixend, suffixstart, se, chromosome_reverse);
             */
            
        }
    }
    
    string Tuple::getNTSeq()
    {
        string s = "";
        getNTSeqFromChromosome();
        if(ntsequence != NULL){
            for(int i = 0; i < length * 3; i++){
                s+=ntsequence[i];
            }
        }
        return s;
    }
    
    /*
     string Tuple::getExcessNTSeq()
     {
     string s = "";
     getNTSeqFromChromosome();
     if(excessntseq != NULL){
     for(int i = 0; i < 2*48; i++){
     s+=excessntseq[i];
     }
     }
     return s;
     }
     */
    
    
    /*
     void Tuple::initializeThSpec(int len){
     //cout << "\nTuple::initializeThSpec len: " << len;
     thspec = new int[len];
     for(int i = 0; i < len; i++){thspec[i]=0;}
     }
     
     void Tuple::deleteThSpec(){
     delete[] thspec;
     thspec = NULL;
     }
     */
    
    long Tuple::getNTLength(){
        return (suffixend - suffixstart + 1) + (prefixend - prefixstart + 1); //getTotal length of nucleotide sequence
    }
    
    bool Tuple::isSpliced(){
        bool rt = false;
        if(getPrefixLength() > 0 && getSuffixLength() > 0) rt = true;
        return rt;
    }
    
    bool Tuple::isDNA(){
        bool rt = false;
        if(chromosome != NULL) rt = true;
        return rt;
    }
    
    int Tuple::getPrefixLength(){
        int pdiff = prefixend - prefixstart;
        if(suffixend >= prefixstart){
            return pdiff + 1;
        }else{
            return (-pdiff) + 1;
        }
    }
    
    int Tuple::getSuffixLength(){
        int sdiff = suffixend - suffixstart;
        if(suffixend >= prefixstart){
            return sdiff + 1;
        }else{
            return (-sdiff) + 1;
        }
    }
    
    int Tuple::getPS(){
        int ps = prefixstart;
        if(chromosome_reverse == true && chromosome != NULL) ps = chromosome->getFullChromSize() - prefixstart;
        return ps;
    }
    
    int Tuple::getPE(){
        int pe = prefixend;
        if(chromosome_reverse == true && chromosome != NULL) pe = chromosome->getFullChromSize() - prefixend;
        return pe;
    }
    
    int Tuple::getSS(){
        int ss = suffixstart;
        if(chromosome_reverse == true && chromosome != NULL) ss = chromosome->getFullChromSize() - suffixstart;
        return ss;
    }
    
    int Tuple::getSE(){
        int se = suffixend;
        if(chromosome_reverse == true && chromosome != NULL) se = chromosome->getFullChromSize() - suffixend;
        return se;
    }
    
    void Tuple::checkSequence(){
        //check sequence for illegal characters
        for(int i = 0; i < length; i++){
            if(aasequence[i] < 65 || aasequence[i] > 89) cout << "\nPM:" << thparentmassMH << " length:" << length << " aasequence:" << getAASeq();
        }
    }
    
    int Tuple::getMissedCleavages(){
        int missedcleavages = 0;
        //do not count last one, this is not a missed one
        for(int i = 0; i < length - 1; i++){
            if(aasequence[i]=='R' || aasequence[i]=='K') missedcleavages++;
        }
        return missedcleavages;
    }
    
    bool Tuple::isSNPok(string label){
        bool r = true;
        if(isSNP && random==false && aasequence[SNPpos]!=SNPnewAA && aasequence[SNPpos] <= 127) r = false;
        if(r == false){
            //start changes BYR
			//cout << "\nTuple::isSNPok():" << label << " rev" << (bool)random << " this:" << (int)this  << " SNPpos:" << SNPpos << " " << SNPoldAA << " " << SNPnewAA << " " << getAASeq();
			cout << "TROUBLE IN TUPLE";
			//for(int i = 0; i < length; i++){cout << " i:" << i << " c:" << (int)aasequence[i];}
        }
        return r;
    }
    
    //unsigned int Tuple::getID(){
    //
    //    return this-> fastaid;
    //
    //}
}