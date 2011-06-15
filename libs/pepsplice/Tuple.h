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

#ifndef TUPLE_H_
#define TUPLE_H_

#include <string>
#include <iostream>
#include <vector>
#include "Chromosome.h"
#include "Protein.h"
#include "pool.h"
#include "definitions.h"

using namespace std;
namespace Pepsplice{



/**
 * @author Franz Roos, ETH Zurich
 * @author Buote Xu (Buote.Xu@stud.uni-heidelberg.de)
 * @date   January, 2010
 * @brief  A tuple represents an AA sequence which will be modified (and penalized accordingly) in different ways before getting scored against a given spectrum
 *
 */
class Tuple
{
private:


public:
	//start changes BYR
	///deprecated constructor, still inside for compatibility reasons
	Tuple(double thpm, long len, long ps, long pe, long ss, long se, unsigned char *aaseq, int aaseqorig, Chromosome *c, Protein *p, bool rev, bool fp, bool trypstart, bool trypend);
    ///standard constructor, it is used to initialize a tuple from an entry in the fasta database.
    Tuple(double thpm, long len, long ps, long pe, long ss, long se, string & aaseq, int aaseqorig, Chromosome *c, Protein *p, bool rev, bool fp, bool trypstart, bool trypend);

    Tuple(const Tuple *tup); //copy constructor
	Tuple(const Tuple& tup); //copy constructor
	//~Tuple();
    ///overloaded new operator, needed to use the customized memory pool
    void * operator new(size_t);
    ///overloaded delete operator, needed to use the customized memory pool
    void operator delete(void * oldpointer);



    ///garbage collection will be done via reference counting
	long pointercount;
	long length;
	double thparentmassMH;
	char *ntsequence;
    long prefixstart;
	long prefixend;
	long suffixstart;
	long suffixend;
    Chromosome *chromosome;
	Protein *protein;
    bool chromosome_reverse;
	bool random;
	unsigned char modifications;
	unsigned char oxidizedmethionines;
    float penalty;
    ///index of the original, unmodified AA-sequence in a database which is logged later
	int aasequenceorig;
    bool trypticstart;
	bool trypticend;
	bool protNterm;
	bool protCterm;
	unsigned char charNterm;
	unsigned char charCterm;
	bool isSNP;
	//added BYR
	bool isSNP2;
    unsigned char SNPoldAA;
	unsigned char SNPnewAA;
    int SNPpos;
    int placeholder;
    ///originally a unsigned char *, the sequences will now be saved inside the class for a 2x performance improvement
    unsigned char aasequence[MAXSEQUENCESIZE]; //BX: Assuming that the biggest resulting Sequence has to fit into these 50 chars.

public:
    inline void increasePointerCount();
    inline void isBeingScored();
    inline void decreasePointerCount();
    bool decreasePointerCount2();
	inline void scoringFinished();
	void setChromosome(Chromosome *c, bool rev);
	void setProtein(Protein *p);
	void setReverse();
	void getNTSeqFromChromosome();
	string getNTSeq();
	inline string getAASeq();
	//string getAASeqOriginal();
	int getPS();
	int getPE();
	int getSS();
	int getSE();
	long getNTLength();
	void initializeThSpec(int len);
	void deleteThSpec();
	bool isSpliced();
	bool isDNA();
	int getPrefixLength();
	int getSuffixLength();
	void checkSequence();
	int getMissedCleavages();
	bool isSNPok(string label);
    //end changes
};


    void Tuple::increasePointerCount(){++pointercount;}

    void Tuple::isBeingScored(){++pointercount;}

    void Tuple::decreasePointerCount()
    {
        pointercount--;
        //if(pointercountreached2 == true) cout << "2" << flush;
        if(pointercount <= 0) delete this;
    }

    void Tuple::scoringFinished()
    {
        pointercount--;
        //if(pointercountreached2 == true) cout << "2" << flush;
        if(pointercount <= 0) delete this;
    }
    string Tuple::getAASeq(){
        return string(reinterpret_cast<const char*>(aasequence),length);
    }
}



#endif /*TUPLE_H_*/

