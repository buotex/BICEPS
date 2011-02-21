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

#ifndef DNAAA_H_
#define DNAAA_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <string>
#include <vector>

using namespace std;
namespace Pepsplice{
/* This class converts DNA triplets to amino acids and their masses
 * 
 * Content:
 * list of 64 triplets and corresponding AA (char) and corresponding masses (double, 4 decimal digits before, 6 digits after comma, exceeds 32 bits)
 * converter {A, C, T, G} to {0, 1, 2, 3}
 * 4*4*4 array containing AA for fast access with DNA triplets -> aa_char[1][3][1] gives AGA gives 
 * 4*4*4 array containing AA masses
 * return AA and mass
 * instantiate object only once, fill in array only once
 * use object many times by calling a function within: dnaToAA
 * 
 */

class DnaAA
{
	
private:
	void initialize();
	void initializeSNPs();
	double aa_monoisotopic(unsigned char aa);
	bool tripletis(char tripletA[3], const char tripletB[]);

public:
	DnaAA();
	double scaling_factor;
	double discretization_factor;
	char dna_triplets[4][4][4][3];
	unsigned char aa_chars[4][4][4];
	double aa_monomasses444[4][4][4];
	double aa_monomasses256[256];
	float  aa_monomasses256float[256];
	int    aa_integermasses256[256];
	double monomassH, monomassO, monomassN, monomassHO, monomassHOH, monomassHOHH, monomassNH3, monomassC13minusC12, monomassPhosphorylation;
	int    integermassH, integermassHOHH;
	
	//SNPs
	bool tripletSNPaa[4][4][4][256]; //given the original triplet, which AAs are possible doing one point mutation?
	bool aaSNPaa[256][256]; //given the original amino acid, which AAs are possible doing one point mutation? enumerates the (unknown) original triplets that are possible
	
	//AA modifications
	double aamod_monomasses[256]; //for each of 256 ASCII codes, the corresponding monoisotopic mass is specified here
	unsigned char aamod_aa[256]; //for each of 256 ASCII codes, the corresponding AA is specified here
	unsigned char aamod_type[256]; //modification type (internal, pepNterm, pepCterm, protNterm, protCterm)
	int aamod_startascii[256]; //for modifiable amino acids, the first ASCII entry of the modifications is given here
	int aamod_nbmodperaa[256]; //for each of 20 amino acids, the number of modifications is specified here
		
	//functions
	char num_nt(char n);
	char nt_num(char n);
	unsigned char nt_aa(char triplet[3]);
	char nt_revnt(char n);
	char nt_nt_clean(char n);
	void checkInitialization();
	string translate(string ntseq);
	double getPRMfromAA(string aaseq);
	double getParentMassMH(unsigned char *aaseq, long len);
	double getParentMassMH(string aaseq);
	void loadAAModifications();
	double string_to_double(string s);
	string intToString(int i);
	string charseq_to_modnumseq(string s1, unsigned char charNterm, unsigned char charCterm);
	string charseq_to_lowercaseseq(string s1);
	double aa_pampenalty(unsigned char aa, unsigned char bb);
};
}
#endif /*DNAAA_H_*/
