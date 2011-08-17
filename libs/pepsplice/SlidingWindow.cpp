#include "SlidingWindow.h"
namespace Pepsplice{
SlidingWindow::SlidingWindow(Services *se0)
{

	se1 = se0;
	timer1 = Timer();
	dnaAA1 = DnaAA();
	
	hot = false; //within hot spot, set by Chromosomes
	
	//parameters
	plenmax = se1->totlenmax;	//maximum prefix length
	slenmax = se1->totlenmax;	//maximum suffix length
	circleextra = 10;		//circle buffer extra length for safety
	circlelength = plenmax + se1->gaplenmax + slenmax + circleextra;
	
	ntsequence 		= new char[se1->totlenmax];
	aasequence		= new unsigned char[se1->totlenmax / 3];
	aamasssequence  = new float[se1->totlenmax / 3];
	
	//arrays for circular access, mirror of sliding window on DNA
	//long, double
	aamass 	       = new float[circlelength];  for(int k = 0; k < circlelength; k++){aamass[k]=0;}
	massruler      = new double[circlelength]; for(int k = 0; k < circlelength; k++){massruler[k]=0;}
	globalposition   = new long[circlelength]; for(int k = 0; k < circlelength; k++){globalposition[k]=-1;}
	
	//char
	nt 				= new char[circlelength];	for(int k = 0; k < circlelength; k++){nt[k]='N';}
	ntnum 			= new char[circlelength];	for(int k = 0; k < circlelength; k++){ntnum[k]=-1;}
	aa 				= new unsigned char[circlelength];	for(int k = 0; k < circlelength; k++){aa[k]='X';}
	readingframe 	= new char[circlelength];	for(int k = 0; k < circlelength; k++){readingframe[k]=-1;}

	//local positions
	sincestop        = new long[circlelength]; for(int k = 0; k < circlelength; k++){sincestop[k]=0;}	
	sinceprefixstart = new long[circlelength]; for(int k = 0; k < circlelength; k++){sinceprefixstart[k]=1;}
	sinceprefixend   = new long[circlelength]; for(int k = 0; k < circlelength; k++){sinceprefixend[k]=1;}
	sincesuffixstart = new long[circlelength]; for(int k = 0; k < circlelength; k++){sincesuffixstart[k]=1;}
	sincesuffixend   = new long[circlelength]; for(int k = 0; k < circlelength; k++){sincesuffixend[k]=1;}

	isstop        = new bool[circlelength]; for(int k = 0; k < circlelength; k++){isstop[k]=true;}	
	isprefixstart = new bool[circlelength]; for(int k = 0; k < circlelength; k++){isprefixstart[k]=false;}
	isprefixend   = new bool[circlelength]; for(int k = 0; k < circlelength; k++){isprefixend[k]=false;}
	issuffixstart = new bool[circlelength]; for(int k = 0; k < circlelength; k++){issuffixstart[k]=false;}
	issuffixend   = new bool[circlelength]; for(int k = 0; k < circlelength; k++){issuffixend[k]=false;}


	
	//access arrays in circular way via mod
	i = 1; imodz = i + 1; imoda = i; imodb = i - 1; imodc = i - 2; imodd = i - 3; //navigate in circular arrays
	
	//4 positions: prefix start, prefix end, suffix start, suffix end
	ps = 0;	pe = 0;	ss = 0;	se = 0;
	monoparentmassMH = 0;
	
	//LOGGING
		//distribution1 = new Distribution(-0.5, max_monoparentmass_M-0.5, long(max_monoparentmass_M)); //parent mass distribution
	illegalnucleotides = 0;
	outputlevel = se1->outputlevel;
	distribution1 = new Distribution(0, plenmax, plenmax + 1, 1);
	j = 0; //increment for each parent mass calculation

	nextout = 100;	
}

SlidingWindow::~SlidingWindow()
{
	delete[] aamass;
	delete[] massruler;
	delete[] globalposition;			
	delete[] nt;
	delete[] ntnum;
	delete[] aa;
	delete[] readingframe;
	delete[] ntsequence;
	delete[] aasequence;
	delete[] aamasssequence;
	delete[] sincestop;
	delete[] sinceprefixstart;
	delete[] sinceprefixend;
	delete[] sincesuffixstart;
	delete[] sincesuffixend;
	delete[] isstop;
	delete[] isprefixstart;
	delete[] isprefixend;
	delete[] issuffixstart;
	delete[] issuffixend;
	
	delete distribution1;
}

bool SlidingWindow::addNextNucleotide(char newNucleotide)
{
	bool rt = false;
	
		//calculate local imod derived from i
		imoda = i % circlelength;
		imodb = cheapUpMod(imoda - 1);
		imodc = cheapUpMod(imoda - 2);
		imodd = cheapUpMod(imoda - 3);
		
		//parse letter, update local information
		ntnum[imoda] = dnaAA1.nt_num(newNucleotide);
		
		//convert illegal nucleotides
		if(ntnum[imoda] >= 0){
			nt[imoda] = newNucleotide;
		}else{
			illegalnucleotides++;
			nt[imoda] = 'N';
		}
		
		//if triplet is legal, add to all levels, else treat as stop codon
		if(ntnum[imoda] >= 0 && ntnum[imodb] >= 0 && ntnum[imodc] >= 0){
			aa[imoda]      = dnaAA1.aa_chars[ntnum[imodc]][ntnum[imodb]][ntnum[imoda]];
			
			//mass ruler
			aamass_temp    = dnaAA1.aa_monomasses444[ntnum[imodc]][ntnum[imodb]][ntnum[imoda]];
			aamass[imoda]  = float(aamass_temp);
			if(aamass_temp > 0){
				massruler[imoda] = massruler[imodd] + aamass_temp;
			}else{
				massruler[imoda] = 0;
				if(outputlevel > 3 && aa[imoda] != '-') cout << "\nSlidingWindow::addNextNucleotide: aamass_temp <= 0 " << aamass_temp << " aa" << aa[imoda] << i;
			}
		}else{
			aa[imoda] = 'X';
			aamass[imoda] = 0;
			massruler[imoda] = 0;
		}

			
			//***********PATTERN DETECTION**********
			//    |---------->        >------------|
			// ooo            oo    oo           ooo
			// K/R   prefix   GT    AG   suffix  K/R


			// iiii
			// mmmm
			// oooo
			// dddd
			// dcba


			//stop codons
			if(aa[imoda] == '-' || aa[imoda] == 'X'){
				sincestop[imoda] = 0;
				isstop[imoda] = true;
			}else{
				sincestop[imoda] = sincestop[imodd] + 3; //valid within reading frame!
				isstop[imoda] = false;
			}	

			//prefix start: K or R
			if(aa[imodb] == 'K' || aa[imodb] == 'R'){
				sinceprefixstart[imoda] = 0;
				isprefixstart[imoda] = true;
			}else{
				sinceprefixstart[imoda] = sinceprefixstart[imodd] + 3; //dependent on suffix length, valid within reading frame!
				isprefixstart[imoda] = false;
			}
			
			//prefix end: GT
			if(se1->noprefixend == true || (nt[imodb] == 'G' && nt[imoda] == 'T')){
				sinceprefixend[imodc] = 0;
				isprefixend[imodc] = true;
			}else{
				sinceprefixend[imodc] = sinceprefixend[imodd] + 1;
				isprefixend[imodc] = false;
			}
	
			//suffix start: AG
			if(se1->nosuffixstart == true || (nt[imodc] == 'A' && nt[imodb] == 'G')){
				sincesuffixstart[imoda] = 0;
				issuffixstart[imoda] = true;
			}else{
				sincesuffixstart[imoda] = sincesuffixstart[imodb] + 1;
				issuffixstart[imoda] = false;
			}

			//suffix end: K or R
			if(aa[imoda] == 'K' || aa[imoda] == 'R'){
				sincesuffixend[imoda] = 0;
				issuffixend[imoda] = true;
			}else{
				sincesuffixend[imoda] = sincesuffixend[imodb] + 1;
				issuffixend[imoda] = false;
			}
		
		globalposition[imoda] = i;
		readingframe[imoda] = i % 3;
		
		//cout.precision(8);		
		//cout  << "\n" << nt[imodd] << " " << aa[imodd] << "\t" << cheapMod(i) << "\t" << globalposition[imodd] << "\t" << (int)imodd << "\t" << (double)massruler[imodd]/se1->dnaAA1->scaling_factor << "\t" << (int)readingframe[imodd] << "  " << sincestop[imodd] << "\t" << sinceprefixstart[imodd] << "\t" << sinceprefixend[imodd] << "\t" << sincesuffixstart[imodd] << "\t" << sincesuffixend[imodd];		
		
		//cout  << "\n" << nt[imodd] << " " << aa[imodd] << "\t" << aamass[imodd] << "\t" << cheapMod(i) << " " << globalposition[imodd] << "  " << imodd << "  " << massruler[imodd] << "\t" << readingframe[imodd] << "  " << sincestop[imodd] << "\t" << sinceprefixstart[imodd] << "   " << sinceprefixend[imodd] << "   " << sincesuffixstart[imodd] << "   " << sincesuffixend[imodd];
		//cout << "\n" << nt[imoda] << (int)ntnum[imoda] << " " << aa[imoda] << "\t" << aamass[imoda] << "\t" << i << "\t" << (double)massruler[imoda];
		
		if(sincesuffixend[imodd]==0) rt = true;
		if(rt == true){
			//enumerate unspliced tryptic peptides anyway
			if(se1->doWholeGenome == true) enumerateEntirePeptides();
			
			//spliced peptides only if parameter == true
			if(se1->spliced == true || hot == true) enumerateSplicedPeptides();
		}
		
		i++;
		return rt;
}

void SlidingWindow::setTuples(Tuples *tps)
{
	tuples = tps;
}

void SlidingWindow::setChromosome(Chromosome *c, bool rev)
{
	activechromosome = c;
	chr_reverse = rev;
}

double SlidingWindow::getBrokenTripletMass()
{
	char ntnum1, ntnum2, ntnum3;
	if(pexcessnt == 0){
		return 0;
	}
	
	if(pexcessnt == 1){
		ntnum1 = ntnum[pe];
		ntnum2 = ntnum[ss];
		ntnum3 = ntnum[cheapDownMod(ss + 1)];
		return dnaAA1.aa_monomasses444[ntnum1][ntnum2][ntnum3];
	}
	if(pexcessnt == 2){
		ntnum1 = ntnum[cheapUpMod(pe - 1)];
		ntnum2 = ntnum[pe];
		ntnum3 = ntnum[ss];
		return dnaAA1.aa_monomasses444[ntnum1][ntnum2][ntnum3];
	}
    return 0;
}

unsigned char SlidingWindow::getBrokenTripletAA()
{
	char ntnum1, ntnum2, ntnum3;
	if(pexcessnt == 0){
		return 0;
	}
	
	if(pexcessnt == 1){
		ntnum1 = ntnum[pe];
		ntnum2 = ntnum[ss];
		ntnum3 = ntnum[cheapDownMod(ss + 1)];
		return dnaAA1.aa_chars[ntnum1][ntnum2][ntnum3];
	}
	if(pexcessnt == 2){
		ntnum1 = ntnum[cheapUpMod(pe - 1)];
		ntnum2 = ntnum[pe];
		ntnum3 = ntnum[ss];
		return dnaAA1.aa_chars[ntnum1][ntnum2][ntnum3];
	}
    return 0;
}


string SlidingWindow::getAASequence(){
	//is a bit slower than getAASequenceFast used to be (see PepSplice0009) due to string vs. char array difference
	
	long k, kstop;
	string sequence;

	if(outputlevel > 5 && totlen % 3 != 0){
		cout << "\nSlidingWindow::getAASequence: totlen is not divisible by 3. totlen: " << totlen << " i: " << i;
		return sequence;
	}
	
	//fill in prefix
	if(plen >= 3){
		k = cheapDownMod(ps + 2);
		kstop = cheapUpMod(pe - pexcessnt);
		while(true){
			sequence.push_back(aa[k]);
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}
	
	//if there is a broken triplet, add it
	if(pexcessnt > 0){
		sequence.push_back(getBrokenTripletAA());
	}
	
	//fill in suffix
	if(slen >= 3){
		k = cheapDownMod(ss + 2 + sexcessnt);
		kstop = se;
		while(true){
			sequence.push_back(aa[k]);
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}

	if(outputlevel > 5 && sequence.length() != totlen / 3){
		cout << "\ni:" << i << " SlidingWindow::getAASequence(): sequence too short or too long! " << "sequence.size(): " << sequence.size() << " " << " totlen: " << totlen << " ps: " << ps << " pe: " << pe << " ss: " << ss << "se: " << se; 
		for(int k = 0; k < totlen/3; k++){cout << k << "." << (long)sequence[k] << " ";}	
	}

	return sequence;
}

string SlidingWindow::getAASequenceFast(){

	long k, kstop, s = 0;
	unsigned char *sequence = aasequence;
	string string_sequence;
	
	//fill in prefix
	if(plen >= 3){
		k = cheapDownMod(ps + 2);
		kstop = cheapUpMod(pe - pexcessnt);
		while(true){
			sequence[s] = aa[k];
			s++;
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}
	
	//if there is a broken triplet, add it
	if(pexcessnt > 0){
		sequence[s] = getBrokenTripletAA();
		s++;
	}
	
	//fill in suffix
	if(slen >= 3){
		k = cheapDownMod(ss + 2 + sexcessnt);
		kstop = se;
		while(true){
			sequence[s] = aa[k];
			s++;
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}
	string seq(sequence + 0, sequence + totlen/3);
	return seq;
}

unsigned char* SlidingWindow::getAASequenceFastChar(){

	long k, kstop, s = 0;
	unsigned char *sequence = new unsigned char[totlen/3];
	
	//fill in prefix
	if(plen >= 3){
		k = cheapDownMod(ps + 2);
		kstop = cheapUpMod(pe - pexcessnt);
		while(true){
			sequence[s] = aa[k];
			s++;
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}
	
	//if there is a broken triplet, add it
	if(pexcessnt > 0){
		sequence[s] = getBrokenTripletAA();
		s++;
	}
	
	//fill in suffix
	if(slen >= 3){
		k = cheapDownMod(ss + 2 + sexcessnt);
		kstop = se;
		while(true){
			sequence[s] = aa[k];
			s++;
			if(k == kstop){
				break;
			}else{
				k = cheapDownMod(k + 3);
			}
		}
	}
	
	return sequence;
}

string SlidingWindow::getNTSequence(){
	long k, s = 0;
	
	//prefix
	k = ps;
	if(plen > 0){
		while(true){
			ntsequence[s] = nt[k];
			s++;
			if(k == pe) break;
			k = cheapDownMod(k + 1);
		}
	}
	
	//suffix
	k = ss;
	if(slen > 0){
		while(true){
			ntsequence[s] = nt[k];
			s++;
			if(k == se) break;
			k = cheapDownMod(k + 1);
		}
	}
	
	if(outputlevel > 5 && s != totlen) cout << "\nSlidingWindow::getNTSequence() i:" << i << " sequence too long or too short!" << ntsequence;
	string seq(ntsequence + 0, ntsequence + totlen);
	return seq;
}

void SlidingWindow::showAAMassSequenceSlow(){

	long pos;
	long count;
	unsigned char* sequence =new unsigned char[plen + slen];
	double aam = 0;
	double pepresm = 0;
	
	pos = ps;
	count = 0;
	cout << "\n";
	while(count < plen + slen){
		pos = cheapMod(pos);
		sequence[count] = ntnum[pos];
		cout << (long)sequence[count];
		count++;
		if(pos == pe){
			pos = ss;
		}else{
			pos++;
		}
	}
	for(int k=2; k <= plen + slen; k = k + 3){
		aam = dnaAA1.aa_monomasses444[sequence[k - 2]][sequence[k-1]][sequence[k-0]];
		cout << "\nSLOW " << (long)sequence[k-2] << (long)sequence[k-1] << (long)sequence[k] << " " << aam;
		pepresm = pepresm + aam;
	}
	
	if(true || abs(pepresm - peptideresiduemass) > 0.01){
		cout << "\nSLOW i: " << i << " ps:" << ps << " pe:" << pe << " ss:" << ss << " se:" << se;
		cout << "\nSLOW pepresm:" << pepresm;
		cout << "\n\n";
	}
    delete sequence;
}

long SlidingWindow::cheapMod(long x)
{
	if (x < 0)             x = x + circlelength;
	if (x >= circlelength) x = x - circlelength;
	return x;
}

long SlidingWindow::cheapUpMod(long x)
{
	if (x < 0) x = x + circlelength;
	return x;
}

long SlidingWindow::cheapDownMod(long x)
{
	if (x >= circlelength) x = x - circlelength;
	return x;
}

long SlidingWindow::cheap3Mod(long x)
{
	if (x < 0)  x = x + 3;
	if (x >= 3) x = x - 3;
	return x;
}

void SlidingWindow::showParameters()
{
	cout << "\n\nPARAMETERS:";
	cout << "\nplenmin:" << se1->plenmin;
	cout << "\nplenmax:" << plenmax;
	cout << "\ngaplenmin:" << se1->gaplenmin;
	cout << "\ngaplenmax:" << se1->gaplenmax;
	cout << "\nslenmin:" << se1->slenmin;
	cout << "\nslenmax:" << slenmax;
	cout << "\ntotlenmax:" << se1->totlenmax;
	cout << "\ncirclelength:" << circlelength;
	cout << "\n";
}

void SlidingWindow::writeCircleContent()
{
	cout << "\n\n\n\nSlidingWindow::writeCircleContent() plen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se << " imodd: " << imodd;
	for(int i = 0; i < circlelength; i++){
		cout << "\n";
		cout << "i nt aa pos: " << i << " " << nt[i] << " " << aa[i] << " " << globalposition[i];
		cout << "  since ps pe ss se stop: " << sinceprefixstart[i] << " " << sinceprefixend[i] << " " << sincesuffixstart[i] << " " << sincesuffixend[i] << " " << sincestop[i];
		cout << "\t is ps pe ss se stop: " << isprefixstart[i] << " " << isprefixend[i] << " " << issuffixstart[i] << " " << issuffixend[i] << " " << isstop[i];
	}
}

void SlidingWindow::writeProgress()
{
	//cout << "\nSW ";
	//cout << "i-1: " << i - 1 << "   j:" << j;
	//cout << " " << nt[imoda] << " " << ps << " " << pe << " " << ss << " " << se;
	//cout << "   totsum:" << (*distribution1).totalsum;
	//cout << "   time:" << timer1.timeSinceStart();
	cout << "\nsec/MB:" << timer1.timeSinceStart() / i * 1000000;
	//cout << "   j/sec:" << j / timer1.timeSinceStart();
	//cout << "   illeg nt: " << illegalnucleotides;
	//distribution1->writeDistribution();
}

void SlidingWindow::enumerateSplicedPeptides() //never actually executed in our code
{
	//int k;
	se = cheapMod(imodd); //suffix end is given by position of neighborhood
	
	// -pe ss-se determine reading frame of ps; break makes sense if stop codon is detected
	// pe3 can be calculated from the three
	
	//LOOP TO GET SUFFIX START
	ss_old = se;
	while(true)
	{
		//cout << "\nSS1 plen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se;
		if(sincesuffixstart[ss_old] >= circlelength) break;				// (A) must not hop further than circle length
		ss = cheapUpMod(ss_old - sincesuffixstart[ss_old]);				//(1) get suffix start
		if(globalposition[ss] > globalposition[ss_old]) break;			// (B) must not hop over circle boundary
		if(sincesuffixstart[ss] != 0) break;							// (C) suffix start must be valid pattern
		//cout << "\nSS2 plen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se;
		
		slen = globalposition[se] - globalposition[ss] + 1;				//(2) get suffix length
		if(slen > slenmax) break;										// (D) suffix must not be too long 
		se3 = se;														//(3) get mass ruler position for suffix end
		if(slen > sincestop[se3]){break;}								// (F) suffix must not contain stop codons

		sexcessnt = cheap3Mod(readingframe[se] - readingframe[ss] + 1);	//(4) get excess nucleotides in suffix
		pexcessnt = cheap3Mod(3 - sexcessnt);							//(4) get excess nucleotides in prefix (depends on suffix)
		ss3 = cheapMod(ss + sexcessnt - 1);								//(5) get mass ruler position for suffix start
		  /*ss3******************************************************
		  * ACgTAcGT  ss=2  se=5  slen=4  sexcessnt = 1  ss3=2
		  * 01234567  knownsuffixmass = massruler[5] - massruler[3-1]
		  **********************************************************/

		knownsuffixmass = massruler[se3] - massruler[ss3];				//(6) get known suffix mass
		if(knownsuffixmass > se1->max_monoparentmassMH) break;			// (G) suffix mass must not exceed maximum peptide mass
		if(slen >= se1->slenmin)												// (E) if suffix is too short, skip if, but extend anyway
		{
			
			//LOOP TO GET PREFIX END
			pe_old = ss;
			while(true) // prefix end must be valid; breaks in this loop increment ss
			{

				if(sinceprefixend[pe_old] >= circlelength) break;				// (A) must not hop further than circle length
				pe = cheapUpMod(pe_old - sinceprefixend[pe_old]);				//(1) get prefix end
				if(globalposition[pe] > globalposition[pe_old]) break;			// (B) must not hop over circle boundary
				if(sinceprefixend[pe] != 0){pe_old = cheapUpMod(pe - 1); break;}// (C) prefix end must be valid pattern
				//cout << "\nPEplen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se;
				gaplen = cheapMod(ss - pe - 1);									//(2) get gap length
				if(gaplen > se1->gaplenmax) break;									// (D) gap must not be too long
				if(gaplen >= se1->gaplenmin)		 									// (E) if gap is too short, skip it, but extend anyway
				{

					pe3 = cheapUpMod(pe - pexcessnt);									//(3) get reading frame of prefix end
					
					//LOOP TO GET PREFIX START
					ps_old = cheapDownMod(pe3 + 1);
					//NN3N
					//7890
					
					while(true) //prefix start must be valid; breaks in this loop increment pe
					{
						//cout << "\nPSplen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se;
						if(sinceprefixstart[ps_old] >= circlelength){break;}			// (A) must not hop further than circle length
						ps = cheapUpMod(ps_old - sinceprefixstart[ps_old]);				//(1) get prefix start
						if(globalposition[ps] > globalposition[ps_old]){break;}			// (B) must not hop over circle boundary
						if(sinceprefixstart[ps] != 0){break;}							// (C) prefix start must be valid pattern
						plen = globalposition[pe] - globalposition[ps] + 1;				//(2) get prefix length
						totlen = plen + slen;
						if(plen > plenmax || (slen + plen) > se1->totlenmax){break;}			// (D) prefix must not be too long, else break and increment pe
						if(getBrokenTripletAA() == '-') break;							// (E) broken triplet must not be a stop codon
						if(plen - pexcessnt > sincestop[pe3]){break;}					// (F) prefix must not contain stop codons
						//xNN3NN
						//012345

												
						if(plen >= se1->plenmin)												// (E, F, H) prefix must be long enough
						{
							ps3 = cheapUpMod(ps - 1);										//(5) get mass ruler position for prefix start
							knownprefixmass = massruler[pe3] - massruler[ps3]; 									//(6) get known prefix mass
							peptideresiduemass = knownprefixmass + knownsuffixmass + getBrokenTripletMass();	//(7) get peptide mass, (method adds prefix and suffix and broken triplet)
							
							monoparentmassMH = peptideresiduemass + dnaAA1.monomassHOHH;
							if(monoparentmassMH <= se1->max_monoparentmassMH && monoparentmassMH >= se1->min_monoparentmassMH){
								
								//cout << "\nplen slen ps pe ss se: " << plen << " " << slen << " " << ps << " " << pe << " " << ss << " " << se << " " << getNTSequence();
									//changes BYR: changed call to include original AASequence						
								tuples->addTuple(new Tuple(monoparentmassMH, totlen/3, globalposition[ps], globalposition[pe], globalposition[ss], globalposition[se],getAASequenceFastChar(),0, activechromosome, NULL, chr_reverse, false, true, true));
								
							}
							
						}//end if plen min
				
						ps_old = cheapUpMod(ps - 3); //get next potential prefixstart
					}//end while ps
					
				}//end if gaplen min
				pe_old = cheapUpMod(pe - 1); //get next potential prefixend
			}//end while pe

		}//end if slen min
		ss_old = cheapUpMod(ss - 1); //get next potential suffixstart
	}//end while ss

}


void SlidingWindow::enumerateEntirePeptides()
{
	pexcessnt = 0;			// there is no broken triplet
	se = cheapMod(imodd); 	//suffix end is given by position of neighborhood
					
					//LOOP TO GET PREFIX START
					ps_old = cheapDownMod(se + 1);
					
					while(true) //prefix start must be valid; breaks in this loop increment se
					{
						if(sinceprefixstart[ps_old] >= circlelength){break;}			// (A) must not hop further than circle length
						ps = cheapUpMod(ps_old - sinceprefixstart[ps_old]);				//(1) get prefix start
						if(globalposition[ps] > globalposition[ps_old]){break;}			// (B) must not hop over circle boundary
						if(sinceprefixstart[ps] != 0){break;}							// (C) prefix start must be valid pattern
						
						ss = cheapDownMod(se + 1); slen = 0; pe = se;		// there is no suffix; 1-1=0 but length 1
						ps3 = cheapUpMod(ps - 1);	// prefix start to calculate parent mass
						plen = globalposition[pe] - globalposition[ps] + 1;				//(2) get getTotal length
						totlen = plen;
						
						if(totlen > se1->totlenmax){break;}									// (D) peptide must not be too long, else break and increment pe
						if(totlen > sincestop[se]){break;}								// (F) peptide must not contain stop codons
												
						if(totlen >= se1->plenmin)												// (E, F, H) prefix must be long enough
						{
							peptideresiduemass = massruler[se] - massruler[ps3]; 			//(6) se3 = se
							monoparentmassMH = peptideresiduemass + dnaAA1.monomassHOHH;		
							if(monoparentmassMH < se1->max_monoparentmassMH && monoparentmassMH > se1->min_monoparentmassMH){
									//changes BYR: changed call to add original sequence					
								tuples->addTuple(new Tuple(monoparentmassMH, totlen/3, globalposition[ps], globalposition[pe], globalposition[ss], globalposition[se], getAASequenceFastChar(),0, activechromosome, NULL, chr_reverse, false, true, true));
                                    //BX: Changing the Tuple constructor for optimization reasons leads to a memleak here, but as that part of the code is never called in the ETS process, it isn't that bad...
							}
						}
				
						ps_old = cheapUpMod(ps - 3); //get next potential prefixstart
					}//end while ps

}

}
