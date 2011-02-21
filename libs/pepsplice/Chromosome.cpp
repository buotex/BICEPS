#include "Chromosome.h"

namespace Pepsplice{
Chromosome::Chromosome(string f, Services *se0)
{
	description_line = "";
	se1 = se0;
	file = f;
	fullchromsize = 0;
	NTsequence = new string();
	//dnaAA1 = DnaAA();
}

Chromosome::~Chromosome()
{
	delete NTsequence;
}

int Chromosome::getFullChromSize()
{
	return fullchromsize;
}

char Chromosome::getNucleotide(long pos)
{
	return (*NTsequence)[pos];
}

char Chromosome::getReverseNucleotide(long pos)
{
	return se1->dnaAA1->nt_revnt((*NTsequence)[pos]);
}

void Chromosome::loadFullSequence()
{
	ifstream inFile;
	inFile.open(file.c_str());
	char newNucleotide;
	string test = "";
	
	//parse file
	getline(inFile, description_line);
	if(se1->outputlevel > 2) cout << "\nChromosome parsing starts. Description line in FASTA file: " << description_line;

	//get first nucleotide
	inFile >> newNucleotide;
	while (inFile.good() == true) //continue until end of file
	{
		NTsequence->push_back(newNucleotide);
		inFile >> newNucleotide; 
	}
	inFile.close();
	
	fullchromsize = NTsequence->size();
	if(se1->outputlevel > 2) cout << "\n" << getFullChromSize() << " nucleotides were parsed. Line breaks etc were not counted.";
}

void Chromosome::unloadFullSequence()
{
	//cout << "\nChromosome::unloadFullSequence(): NT sequence is being evicted" << flush;
	delete NTsequence;
	NTsequence = new string();
	(*NTsequence) = "";
	unloaded = true;
}

char* Chromosome::getNTSeqForTuple(long len, long ps, long pe, long ss, long se, bool rev)
{
	int size = getFullChromSize();
	if(unloaded == true){
		cout << "\nChromosome::getNTSeqForTuple: chromosome " << description_line << " was already unloaded" << flush;
		return NULL;
	}
	
	//test length
	int plen = pe - ps + 1;
	int slen = se - ss + 1;
	int tlen = plen + slen;
	if(tlen/3 != len){
		cout << "\nChromosome::getNTSeqForTuple: length wrong, tlen is " << tlen << flush;
		return NULL;
	}

	char *ntseq = new char[len*3];
	string seq = "";
		
	if(rev == false){
		//forward tuple
		//cout << "\n\nPS: " << ps << " PE: " << pe << " SS: " << ss << " SE: " << se << "\n";
		for(int i = ps; i <= pe; i++){
			seq += (*NTsequence)[i - 1];
			//cout << (*NTsequence)[i - 1];
		}
		//cout << "\t";
		for(int i = pe + 1; i < ss; i++){
			//cout << (*NTsequence)[i - 1];
		}
		//cout << "\t";
		for(int i = ss; i <= se; i++){
			seq += (*NTsequence)[i - 1];
			//cout << (*NTsequence)[i - 1];
		}
	}else{
		//chromosome_reverse tuple
		long psrev = size - ps;
		long perev = size - pe;
		long ssrev = size - ss;
		long serev = size - se;
		for(int i = psrev; i >= perev; i--){
			seq += se1->dnaAA1->nt_revnt((*NTsequence)[i]);
		}
		for(int i = ssrev; i >= serev; i--){
			seq += se1->dnaAA1->nt_revnt((*NTsequence)[i]);
		}
	}
		
	//convert string to dynamic character array
	if(seq.size() != len * 3) cout << "\nChromosome::getTupleSeq: sequence length is wrong! seq.size() " << seq.size() << " len: " << len << " ps: " << ps << " pe: " << pe << " ss: " << ss << " se: " << se;
	for(int i = 0; i < len*3; i++){
		ntseq[i] = seq[i];
	}
	return ntseq;
}
}

