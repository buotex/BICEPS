#include "UnreachableBins.h"
namespace Pepsplice{
UnreachableBins::UnreachableBins(Services *se0)
{
	se1 = se0;
	initialized = false;
	lowmassbins_tryptic = new bool[1000];
	for(int i = 0; i < 1000; i++){lowmassbins_tryptic[i] = 1;}
	lowmassbins_nontryptic = new bool[1000];
	for(int i = 0; i < 1000; i++){lowmassbins_nontryptic[i] = 1;}
	highmassbins_tryptic = new bool[1000];
	for(int i = 0; i < 1000; i++){highmassbins_tryptic[i] = 1;}	
	highmassbins_nontryptic = new bool[1000];
	for(int i = 0; i < 1000; i++){highmassbins_nontryptic[i] = 1;}
	Nspectrum = new bool[100000];
	for(int i = 0; i < 100000; i++){Nspectrum[i] = 0;}
}

UnreachableBins::~UnreachableBins()
{
    delete[] lowmassbins_tryptic;
    delete[] lowmassbins_nontryptic;
    delete[] highmassbins_tryptic;
    delete[] highmassbins_nontryptic;
    delete[] Nspectrum;
}


/***************************************************************************
 * REACHABILITY OF BINS, correction of N and K
 * distinction between tryptic and non-tryptic
 * distinction between low and high mass region (high mass region is smeared 
 * due to parent mass tolerance
 * ************************************************************************/

bool UnreachableBins::isBinReachable(int bin, int trypticends){
	if(initialized == false) initialize();
	if(bin > 400) return true;
	if(trypticends < 2){
		return lowmassbins_nontryptic[bin];
	}else{
		return lowmassbins_tryptic[bin];
	}
}

int UnreachableBins::getTotal(int trypticends){
	if(initialized == false) initialize();	
	if(trypticends < 2){
		return total_nontryptic;
	}else{
		return total_tryptic;
	}
}

int UnreachableBins::estimateK(bool *binaryB, int parentmassMHdisc, int trypticends){
	if(initialized == false) initialize();
	
	int K = 0;

	//choose tryptic or non-tryptic array for lookup
	bool *lowmassbins;
	bool *highmassbins;
	if(trypticends < 2){
		lowmassbins = lowmassbins_nontryptic;
		highmassbins = highmassbins_nontryptic;
	}else{
		lowmassbins = lowmassbins_tryptic;
		highmassbins = highmassbins_tryptic;
	}
	
	for(int i = 0; i < parentmassMHdisc; i++){
		if(binaryB[i] > 0){
			bool isreachable = true;
			if(i < 300 && lowmassbins[i] == 1) isreachable = false;
			if(parentmassMHdisc + 1 - i < 300 && highmassbins[parentmassMHdisc - i] == 1) isreachable = false;
			if(isreachable == true) K++;
		}
	}	
	
	return K;	
}

void UnreachableBins::estimateNK(bool *binaryB, int parentmassMHdisc, int trypticends, int *N, int *K){
	if(initialized == false) initialize();
	
	(*K) = 0;
	(*N) = 0;

	//choose tryptic or non-tryptic array for lookup
	bool *lowmassbins;
	bool *highmassbins;
	if(trypticends < 2){
		lowmassbins = lowmassbins_nontryptic;
		highmassbins = highmassbins_nontryptic;
	}else{
		lowmassbins = lowmassbins_tryptic;
		highmassbins = highmassbins_tryptic;
	}

	//clean Nspectrum as far as needed
	for(int i = 0; i < parentmassMHdisc + 10; i++){Nspectrum[i] = 0;}
	
	//label low empty region
	for(int i = 0; i < parentmassMHdisc; i++){
		if(binaryB[i] == 1) break;
		Nspectrum[i] = 1;
	}
	
	//label high empty region
	for(int i = parentmassMHdisc - 1; i >= 0; i--){
		if(binaryB[i] == 1) break;
		Nspectrum[i] = 1;
	}
	
	//label combinatorially impossible bins
	for(int i = 0; i < 300; i++){
		if(lowmassbins[i] == 1) Nspectrum[i] = 1;
		if(highmassbins[i] == 1) Nspectrum[parentmassMHdisc - i] = 1;		
	}
	
	for(int i = 0; i < parentmassMHdisc; i++){
		//if(i % 50 == 0) cout << "  " << i << "\n";
		if(Nspectrum[i] == 0){
			if(binaryB[i] == 0){
				(*N)++;
				//cout << "'";
			}
			if(binaryB[i] == 1){
				(*K)++;
				//cout << ":";
			}
		}
		//if(Nspectrum[i] == 1 && binaryB[i] == 0) cout << " ";
		//if(Nspectrum[i] == 1 && binaryB[i] == 1) cout << ".";
	}
	//cout << parentmassMHdisc << " N: " << (*N) << " K: " << (*K);
}


void UnreachableBins::initialize(){
	//cout << "\n\nUnreachableBins is being initialized\n";
	
	int lt = fillBins2(true, lowmassbins_tryptic);
	int ht = fillBins2(true, highmassbins_tryptic);
	//int ht = smearBins(lowmassbins_tryptic, highmassbins_tryptic, 500, (int)se1->masstol_belowDa, (int)se1->masstol_aboveDa);
	
	int ln = fillBins2(false, lowmassbins_nontryptic);
	int hn = fillBins2(false, lowmassbins_nontryptic);
	//int hn = smearBins(lowmassbins_nontryptic, highmassbins_nontryptic, 500, (int)se1->masstol_belowDa, (int)se1->masstol_aboveDa);
	
	total_tryptic = lt + ht;
	total_nontryptic = ln + hn;
	//cout << "\n\nTotal unreachable bins tryptic: " << total_tryptic << " non-tryptic: " << total_nontryptic;
	
	initialized = true;
}

int UnreachableBins::fillBins(bool tryptic, bool* unreachablebins)
{
	for(int i = 0; i < 1000; i++){unreachablebins[i]=1;}
	
	int residuemasses[20] = {0, 57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186};
	for(int i = 0; i < 19; i++){
		//differentiate for tryptic peptides and semi- or non-tryptic peptides

		for(int j = 0; j < 19; j++){
			for(int k = 0; k < 19; k++){
				for(int l = 0; l < 19; l++){
					for(int m = 0; m < 19; m++){
						int rm = residuemasses[i] + residuemasses[j] + residuemasses[k] + residuemasses[l] + residuemasses[m];
						if(rm < 57) rm = 57;
						unreachablebins[rm+1] = 0;
						if(tryptic == 0 || residuemasses[i] == 128 || residuemasses[i] == 156){
							unreachablebins[rm+19] = 0;
						}
					}
				}
			}
		}

	}
	
	//count
	int count=0;
	//cout << "\n";
	for(int i = 1; i <= 350; i++){
		//cout << (int)unreachablebins[i];
		if(unreachablebins[i]==1) count++;
		//if(i%50 == 0) cout << "\t" << i << "\n";
	}
	//cout << "\nTotal unreachable bins in the lowest mass region if tryptic is " << tryptic << ": " << count << "\n";
	return count;
}

int UnreachableBins::fillBins2(bool tryptic, bool* unreachablebins)
{
	for(int i = 0; i < 1000; i++){unreachablebins[i]=1;}
	
	int residuemasses[18] = {57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186};
	
	//set start positions
	if(tryptic){
		for(int j = 0; j < 18; j++){
			unreachablebins[residuemasses[j] + 1] = 0;
		}
		unreachablebins[128 + 19] = 0;
		unreachablebins[156 + 19] = 0;
	}else{
		for(int j = 0; j < 18; j++){
			unreachablebins[residuemasses[j] + 1] = 0;
			unreachablebins[residuemasses[j] + 19] = 0;
		}			
	}
	
	for(int i = 0; i < 500; i++){
		if(unreachablebins[i]==0){
			for(int j = 0; j < 18; j++){
				unreachablebins[i+residuemasses[j]]=0;
			}
		}
	}
	
	//count
	int count=0;
	//cout << "\n";
	for(int i = 1; i <= 350; i++){
		//cout << (int)unreachablebins[i];
		if(unreachablebins[i]==1) count++;
		//if(i%50 == 0) cout << "\t" << i << "\n";
	}
	//cout << "\nTotal unreachable bins in the lowest mass region if tryptic is " << tryptic << ": " << count << "\n";
	return count;
}


int UnreachableBins::smearBins(bool *sourcebins, bool *targetbins, int length, int negtol, int postol){
	int count = 0;
	if(postol + negtol < 1) postol = -negtol + 1;
	cout << "\n negtol: " << negtol << " postol: " << postol;
	for(int i = 0; i < length; i++){
		if(sourcebins[i] == 0){
			for(int j = -negtol; j < postol; j++){
				targetbins[i+j]=0;
			}
		}
	}
	
	//count
	cout << "\n";
	for(int i = 1; i <= 350; i++){
		cout << (int)targetbins[i];
		if(targetbins[i]==1) count++;
		if(i%50 == 0) cout << "\t" << i << "\n";
	}
	cout << "\nTotal unreachable bins (smeared): " << count << "\n";
	return count;
}
}