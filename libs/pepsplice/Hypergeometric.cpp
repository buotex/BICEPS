#include "Hypergeometric.h"
namespace Pepsplice{
Hypergeometric::Hypergeometric()
{
	initLnFactorials();
}

Hypergeometric::~Hypergeometric()
{
	delete[] log10factorials;
	delete[] log10factorialsFloat;
}

double Hypergeometric::log10HypergeometricPDF(int N, int K, int n, int k) const { //Changes by BX
	
	/*
	 * pdf(N, K, n, k) = C(K, k)*C(N-K, n-k)/C(N-n)
	 * 
	 * C(K, k) = K choose k
	 * 
	 */
		
	double lnProbability = 0;
	lnProbability += (+ log10factorials[K  ] - log10factorials[k  ] - log10factorials[K-k]);
	lnProbability += (+ log10factorials[N-K] - log10factorials[n-k] - log10factorials[N-K-n+k]);
	lnProbability -= (+ log10factorials[N  ] - log10factorials[n  ] - log10factorials[N-n]);
	return lnProbability;
}

float Hypergeometric::log10HypergeometricPDFfloat(int N, int K, int n, int k) const {//Changes by BX
	
	/*
	 * pdf(N, K, n, k) = C(K, k)*C(N-K, n-k)/C(N-n)
	 * 
	 * C(K, k) = K choose k
	 * 
	 */
	
	float lnProbability = 0;
	lnProbability += (+ log10factorialsFloat[K  ] - log10factorialsFloat[k  ] - log10factorialsFloat[K-k]);
	lnProbability += (+ log10factorialsFloat[N-K] - log10factorialsFloat[n-k] - log10factorialsFloat[N-K-n+k]);
	lnProbability -= (+ log10factorialsFloat[N  ] - log10factorialsFloat[n  ] - log10factorialsFloat[N-n]);
	return lnProbability;
}

double Hypergeometric::log10HypergeometricPvalue(int N, int K, int n, int k) const{ //Changes by BX
	
	double pvalue = 0;
	//cout << "\n\n N K n k: " << N << " " << K << " " << n << " " << k;
	for(int i = n; i >= k; i--){
		pvalue += pow(10, log10HypergeometricPDF(N, K, n, i));
			//cout << "\n" << log10HypergeometricPDF(N, K, n, i);
			//cout << "\t" << pow(10, log10HypergeometricPDF(N, K, n, i));
			//cout << "\t" << pvalue;
	}
	return log10(pvalue);
}


void Hypergeometric::initLnFactorials(){ //Changes by BX
	int maxparentmass = 100000;
	double lni = 0;
	log10factorials = new double[maxparentmass];
	log10factorialsFloat = new float[maxparentmass];
    log10factorials[0] = 0.;
    log10factorialsFloat[0] = 0.0f;
	for(int i = 1; i < maxparentmass; i++){
		//if(i > 0) lni = lni + log10((double)i); //not needed
        lni += log10((double)i);
		log10factorials[i] = lni;
		log10factorialsFloat[i] = lni;
	}
}
    
}
