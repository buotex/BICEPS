#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
void gmm_bic(int numSegments, int numScores, const char* filename, 
             std::vector<double> & bicMu, std::vector<double> & bicSigma, std::vector<int> & labels);

