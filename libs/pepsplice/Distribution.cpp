#include "Distribution.h"
namespace Pepsplice{
Distribution::Distribution(double min, double max, int wid, long len)
{
	minimum = min;
	maximum = max;
	length = len;
	if(length < 1) length = 1;
	width = wid;
	bins = new double*[width];
	
	//initialize with zeroes
	for(int i=0; i<width; i++){
		bins[i] = new double[length];
		for(int j=0; j<length; j++){bins[i][j]=0;}
	}
	granularity = (maximum - minimum) / length;
}

Distribution::~Distribution()
{
	
	for(int i = 0; i < width; i++){	
		delete[] bins[i];
	}
	delete[] bins;

}

int Distribution::addElement(int series, double element)
{
	int binpos = getBinPos(series, element);
	if((series < 0) || (series >= width)) cout << "\nDistribution:: series out of bounds";
	bins[series][binpos] = bins[series][binpos] + 1;
	return binpos;
}

int Distribution::addElement(int series, double element, double increase)
{	
	int binpos = getBinPos(series, element);
	if((series < 0) || (series >= width)) cout << "\nDistribution::addElement series out of bounds";
	bins[series][binpos] = bins[series][binpos] + increase;
	return binpos;
}

int Distribution::getBinPos(int series, double element)
{
	long binpos = (long)((element - minimum) / granularity);
	if(binpos >= length) binpos = length - 1;
	if(binpos < 0) binpos = 0;
	return binpos;
}

void Distribution::setBinValue(int series, long binpos, double count)
{
	if(binpos >= length) binpos = length - 1;
	if(binpos < 0) binpos = 0;
	if((series < 0) || (series >= width)) cout << "\nDistribution::setBin series out of bounds";
	bins[series][binpos] = count;
}

double Distribution::getBinValue(int series, long binpos)
{
	if(binpos >= length) binpos = length - 1;
	if(binpos < 0) binpos = 0;
	if((series < 0) || (series >= width)) cout << "\nDistribution::getBin series out of bounds";
	return bins[series][binpos];
}

double Distribution::getBinMiddle(long binpos)
{
	if(binpos >= length) binpos = length - 1;
	if(binpos < 0) binpos = 0;
	return minimum + (binpos + 0.5) * granularity;
}

double Distribution::sumSeries(int series)
{
	if((series < 0) || (series >= width)) cout << "\nDistribution::sumSeries series out of bounds";
	double sum = 0;
	for(int i = 0; i < length; i++){
		sum += bins[series][i];
	}
	return sum;
}

void Distribution::writeDistribution(string filename)
{
	ofstream outFile;
	outFile.open(filename.c_str());
	for(int j = 0; j < length; j++){
		outFile << "\n";
		outFile << j << "\t";
		outFile << getBinMiddle(j) << "\t";
		for(int i = 0; i < width; i++){	
			outFile << bins[i][j] << "\t";
		}
	}
	outFile.close();
}
    
}
