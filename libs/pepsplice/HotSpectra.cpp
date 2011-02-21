#include "HotSpectra.h"
namespace Pepsplice{
HotSpectra::HotSpectra(Services *se0)
{
	se1 = se0;
	scanned_dtas = 0;
	d = new Distribution(-10, 10, 2, 2000);
}

HotSpectra::~HotSpectra()
{
    delete d;
}

void HotSpectra::parseHotSpectrumFiles(string filename)
{
	ifstream inFile;
	inFile.open(filename.c_str());
	string line;
	while(  getline( inFile, line ) )
	{
		parseHotSpectrumFile(line);		
	}	
	inFile.close();
	d->writeDistribution("distribution_hotspectra.txt");	
}

void HotSpectra::parseHotSpectrumFile(string filename){
	cout << "\n\nParsing hot spectrum file: " << filename << "\t";

	ifstream inFile;
	inFile.open(filename.c_str());
	string line;
	while(  getline( inFile, line ) )
	{
		scanned_dtas++;
		
		istringstream iss( line );
		string field;
		string dtaname = "";
		double qualityscore = 0;
		int i = 0;
		while( getline(iss, field, '\t') ){
			i++;
			if(i == 1) dtaname = field;
			if(i == 2) qualityscore = se1->string_to_double(field);
		}
		
		d->addElement(0, qualityscore);
		if(true) hotspecdtas.push_back(dtaname);
		//if(qualityscore >= 1) hotspecdtas.push_back(dtaname);
		
	}	
	sort(hotspecdtas.begin(), hotspecdtas.end());
	inFile.close();
	

	
	cout << " ... " << scanned_dtas << " dta file names scanned and " << hotspecdtas.size() << " above cutoff parsed.\n";
	for(int i = 0; i < 10 && i < hotspecdtas.size(); i++){
		cout << "\n" << hotspecdtas[i];
	}
}

bool HotSpectra::checkIfHot(string specname){
	bool wasfound = false;
	int hssize = hotspecdtas.size();
	int a = distance(hotspecdtas.begin(), lower_bound(hotspecdtas.begin(), hotspecdtas.end(), specname) );
	if(a < hssize){
		if(hotspecdtas[a] == specname) wasfound = true;
		//cout << "\n" << a << "\t" << specname << "\t" << hotspecdtas[a] << "\t" << wasfound;
	}
	return wasfound;
}
}