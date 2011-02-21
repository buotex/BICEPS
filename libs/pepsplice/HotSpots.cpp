#include "HotSpots.h"
namespace Pepsplice{
HotSpots::HotSpots(Services *se0)
{
	se1 = se0;
}

HotSpots::~HotSpots()
{
}

void HotSpots::parseHotSpots(string file)
{
	//cout << "\nparsing hot spots: " << file << flush;

	string line;
	ifstream inFile;
	inFile.open(file.c_str());

	//map< size_t,  vector< string> > Records;
	//size_t PrimaryKey( 0 );
	
	string chrname;
	bool rev;
	long ps;
	long pe;
	
	while(  getline( inFile, line ) )
	{
		hotspotstrings.push_back(line);
		//++PrimaryKey;
	}	
	
	cout << "\n" << hotspotstrings.size() << " hot spots parsed.";

	inFile.close();
	//getHotSpots(">CHR1v01212004", 1);

}


vector<int> HotSpots::getHotSpots(string chromname, bool rev1)
{
	string line;
	string field;
	vector<int> hotspotints;
	
	cout << "\n";
	for(int i = 0; i < hotspotstrings.size(); i++){
		istringstream iss( hotspotstrings[i] );
		//cout << "\n" << i << "---" << line << "---";
		int f = 0;
		bool hit = true;
		for(int j = 1; j <= 4; j++){
			getline( iss, field, '\t' );
			if(j == 1){
				if(field != chromname) hit = false;
				//cout << "CHR:" << field << "." << hit << "\t";
			}else if(j == 2){
				bool rev2 = false;
				if(field == "fwd"){rev2 = false;}else{rev2 = true;}
				//IGNORE DIRECTION ON GENOME
				//if(rev2 != rev1) hit = false;
			}else if(j == 3){
				if(hit == true){
					hotspotints.push_back((int)se1->string_to_double(field));
					//if(i < 3) cout << "\n" <<  field << " hotspotints.size(): " << hotspotints.size();
				}
				//cout << "PS:" << field << ".\t";
			}else if(j == 4){
				if(hit == true) hotspotints.push_back((int)se1->string_to_double(field));				
				//cout << "PE:" << field << ".\t";
			}
		}
	}
	
	sort(hotspotints.begin(), hotspotints.end());

	return hotspotints;
}
}
