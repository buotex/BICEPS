#include "RegularRankModel.h"


void RegularRankModel::init_default_thresholds()
{
//	const int levels[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,35,55,80,999};
//	const int levels[]={0,1,2,3,4,5,10,15,20,50,100,9999};
	const int levels[]={0,3,7,12,19,30,55,9999};
	
	const int num_levels = sizeof(levels)/sizeof(int);

	level_thresholds.clear();
	int i;
	for (i=0; i<num_levels; i++)
		level_thresholds.push_back(levels[i]);

	num_peak_levels = num_levels;
}



void RegularRankModel::write_model_specifics(ostream &os) const
{
	os << "#REGULAR_RANK_MODEL" <<endl;
	os << level_thresholds.size();
	int i;

	for (i=0; i<level_thresholds.size(); i++)
		os << " " << level_thresholds[i];
	os << endl;
}


void RegularRankModel::read_model_specifics(istream &is)
{
	char buff[256];
	int i;

	is.getline(buff,256);

	if (strncmp("#REGULAR_RANK_MODEL",buff,18))
	{
		cout << "Error: expecting to find \"#REGULAR_RANK_MODEL\", not: " << buff << endl;
		exit(1);
	}

	
	is.getline(buff,1024);
	istringstream iss(buff);
	iss >> num_peak_levels;
	level_thresholds.resize(num_peak_levels);
	for (i=0; i<num_peak_levels; i++)
		iss >> level_thresholds[i];
}





void RegularRankModel::print_level_legend( ostream& os) const
{
	os << "Rank Levels: " << endl;
	os << "0 - no peak" << endl;
	int r;
	for (r=1; r<this->level_thresholds.size(); r++)
	{
		if (level_thresholds[r]-1 == level_thresholds[r-1])
		{
			os << r << " - rank " << level_thresholds[r] << endl;
		}
		else
			os << r << " - ranks " << level_thresholds[r-1]+1 << " to " <<
				level_thresholds[r] << endl;
	}
	os << endl;
}

