#include "Config.h"
#include "auxfun.h"
#include "bicepsdefinitions.h"

Config::Config() : resource_dir(biceps::bicepsconfigpath.append("/Models").c_str()), min_exclude_range(9999999), max_exclude_range(NEG_INF) {};


void Config::init_with_defaults()
{
	int i;

	mass_spec_type=ESI_MASS_SPEC; // default type
	config_file = "";
	model_name = "";

	ind_read_PTM_file = false;

	max_n_term_mod = 0;
	max_c_term_mod = 0;
	min_n_term_mod = 0;
	min_c_term_mod = 0; 

	tolerance = 0.5;
	pm_tolerance = 2.5;

	local_window_size=200;
	max_number_peaks_per_local_window=15;
	number_of_strong_peaks_per_local_window=10;
//	random_prob = 0.1;

	max_combo_mass = 0;

	max_edge_length = 2;

	max_charge_for_size=10;       // the size  of size_thresholds

	set_digest_type(0);

	need_to_estimate_pm = 1;

	use_spectrum_charge = 0;

	use_spectrum_mz = 0;

	filter_flag = 1;

	need_to_normalize = 1;

	itraq_mode = 0;

	terminal_score = 10;

	digest_score = 10;

	forbidden_pair_penalty = 25;
	strong_type1_idx =-1;
	strong_type2_idx =-1;

	
	min_ranges.clear();
	max_ranges.clear();
	min_exclude_range=999999;
	max_exclude_range=NEG_INF;

	init_standard_aas();

	session_aas = standard_aas;

	init_model_size_and_region_thresholds();

	int c;
	for (c=max_charge_for_size; c>=0; c--)
		init_regional_fragment_set_defaults(0,c);

	// these all operate on the original aas
	session_tables.init_for_standard_aas();

	// insert labels of original aas
	label2aa.clear();
	const vector<string>& aa2label = get_aa2label();
	for (i=0; i<aa2label.size(); i++)
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i],i));

	calc_aa_combo_masses();
	set_aa_variants();
	fill_allowed_double_edges();

	init_allowed_node_masses();

	char PTM_file[256];
	strcpy(PTM_file,get_resource_dir().c_str());
	strcat(PTM_file,"/");
	strcat(PTM_file,"PepNovo_PTMs.txt"); 

	read_PTM_file(PTM_file);

//	all_fragments.print();
}





void Config::init_standard_aas()
{
	int i;
	standard_aas.clear();
	for (i=Ala; i<=Val; i++)
		standard_aas.push_back(i);
}


void Config::init_model_size_and_region_thresholds()
{


	size_thresholds.clear();
	size_thresholds.resize(max_charge_for_size+1);

	// region thresholds
	region_thresholds.resize(2);
	region_thresholds[0]=0.225;
	region_thresholds[1]=0.775;
}



void Config::add_exclude_range(mass_t min_range, mass_t max_range)
{
	if (min_range>=max_range)
		return;

	this->min_ranges.push_back(min_range);
	this->max_ranges.push_back(max_range);
	if (this->min_exclude_range>min_range)
		min_exclude_range = min_range;

	if (this->max_exclude_range<max_range)
		max_exclude_range = max_range;
}


/********************************************************
Sets the diegest type.
Assumes that each digest type has a set of N or C terminal
preferred amino acids
*********************************************************/
void Config::set_digest_type(int type)
{
	digest_type = type;
	n_term_digest_aas.clear();
	c_term_digest_aas.clear();

	if (digest_type == NON_SPECIFIC_DIGEST)
		return;

	if (digest_type == TRYPSIN_DIGEST)
	{
		c_term_digest_aas.push_back(Lys);
		c_term_digest_aas.push_back(Arg);
		return;
	}

	cout << "Error: digest type not supported: " << type << endl;
	exit(1);
}

/*********************************************************
Inits the regional fragment sets (resizes according to the
charges, size_idxs, and region_idxs). If set_type = 0, 
each region is given all fragments that are relevant.
If set type = 0, uses all
**********************************************************/
void Config::init_regional_fragment_set_defaults(int set_type, int charge)
{
	const int num_regions = region_thresholds.size()+1;

	if (max_charge_for_size <= charge)
	{
		max_charge_for_size = charge;
		regional_fragment_sets.resize(max_charge_for_size+1);
	}

	regional_fragment_sets[charge].resize(size_thresholds[charge].size());
	int j;
	for (j=0; j<size_thresholds[charge].size(); j++)
	{
		int k;
		regional_fragment_sets[charge][j].resize(num_regions);
		for (k=0; k<num_regions; k++)
		{
			if (set_type == 0)
			{
				regional_fragment_sets[charge][j][k].init_with_all_types(charge,this);
			}
			else
			{
				cout << "Error: bad option for init!"<< endl;
				exit(1);
			}
		}
	}
}



// selects the fragments to be used, uses a cutoff that is X times random prob)
void Config::select_fragments_in_sets(score_t min_prob_coef, int max_num_frags)
{
	int c;
	for (c=1; c<regional_fragment_sets.size(); c++)
	{
		int s;
		for (s=0; s<regional_fragment_sets[c].size(); s++)
		{
			int r;
			
			for (r=0; r<regional_fragment_sets[c][s].size(); r++)
			{
				float min_prob = min_prob_coef * get_regional_random_probability(c,s,r);
				regional_fragment_sets[c][s][r].select_fragments_with_minimum_prob(min_prob,
					max_num_frags);
			}
		}
	}
}





// For each regional fragments selects all fragments that have a minimal probability
// to be strong.
// also look for the two strongest fragment types and store them in the
// strong_type1_idx and strong_type2_idx variables
void Config::select_strong_fragments(int charge,
						score_t min_prob, int max_num_strong, bool verbose)
{
	const vector<FragmentType>& all_fragment_types = get_all_fragments();
	
	int s;
	for (s=0; s<regional_fragment_sets[charge].size(); s++)
	{
		int r;
		for (r=0; r<regional_fragment_sets[charge][s].size(); r++)
		{
			regional_fragment_sets[charge][s][r].select_strong_fragments(this,min_prob,max_num_strong);
			const vector<int>& strong_type_idxs = regional_fragment_sets[charge][s][r].get_strong_frag_type_idxs();

			// see if these frags should be put in the strong two fragments
			int i;
			for (i=0; i<strong_type_idxs.size(); i++)
			{
				int frag_idx = strong_type_idxs[i];

				if (strong_type1_idx == -1)
				{
					strong_type1_idx = frag_idx;
					continue;
				}

				if (strong_type2_idx == -1 &&
					frag_idx != strong_type1_idx)
				{
					strong_type2_idx = frag_idx;
					continue;
				}

				// make sure that type2 has the lower probability
				if (all_fragment_types[strong_type1_idx].prob<
					all_fragment_types[strong_type2_idx].prob)
				{
					int tmp;
					tmp=strong_type1_idx;
					strong_type1_idx=strong_type1_idx;
					strong_type1_idx=tmp;
				}
				if (all_fragment_types[frag_idx].prob>
					all_fragment_types[strong_type2_idx].prob)
				{
					strong_type2_idx = frag_idx;
				}
			}
		}
	}	
}

// returns the region idx for a (breakge) mass
int	 Config::calc_region_idx(mass_t break_mass, mass_t pm_with_19, int charge,
					 mass_t min_peak_mass, mass_t max_peak_mass) const
{
//	return 0;

	const int n_threshes = region_thresholds.size();
	int i;
	vector<mass_t> threshes;
	threshes.resize(region_thresholds.size());
	for (i=0; i<region_thresholds.size(); i++)
		threshes[i] = pm_with_19 * region_thresholds[i];

	if (charge<=2)
	{
		mass_t alt_first_thresh = (min_peak_mass>0) ? min_peak_mass + 150 : threshes[0];
		mass_t alt_last_thresh =  (max_peak_mass>0) ? max_peak_mass - 150 : threshes[n_threshes-1];
		
		if (break_mass<threshes[0] || break_mass<alt_first_thresh)
			return 0;

		if (break_mass>threshes[n_threshes-1] || break_mass>alt_last_thresh)
			return n_threshes;


		for (i=1; i<n_threshes; i++)
			if (break_mass<threshes[i])
				break;
		return i;
	}
	else
	{
		if (break_mass<threshes[0])
			return 0;

		if (break_mass>threshes[n_threshes-1])
			return n_threshes;


		for (i=1; i<n_threshes; i++)
			if (break_mass<threshes[i])
				break;
		return i;
	}
}


void Config::set_size_thresholds_according_to_set_of_masses(int charge,
											vector<mass_t>& spectra_masses)
{
	const int num_spectra = spectra_masses.size();

	if (size_thresholds.size()<=charge)
		size_thresholds.resize(charge+1);

	size_thresholds[charge].clear();

	if (num_spectra<1500)
		return;

	sort(spectra_masses.begin(),spectra_masses.end());

	if (num_spectra<3000)
	{
		int idx = num_spectra/2;
		size_thresholds[charge].push_back(spectra_masses[idx]);
		return;
	}

	if (num_spectra<10000)
	{
		// use 3 sizes
		
		int idx1 = num_spectra/3;
		int idx2 = idx1+idx1;

		size_thresholds[charge].push_back(spectra_masses[idx1]);
		size_thresholds[charge].push_back(spectra_masses[idx2]);
	}

	// for large datasets of spectra, use the following rule
	const mass_t min_pm_diff = 200.0;
	const int    min_num_spectra_per_model = 25000;
	const int	 num_models = spectra_masses.size()/min_num_spectra_per_model;

	if (num_models <=1)
	{
		size_thresholds[charge].clear();
		size_thresholds[charge].push_back(POS_INF);
		return;
	}

	int i=0;
	int prev_i=0;
	int next_idx = min_num_spectra_per_model;
	mass_t last_mass =0;

	size_thresholds[charge].clear();

	while (i<spectra_masses.size())
	{
		if (spectra_masses[i]-last_mass>min_pm_diff && 
			i>= next_idx)
		{
			size_thresholds[charge].push_back(spectra_masses[i]);
			cout << "charge " << charge << "\t size " << size_thresholds[charge].size() << " \tmasss " << 
				setprecision(1) << fixed << size_thresholds[charge][size_thresholds[charge].size()-1] << 
				"\t  (spectra " << i-prev_i << ")" << endl;

			last_mass = spectra_masses[i];
			prev_i = i;

			int md_count = size_thresholds[charge].size();
			if (md_count == num_models -1 || spectra_masses.size() - i < 1.5*min_num_spectra_per_model)
				break;
			
			next_idx = i + (spectra_masses.size() - i)/(num_models-md_count);
	
		}
		i++;
	}

	size_thresholds[charge].push_back(POS_INF);
	
	cout << "charge " << charge << "\t size " << size_thresholds[charge].size() << 
		" \tmasss " << setprecision(1) << fixed << size_thresholds[charge][size_thresholds[charge].size()-1] <<
		"\t  (spectra " << spectra_masses.size()-prev_i << ")" << endl;
	
}

void Config::print_size_thresholds() const
{
	int charge;
	for (charge=0; charge<size_thresholds.size(); charge++)
	{
		cout << charge << "\t" << size_thresholds[charge].size();
		int t;
		for (t=0; t<size_thresholds[charge].size(); t++)
			cout << "\t" << size_thresholds[charge][t];
		cout << endl;
	}
}


int Config::determine_charge(Peptide& pep, mass_t m_over_z) const
{
	int c;
	pep.calc_mass(this);
	mass_t pep_mass_with_19 = pep.get_mass()+MASS_OHHH;

	for (c=1; c<10; c++)
	{
		mass_t calc_mass_with_19 = m_over_z * c - (1.0023 * (c-1));
		if (fabs(calc_mass_with_19- pep_mass_with_19)<10)
			return c;
	}
	return -1;
}


int	 Config::get_max_session_aa_idx() const
{
	const vector<int>& aas = get_session_aas();
	int max_aa=-1;
	int i;
	for (i=0; i<aas.size(); i++)
		if (aas[i]>max_aa)
			max_aa=aas[i];

	return max_aa;
}



/***********************************************************
// returns the idx of an aa from its label
// -1 if label is not found
***********************************************************/
int Config::get_aa_from_label(const string& label) const
{
	STRING2INT_MAP::const_iterator iter = label2aa.find(label);

	if (iter == label2aa.end())
		return -1;

	return (*iter).second;
}

int Config::get_aa_with_position_from_label(const string& label, int position) const
{
	STRING2INT_MAP::const_iterator iter;
	for (iter = label2aa.begin();iter != label2aa.end(); iter++)
	{
		int aa_idx = (*iter).second;
		if (session_tables.get_aa_position(aa_idx) != position ||
			label[0] != (*iter).first[0])
			continue;

		if (label == session_tables.get_aa2label(aa_idx))
			return aa_idx;
	}


	return -1;
}


void Config::print_fragments(ostream &os) const
{
	int i;
	os << "#ALL FRAGMENTS " << all_fragments.fragments.size() << endl;
	for (i=0; i<all_fragments.fragments.size(); i++)
		all_fragments.fragments[i].write_fragment(os);
}


void Config::print_regional_fragment_sets(ostream &os) const
{
	os << "#FRAGMENT SETS" << endl;

	int c;
	for (c=regional_fragment_sets.size()-1; c>=0; c--)
	{
		int s;
		for (s=regional_fragment_sets[c].size()-1; s>=0; s--)
		{
			int r;
			for (r=regional_fragment_sets[c][s].size()-1; r>=0; r--)
			{
				int f;
				const vector<int>& frag_idxs = regional_fragment_sets[c][s][r].get_frag_type_idxs();

				if (frag_idxs.size()>0)
				{
					os << c << " " << s << " " << r << " " << frag_idxs.size() << " " <<
						setprecision(5) << regional_fragment_sets[c][s][r].get_rand_prob() << endl;
					int i;

					// write fragments
					for (f=0; f<frag_idxs.size(); f++)
					{
						os << left << setw(9) << all_fragments.get_fragment(frag_idxs[f]).label 
							<< " " << setprecision(4) << regional_fragment_sets[c][s][r].get_frag_prob(f)<< "  " << endl;
						//	setprecision(5) << all_fragments.get_fragment(frag_idxs[f]).offset << endl;
					}

					// write strong
					const vector<int>& strong = regional_fragment_sets[c][s][r].get_strong_frag_type_idxs();
					os << "strong " << strong.size();
					for (i=0; i<strong.size(); i++)
						os << " " << get_fragment(strong[i]).label;
					os << endl;

					// write combo
					const vector<FragmentCombo>& combos = regional_fragment_sets[c][s][r].get_frag_type_combos();
					os << "combos " << combos.size() << endl;
					for (i=0; i<combos.size(); i++)
						combos[i].print_combo(this,os);

					os << endl;
				}
			}
		}
	}
}





/************************************************************
Reads the fragments from the stream
*************************************************************/
void Config::read_fragments(istream& is)
{
	int i,num_frags=-1;
	char buff[128];
	is.getline(buff,128);

	if (sscanf(buff,"#ALL FRAGMENTS %d",&num_frags) != 1)
	{
		cout << "Error: bad line in fragments file: " << buff << endl;

		if (strlen(buff)<4)
		{
			cout << "This error can possibly be due to Wnix/Windows issues.";
			cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
		}
		exit(1);
	}

	all_fragments.clear_set();
	all_strong_fragment_type_idxs.clear();
	
	for (i=0; i<num_frags; i++)
	{
		FragmentType ft;
		ft.read_fragment(is);
		all_fragments.add_fragment_type(ft);
	}
}

/************************************************************
reads the fragment sets from a stream
*************************************************************/
void Config::read_regional_fragment_sets(istream& is)
{
	int i,num_frags=-1;
	char buff[128];
	is.getline(buff,128);

	
	if (strncmp(buff,"#FRAGMENT SETS",12))
	{
		cout << "Error: bad line in fragments file: " << buff << endl;
		if (strlen(buff)<4)
		{
			cout << "This error can possibly be due to Wnix/Windows issues.";
			cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
		}

		exit(1);
	}

	regional_fragment_sets.clear();

	while ( ! is.eof())
	{
		int charge,size_idx,region_idx,num_fragments;
		is.getline(buff,128);
		if (is.gcount()<2)
			continue;

		// some other parameter set is in this file, put back line and return
		if (buff[0]=='#')
		{
			int i;
			for (i=strlen(buff)-1; i>=0; i--)
				is.putback(buff[i]);
			return;
		}

		float rand_prob=-1;
		if (sscanf(buff,"%d %d %d %d %f",&charge,&size_idx,&region_idx,&num_fragments,
			&rand_prob) != 5)
		{
			cout << "Error: bad line in fragments file: " << buff << endl;
			if (strlen(buff)<4)
			{
				cout << "This error can possibly be due to Wnix/Windows issues.";
				cout << "Try running dos2unix on all \".txt\" files in the Models directory (run \"dos2unix Models/*.*\" and \"dos2unix Models/*/*.*\")." << endl;
			}

			exit(1);
		}
		
		// make sure fragment_sets vectors are large enough
		if (regional_fragment_sets.size()<charge+1)
			regional_fragment_sets.resize(charge+1);
		if (regional_fragment_sets[charge].size()<size_idx+1)
			regional_fragment_sets[charge].resize(size_idx+1);
		if (regional_fragment_sets[charge][size_idx].size()<region_idx+1)
			regional_fragment_sets[charge][size_idx].resize(region_idx+1);
		
		RegionalFragments& rf = regional_fragment_sets[charge][size_idx][region_idx];

		rf.set_rand_prob(rand_prob);

		rf.frag_type_idxs.clear();

		// read fragments
		for (i=0; i<num_fragments; i++)
		{
			string  label;
			score_t prob = 0;
			is.getline(buff,128);
			istringstream iss(buff);

			iss >> label >> prob;
			
			const int frag_idx = this->get_frag_idx_from_label(label);
			if (frag_idx<0)
			{
				cout << "Error: unrecognized fragment: " << label << endl;
				exit(1);
			}

			all_fragments.fragments[frag_idx].prob += prob;

			rf.frag_type_idxs.push_back(frag_idx);
			rf.frag_probs.push_back(prob);
		}
		
		// read strong fragments
		int i;
		is.getline(buff,128);
		if (strncmp(buff,"strong",6))
		{
			cout << "Error: expected to see list of strong frags! : " << buff << endl;
			exit(1);
		}
		int num_strong=0;
		istringstream iss(buff+7);
		iss >> num_strong;
		rf.strong_frag_type_idxs.clear();
		for (i=0; i<num_strong; i++)
		{
			string label;
			iss >> label;
			int f_idx=get_frag_idx_from_label(label);
			if (f_idx>=0)
			{
				rf.strong_frag_type_idxs.push_back(f_idx);
				int j;
				for (j=0; j<all_strong_fragment_type_idxs.size(); j++)
					if (all_strong_fragment_type_idxs[j] == f_idx)
						break;
				if (j==all_strong_fragment_type_idxs.size())
					all_strong_fragment_type_idxs.push_back(f_idx);
			}
			else
				break;	
		}

		// read combos
		is.getline(buff,128);
		if (strncmp(buff,"combos",6))
		{
			cout << "Error: expected to see list of combos! : " << buff << endl;
			exit(1);
		}
		int num_combos;
		istringstream iss2(buff+7);
		iss2 >> num_combos;
		rf.frag_type_combos.clear();
		for (i=0; i<num_combos; i++)
		{
			FragmentCombo fc;
			fc.read_combo(this,is);
			int j;
			int  strong_frag_pos=-1;
			for (j=0; j<fc.frag_inten_idxs.size(); j++)
			{
				int k;
				for (k=0; k<rf.strong_frag_type_idxs.size(); k++)
				{
					if (fc.frag_inten_idxs[j] == rf.strong_frag_type_idxs[k])
					{
						strong_frag_pos = j;
						break;
					}
				}
			}

			if (strong_frag_pos<0)
			{
				cout << "Error: combo must have at least one strong fragment type with intensity!" << endl;
				exit(1);
			}

			// swap to make the strong frag first
			if (strong_frag_pos>0)
			{
				int tmp;
				tmp=fc.frag_inten_idxs[strong_frag_pos];
				fc.frag_inten_idxs[strong_frag_pos] = fc.frag_inten_idxs[0];
				fc.frag_inten_idxs[0] = tmp;
			}
			rf.frag_type_combos.push_back(fc);
		}
	}

	// clone regional fragment sets if needed
	if (regional_fragment_sets[2].size()>0)
	{
		if (regional_fragment_sets[1].size() == 0)
			clone_regional_fragment_sets(2,1);

		if (regional_fragment_sets.size()<=4)
			regional_fragment_sets.resize(5);

		if (regional_fragment_sets[3].size() == 0)
			clone_regional_fragment_sets(2,3);

		if (regional_fragment_sets[4].size() == 0)
			clone_regional_fragment_sets(3,4);
	}



}


void Config::clone_regional_fragment_sets(int source_charge, int target_charge)
{
	if (regional_fragment_sets.size()<=target_charge)
		regional_fragment_sets.resize(target_charge+1);

	regional_fragment_sets[target_charge] = regional_fragment_sets[source_charge];
}



bool read_mass_type_line(const char* prefix, char *line, mass_t& val)
{
	int len = strlen(prefix);
	if (strncmp(prefix,line,len))
		return false;
	istringstream is(line+len);
	is >> val;
	return true;
}

bool read_score_type_line(const char* prefix, char *line, score_t& val)
{
	int len = strlen(prefix);
	if (strncmp(prefix,line,len))
		return false;
	istringstream is(line+len);
	is >> val;
	return true;
}

/**********************************************************************
// parses a line that is assumed to be from a config file
// all parameters are assumed to start with
// #CONF <PARAMETER_NAME> VALUES
***********************************************************************/
void Config::parse_config_parameter(char *current_line)
{
	char buff[256];
	int number;

	// fragments file
	if ( sscanf(current_line,"#CONF FRAGMENTS_FILE %s",buff) == 1)
	{
		string path = resource_dir + "/" + string(buff);

		fstream fs(path.c_str(),ios::in);
		if (! fs.good())
		{
			//cout << "Error opening fragment sets: " << path << endl;
			throw string("Error opening fragment sets: " + path);
            //exit(1);
		}
		this->read_fragments(fs);
		fragments_file=buff;
		return;
	}

	// regional fragment sets
	if ( sscanf(current_line,"#CONF REGIONAL_FRAGMENT_SETS_FILE %s",buff) == 1)
	{
		string path = resource_dir + "/" + string(buff);

		fstream fs(path.c_str(),ios::in);
		if (! fs.good())
		{
			//cout << "Error opening fragment sets: " << buff << endl;
			throw string("Error opening fragment sets: ") + string(buff); 
		}
		this->read_regional_fragment_sets(fs);
		regional_fragment_sets_file=buff;
		set_all_regional_fragment_relationships();
		return;
	}


	// mass spec type
	if ( sscanf(current_line,"#CONF MASS_SPEC_TYPE %d",&number) == 1)
	{
		mass_spec_type = number;
		return;
	}

	// digest
	if ( sscanf(current_line,"#CONF DIGEST_TYPE %d",&number) == 1)
	{
		set_digest_type(number);
		return;
	}


	// resource dir
	if ( sscanf(current_line,"#CONF RESOURCE_DIR %s",buff) == 1)
	{
		resource_dir = string(buff);
		return;
	}

	if ( sscanf(current_line,"#CONF MAX_NUMBER_OF_PEAKS_PER_LOCAL_WINDOW %d",&number) == 1)
	{
		max_number_peaks_per_local_window = number;
		return;
	}

	if ( sscanf(current_line,"#CONF NUMBER_OF_STRONG_PEAKS_PER_LOCAL_WINDOW %d",&number) == 1)
	{
		number_of_strong_peaks_per_local_window = number;
		return;
	}

	if (read_mass_type_line("#CONF LOCAL_WINDOW_SIZE",current_line,local_window_size))
		return;

	if (read_mass_type_line("#CONF TOLERANCE",current_line,tolerance))
		return;

	if (read_mass_type_line("#CONF PM_TOLERANCE",current_line,pm_tolerance))
		return;

	if ( sscanf(current_line,"#CONF MAX_EDGE_LENGTH %d",&number) == 1)
	{
		max_edge_length = number;
		return;
	}

	
	if (read_score_type_line("#CONF TERMINAL_SCORE",current_line,terminal_score))
		return;

	if (read_score_type_line("#CONF DIGEST_SCORE",current_line,digest_score))
		return;

	if (read_score_type_line("#CONF FORBIDDEN_PAIR_PENALTY",current_line,forbidden_pair_penalty))
		return;

	if (sscanf(current_line,"#CONF NEED_TO_CORRECT_PM %d",&number) == 1)
	{
		need_to_estimate_pm = number;
		return;
	}


	if (!strncmp("#CONF SIZE_THRESHOLDS",current_line,21))
	{
		istringstream is(current_line+22);
		int i,num_sizes;

		is >> num_sizes;
		size_thresholds.resize(num_sizes);
		for (i=0; i<num_sizes; i++)
		{
			int j,size;
			is >> size;
			size_thresholds[i].resize(size);
			for (j=0; j<size; j++)
				is >> size_thresholds[i][j];
		}
		return;
	}

	if (!strncmp("#CONF REGION_THRESHOLDS",current_line,22))
	{
		istringstream is(current_line+23);
		int i,num_sizes;

		is >> num_sizes;
		region_thresholds.resize(num_sizes); 
		for (i=0; i<num_sizes; i++)
			is >> region_thresholds[i];

		return;
	}

	if ( !strncmp(current_line,"#CONF MASS_EXCLUDE_RANGE",23))
	{
		istringstream is(current_line+24);
		mass_t min_inp_range=-1.0;
		mass_t max_inp_range=-2.0;

		is >> min_inp_range >> max_inp_range;
		if (max_inp_range<min_inp_range)
		{
			cout << "Error: paramter \"#CONF MASS_EXCLUDE_RANGE\" should be followed a pair of numbers: minimal mass range and maximal mass range to exclude!" << endl;
			cout << "BAD Line: " << current_line << endl;
		}

		add_exclude_range(min_inp_range,max_inp_range);
	
		return;
	}

	if (current_line[0] == '#' || strlen(current_line)<3)
		return;

	cout << "Warning: bad line in config file: " << current_line << endl;
}




void Config::read_config(const char* file_name)
{
	config_file = file_name;
	
	init_with_defaults();
	
	string path = resource_dir + "/" + string(file_name);

	fstream fs(path.c_str(),ios::in);

	while (! fs.eof() && fs.good())
	{
		char buff[1024];

		fs.getline(buff,1024);
		if (fs.gcount()<4)
			continue;

		parse_config_parameter(buff);
	}
}



void Config::write_config()
{

	if (fragments_file.length()>0)
	{
		fragments_file = model_name + "_fragments.txt";
		string path = resource_dir + "/" + fragments_file;
		ofstream fs(path.c_str(),ios::out);
		print_fragments(fs);
	}

	if (regional_fragment_sets_file.length()>0)
	{
		regional_fragment_sets_file = model_name + "_fragment_sets.txt";
		string path = resource_dir + "/" + regional_fragment_sets_file;
		ofstream rfs(path.c_str(),ios::out);
		print_regional_fragment_sets(rfs);
	}

	ofstream os(config_file.c_str(), ios::out);
	print_config_parameters(os);
}



/****************************************************************
	calclates the masses of the different aa_combos
	combos are sorted aa lists.
	lazy version, hardcoded, more elegant to do recursive way 
	(without real recursion).

  And also finds the maximum combo mass

  Fills in the aa_edge_combos vector which holds all lengths together
*****************************************************************/
void Config::calc_aa_combo_masses()
{
	const vector<int>& session_aas = get_session_aas();
	const vector<mass_t>& aa2mass = get_aa2mass();
	const int num_aas = session_aas.size();
	const int last_aa = session_aas[num_aas-1];


	vector< vector<AA_combo> > aa_combo_by_length; // first dim is edge length
	aa_combo_by_length.resize(max_edge_length+1);

	int i;

	if (max_edge_length > MAX_EDGE_SIZE)
	{
		cout << "Error: code doesn't support edges larger than "<< MAX_EDGE_SIZE << endl;
		exit(0);
	}

	if (max_edge_length>5)
	{
		cout << "Error: code doesn't support such large edges!"<< endl;
		exit(0);
	}

	aa_combo_by_length[1].clear();
	for (i=0; i<num_aas; i++)
	{
		const int aa1 = session_aas[i];
		if (aa1 == Ile)
			continue;

		AA_combo ac;
		ac.amino_acids[0]=aa1;
		ac.num_aa=1;
		ac.total_mass=aa2mass[aa1];
		
		aa_combo_by_length[1].push_back(ac);
	}
	sort(aa_combo_by_length[1].begin(),aa_combo_by_length[1].end());
	

	if (max_edge_length>=2)
	{
		aa_combo_by_length[2].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;
		
				AA_combo ac;
				ac.amino_acids[0]=aa1;
				ac.amino_acids[1]=aa2;
				ac.num_aa=2;
				ac.total_mass+=mass1 + aa2mass[aa2];
				
				aa_combo_by_length[2].push_back(ac);
			}
		}
		sort(aa_combo_by_length[2].begin(),aa_combo_by_length[2].end());
	}

	if (max_edge_length>=3)
	{
		aa_combo_by_length[3].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					AA_combo ac;
					ac.amino_acids[0]=aa1;
					ac.amino_acids[1]=aa2;
					ac.amino_acids[2]=aa3;
					ac.num_aa=3;
					ac.total_mass+=mass2 + aa2mass[aa3];
					
					aa_combo_by_length[3].push_back(ac);
				}
			}
		}
		sort(aa_combo_by_length[3].begin(),aa_combo_by_length[3].end());
	}


	if (max_edge_length>=4)
	{
		aa_combo_by_length[4].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					const mass_t mass3 = mass2 + aa2mass[aa3];
					int l;

					for (l=k; l<num_aas; l++)
					{
						const int aa4 = session_aas[l];
						if (aa4 == Ile)
							continue;

						AA_combo ac;
						ac.amino_acids[0]=aa1;
						ac.amino_acids[1]=aa2;
						ac.amino_acids[2]=aa3;
						ac.amino_acids[3]=aa4;
						ac.num_aa=4;
						ac.total_mass+=mass3 + aa2mass[aa4];
						
						aa_combo_by_length[4].push_back(ac);
					}
				}
			}
		}
		sort(aa_combo_by_length[4].begin(),aa_combo_by_length[4].end());
	}

	if (max_edge_length>=5)
	{
		aa_combo_by_length[5].clear();
		for (i=0; i<num_aas; i++)
		{
			int j;
			const int aa1 = session_aas[i];
			const mass_t mass1 = aa2mass[aa1];
			if (aa1 == Ile)
				continue;

			for (j=i; j<num_aas; j++)
			{
				const int aa2 = session_aas[j];
				if (aa2 == Ile)
					continue;

				const mass_t mass2 = mass1 + aa2mass[aa2];
				int k;

				for (k=j; k<num_aas; k++)
				{
					const int aa3 = session_aas[k];
					if (aa3 == Ile)
						continue;

					const mass_t mass3 = mass2 + aa2mass[aa3];
					int l;

					for (l=k; l<num_aas; l++)
					{
						const int aa4 = session_aas[l];
						if (aa4 == Ile)
							continue;

						const mass_t mass4 = mass3 + aa2mass[aa4];
						int p;

						for (p=l; p<num_aas; p++)
						{
							const int aa5 = session_aas[p];
							if (aa5 == Ile)
								continue;

							AA_combo ac;
							ac.amino_acids[0]=aa1;
							ac.amino_acids[1]=aa2;
							ac.amino_acids[2]=aa3;
							ac.amino_acids[3]=aa4;
							ac.amino_acids[4]=aa5;
							ac.num_aa=5;
							ac.total_mass+=mass4 + aa2mass[aa5];
						
							aa_combo_by_length[5].push_back(ac);
						}
					}
				}
			}
		}
		sort(aa_combo_by_length[5].begin(),aa_combo_by_length[5].end());
	}
	
	max_combo_mass = aa_combo_by_length[max_edge_length][aa_combo_by_length[max_edge_length].size()-1].total_mass + 1.0;

	// fill in aa_edge_combos
	aa_edge_combos.clear();
	aa_edge_combos.reserve((int)(1.5 * aa_combo_by_length[max_edge_length].size()));
	for (i=1; i<aa_combo_by_length.size(); i++)
	{
		int j;
		for (j=0; j<aa_combo_by_length[i].size(); j++)
			aa_edge_combos.push_back(aa_combo_by_length[i][j]);	
	}

	sort(aa_edge_combos.begin(),aa_edge_combos.end());

	// make edge variant vector and fill in the variant pointers
	// and make combo_idxs_by_length vectors
	combo_idxs_by_length.resize(max_edge_length+1);
	for (i=0; i<=max_edge_length; i++)
		combo_idxs_by_length[i].clear();

	variant_vector.clear();
	variant_vector.reserve(aa_edge_combos.size()*5);

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		AA_combo& combo = aa_edge_combos[i];
		combo_idxs_by_length[combo.num_aa].push_back(i);

		if (combo.num_aa == 1)
		{
			combo.num_variants = 1;
			combo.variant_start_idx = variant_vector.size();

			variant_vector.push_back(1);
			variant_vector.push_back(combo.amino_acids[0]);
			continue;
		}


		if (combo.num_aa == 2)
		{
			if (combo.amino_acids[0]==combo.amino_acids[1])
			{
				combo.num_variants = 1;
				combo.variant_start_idx = variant_vector.size();

				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[0]);
				variant_vector.push_back(combo.amino_acids[1]);
				continue;
			}
			else
			{
				combo.num_variants = 2;
				combo.variant_start_idx = variant_vector.size();

				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[0]);
				variant_vector.push_back(combo.amino_acids[1]);
				variant_vector.push_back(2);
				variant_vector.push_back(combo.amino_acids[1]);
				variant_vector.push_back(combo.amino_acids[0]);
				continue;
			}
		}

		// generate variant using permutation generating function
		vector<int> org_perm;
		vector< vector<int> > all_perms;

		int i;
		org_perm.clear();
		for (i=0; i<combo.num_aa; i++)
			org_perm.push_back(combo.amino_acids[i]);

		generate_all_permutations(org_perm,all_perms);
		combo.num_variants = all_perms.size();
		combo.variant_start_idx = variant_vector.size();

		for (i=0; i<all_perms.size(); i++)
		{
			int j;
			variant_vector.push_back(combo.num_aa);
			for (j=0; j<combo.num_aa; j++)
				variant_vector.push_back(all_perms[i][j]);
		}
		continue;
	}

	// fill the combo_start_idxs
	int last_combo_idx = aa_edge_combos.size()-1;
	mass_t largest_combo_mass = aa_edge_combos[last_combo_idx].total_mass;
	int last_mass_idx = (int)(largest_combo_mass + 1.0);

	combo_start_idxs.resize(last_mass_idx+5,-1);

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		int mass_idx = (int)(aa_edge_combos[i].total_mass);
		if (combo_start_idxs[mass_idx]<0)
			combo_start_idxs[mass_idx]=i;
	}

	// fill in combo idxs for entries with -1
	int previous=-1;
	for (i=0; i<combo_start_idxs.size(); i++)
	{
		if (combo_start_idxs[i]<0)
		{
			combo_start_idxs[i]=previous;
		}
		else
			previous = combo_start_idxs[i];
	}

	
}


const int * Config::get_first_variant_ptr(int combo_idx) const
{ 
	return &variant_vector[aa_edge_combos[combo_idx].variant_start_idx]; 
}


// returns a pointer and number of variants which have masses in the given mass ranges
int Config::get_ptrs_for_combos_in_mass_range(mass_t min_mass, mass_t max_mass, 
												   int& num_combos) const
{
	num_combos = 0;
	int min_idx = (int)(min_mass);

	if (min_idx>=combo_start_idxs.size())
		return -1;

	int combo_idx = combo_start_idxs[min_idx];
	while (combo_idx<aa_edge_combos.size())
	{
		if (aa_edge_combos[combo_idx].total_mass>max_mass)
			return -1;

		if (aa_edge_combos[combo_idx].total_mass>=min_mass)
			break;

		combo_idx++;
	}

	if (combo_idx == aa_edge_combos.size())
		return -1;

	int combos_start_idx = combo_idx;
	while (combo_idx<aa_edge_combos.size())
	{
		if (aa_edge_combos[combo_idx].total_mass>max_mass)
			break;

		combo_idx++;
	}

	num_combos = combo_idx-combos_start_idx;
	
	return combos_start_idx;
}

// returns true if there is a combo that contains the ordered variant of the given aas
bool Config::combos_have_variant(const vector<int>& combos, int num_aa, int *var_aas) const
{
	int c;
	for (c=0; c<combos.size(); c++)
	{
		const AA_combo& combo =  aa_edge_combos[combos[c]];
		int pos = combo.variant_start_idx;
		int v;

		for (v=0; v<combo.num_variants; v++)
			if (variant_vector[pos++]==num_aa)
			{
				int a;
				for (a=0; a<num_aa; a++)
					if (variant_vector[pos+a] != var_aas[a])
						break;
				
				if (a==num_aa)
					return true;

				pos += num_aa;
			}
	}
	return false;
}






/**********************************************************
	calculates the aa_variants vectors (for terminalss and
	amino acids Ala-Val
***********************************************************/
void Config::set_aa_variants()
{
	const vector<int>& org_aas = this->get_org_aa();
	int a;
	aa_variants.clear();
	aa_variants.resize(Val+1);

	aa_variants[N_TERM].push_back(N_TERM);
	aa_variants[C_TERM].push_back(C_TERM);
	aa_variants[Xle].push_back(Xle);

	for (a=0; a<session_aas.size(); a++)
	{
		int aa = session_aas[a];
		int org_aa=org_aas[aa];
		aa_variants[org_aa].push_back(aa);
	}
}


/**********************************************************
Fills the array of amino acid combinations that can be
double edeges. Based on Yingying Huang et al  2003
***********************************************************/
void Config::fill_allowed_double_edges(bool allow_all)
{
	int i;
	int max_aa = Val;
	for (i=0; i<session_aas.size(); i++)
		if (session_aas[i]>max_aa)
			max_aa= session_aas[i];

	int num_aa = max_aa+1;

	allowed_double_edge.clear();
	double_edge_with_same_mass_as_single.clear();

	allowed_double_edge.resize(num_aa);
	double_edge_with_same_mass_as_single.resize(num_aa);

	for (i=0; i<num_aa; i++)
	{
		allowed_double_edge[i].resize(num_aa,allow_all);
		double_edge_with_same_mass_as_single[i].resize(num_aa,false);
	}

	if (tolerance<0.05)
		allow_all=true;


	if (! allow_all)
	{
		for (i=Ala; i<=Val; i++)
		{
			allowed_double_edge[Pro][i]=true;
			allowed_double_edge[Gly][i]=true;
			allowed_double_edge[Ser][i]=true;
		}

	//	allowed_double_edge[Gly][Ala]=false;
		allowed_double_edge[Ser][Ser]=false;
		allowed_double_edge[Ser][Tyr]=false;
		allowed_double_edge[Ser][Pro]=false;
		allowed_double_edge[Ser][His]=false;
		allowed_double_edge[Ser][Gly]=false;
		allowed_double_edge[Ser][Val]=false;



		allowed_double_edge[Thr][Val]=true;
		allowed_double_edge[Thr][Gln]=true;
		allowed_double_edge[Thr][His]=true;
		allowed_double_edge[Thr][Glu]=true;
		allowed_double_edge[Thr][Asp]=true;

		allowed_double_edge[Asn][Asp]=true;
		allowed_double_edge[Asn][Glu]=true;
		allowed_double_edge[Asn][His]=true;
		allowed_double_edge[Asn][Ile]=true;
		allowed_double_edge[Asn][Leu]=true;
		allowed_double_edge[Asn][Met]=true;
		allowed_double_edge[Asn][Asn]=true;
		allowed_double_edge[Asn][Gln]=true;
		allowed_double_edge[Asn][Val]=true;
		allowed_double_edge[Asn][Tyr]=true;

		allowed_double_edge[Tyr][Asp]=true;
		allowed_double_edge[Tyr][Glu]=true;
		allowed_double_edge[Tyr][Thr]=true;

		for (i=Ala; i<=Val; i++)
		{
			allowed_double_edge[i][Lys]=true;
			allowed_double_edge[i][Arg]=true;
		}

		allowed_double_edge[Phe][Glu]=true;

		allowed_double_edge[Ala][Leu]=true;
		allowed_double_edge[Leu][Ala]=true;
		allowed_double_edge[Ala][Ala]=true;
		allowed_double_edge[Trp][Glu]=true;

		// add these double edges, give penalty if they are used as a double edge
		allowed_double_edge[Ala][Gly]=true;
		allowed_double_edge[Gly][Ala]=true;
		allowed_double_edge[Gly][Gly]=true;
		allowed_double_edge[Ala][Asp]=true;
		allowed_double_edge[Asp][Ala]=true;
		allowed_double_edge[Val][Ser]=true;
		allowed_double_edge[Ser][Val]=true;
		allowed_double_edge[Gly][Glu]=true;
		allowed_double_edge[Glu][Gly]=true;

	}	

	double_edge_with_same_mass_as_single[Ala][Gly]=true;
	double_edge_with_same_mass_as_single[Gly][Ala]=true;
	double_edge_with_same_mass_as_single[Gly][Gly]=true;
	double_edge_with_same_mass_as_single[Ala][Asp]=true;
	double_edge_with_same_mass_as_single[Asp][Ala]=true;
	double_edge_with_same_mass_as_single[Val][Ser]=true;
	double_edge_with_same_mass_as_single[Ser][Val]=true;
	double_edge_with_same_mass_as_single[Gly][Glu]=true;
	double_edge_with_same_mass_as_single[Glu][Gly]=true;

	// add fix for PTMs
	const vector<int>& org_aas = session_tables.get_org_aa();
	for (i=0; i<session_aas.size(); i++)
	{
		int j;
		for (j=0; j<session_aas.size(); j++)
		{
			int aa1 = session_aas[i];
			int aa2 = session_aas[j];

			if (org_aas[aa1]>=Ala && org_aas[aa2]>=Ala &&
				allowed_double_edge[org_aas[aa1]][org_aas[aa2]])
				allowed_double_edge[aa1][aa2]=true;
		}
	}
}


/************************************************************
Add new fragments if they do not appear in list
*************************************************************/
void Config::add_fragment_types(const FragmentTypeSet& fts)
{
	int i;
	for (i=0; i<fts.get_fragments().size(); i++)
	{
		all_fragments.add_fragment_type(fts.get_fragment(i));
	}
}



void Config::print_config_parameters(ostream& os) const
{
	int i;

	if (fragments_file.length()>0)
		os << "#CONF FRAGMENTS_FILE " << fragments_file << endl;

	if (regional_fragment_sets_file.length()>0)
		os << "#CONF REGIONAL_FRAGMENT_SETS_FILE " << regional_fragment_sets_file << endl;

	if (aa_combo_file.length()>0)
		os << "#CONF AA_COMBO_FILE " << aa_combo_file << endl;

//	if (resource_dir.length()>0)
//		os << "#CONF RESOURCE_DIR " << resource_dir << endl;


	os << "#CONF MASS_SPEC_TYPE " << mass_spec_type << "(currently only option 1) " << endl;

	os << "#CONF DIGEST_TYPE " << digest_type << " (" << NON_SPECIFIC_DIGEST << " - non specific, " <<
		 TRYPSIN_DIGEST << " - trypsin)" << endl;

	os << "#CONF MAX_NUMBER_OF_PEAKS_PER_LOCAL_WINDOW " << 
		max_number_peaks_per_local_window << endl;

	os << "#CONF NUMBER_OF_STRONG_PEAKS_PER_LOCAL_WINDOW " <<
		number_of_strong_peaks_per_local_window << endl;


	os << "#CONF LOCAL_WINDOW_SIZE " << setprecision(4) << local_window_size << endl;

	os << "#CONF TOLERANCE " << setprecision(4) << tolerance << endl;

	os << "#CONF PM_TOLERANCE " << setprecision(4) << pm_tolerance << endl;

	os << "#CONF MAX_EDGE_LENGTH " << max_edge_length << endl;

	os << "#CONF TERMINAL_SCORE " << terminal_score << endl;

	os << "#CONF DIGEST_SCORE " << digest_score << endl;

	os << "#CONF FORBIDDEN_PAIR_PENALTY " << forbidden_pair_penalty << endl;

	os << "#CONF NEED_TO_CORRECT_PM " << need_to_estimate_pm << " (0 - no, 1 - yes)" << endl;

	os << "#CONF SIZE_THRESHOLDS " << size_thresholds.size() << " ";
	for (i=0; i<size_thresholds.size(); i++)
	{
		int j;
		os << size_thresholds[i].size() << " ";
		for (j=0; j<size_thresholds[i].size(); j++)
			os << fixed << setprecision(4) << size_thresholds[i][j] << " ";
	}
	os << endl;

	os << "#CONF REGION_THRESHOLDS " << region_thresholds.size();
	for (i=0; i<region_thresholds.size(); i++)
		os << " " << setprecision(4) << region_thresholds[i];

	os << endl;

	for (i=0; i<min_ranges.size(); i++)
	{
		os << "#CONF MASS_EXCLUDE_RANGE " << setprecision(3) << min_ranges[i] << " " << max_ranges[i] << endl;
	}
}


/************************************************************
outputs the selected aas for the different regions
*************************************************************/
void Config::print_table_aas(const ConversionTables& table, 
							 const vector<int>& aas) const
{
	int i;
	cout << aas.size() << " amino acids:" << endl;

	for (i=0; i<aas.size(); i++)
		cout << setw(6) << left << aas[i] << setprecision(4) << right << fixed << setw(8) 
			  << table.get_aa2mass(aas[i]) << "   " << left << 
			  table.get_aa2label(aas[i]) << endl;
}

void Config::print_session_aas() const
{
	cout << endl << "AMINO ACIDS" << endl;
	print_table_aas(session_tables,session_aas);
	cout << "N_TERM " << session_tables.get_aa2mass(N_TERM) << endl;
	cout << "C_TERM " << session_tables.get_aa2mass(C_TERM) << endl;
}


void Config::print_aa_variants() const
{
	int i;
	const vector<string>& aa2label = get_aa2label();

	for (i=0; i<aa_edge_combos.size(); i++)
	{
		cout << left << i << " " << aa_edge_combos[i].num_variants << " , " <<
			aa_edge_combos[i].variant_start_idx << "  ";
		int j;

		int pos = aa_edge_combos[i].variant_start_idx;
		cout << "(" << pos << ") ";
		for (j=0; j<aa_edge_combos[i].num_variants; j++)
		{
			int k;
			int num_aa = variant_vector[pos++];
			cout << " ";
			for (k=0; k<num_aa; k++)
				cout << aa2label[variant_vector[pos++]];
		}
		cout << endl;
	}
}



