#include "DiscretePeakModel.h"

void DiscretePeakModel::clone_charge_model(int source_charge, int target_charge)
{
	config.clone_regional_fragment_sets(source_charge, target_charge);

	if (regional_models.size()<=target_charge)
		regional_models.resize(target_charge+1);

	regional_models[target_charge] = regional_models[source_charge];
}


/********************************************************************
Learns all different table models, and writes them to the resource
directory.
*********************************************************************/
void DiscretePeakModel::train_score_model(const char *name, const FileManager& fm, 
		 int charge, int size_idx, int region_idx)
{
	int c;

	for (c=0; c<=fm.get_max_charge(); c++)
	{
		if (charge>0 && charge != c)
			continue;

		if (fm.get_num_spectra(c) < 1)
			continue;

		learn_parent_dependencies(fm,c);

		learn_q_rand(fm,c);

		print_report(c,0,1);
	}


	convert_probs_to_scores();
	// write score model
	
	write_score_model(name);

	// trains and writes edge models
//	edge_model.train_all_edge_models(fm,this);

//	this->write_model();
//	this->read_model(name);

//	edge_model.read_edge_models(&config,(char *)config.get_model_name().c_str());
	// trains and writes tag models 
//	tag_model.train_models(fm,this);
}







// learns the probabilities of randomly observing peaks
void DiscretePeakModel::learn_q_rand(const FileManager& fm, int charge)
{
	int i;
	FileSet fs;

	const vector<FragmentType>& all_fragments = config.get_all_fragments();
	const int max_frag_type_idx = config.get_all_fragments().size();
	
	fs.select_files(fm,0,2500,-1,-1,charge);

	vector<double> noise_probs;
	vector<double> num_files_with_rank_level; // to nake unbaised calculation of q_rand


	noise_probs.resize(level_thresholds.size(),0.001);
	num_files_with_rank_level.resize(level_thresholds.size()+1,0);

	for (i=0; i<fs.get_total_spectra(); i++)
	{
		AnnotatedSpectrum as;

		fs.get_next_spectrum(fm,&config,&as);
		this->init_model_for_scoring_spectrum(&as);

		int spec_rank = this->get_lowest_level_in_spectrum(&as);
		int j;
		for (j=0; j<=spec_rank; j++)
			num_files_with_rank_level[j]++;

		as.annotate_spectrum(as.get_true_mass_with_19());

		mass_t max_peak_mass = as.get_max_peak_mass();
		if (max_peak_mass > as.get_org_pm_with_19())
			max_peak_mass = as.get_org_pm_with_19();

		mass_t peak_area = max_peak_mass - as.get_min_peak_mass();

		// add annotations to the peak statistics only from the breakages 1 to n-1
		double bin_prob = config.get_pm_tolerance()*2 / peak_area;
		int b;
		vector<Breakage>& breakages = as.get_non_const_breakages();
		const vector<Peak>& peaks = as.get_peaks();
		vector<bool> used_peaks;
		used_peaks.resize(peaks.size(),false);
		for (b=1; b<breakages.size()-1; b++)
		{
			this->set_breakage_peak_levels(&breakages[b]);
			int f;
			for (f=0; f<breakages[b].fragments.size(); f++)
			{
				const int f_idx = breakages[b].fragments[f].frag_type_idx;
				const int p_idx = breakages[b].fragments[f].peak_idx;
				used_peaks[p_idx]=true;
			}
		}

		// add unused peaks to noise statistics
		const vector< vector<PeakAnnotation> >& peak_annotations = as.get_peak_annotations();

		int p_idx;
		for (p_idx =0; p_idx<peaks.size(); p_idx++)
		{
			if (! used_peaks[p_idx] && peak_annotations[p_idx].size() == 0)
			{
				noise_probs[get_peak_level(p_idx)]+=bin_prob;
			}
		}
	}

	double inten_prob = 0;
	for (i=1; i<noise_probs.size(); i++)
	{
		if (num_files_with_rank_level[i]>0)
		{
			noise_probs[i] /= num_files_with_rank_level[i];
			inten_prob += noise_probs[i];
		}
	}
	noise_probs[0] = 1.0 - inten_prob;

	// correct frag_probs to acocunt for noise (some of the annotated peaks might be noise)

	// calculate q values
	q_rand = noise_probs;
}




void DiscretePeakModel::print_frag_probs(const vector<int>& frag_type_idxs, 
										 int charge, int size, int region, ostream& os) const
{
	const vector<FragmentType>& all_fragments = config.get_all_fragments();

	int f;
	int len=15;
	os << "Rank     " << "q_rand   ";
	for (f=0; f<frag_type_idxs.size(); f++)
	{
		int f_idx = frag_type_idxs[f];
		int sl = (len - all_fragments[f_idx].label.length())/2;
		os << setw(sl) << left << " " << setw(len -sl) << left << all_fragments[f_idx].label;
	}
	os << endl;
	os << "                  ";
	for (f=0; f<frag_type_idxs.size(); f++)
		os << " prob   score  ";

	os<< endl;


//	regional_models[2][0][1].independent_frag_tables[0].

//	const vector< vector<double> >& f_probs = regional_models[charge][size][region].independent_frag_tables;

	vector<double> totals;
	totals.resize(frag_type_idxs.size(),0);
	double total_rand=0;
	int r;
	for (r=0; r<level_thresholds.size(); r++)
	{
		if (r==0)
		{
			os << setw(8) << left << "no peak";
		}
		else if (level_thresholds[r]-level_thresholds[r-1]==1)
		{
			os << setw(8) << left << level_thresholds[r];
		}
		else
		{
			os << setw(3) << left << level_thresholds[r-1]+1 << "-" << setw(3) << right << level_thresholds[r] << " ";
		}
		os << " " << setprecision(5) << left << q_rand[r] << "   " ;
		total_rand+=q_rand[r];
		for (f=0; f<frag_type_idxs.size(); f++)
		{
			int f_idx = frag_type_idxs[f];
			score_t frag_prob = regional_models[charge][size][region].independent_frag_tables[f_idx].get_score_prob_from_idx(r); 
			os << setprecision(4) << frag_prob << "  ";
			os << setw(5) << setprecision(2) << log(frag_prob / q_rand[r]) << "  ";
			totals[f]+= frag_prob;
		}
		os << endl;
		if (r == 0)
			os << endl;
	}

	os << endl;
	os << "Toatal:  ";
	os << setw(10) << left << setprecision(4) << total_rand;
	for (f=0; f<frag_type_idxs.size(); f++)
	{
		os << setw(len) << left << setprecision(4) << totals[f];
	}
	os << endl;
}


void DiscretePeakModel::print_report(int charge, int size, int region,
									 ostream& os) const
{
	const vector<int> frag_type_idxs = config.get_regional_fragment_type_idxs(charge,size,region);

	int i;
	for (i=0; i<frag_type_idxs.size(); i+=4)
	{
		vector<int> f_idxs;
		f_idxs.clear();
		int f;
		for (f=i; f<i+4 && f<frag_type_idxs.size(); f++)
			f_idxs.push_back(frag_type_idxs[f]);

//		this->print_frag_probs(f_idxs, charge, size, region, os);
		os << endl << endl;
	}

	
}



struct idx_dkl {
	bool operator< ( const idx_dkl& other) const
	{
		return dkl_sum > other.dkl_sum;
	}

	int idx;
	double dkl_sum;
};

/*************************************************************************
// creates models allowing each fragment to have upto two parents
// (chooses the ones that give the best difference in probability
// in terms of DKL). Parent fragments must appear before the current
// fragment in the order of the regional fragment set
*************************************************************************/
void DiscretePeakModel::learn_parent_dependencies(const FileManager& fm, int charge)
{
	FileSet fs;
	const vector< vector< vector< RegionalFragments > > >& regional_fragments = config.get_regional_fragment_sets();

	fs.select_files(fm,0,1E6,-1,-1,charge);


	// init regional models
	if (regional_models.size() < charge+1)
		regional_models.resize(charge+1);

	// This holds all the possible tables (all combination of up to 2 parents)
	// for each fragments in each region
	vector< vector< vector< vector< FragProbTable > > > > all_tables; // size,region,frag,table

	int size_idx;
	regional_models[charge].resize(regional_fragments[charge].size());
	all_tables.resize(regional_fragments[charge].size());
	for (size_idx=0; size_idx<regional_fragments[charge].size(); size_idx++)
	{
		int region_idx;
		regional_models[charge][size_idx].resize(regional_fragments[charge][size_idx].size());
		all_tables[size_idx].resize(regional_fragments[charge][size_idx].size());

		cout << "SIZE " << size_idx << endl;

		for (region_idx=0; region_idx<regional_fragments[charge][size_idx].size(); region_idx++)
		{
			const vector<int>& frag_type_idxs = regional_fragments[charge][size_idx][region_idx].get_frag_type_idxs();
			regional_models[charge][size_idx][region_idx].frag_type_idxs = frag_type_idxs;
			regional_models[charge][size_idx][region_idx].independent_frag_tables.resize(frag_type_idxs.size());
			regional_models[charge][size_idx][region_idx].single_breakage_tables.resize(frag_type_idxs.size());

			all_tables[size_idx][region_idx].resize(frag_type_idxs.size());

			cout << "  REGION " << region_idx << endl;

			// create all possible tables for a given fragments (set parents according
			// to order in frag_type_idxs)
			int f;
			for (f=0; f<frag_type_idxs.size(); f++)
			{
				all_tables[size_idx][region_idx][f].clear();

				cout << "     FRAG " << config.get_fragment(frag_type_idxs[f]).label << endl;

				// first add table with no parents
				FragProbTable fpt;
				vector<int> fields,num_vals;
				fields.resize(5,-1);
				fields[0]=frag_type_idxs[f];
				num_vals.resize(5,0);
				num_vals[0]= num_peak_levels;

				fpt.init_fields(&config,fields,num_vals);
				fpt.init_counts();
				 
				all_tables[size_idx][region_idx][f].push_back(fpt);

				// add all tables with one parent
				if (f<1)
					continue;

				int p;
				for (p=0; p<f; p++)
				{
					FragProbTable fpt;
					vector<int> fields,num_vals;
					fields.resize(5,-1);
					fields[0]=frag_type_idxs[f];
					fields[1]=frag_type_idxs[p];
					num_vals.resize(5,0);
					num_vals[0]= num_peak_levels;
					num_vals[1] = 2 ; // binary values
					num_vals[1] = num_peak_levels;

					fpt.init_fields(&config,fields,num_vals);
					fpt.init_counts();
					 
					all_tables[size_idx][region_idx][f].push_back(fpt);
				}

				// add all tables with two parents
			//	if (f<2)
					continue;

				int p1,p2;
				for (p1=0; p1<f; p1++)
					for (p2=0; p2<f; p2++)
					{
						if (p1 == p2)
							continue;

						FragProbTable fpt;
						vector<int> fields,num_vals;
						fields.resize(5,-1);
						fields[0]=frag_type_idxs[f];
						fields[1]=frag_type_idxs[p1];
						fields[2]=frag_type_idxs[p2];
						num_vals.resize(5,0);
						num_vals[0]= num_peak_levels;
						num_vals[1] = 2 ; // binary values
						num_vals[1] = num_peak_levels;
						num_vals[2] = 2 ; // binary values

						if (p1>p2 && num_vals[1] == num_vals[2])
							continue;

						fpt.init_fields(&config,fields,num_vals);
						fpt.init_counts();
						 
						all_tables[size_idx][region_idx][f].push_back(fpt);
					}
			}
		}

	}


	// Add instances to all tables
	int counter=0;
	while(1)
	{
		AnnotatedSpectrum as;
		if (! fs.get_next_spectrum(fm,&config,&as))
			break;

	//	if (counter>=300)
	//		break;

		as.annotate_spectrum(as.get_true_mass_with_19());
		const int size_idx = as.get_size_idx();
		init_model_for_scoring_spectrum(&as);

		cout << counter++ << " " << as.get_file_name() << endl;
		

		int b;
		vector<Breakage>& breakages = as.get_non_const_breakages();
		for (b=0; b<breakages.size(); b++)
		{
			set_breakage_peak_levels(&breakages[b]);
			int region_idx = config.calc_region_idx(breakages[b].mass,
				as.get_true_mass_with_19(), as.get_charge(), as.get_min_peak_mass(),as.get_max_peak_mass());
			
			// add an instance to each table
			cout << region_idx;
			const vector<int>& frag_type_idxs = regional_models[charge][size_idx][region_idx].frag_type_idxs;
			int f,t;

			// only add instances if the expected peak position is within the min/max range
			for (f=0; f<frag_type_idxs.size(); f++)
			{
				mass_t exp_frag_mass = config.get_fragment(frag_type_idxs[f]).calc_expected_mass(breakages[b].mass,as.get_true_mass_with_19());
				if (exp_frag_mass < as.get_min_peak_mass() || exp_frag_mass > as.get_max_peak_mass())
					continue;

				for (t=0; t<all_tables[size_idx][region_idx][f].size(); t++)
				{
					if ( ! all_tables[size_idx][region_idx][f][t].are_relevant_fragments_visible(&breakages[b]) )
						 continue;

					all_tables[size_idx][region_idx][f][t].add_instance(&breakages[b]);
				}
			}
		}
		cout << endl;
	}

	// calculate probabilities and choose the best table for each fragment
	// (the one that has maximal dkl).
	// also assign for each fragment a table of independent probabilities
	// (what you get if the fragment has no parents).
	for (size_idx =0; size_idx < all_tables.size(); size_idx++)
	{
		int region_idx;
		for (region_idx=0; region_idx<all_tables[size_idx].size(); region_idx++)
		{
			const vector<int>& frag_type_idxs = regional_fragments[charge][size_idx][region_idx].get_frag_type_idxs();

			int f;
			for (f=0; f<all_tables[size_idx][region_idx].size(); f++)
			{
				int t;
				vector<double> dkl_sums;
				dkl_sums.resize(all_tables[size_idx][region_idx][f].size(),0);

				for (t=0; t<all_tables[size_idx][region_idx][f].size(); t++)
					all_tables[size_idx][region_idx][f][t].calc_probs();

				vector<score_t> ind_probs;
				ind_probs.resize(num_peak_levels);
				for (t=0; t<num_peak_levels; t++)
					ind_probs[t] = all_tables[size_idx][region_idx][f][0].score_probs[t];

				cout << "Tables for " << config.get_fragment(frag_type_idxs[f]).label << " ( charge: " <<
					charge << "  size: " << size_idx << "  region: " << region_idx << " )" << endl;

				all_tables[size_idx][region_idx][f][0].print_pretty(&config);
				cout<< endl;

				regional_models[charge][size_idx][region_idx].charge = charge;
				regional_models[charge][size_idx][region_idx].size_idx = size_idx;
				regional_models[charge][size_idx][region_idx].region_idx = region_idx;

				// set the independent probs 
				regional_models[charge][size_idx][region_idx].independent_frag_tables[f] =
					all_tables[size_idx][region_idx][f][0];

				// use the independent model if this fragment has no parents
				if (all_tables[size_idx][region_idx][f].size() == 1)
				{
					regional_models[charge][size_idx][region_idx].single_breakage_tables[f] =
						all_tables[size_idx][region_idx][f][0];

					continue;
				}

				vector<idx_dkl> pairs;
				pairs.clear();
				for (t=1; t<all_tables[size_idx][region_idx][f].size(); t++)
				{
					idx_dkl p;

					p.idx = t;
					p.dkl_sum = all_tables[size_idx][region_idx][f][t].calc_dkl_sum(ind_probs);

				//	all_tables[size_idx][region_idx][f][t].print_pretty(&config);	
				//	cout << "DKL SUM: " << all_tables[size_idx][region_idx][f][t].calc_dkl_sum(ind_probs) << endl << endl;
					pairs.push_back(p);
				}

				sort(pairs.begin(),pairs.end());
				int i;
				for (i=0; i<pairs.size() && i<7; i++)
				{
					string name;
					all_tables[size_idx][region_idx][f][pairs[i].idx].make_table_name(&config,name);
					cout << setw(6) << setprecision(4) << pairs[i].dkl_sum << " " << name << endl;
				}
				cout << endl;

				// find combo with highest dkl, might want to give a bit preference to 
				// tables with 1 parent, or strong parents...

				int two_frag_idx=-1, one_frag_idx =-1;
				double max_two_frag_dkl=0, max_one_frag_dkl=0;

				for (i=0; i<pairs.size(); i++)
				{
					const vector<int>& fields = all_tables[size_idx][region_idx][f][pairs[i].idx].get_fields();
					if (fields[1]>=0 && fields[2]>=0)
					{
						if (pairs[i].dkl_sum>max_two_frag_dkl)
						{
							max_two_frag_dkl = pairs[i].dkl_sum;
							two_frag_idx = pairs[i].idx;
						}
					}
					else
					{
						if (pairs[i].dkl_sum>max_one_frag_dkl)
						{
							max_one_frag_dkl = pairs[i].dkl_sum;
							one_frag_idx = pairs[i].idx;
						}
					}
				}

				int selected_table_idx=-1;
				if (max_two_frag_dkl>max_one_frag_dkl && (max_two_frag_dkl/max_one_frag_dkl > 1.1))
				{
					selected_table_idx = two_frag_idx;
				}
				else
					selected_table_idx = one_frag_idx;

				regional_models[charge][size_idx][region_idx].single_breakage_tables[f] =
						all_tables[size_idx][region_idx][f][selected_table_idx];

				string name;
				all_tables[size_idx][region_idx][f][selected_table_idx].make_table_name(&config,name);
				cout << "Selected: " << name << endl << endl; 
			}
		}
	}

	

}



// writes all relevant info:
// thresholds and tables
void DiscretePeakModel::write_score_model(const char *name) const
{
	string file_name = config.get_resource_dir() + "/" + string(name) + "_break_score.txt";
	ofstream os(file_name.c_str());

	if (! os.is_open() || ! os.good())
	{
		cout << "Error: couldn't open score model for writing: " << file_name << endl;
		exit(1);
	}
	write_model_specifics(os);

	int c;
	int max_charge = 0;
	for (c=0; c<regional_models.size(); c++)
		if (regional_models[c].size()>0)
			max_charge = c;

	os << "#MAX_CHARGE " << max_charge << endl;

	for (c=1; c<=max_charge; c++)
	{
		if (regional_models[c].size() == 0)
			continue;

		os << c << " " << regional_models[c].size() << endl; // num sizes
		int size_idx;
		for (size_idx=0; size_idx< regional_models[c].size(); size_idx++)
		{
			os << regional_models[c][size_idx].size() << endl; // num regions
			int region_idx;
			for (region_idx=0; region_idx<regional_models[c][size_idx].size(); region_idx++)
			{
				regional_models[c][size_idx][region_idx].write_regional_model(os);
			}
		}
	}

	os.close();
}



// converts all the probabilities to socres in the tables
void DiscretePeakModel::convert_probs_to_scores()
{
	int c,s,r;

	for (c=0; c<regional_models.size(); c++)
		for (s=0; s<regional_models[c].size(); s++)
			for (r=0; r<regional_models[c][s].size(); r++)
				regional_models[c][s][r].convert_to_scores(q_rand);

}


void DiscretePeakModel::read_score_model(const char *name, bool silent_ind)
{
	char buff[128];
	string file_name = config.get_resource_dir() + "/" + string(name) + "_break_score.txt";
	ifstream is(file_name.c_str());
	if (! is.good() || ! is.is_open())
	{
		cout << "Error: couldn't open break score model for reading: " << file_name << endl;
		exit(1);
	}

	read_model_specifics(is);
	
	is.getline(buff,128);

	if (sscanf(buff,"#MAX_CHARGE %d",&max_score_model_charge) != 1)
	{
		cout << "Error: bad line in score model: " << buff << endl;
		exit(1);
	}

	this->config.set_max_charge_for_size(max_score_model_charge);

	regional_models.resize(max_score_model_charge+3);

	int charge = 0;
	while (charge<max_score_model_charge)
	{
		int num_sizes,size_idx;
		is.getline(buff,128);
		if (sscanf(buff,"%d %d",&charge,&num_sizes) != 2)
		{
			cout << "Error: bad line in score model: " << buff << endl;
			exit(1);
		}

		regional_models[charge].resize(num_sizes);
		for (size_idx=0; size_idx<num_sizes; size_idx++)
		{
			int num_regions, region_idx;
			is.getline(buff,128);
			if (sscanf(buff,"%d",&num_regions) != 1)
			{
				cout << "Error: bad line in score model: " << buff << endl;
				exit(1);
			}
			regional_models[charge][size_idx].resize(num_regions);
			for (region_idx=0; region_idx<num_regions; region_idx++)
			{
				regional_models[charge][size_idx][region_idx].read_regional_model(&config,is);
				int f;
				for (f=0; f<regional_models[charge][size_idx][region_idx].single_breakage_tables.size(); f++)
				{
					regional_models[charge][size_idx][region_idx].single_breakage_tables[f].set_model_tolerance(config.get_tolerance());
					regional_models[charge][size_idx][region_idx].single_breakage_tables[f].config=&config;
				}
			}
		}
	}

	// copy models for other charges. Prefer charge 2 model
	int min_model_idx = -1;
	int max_model_idx = -1;

	if (max_score_model_charge<4)
		max_score_model_charge=4;

	int c;
	for (c=0; c<=max_score_model_charge; c++)
		if (regional_models[c].size() >0)
		{
			min_model_idx=c;
			break;
		}
	
	for (c=max_score_model_charge; c>=0; c--)
		if (regional_models[c].size() >0)
		{
			max_model_idx=c;
			break;
		}

	if (min_model_idx<0)
	{
		cout << "Error: no score models were read!" << endl;
		exit(1);
	}

	bool has_charge_2 = (regional_models[2].size()>0);
	if (regional_models[1].size()==0)
	{
		if (has_charge_2)
		{
			regional_models[1]=regional_models[2];
		}
		else
			regional_models[1]=regional_models[min_model_idx];
	}

	for (c=2; c<=max_score_model_charge; c++)
	{
		if (regional_models[c].size() == 0)
			regional_models[c]=regional_models[max_model_idx];
	}

	// read edge models..
	edge_model.read_edge_models(&config,(char *)config.get_model_name().c_str());


	is.close();

}





// calcs the score of participating fragments that have istopic peak scores
void DiscretePeakModel::add_breakage_iso_score(Spectrum *spec, 
												Breakage *breakage) const
{
	if (! breakage)
		return;

	const vector<int>& strong_frags = config.get_regional_fragment_sets()
		[breakage->parent_charge][breakage->parent_size_idx][breakage->region_idx].get_strong_frag_type_idxs();
	

	int i;
	for (i=0; i<strong_frags.size(); i++)
	{
		int pos = breakage->get_position_of_frag_idx(strong_frags[i]);
		if (pos>=0)
		{
			float iso_level = spec->get_peak_iso_level(breakage->fragments[pos].peak_idx);
			if (iso_level>0)
				breakage->score -= (2 + iso_level * 3);
		}
	}
}

void DiscretePeakModel::print_all_table_names(ostream& os) const
{
	int i,j,k;

	for (i=0; i<regional_models.size(); i++)
		for (j=0; j<regional_models[i].size(); j++)
			for (k=0; k<regional_models[i][j].size(); k++)
				regional_models[i][j][k].print_table_names(&config,os);
}


// prints joint scores of two top  fragments
void DiscretePeakModel::print_joint_scores(ostream& os) const
{
	print_level_legend(os);
	RegionalPepNovoModel rpm = regional_models[2][0][1];

	int f1_idx = rpm.frag_type_idxs[0];
	int f2_idx = rpm.frag_type_idxs[1];

	os << "joint score for fragments " << config.get_fragment(f1_idx).label << " and " <<
		config.get_fragment(f2_idx).label << endl;

	vector< vector<score_t> > joint_scores;
	joint_scores.resize(num_peak_levels);
	int i;
	for (i=0; i<num_peak_levels; i++)
		joint_scores[i].resize(num_peak_levels,NEG_INF);

	for (i=0; i<num_peak_levels; i++)
	{
		int j;
		for (j=0; j<num_peak_levels; j++)
		{
			table_entry e;

			e[0]=i;
			e[1]=j;
			e[2]=0;
			e[3]=0;
			e[4]=0;
					  
			int idx1= rpm.single_breakage_tables[0].calc_table_idx(e);
			int idx1b= rpm.independent_frag_tables[0].calc_table_idx(e);

			e[0]=j;
			e[1]=i;

			int idx2= rpm.single_breakage_tables[1].calc_table_idx(e);
			int idx2b= rpm.independent_frag_tables[1].calc_table_idx(e);

			joint_scores[i][j]= rpm.single_breakage_tables[0].score_probs[idx1] +
								rpm.single_breakage_tables[1].score_probs[idx2];
			joint_scores[i][j]-= rpm.independent_frag_tables[0].score_probs[idx1b] +
								 rpm.independent_frag_tables[1].score_probs[idx2b];

		}
	}

	os << "   ";
	for (i=0; i<num_peak_levels; i++)
		os << setw(3) << right << i << "    ";
	os << endl << endl;

	for (i=0; i<num_peak_levels; i++)
	{
		int j;
		os << left << setw(3) << i;
		for (j=0; j<num_peak_levels; j++)
		{
			os << setw(6) << right << joint_scores[i][j] << " ";
		}
		os << endl << endl;
	}


}

