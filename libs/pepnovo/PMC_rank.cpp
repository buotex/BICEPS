#include "PMCSQS.h"
#include "auxfun.h"

PMCSQS_Scorer::~PMCSQS_Scorer(){ //BX

for (int i = 0; i < pmc_rank_models.size() ; i++)
	for (int j = 0; j < pmc_rank_models[i].size(); j++)
		if (pmc_rank_models[i][j])
            delete pmc_rank_models[i][j];

for (int i = 0; i < sqs_models.size(); i++)
	for (int j = 0; j < sqs_models[i].size(); j++)
		for (int k =0; k < sqs_models[i][k].size(); k++)
			if (sqs_models[i][j][k])
                delete sqs_models[i][j][k];

}

bool PMCSQS_Scorer::read_pmc_rank_models(Config *_config, char *file_name)
{
	config = _config;

	string path = config->get_resource_dir() + "/" + file_name;
	ifstream in_stream(path.c_str(),ios::in);
	if (! in_stream.good())
	{
		cout << "Warning: couldn't open pmc rank model for reading: " << path << endl;
		return false;
	}


	char buff[512];
	int num_charges=-1;

	in_stream.getline(buff,256);
	istringstream iss1(buff);

	frag_pair_sum_offset=NEG_INF;
	bin_increment=NEG_INF;
	iss1 >> bin_increment >> this->frag_pair_sum_offset;
	if (frag_pair_sum_offset==NEG_INF || bin_increment == NEG_INF)
	{
		cout << "Error in pmc model file!" << endl;
		exit(1);
	}

	in_stream.getline(buff,256);
	istringstream iss(buff);

	iss >> num_charges;
	max_model_charge=num_charges-1;
	
	pmc_rank_models.resize(num_charges);
	pmc_rank_mass_thresholds.resize(num_charges);
	pmc_charge_mz_biases.resize(num_charges);


	int i;
	for (i=0; i<num_charges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_threshes=0;
		iss >> num_threshes;
		
		pmc_rank_mass_thresholds[i].resize(num_threshes,NEG_INF);
		int j;
		for (j=0; j<num_threshes; j++)
			iss >> pmc_rank_mass_thresholds[i][j];
	}

	for (i=0; i<num_charges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_biases=0;
		iss >> num_biases;
		
		pmc_charge_mz_biases[i].resize(num_biases,NEG_INF);
		int j;
		for (j=0; j<num_biases; j++)
			iss >> pmc_charge_mz_biases[i][j];
	}
	
	// read Boost models
	for (i=0; i<num_charges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);

		int num_models=-1;
		iss >> num_models;

		if (num_models<0)
		{
			cout << "Error: bad parsing of PMCR model file!" << endl;
			exit(0);
		}
		pmc_rank_models[i].resize(num_models,NULL);

		int j;
		for (j=0; j<num_models; j++)
		{
			pmc_rank_models[i][j]=new RankBoostModel;
			pmc_rank_models[i][j]->read_rankboost_model(in_stream);
		}
		
	}
	in_stream.close();

//	this->write_pmc_rank_models("XXX.txt");

	this->ind_initialized_pmcr = true;
	return true;
}


void PMCSQS_Scorer::write_pmc_rank_models(const char *path) const
{
	ofstream out_stream(path,ios::out);
	if (! out_stream.good())
	{
		cout << "Error: couldn't open pmc model for writing: " << path << endl;
		exit(1);
	}

	out_stream << this->bin_increment << " " << this->frag_pair_sum_offset << endl;
	out_stream << this->pmc_rank_mass_thresholds.size() << endl;
	out_stream << setprecision(3) << fixed;

	int i;
	for (i=0; i<this->pmc_rank_mass_thresholds.size(); i++)
	{
		out_stream << pmc_rank_mass_thresholds[i].size();
		int j;
		for (j=0; j<pmc_rank_mass_thresholds[i].size(); j++)
			out_stream << " " << pmc_rank_mass_thresholds[i][j];
		out_stream << endl;
	}

	
	for (i=0; i<this->pmc_charge_mz_biases.size(); i++)
	{
		out_stream << pmc_charge_mz_biases[i].size();
		int j;
		for (j=0; j<pmc_charge_mz_biases[i].size(); j++)
			out_stream << " " << pmc_charge_mz_biases[i][j];
		out_stream << endl;
	}

	// write RankBoost models
	for (i=0; i<pmc_rank_models.size(); i++)
	{
		int j;
		
		if (pmc_rank_models[i].size()==1 && ! pmc_rank_models[i][0])
		{
			out_stream << 0 << endl;
			continue;
		}

		out_stream << pmc_rank_models[i].size() << endl;
		for (j=0; j<pmc_rank_models[i].size(); j++)
		{
			if (pmc_rank_models[i][j])
			{
				pmc_rank_models[i][j]->write_rankboost_model(out_stream,true);
			}
			else
			{
				cout << "Error: non intialized rank pmc model!" << endl;
				exit(1);
			}
		}
	}
	
	out_stream.close();
}


void PMCSQS_Scorer::set_pmc_mass_thresholds(int option)
{
	if (option==0)
	{
		pmc_rank_mass_thresholds = config->get_size_thresholds();
		int i;
		for (i=0; i<pmc_rank_mass_thresholds.size(); i++)
		{
			if (pmc_rank_mass_thresholds[i].size()>0)
			{
				if (pmc_rank_mass_thresholds[i][pmc_rank_mass_thresholds[i].size()-1]>10000)
					pmc_rank_mass_thresholds[i].pop_back();
			}
		}
	}

	if (option==1)
	{
		pmc_rank_mass_thresholds.resize(4);
		pmc_rank_mass_thresholds[1].push_back(1150.0);
		pmc_rank_mass_thresholds[1].push_back(1400.0);
 
		pmc_rank_mass_thresholds[2].push_back(1100.0);
		pmc_rank_mass_thresholds[2].push_back(1300.0);
		pmc_rank_mass_thresholds[2].push_back(1600.0);
		pmc_rank_mass_thresholds[2].push_back(1900.0);
		pmc_rank_mass_thresholds[2].push_back(2400.0);

		pmc_rank_mass_thresholds[3].push_back(1950.0);
		pmc_rank_mass_thresholds[3].push_back(2450.0);
		pmc_rank_mass_thresholds[3].push_back(3000.0);
	}

	if (option==2)
	{
		pmc_rank_mass_thresholds.resize(4);
		
		pmc_rank_mass_thresholds[2].push_back(1300.0);
	
		pmc_rank_mass_thresholds[2].push_back(1900.0);
	}
}



void convert_ME_to_RankBoostSample(const ME_Regression_Sample& me,
								   RankBoostSample& rbs)
{
	rbs.clear();
	int i;
	for (i=0; i<me.f_vals.size(); i++)
		rbs.add_real_feature(me.f_vals[i].f_idx,me.f_vals[i].val);
}


/******************************************************************************
Train PMC models from positive example files
*******************************************************************************/
void PMCSQS_Scorer::train_pmc_rank_models(Config *config, const FileManager& fm, 
										  int sel_charge, bool overwrite)
{	
	const bool sample_diagnostic = false;
	const vector<int>& spectra_counts = fm.get_spectra_counts();
	
	max_model_charge=0;

	int charge;
	for (charge=1; charge<spectra_counts.size(); charge++)
	{
		if (spectra_counts[charge]>=MIN_SPECTRA_FOR_PMCSQS_MODEL)
			max_model_charge=charge;
	}

	const int max_to_read_per_file = 40000;
	
	vector<string> real_names;
	init_PMC_feature_names(real_names);


	// try and read existing pmc model, otherwise init a new one
	string pmc_path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	ifstream model_stream(pmc_path.c_str());
	if (model_stream.is_open() && model_stream.good())
	{
		model_stream.close();
		string pmcr_name = config->get_model_name() + "_PMCR.txt";
		const char *path = pmc_path.c_str();
		this->read_pmc_rank_models(config,(char *)pmcr_name.c_str());
	}
	else
	{
		set_pmc_mass_thresholds();
	
		this->set_frag_pair_sum_offset(MASS_PROTON); // b+y - PM+19
		this->set_bin_increment(0.1);
		pmc_rank_models.resize(pmc_rank_mass_thresholds.size());
		pmc_charge_mz_biases.resize(pmc_rank_mass_thresholds.size());
	}
	
	const double prop_train = 0.5;


	// It is assumed that the mass thresholds were set according to the training data
	// (this is done manually with values encoded in the set_mass_threhsolds function)
	for (charge=1; charge<=max_model_charge; charge++)
	{
		if (sel_charge>0 && charge != sel_charge)
			continue;

		const int num_sizes = pmc_rank_mass_thresholds[charge].size();
		pmc_rank_models[charge].resize(num_sizes+1,NULL);
		pmc_charge_mz_biases[charge].resize(num_sizes+1,0);

		
		int size_idx;
		for (size_idx=0; size_idx<=num_sizes; size_idx++)
		{
			if (pmc_rank_models[charge][size_idx] && ! overwrite)
				continue;

			vector<SingleSpectrumFile *> test_ssfs;
			BasicSpecReader bsr;
			static QCPeak peaks[5000];
			RankBoostDataset train_ds, test_ds, pos_ds, neg_ds;

			mass_t min_mass =0;
			mass_t max_mass = POS_INF;

			if (size_idx>0)
				min_mass = pmc_rank_mass_thresholds[charge][size_idx-1];

			if (size_idx<num_sizes)
				max_mass = pmc_rank_mass_thresholds[charge][size_idx];

			// these ranges are given according to pm_with_19
			// so files should be selected through select_files and not
			// select_file_in_mz_range
			FileSet fs;		
			fs.select_files(fm,min_mass,max_mass,-1,-1,charge);

			if (fs.get_total_spectra()<500)
				continue;

			
			int num_groups_in_train=0;
			int num_groups_in_test=0;

			cout << "TRAINING charge " << charge << " size " << size_idx << "  (" <<
				min_mass << "-" << max_mass << ")" << endl;

			fs.randomly_reduce_ssfs(max_to_read_per_file);
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			const int num_samples = all_ssf.size();
			
			// first find the bias in number of bins between the true m/z bin and
			// the optimal m/z bin
			vector<bool> skipped_idxs;
			skipped_idxs.resize(num_samples,false);
			int skipped_bad_mz=0;
			mass_t total_bias=0;
			int i;
			for (i=0; i<num_samples; i++)
			{
				SingleSpectrumFile* ssf = all_ssf[i];
				BasicSpectrum bs;
			
				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
				bs.peaks = peaks;
				bs.ssf = ssf;

				ssf->peptide.calc_mass(config);
				
				const mass_t true_mz = (ssf->peptide.get_mass()+MASS_H2O+(mass_t)charge)/(mass_t)charge;

				if (fabs(true_mz - bs.ssf->m_over_z)>2.5)
				{
					//cout << setprecision(2) << true_mz << " <---> " << bs.ssf->m_over_z << " skipping" << endl;
					skipped_bad_mz++;
					skipped_idxs[i]=true;
					continue;
				} 

				init_for_current_spec(config,bs);
				calculate_curr_spec_pmc_values(bs, bin_increment);

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = curr_spec_rank_pmc_tables[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<true_mz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-true_mz>true_mz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				if (opt_bin_idx <=0 || opt_bin_idx == pmc_stats.size()-1)
				{
					skipped_bad_mz++;
					skipped_idxs[i]=true;
					continue;
				}

				total_bias += (pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z);

				if (fabs(pmc_stats[opt_bin_idx].m_over_z - pmc_stats[true_mz_bin_idx].m_over_z)>4.0)
				{
					cout << "opt bin: " << opt_bin_idx << " (" << pmc_stats[opt_bin_idx].m_over_z << ")  ";
					cout << "tru bin: " << true_mz_bin_idx << " ("<< pmc_stats[true_mz_bin_idx].m_over_z << ")" << endl;
				}
			} 

			mass_t mz_bias = total_bias / (mass_t)(num_samples-skipped_bad_mz);
			pmc_charge_mz_biases[charge][size_idx]=mz_bias;

			cout << "m/z bias: " << setprecision(4) << mz_bias << endl;
			cout << "skipped " << skipped_bad_mz << "/" << num_samples <<
				"  because of m/z more than 2.5 away from observed..." << endl; 

		//	pmc_charge_mz_biases[charge][size_idx] = 0;

			for (i=0; i<num_samples; i++)
			{
				if (skipped_idxs[i])
					continue;

				SingleSpectrumFile* ssf = all_ssf[i];
				BasicSpectrum bs;
			
				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
				bs.peaks = peaks;
				bs.ssf = ssf;
				const mass_t true_mz = (ssf->peptide.get_mass()+MASS_H2O+(mass_t)charge)/(mass_t)charge;

				init_for_current_spec(config,bs);
				calculate_curr_spec_pmc_values(bs, bin_increment);

				// find the true_mz_bin_idx
				
				const vector<PMCRankStats>& pmc_stats = curr_spec_rank_pmc_tables[charge];
				int true_mz_bin_idx=0;
				while (true_mz_bin_idx<pmc_stats.size() && pmc_stats[true_mz_bin_idx].m_over_z<true_mz)
					true_mz_bin_idx++;

				if (true_mz_bin_idx == pmc_stats.size())
					true_mz_bin_idx--;

				if (true_mz_bin_idx>0 && pmc_stats[true_mz_bin_idx].m_over_z-true_mz>true_mz-pmc_stats[true_mz_bin_idx-1].m_over_z)
					true_mz_bin_idx--;

				int opt_bin_idx = get_optimal_bin(true_mz_bin_idx, charge);

				
				static vector<RankBoostSample> spec_samples;
				fill_RankBoost_smaples_with_PMC(bs, charge, spec_samples);

				// select samples and add them to pmc_ds
				int good_idx;
				vector<int> bad_idxs;
				select_training_sample_idxs(charge,spec_samples,bs,good_idx,bad_idxs);

				const bool ind_add_to_train = (my_random()<prop_train);
				int group_idx;
				
				if (ind_add_to_train)
				{
					group_idx= num_groups_in_train++;	
				}
				else
				{
					group_idx= num_groups_in_test++;
					test_ssfs.push_back(ssf);
				}
				
				
				RankBoostDataset& ds = (ind_add_to_train ? train_ds : test_ds);

				const int pos_index  = ds.get_num_samples();
				spec_samples[good_idx].group_idx = group_idx;
				spec_samples[good_idx].rank_in_group=0;

				ds.add_sample(spec_samples[good_idx]);
				if (sample_diagnostic)
					pos_ds.add_sample(spec_samples[good_idx]);

				int j;
				for (j=0; j<bad_idxs.size(); j++)
				{
					const int bad_idx = bad_idxs[j];
					if (bad_idx < 0 || bad_idx>= spec_samples.size())
						continue;
		
					spec_samples[bad_idx].group_idx=group_idx;
					spec_samples[bad_idx].rank_in_group=1;

					ds.add_to_phi_vector(ds.get_num_samples(),pos_index);
					ds.add_sample(spec_samples[bad_idx]);

					if (sample_diagnostic)
						neg_ds.add_sample(spec_samples[bad_idx]);
				}						   
			}

			train_ds.set_num_groups(num_groups_in_train);
			test_ds.set_num_groups(num_groups_in_test);
			
			train_ds.compute_total_phi_weight();
			train_ds.initialize_potenital_lists();
			train_ds.initialzie_real_feature_table(real_names.size());

			test_ds.compute_total_phi_weight();

			if (pmc_rank_models[charge][size_idx])
				delete pmc_rank_models[charge][size_idx];
			
			pmc_rank_models[charge][size_idx] = new RankBoostModel;
		

			RankBoostModel* boost = pmc_rank_models[charge][size_idx];

			vector<string> empty;
			empty.clear();
			boost->init_rankboost_model_feature_names(empty,real_names);
			boost->init_rankboost_model_for_training(train_ds,100,25);

			train_ds.initialize_real_vote_lists(*boost);

			if (sample_diagnostic)
			{
				boost->summarize_features_pos_neg(pos_ds.get_samples(),neg_ds.get_samples());
			}
			else
				boost->summarize_features(train_ds.get_samples());

			boost->train_rankboost_model(train_ds,4000,NULL,&test_ds);
			
			boost->ouput_ranked_feature_list();

		//	output_pmc_rank_results(fm,charge,test_ssfs);

		//	exit(0);

			ind_initialized_pmcr = true;
		//	string path;
		//	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCRtt.txt";
		//	this->write_pmc_rank_models(path.c_str());
			
		}
	}

	string path;
	path = config->get_resource_dir() + "/" + config->get_model_name() + "_PMCR.txt";
	this->write_pmc_rank_models(path.c_str());
	ind_initialized_pmcr = true;
}


struct offset_pair {
	offset_pair() : offset(POS_INF), inten_sum(0) {};
	offset_pair(mass_t off,float inten) : offset(off), inten_sum(inten) {};
	mass_t offset;
	float inten_sum;
};

bool cmp_offset_pair_offset (const offset_pair& a, const offset_pair& b)
{
	return (a.offset<b.offset);
}

bool cmp_offset_pair_inten (const offset_pair& a, const offset_pair& b)
{
	return (a.inten_sum>b.inten_sum);
}


float calc_mean_abs_offset(const vector<float>& offsets_by_inten)
{
	const float missing_pair_offset = 0.5;
	const int   num_offsets         = 3;

	if (offsets_by_inten.size()==0)
		return 1000;

	float abs_off=0;
	int i;
	for (i=0; i<num_offsets && i<offsets_by_inten.size(); i++)
		abs_off+=fabs(offsets_by_inten[i]);

	abs_off += (3-i)*missing_pair_offset;
	
	return (abs_off/num_offsets);
}


void calc_pmc_rank_stats_for_mass(const QCPeak *peaks, 
										  int num_peaks, 
										  mass_t single_charge_pair_sum,
										  mass_t tolerance, 
										  const vector<float>& iso_levels,
										  const vector<bool>& strong_inds,
										  const vector<bool>& strict_iso_inds,
										  PMCRankStats& stats)
{
	const mass_t min_single_sum = single_charge_pair_sum - tolerance;
	const mass_t max_single_sum = single_charge_pair_sum + tolerance;

	const mass_t min_double_sum = min_single_sum + 1.0;
	const mass_t max_double_sum = max_single_sum + 1.0;
	const mass_t double_charge_pair_sum = single_charge_pair_sum +1.0;

	const mass_t min_single_h2o_sum = min_single_sum - MASS_H2O;
	const mass_t max_single_h2o_sum = max_single_sum - MASS_H2O;
	const mass_t single_charge_pair_h2o_sum = single_charge_pair_sum - MASS_H2O;

	const mass_t min_double_h2o_sum = min_double_sum - MASS_H2O;
	const mass_t max_double_h2o_sum = max_double_sum - MASS_H2O;
	const mass_t double_charge_pair_h2o_sum = double_charge_pair_sum - MASS_H2O;

	static vector<offset_pair> by_pairs,  strong_pairs;
	static vector<offset_pair> c2_pairs,  strong_c2_pairs;
	static vector<offset_pair> h2o_pairs, c2_h2o_pairs;

	by_pairs.clear();
	strong_pairs.clear();
	c2_pairs.clear();
	strong_c2_pairs.clear();
	h2o_pairs.clear();
	c2_h2o_pairs.clear();

	stats.clear();

	int forward_idx = -1;
	int back_idx = num_peaks-1;

	// find pairs of b/y
	while (forward_idx<back_idx)
	{
		forward_idx++;
		if (iso_levels[forward_idx]>0)
		{
			continue;
		}

		while (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass>max_single_sum)
			back_idx--;

		if (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass > min_single_sum)
		{
			if (iso_levels[back_idx]>0)
				continue;

			const mass_t offset = fabs(peaks[forward_idx].mass + peaks[back_idx].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
					
			by_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_frag_pairs += inten_sum;

			if (strong_inds[forward_idx] || strong_inds[back_idx])
			{
				strong_pairs.push_back(offset_pair(offset,inten_sum));
				stats.inten_strong_pairs += inten_sum;
			}
		}
	}

	// find pairs b/y2
	forward_idx = -1;
	back_idx = num_peaks-1;

	const int last_idx =num_peaks-1;
	while (forward_idx<last_idx && back_idx>=0)
	{
		forward_idx++;
		if (iso_levels[forward_idx]>0)
			continue;
			
		mass_t sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		while (sum>max_double_sum)
		{
			back_idx--;
			if (back_idx<0)
				break;
			sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		}

		if (back_idx>=0 && sum > min_double_sum)
		{
			if (iso_levels[back_idx]>0)
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
			
			c2_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_c2_pairs += inten_sum;

			if (strong_inds[forward_idx] || strong_inds[back_idx])
			{
				strong_c2_pairs.push_back(offset_pair(offset,inten_sum));
				stats.inten_c2_strong_pairs = inten_sum;
			}
		}
	}

	// find pairs of b/y-H2O
	forward_idx = -1;
	back_idx = num_peaks-1;

	while (forward_idx<back_idx)
	{
		forward_idx++;
		if (iso_levels[forward_idx]>0)
			continue;

		while (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass>max_single_h2o_sum)
			back_idx--;

		if (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass > min_single_h2o_sum)
		{
			if (iso_levels[back_idx]>0)
				continue;

			const mass_t offset = fabs(peaks[forward_idx].mass + peaks[back_idx].mass - single_charge_pair_h2o_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
					
			h2o_pairs.push_back(offset_pair(offset,inten_sum));
			stats.inten_h2o_loss_frag_pairs += inten_sum;
		}
	}

	// find pairs b/y2 - H2O
	forward_idx = -1;
	back_idx = num_peaks-1;

	while (forward_idx<last_idx && back_idx>=0)
	{
		forward_idx++;
		if (iso_levels[forward_idx]>0)
			continue;
			
		mass_t sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		while (sum>max_double_h2o_sum)
		{
			back_idx--;
			if (back_idx<0)
				break;
			sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		}

		if (back_idx>=0 && sum > min_double_h2o_sum)
		{
			if (iso_levels[back_idx]>0)
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_h2o_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
			
			c2_h2o_pairs.push_back(offset_pair(offset,inten_sum));
			stats.itnen_h2o_loss_c2_frag_pairs += inten_sum;
		}
	}

	stats.num_frag_pairs = by_pairs.size();
	stats.num_strong_frag_pairs = strong_pairs.size();
	stats.num_c2_frag_pairs = c2_pairs.size();
	stats.num_strong_c2_frag_pairs = strong_c2_pairs.size();
	stats.num_h2o_loss_frag_pairs = h2o_pairs.size();
	stats.num_h2o_loss_c2_frag_pairs = c2_h2o_pairs.size();

	int i;

	vector<float>& offset_pairs_ordered_by_inten = stats.offset_pairs_ordered_by_inten;
	sort(by_pairs.begin(),by_pairs.end(),cmp_offset_pair_inten);
	offset_pairs_ordered_by_inten.resize(by_pairs.size());
	for (i=0; i<by_pairs.size(); i++)
		offset_pairs_ordered_by_inten[i]=by_pairs[i].offset;
	stats.mean_offset_pairs=calc_mean_abs_offset(offset_pairs_ordered_by_inten);

	vector<float>& strong_offset_pairs_ordered_by_inten = stats.strong_offset_pairs_ordered_by_inten;
	sort(strong_pairs.begin(),strong_pairs.end(),cmp_offset_pair_inten);
	strong_offset_pairs_ordered_by_inten.resize(strong_pairs.size());
	for (i=0; i<strong_pairs.size(); i++)
		strong_offset_pairs_ordered_by_inten[i]=strong_pairs[i].offset;
	stats.mean_offset_strong_pairs=calc_mean_abs_offset(strong_offset_pairs_ordered_by_inten);

	vector<float>& c2_offset_pairs_ordered_by_inten = stats.c2_offset_pairs_ordered_by_inten;
	sort(c2_pairs.begin(),c2_pairs.end(),cmp_offset_pair_inten);
	c2_offset_pairs_ordered_by_inten.resize(c2_pairs.size());
	for (i=0; i<c2_pairs.size(); i++)
		c2_offset_pairs_ordered_by_inten[i]=c2_pairs[i].offset;
	stats.mean_offset_c2_pairs=calc_mean_abs_offset(c2_offset_pairs_ordered_by_inten);



	// fill in additional iso sum features (look at pairs that sum to expected, expected+1 expected+2)
	
	// find pairs of b/y

	static vector<offset_pair> pairs0,  pairs1, pairs2;
	static vector<offset_pair> c2_pairs0,  c2_pairs1, c2_pairs2;
	
	pairs0.clear();
	forward_idx = -1;
	back_idx = num_peaks-1;
	while (forward_idx<back_idx)
	{
		forward_idx++;
		if (strict_iso_inds[forward_idx])
			continue;

		while (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass>max_single_sum)
			back_idx--;

		if (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass > min_single_sum)
		{
			if (strict_iso_inds[back_idx])
				continue;

			const mass_t offset = fabs(peaks[forward_idx].mass + peaks[back_idx].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
					
			pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs1.clear();
	forward_idx = -1;
	back_idx = num_peaks-1;
	const mass_t max1 = max_single_sum+1.0;
	const mass_t min1 = min_single_sum+1.0;
	while (forward_idx<back_idx)
	{
		forward_idx++;
	
		while (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass>max1)
			back_idx--;

		if (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass > min1)
		{
			if (! (strict_iso_inds[back_idx] || strict_iso_inds[forward_idx]))
				continue;

			const mass_t offset = fabs(peaks[forward_idx].mass + peaks[back_idx].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
					
			pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs2.clear();
	forward_idx = -1;
	back_idx = num_peaks-1;
	const mass_t max2 = max_single_sum+2.0;
	const mass_t min2 = min_single_sum+2.0;
	while (forward_idx<back_idx)
	{
		forward_idx++;
	
		while (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass>max2)
			back_idx--;

		if (back_idx>=0 && peaks[forward_idx].mass + peaks[back_idx].mass > min2)
		{
			if (! (strict_iso_inds[back_idx] || strict_iso_inds[forward_idx]))
				continue;

			const mass_t offset = fabs(peaks[forward_idx].mass + peaks[back_idx].mass - single_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
					
			pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs0.clear();
	forward_idx = -1;
	back_idx = num_peaks-1;
	while (forward_idx<back_idx)
	{
		forward_idx++;
		if (strict_iso_inds[forward_idx])
			continue;
			
		mass_t sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		while (back_idx>=0 && sum>max_double_sum)
		{
			back_idx--;
			if (back_idx<0)
				break;
			sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		}

		if (back_idx>=0 && sum > min_double_sum)
		{
			if (strict_iso_inds[back_idx])
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
			
			c2_pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	c2_pairs1.clear();
	const mass_t maxc21 = max_double_sum + 1.0;
	const mass_t minc21 = min_double_sum + 1.0;
	forward_idx = -1;
	back_idx = num_peaks-1;
	while (forward_idx<back_idx)
	{
		forward_idx++;
	
		mass_t sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		while (back_idx>=0 && sum>maxc21)
		{
			back_idx--;
			if (back_idx<0)
				break;
			sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		}

		if (back_idx>=0 && sum > minc21)
		{
			if (! (strict_iso_inds[back_idx] || strict_iso_inds[forward_idx]) )
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
			
			c2_pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs2.clear();
	const mass_t maxc22 = max_double_sum + 2.0;
	const mass_t minc22 = min_double_sum + 2.0;
	forward_idx = -1;
	back_idx = num_peaks-1;
	while (forward_idx<back_idx)
	{
		forward_idx++;
	
		mass_t sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		while (back_idx>=0 && sum>maxc22)
		{
			back_idx--;
			if (back_idx<0)
				break;
			sum = 2*peaks[forward_idx].mass + peaks[back_idx].mass;
		}

		if (back_idx>=0 && sum > minc22)
		{
			if (! (strict_iso_inds[back_idx] || strict_iso_inds[forward_idx]) )
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_sum);
			const float inten_sum = peaks[forward_idx].intensity + peaks[back_idx].intensity;
			
			c2_pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}

	// use the first 4 peaks
	stats.inten_strict_pairs0=0;
	stats.num_strict_pairs0 = pairs0.size();
	sort(pairs0.begin(),pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs0.size(); i++)
	{
//		stats.offset_strict_pairs0.push_back(pairs0[i].offset);
		stats.inten_strict_pairs0+=pairs0[i].inten_sum;
	}


	stats.inten_strict_pairs1=0;
	stats.num_strict_pairs1 = pairs1.size();
	sort(pairs1.begin(),pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs1.size(); i++)
	{
//		stats.offset_strict_pairs1.push_back(pairs1[i].offset);
		stats.inten_strict_pairs1+=pairs1[i].inten_sum;
	}

	stats.inten_strict_pairs2=0;
	stats.num_strict_pairs2 = pairs2.size();
	sort(pairs2.begin(),pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<pairs2.size(); i++)
	{
//		stats.offset_strict_pairs2.push_back(pairs2[i].offset);
		stats.inten_strict_pairs2+=pairs2[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs0=0;
	stats.c2_num_strict_pairs0 = c2_pairs0.size();
	sort(c2_pairs0.begin(),c2_pairs0.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs0.size(); i++)
	{
//		stats.c2_offset_strict_pairs0.push_back(c2_pairs0[i].offset);
		stats.c2_inten_strict_pairs0+=c2_pairs0[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs1=0;
	stats.c2_num_strict_pairs1 = c2_pairs1.size();
	sort(c2_pairs1.begin(),c2_pairs1.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs1.size(); i++)
	{
//		stats.c2_offset_strict_pairs1.push_back(c2_pairs1[i].offset);
		stats.c2_inten_strict_pairs1+=c2_pairs1[i].inten_sum;
	}
	
	stats.c2_inten_strict_pairs2=0;
	stats.c2_num_strict_pairs2 = c2_pairs2.size();
	sort(c2_pairs2.begin(), c2_pairs2.end(),cmp_offset_pair_inten);
	for (i=0; i<4 && i<c2_pairs2.size(); i++)
	{
//		stats.c2_offset_strict_pairs2.push_back(c2_pairs2[i].offset);
		stats.c2_inten_strict_pairs2+=c2_pairs2[i].inten_sum;
	}

}



void fill_rank_PMC_stats(int charge,
						  const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  mass_t minus_range, 
						  mass_t plus_range,
						  mass_t increment,
						  Config *config,
						  const BasicSpectrum& bs,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& iso_inds,
						  vector<PMCRankStats>& pmc_stats_vec)
{

	const mass_t tolerance = config->get_tolerance()*0.55;
	const int num_bins_per_Da = (int)(1.0/increment);
	const int num_minus_bins = (int)((-minus_range)*num_bins_per_Da)+1;
	const int num_plus_bins = (int)(plus_range*num_bins_per_Da)+1;
	const mass_t one_over_charge = 1.0/(mass_t)charge;

	const int total_num_bins = num_minus_bins + num_plus_bins;
	if (pmc_stats_vec.size() != total_num_bins)
		pmc_stats_vec.resize(total_num_bins);

	int i;
	for (i=0; i<num_minus_bins; i++)
	{	
		const mass_t delta = (i - num_minus_bins)*increment;

		calc_pmc_rank_stats_for_mass(bs.peaks,bs.num_peaks, single_charge_pair_sum+delta,
			tolerance, iso_levels, strong_inds, iso_inds, pmc_stats_vec[i]);

		pmc_stats_vec[i].m_over_z = (single_charge_pair_sum+delta+(charge-2)*MASS_PROTON)/charge;
	}

	const int start_bin_idx = i;

	for (i=0; i<num_plus_bins; i++)
	{
		const mass_t delta = i*increment;

		calc_pmc_rank_stats_for_mass(bs.peaks,bs.num_peaks, single_charge_pair_sum+delta,
			tolerance, iso_levels, strong_inds, iso_inds, pmc_stats_vec[start_bin_idx+i]);

		pmc_stats_vec[start_bin_idx+i].m_over_z = (single_charge_pair_sum+delta+(charge-2)*MASS_PROTON)/charge;
	}
}




void PMCRankStats::clear()
{
	m_over_z=0;

	rank_score = NEG_INF;

	num_frag_pairs=0;
	num_strong_frag_pairs=0;
	num_c2_frag_pairs=0;
	num_strong_c2_frag_pairs=0;
	num_h2o_loss_frag_pairs=0;

	inten_frag_pairs=0;
	inten_strong_pairs=0;
	inten_c2_pairs=0;
	inten_c2_strong_pairs=0;
	inten_h2o_loss_frag_pairs=0;
	itnen_h2o_loss_c2_frag_pairs=0;

	mean_offset_pairs=0;
	mean_offset_strong_pairs=0;
	mean_offset_c2_pairs=0;
	mean_offset_c2_strong_pairs=0;
	mean_offset_h2o_pairs=0;
	mean_offset_c2_h2o_pairs=0;

	ind_pairs_with_min_tol=false;			 
	ind_strong_pairs_with_min_tol=false;
	ind_c2_pairs_with_min_tol=false;
	ind_c2_strong_pairs_with_min_tol=false;
	log_dis_from_pairs_min_tol=0;			 
	log_dis_from_strong_pairs_min_tol=0;
	log_dis_from_c2_pairs_min_tol=0;		 
	log_dis_from_c2_strong_pairs_min_tol=0;

	offset_pairs_ordered_by_inten.clear();
	strong_offset_pairs_ordered_by_inten.clear();
	c2_offset_pairs_ordered_by_inten.clear();


	num_strict_pairs0=0; inten_strict_pairs0=0;
	num_strict_pairs1=0; inten_strict_pairs1=0;
	num_strict_pairs2=0; inten_strict_pairs2=0;
}




/**************************************************************************

  Fills in the RankBoost feature data
***************************************************************************/
void PMCSQS_Scorer::fill_RankBoost_smaples_with_PMC(
									const BasicSpectrum& bs,
									int charge,
									vector<RankBoostSample>& samples) const
{

	const int num_samples = curr_spec_rank_pmc_tables[charge].size();
	const int idx_skip = int((1.0/bin_increment)+0.00001);
	vector<int> idx_offsets;
	int i;

	idx_offsets.clear();
	idx_offsets.push_back(-2*idx_skip);
	idx_offsets.push_back(-1*idx_skip);
	idx_offsets.push_back(idx_skip);
	idx_offsets.push_back(2*idx_skip);

	if (samples.size() != num_samples)
		samples.resize(num_samples);

	for (i=0; i<num_samples; i++)
	{
		const PMCRankStats& stats = curr_spec_rank_pmc_tables[charge][i];
		RankBoostSample& sam = samples[i];

		const float inten_norm = 1.0/(curr_spec_total_intensity+1.0);
		int r_idx=0;
		const mass_t mz_offset = (stats.m_over_z - bs.ssf->m_over_z);

		sam.clear();
		sam.add_real_feature(r_idx++,mz_offset);

		if (stats.num_frag_pairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.num_frag_pairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.num_strong_frag_pairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

		if (stats.num_c2_frag_pairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.num_c2_frag_pairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.num_strong_c2_frag_pairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

			
	/*	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >4");

		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >4");*/

		sam.add_real_feature(r_idx++,stats.num_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_strong_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_c2_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_strong_c2_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs);
		sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs);

		sam.add_real_feature(r_idx++,stats.inten_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_c2_strong_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.inten_h2o_loss_frag_pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.itnen_h2o_loss_c2_frag_pairs * inten_norm);

		// averages of top k offsets

		float avg=0;
		int j;
		for (j =0; j<7 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;

		avg=0;
		for (j =0; j<7 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
		{
			avg += fabs(stats.c2_offset_pairs_ordered_by_inten[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;


		// offset data
	
		if (stats.mean_offset_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_pairs/(1.0+stats.num_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_strong_pairs/(1.0+stats.num_strong_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_pairs/(1.0+stats.num_c2_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_strong_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_strong_pairs/(1.0+stats.num_strong_c2_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_h2o_pairs/(1.0+stats.num_h2o_loss_frag_pairs));
		}
		else
			r_idx+=2;

		if (stats.mean_offset_c2_h2o_pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs);
			sam.add_real_feature(r_idx++,stats.mean_offset_c2_h2o_pairs/(1.0+stats.num_h2o_loss_c2_frag_pairs));
		}
		else
			r_idx+=2;

		// individual offsets
		for (j=0; j<5 && j<stats.offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

		for (j=0; j<5 && j<stats.c2_offset_pairs_ordered_by_inten.size(); j++)
			sam.add_real_feature(r_idx+j,stats.c2_offset_pairs_ordered_by_inten[j]);
		r_idx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.add_real_feature(r_idx++,stats.num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.inten_strict_pairs2 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs0);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs1);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.c2_num_strict_pairs2);
		sam.add_real_feature(r_idx++,stats.c2_inten_strict_pairs2 * inten_norm);

		// add comparative features to -2 -1 +1 +2 Da away
		for (j=0; j<idx_offsets.size(); j++)
		{
			const int other_idx = i + idx_offsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				r_idx+=12;
				continue;
			}

			const PMCRankStats& other = curr_spec_rank_pmc_tables[charge][other_idx];

			sam.add_real_feature(r_idx++,stats.num_frag_pairs - other.num_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_strong_frag_pairs - other.num_strong_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_c2_frag_pairs - other.num_c2_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_strong_c2_frag_pairs - other.num_strong_c2_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_frag_pairs - other.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,stats.num_h2o_loss_c2_frag_pairs - other.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(stats.inten_frag_pairs - other.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_strong_pairs - other.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_pairs - other.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_c2_strong_pairs - other.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.inten_h2o_loss_frag_pairs - other.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.itnen_h2o_loss_c2_frag_pairs - other.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}

		const int plus_idx = i + idx_skip;
		const int minus_idx = i-idx_skip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = curr_spec_rank_pmc_tables[charge][plus_idx];
			const PMCRankStats& minus = curr_spec_rank_pmc_tables[charge][minus_idx];

			sam.add_real_feature(r_idx++,plus.num_frag_pairs - minus.num_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_strong_frag_pairs - minus.num_strong_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_c2_frag_pairs - minus.num_c2_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_strong_c2_frag_pairs - minus.num_strong_c2_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_frag_pairs - minus.num_h2o_loss_frag_pairs);
			sam.add_real_feature(r_idx++,plus.num_h2o_loss_c2_frag_pairs - minus.num_h2o_loss_c2_frag_pairs);

			sam.add_real_feature(r_idx++,(plus.inten_frag_pairs - minus.inten_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_strong_pairs - minus.inten_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_pairs - minus.inten_c2_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_c2_strong_pairs - minus.inten_c2_strong_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.inten_h2o_loss_frag_pairs - minus.inten_h2o_loss_frag_pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.itnen_h2o_loss_c2_frag_pairs - minus.itnen_h2o_loss_c2_frag_pairs) * inten_norm);
		}
	}
}


void init_PMC_feature_names(vector<string>& names)
{
	names.clear();
	int i;

	names.push_back("OFFSET FROM MEASURED M/Z");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >=3");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >=3");
	
	names.push_back("# PAIRS");
	names.push_back("# STRONG PAIRS");
	names.push_back("# C2 PAIRS");
	names.push_back("# STRONG C2 PAIRS");
	names.push_back("# H2O PAIRS");
	names.push_back("# C2 H2O PAIRS");

	names.push_back("INTEN PAIRS");
	names.push_back("INTEN STRONG PAIRS");
	names.push_back("INTEN C2 PAIRS");
	names.push_back("INTEN STRONG C2 PAIRS");
	names.push_back("INTEN H2O PAIRS");
	names.push_back("INTEN C2 H2O PAIRS");

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP (STRONG %d)",i);
		names.push_back(name);
	}

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP C2 (STRONG %d)",i);
		names.push_back(name);
	}

	names.push_back("MEAN OFFSET PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET PAIRS");

	names.push_back("MEAN OFFSET STRONG PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG PAIRS");

	names.push_back("MEAN OFFSET C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 PAIRS");

	names.push_back("MEAN OFFSET STRONG C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG C2 PAIRS");

	names.push_back("MEAN OFFSET H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET H2O PAIRS");

	names.push_back("MEAN OFFSET C2 H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 H2O PAIRS");

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"C2 PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}



	names.push_back("NUM STRICT 0");
	names.push_back("INTEN STRICT 0");

	names.push_back("NUM STRICT 1");
	names.push_back("INTEN STRICT 1");

	names.push_back("NUM STRICT 2");
	names.push_back("INTEN STRICT 2");

	names.push_back("NUM C2 STRICT 0");
	names.push_back("INTEN C2 STRICT 0");

	names.push_back("NUM C2 STRICT 1");
	names.push_back("INTEN C2 STRICT 1");

	names.push_back("NUM C2 STRICT 2");
	names.push_back("INTEN C2 STRICT 2");

	// diff features with -2 -1 +1 +2
	const string dis_labels[]={"-2","-1","+1","+2"};
	for (i=0; i<4; i++)
	{
		const string prefix = "DIFF WITH "+dis_labels[i]+" ";

		names.push_back(prefix+"# PAIRS");
		names.push_back(prefix+"# STRONG PAIRS");
		names.push_back(prefix+"# C2 PAIRS");
		names.push_back(prefix+"# STRONG C2 PAIRS");
		names.push_back(prefix+"# H2O PAIRS");
		names.push_back(prefix+"# C2 H2O PAIRS");

		names.push_back(prefix+"INTEN PAIRS");
		names.push_back(prefix+"INTEN STRONG PAIRS");
		names.push_back(prefix+"INTEN C2 PAIRS");
		names.push_back(prefix+"INTEN STRONG C2 PAIRS");
		names.push_back(prefix+"INTEN H2O PAIRS");
		names.push_back(prefix+"INTEN C2 H2O PAIRS");
	}

	names.push_back("DIFF +1/-1 # PAIRS");
	names.push_back("DIFF +1/-1 # STRONG PAIRS");
	names.push_back("DIFF +1/-1 # C2 PAIRS");
	names.push_back("DIFF +1/-1 # STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 # H2O PAIRS");
	names.push_back("DIFF +1/-1 # C2 H2O PAIRS");

	names.push_back("DIFF +1/-1 INTEN PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN H2O PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 H2O PAIRS");
	cout << "Initialized: " << names.size() << " real feature names..." << endl;
}




void PMCSQS_Scorer::output_pmc_rank_results(const FileManager& fm,
											int charge,
											const vector<SingleSpectrumFile *>& test_ssfs) 
{
	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	vector<int> org_offset_counts, new_offset_counts;
	org_offset_counts.resize(201,0);
	new_offset_counts.resize(201,0);

	vector<mass_t> org_offsets;
	vector<mass_t> corr_offsets;

	org_offsets.clear();
	corr_offsets.clear();

	int i;
	for (i=0; i<test_ssfs.size(); i++)
	{
		SingleSpectrumFile* ssf = test_ssfs[i];
		BasicSpectrum bs;
	
		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

		init_for_current_spec(config,bs);
		calculate_curr_spec_pmc_values(bs, bin_increment);

		PmcSqsChargeRes res;
		find_best_mz_values_from_rank_model(bs, charge, config->get_pm_tolerance(),res);

		ssf->peptide.calc_mass(config);
		mass_t true_mz = (ssf->peptide.get_mass() + 18.01 + charge)/charge;

		org_offsets.push_back(true_mz - ssf->m_over_z);
		corr_offsets.push_back(true_mz - res.mz1);
	}

	mass_t m_org,sd_org,m_corr,sd_corr;
	calc_mean_sd(org_offsets,&m_org, &sd_org);
	calc_mean_sd(corr_offsets,&m_corr,&sd_corr);

	cout << "CHARGE: " << charge << endl;
	cout << "ORG:  mean " << m_org << " " << sd_org << endl;
	cout << "CORR: mean " << m_corr << " " << sd_corr << endl;

	for (i=0; i<org_offsets.size(); i++)
	{
		int org_idx = 100 + int(org_offsets[i] * 20);
		if (org_idx<0)
			org_idx = 0;
		if (org_idx>200)
			org_idx=200;
		org_offset_counts[org_idx]++;

		int new_idx = 100 + int(corr_offsets[i] * 20);
		if (new_idx<0)
			new_idx = 0;
		if (new_idx>200)
			new_idx=200;
		new_offset_counts[new_idx]++;
	}

	int cum_org=0;
	int cum_new=0;
	for (i=0; i<=200; i++)
	{

		if (org_offset_counts[i]==0 && new_offset_counts[i]==0)
			continue;
		
		cum_org+=org_offset_counts[i];
		cum_new+=new_offset_counts[i];
		cout << fixed << setprecision(3) << i*0.05 - 5.0 << "\t" <<
			org_offset_counts[i]/(float)org_offsets.size() << "\t" <<
			new_offset_counts[i]/(float)corr_offsets.size() << "\t" <<
			cum_org/(float)org_offsets.size() << "\t"<<
			cum_new/(float)corr_offsets.size() << endl;

	}


}


/**********************************************************************************
Finds the best m/z values for a given charge
If the pm_tolerance is low (less than 0.5) then the m/z vlaues need to reflect +-X Da
from the recorded m/z
***********************************************************************************/
void PMCSQS_Scorer::find_best_mz_values_from_rank_model(
										const BasicSpectrum& bs, 
										int charge,
										mass_t pm_tolerance,
										PmcSqsChargeRes& res)
{
	static vector<RankBoostSample> spec_samples;
	static vector<float> rank_scores;

	const mass_t spec_mz = bs.ssf->m_over_z;
	const mass_t allowed_mz_diff = pm_tolerance / charge;

	fill_RankBoost_smaples_with_PMC(bs, charge, spec_samples);

	if (rank_scores.size() != spec_samples.size())
		rank_scores.resize(spec_samples.size(),NEG_INF);

	const mass_t pm_with_19 = bs.ssf->m_over_z * charge - (charge + 1);
	const int size_idx = get_rank_model_size_idx(charge, pm_with_19);
	
	int best_idx=-1;
	float best_score=NEG_INF;

	if (charge>= pmc_rank_models.size() ||
		size_idx>= pmc_rank_models[charge].size() ||
		! pmc_rank_models[charge][size_idx])
	{
		//
	}
	else
	{
	//	cout << "spec: " << spec_mz << "  charge: " << charge << endl;
		int i;
		for (i=0; i<spec_samples.size(); i++)
		{
			// make sure the m/z is in an allowed range (i.e., it has a mass that can possibly be a +-X Da
			// shift from the true pm). Assume shift can be at most
			if (pm_tolerance<0.5)
			{
				const mass_t table_mz = curr_spec_rank_pmc_tables[charge][i].m_over_z;
				const mass_t one_over_charge = 1.0033 / (mass_t)charge; // ~mass of isotopic peak difference
				int d;
				for (d=-4; d<=4; d++)
				{
					const mass_t mz_diff = fabs(table_mz + d*one_over_charge - spec_mz);
					if (mz_diff<allowed_mz_diff)
						break;
				}

				if (d>4)
				{
					rank_scores[i] = NEG_INF;
					continue;
				}
			//	cout << "ok : " << table_mz << endl;
			}

			rank_scores[i]=pmc_rank_models[charge][size_idx]->calc_rank_score(spec_samples[i]);
			curr_spec_rank_pmc_tables[charge][i].rank_score = rank_scores[i];
			if (rank_scores[i]>best_score)
			{
				best_score=rank_scores[i];
				best_idx = i;
			}
		}
	}

	// no suitable models were found for this spectrum
	if (best_idx<0)
	{
		res.mz1 = bs.ssf->m_over_z;
		res.score1 = 10.0;
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
		return;
	}

	
	res.mz1 = curr_spec_rank_pmc_tables[charge][best_idx].m_over_z;
	res.score1 = best_score;

	// look for additional m/z
	int second_best_idx=-1;
	float second_best_score=NEG_INF;

	const mass_t mz_diff = curr_spec_rank_pmc_tables[charge][1].m_over_z - 
						   curr_spec_rank_pmc_tables[charge][0].m_over_z;

	const int idx_diff = (int)(0.45/(charge * mz_diff));

	int i;
	for (i=0; i<spec_samples.size(); i++)
	{
		if (rank_scores[i]>NEG_INF && fabs(float(i-best_idx))<idx_diff)
			continue;

		if (rank_scores[i]>second_best_score)
		{
			second_best_score=rank_scores[i];
			second_best_idx = i;
		}

	}
 
	if (second_best_idx>=0)
	{
		res.mz2 = curr_spec_rank_pmc_tables[charge][second_best_idx].m_over_z;
		res.score2 = second_best_score;
	} 
	else
	{
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
	}

//	const int size_idx = this->get_rank_model_size_idx(charge,res.mz1*charge-charge+1);
//	res.mz1 -= pmc_charge_mz_biases[charge][size_idx];
//	res.mz2 -= pmc_charge_mz_biases[charge][size_idx];

//	cout << charge << " ]\t" << res.mz1 << "\t" << res.prob1 << "\t" << res.mz2 << "\t" << res.prob2 << endl;


}




