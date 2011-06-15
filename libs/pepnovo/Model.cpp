#include "Model.h"
#include "FragmentSelection.h"


// reads a model and all relevant files
// the model files are assumed to be in the resource_dir
// all this model's files are assumed to have a name <model_name>_XXXXX.txt
// the main model file is <model_name>.txt
void Model::read_model(const char* name, bool silent_ind)
{
	char file[256];

	model_name = name;
	if (config.get_resource_dir().length()<2)
	{
		config.set_resource_dir("Models");
	}

	config.set_model_name(string(name));

	strcpy(file,config.get_resource_dir().c_str());
	strcat(file,"/");
	strcat(file,name); 
	strcat(file,".txt");   

	fstream fs(file,ios::in);
	if (! fs.good() )  
	{
		cout << "Error: couldn't open model file: " << file << endl;
		exit(1);
	}

	while (! fs.eof())
	{
		char buff[1024];
		fs.getline(buff,1024);
		if (fs.gcount()<4)
			continue;

		char arg[128];
		if (sscanf(buff,"#CONFIG_FILE %s",arg) == 1)
		{
			config.read_config(arg);
			config.set_model_name(string(model_name));
			continue;
		}

		if (! strncmp("#CONF",buff,5))
		{
			string path = config.get_resource_dir() + "/" + string(buff);
			config.parse_config_parameter((char *)path.c_str());
			continue;
		}

		if (sscanf(buff,"#BREAK_SCORE_MODEL %s",arg) ==1)
		{
			read_score_model(arg,silent_ind);
			continue;
		}

		if (sscanf(buff,"#EDGE_MODEL %s",arg) ==1)
		{
			edge_model.read_edge_models(&config,arg,silent_ind);
			continue;
		}

		if (sscanf(buff,"#SQS_MODEL %s",arg) == 1)
		{
			pmcsqs.read_sqs_models(&config,arg);
			continue;
		}
		
		if (sscanf(buff,"#PMCR_MODEL %s",arg) == 1)
		{
			pmcsqs.read_pmc_rank_models(&config,arg);
			continue;
		}

		if (sscanf(buff,"#COMP_ASSIGNER %s",arg) == 1)
		{
			comp_assigner.read_and_init_from_tables(&config,arg);
			continue;
		}

		if (sscanf(buff,"#AAP_MODEL %s",arg) == 1)
		{
			amino_acid_probs.read_amino_acid_prob_models(&config,arg);
			continue;
		}
	}

	// check if some of the defaults need to be changed
	if (config.get_max_edge_length() != 2)
		config.calc_aa_combo_masses();

}



// writes a model and all relevant files
// the model files are assumed to be in the resource_dir
// all this model's files are assumed to have a name <model_name>_XXXXX.txt
// the main model file is <model_name>.txt
void Model::write_model()
{
	string model_file;

	model_file = config.get_resource_dir() + "/" + model_name + ".txt";

	fstream os(model_file.c_str(),ios::out);
	if ( ! os.good())
	{
		cout << "Error writing model to " << model_file << endl;
		exit(1);
	}


	string config_file = config.get_resource_dir() + "/" + model_name + "_config.txt";
	config.set_config_file(config_file);
	config.set_model_name(model_name);
	os << "#CONFIG_FILE " << model_name + "_config.txt" << endl;
	config.write_config();


	if (pmcsqs.get_ind_initialized_pmcr())
	{
		os << "#PMCR_MODEL " << model_name + "_PMCR.txt" << endl;
		string path = config.get_resource_dir() + "/" + model_name + "_PMCR.txt";
		pmcsqs.write_pmc_rank_models(path.c_str());
	}

	if (pmcsqs.get_ind_initialized_sqs())
	{
		os << "#SQS_MODEL " << model_name + "_SQS.txt" << endl;
		string path = config.get_resource_dir() + "/" + model_name + "_SQS.txt";
		pmcsqs.write_sqs_models(path.c_str());
	}

	if (comp_assigner.get_ind_was_initialized())
	{
		os << "#COMP_ASSIGNER " << comp_assigner.get_model_name() << endl;
	}

	os << "#BREAK_SCORE_MODEL " << model_name << endl;

	write_score_model(model_name.c_str());

	os << "#EDGE_MODEL " << model_name << endl;

	edge_model.write_edge_models(model_name.c_str());

	if (amino_acid_probs.get_ind_initialized())
	{
		os << "#AAP_MODEL " << model_name << endl;
		string path = config.get_resource_dir() + "/" + model_name + "_AAP.txt";
		amino_acid_probs.write_amino_acid_prob_models(path.c_str());
	}
}


/*************************************************************************
This function performs the entire training process of the model
Allows for training in stages, gives better output and checks that
previous stages are intialized
**************************************************************************/
void Model::train_model_in_stages(
			const char *name, 
			const FileManager& fm, 
			mass_t initial_tolerance, 
			int start_stage, 
			int end_stage,
			int specific_charge, 
			int specific_size, 
			int specific_region,
			char *neg_sqs_list)
{
	if (end_stage>1000)
		end_stage = 20;
	stages_intialized.resize(end_stage,false);
	
	model_name = name;
	config.set_model_name(string(name));


	cout << endl << "STAGE 0: Partitioning according to size/charge " << endl;
	cout <<         "**********************************************" <<endl;
	if (start_stage>0)
	{
		cout << endl << "Already done." << endl;	
	}
	else
	{
		cout << endl;
		int charge;
		for (charge = fm.get_min_charge(); charge<= fm.get_max_charge(); charge++)
		{
			vector<mass_t> spectra_masses;
			FileSet fs;
			fs.select_all_files(fm);
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
			int i;
			for (i=0; i<all_ssf.size(); i++)
				if (all_ssf[i]->charge == charge)
					spectra_masses.push_back(all_ssf[i]->org_pm_with_19);

			config.set_size_thresholds_according_to_set_of_masses(charge,spectra_masses);
		}
	}
	cout << endl << "Using following thresholds:" << endl;
	config.print_size_thresholds();


	cout << endl << "STAGE 1: Select Fragment types" << endl;
	cout <<         "******************************" <<endl;
	if (start_stage>1)
	{
		cout << endl << "Already done." << endl;	
	}
	else
	{
		config.set_tolerances(initial_tolerance);
		cout << endl;
		select_fragments(name,fm,15,0.01);
		config.set_all_regional_fragment_relationships();
	}
	cout << endl << "Fragments being used:" << endl;
	config.print_all_fragments();

	cout << endl << "STAGE 2: calculating fragment and PM tolerances" << endl;
	cout <<         "***********************************************" <<endl;
	if (start_stage>2)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		int c;
		for (c=0; c<config.get_max_charge_for_size(); c++)
			if (config.get_size_thresholds()[c].size()>0)
				config.select_strong_fragments(c,0.5,3);
		
		cout << "Calculating precursor mass tolerance..." << endl;
		mass_t pm_tol = calc_parent_mass_tolerance_distribution(this, fm, 0.95);

		cout << "Calculating fragment mass tolerance..." << endl;
		mass_t tol    = calc_tolerance_distribution(this, fm , initial_tolerance*1.2,0.96);

		config.set_pm_tolerance(pm_tol);

		if (pm_tol <0.000001)
		{
			pm_tol = tol;
		}

		if (pm_tol<tol)
		{
			config.set_tolerance(tol+pm_tol);
		}
		else
			config.set_tolerance(tol);
	}
	cout << endl << "PM tolerance " << fixed << setprecision(4) << config.get_pm_tolerance() << endl;
	cout << "Need to correct PM: " << config.get_need_to_estimate_pm() << endl;
	cout << "Fragment tolerance  " << config.get_tolerance() << endl;

	
//	config.print_all_regional_fragment_relationships();
	
	cout << endl << "STAGE 3: Train breakage score models" << endl;
	cout <<         "************************************" <<endl;
	cout << endl;
	if (start_stage>3)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only Specified model  " << specific_charge << " " << 
					specific_size << " " << specific_region << endl << endl;

		this->train_score_model(name,fm,specific_charge, specific_size, specific_region);	
	}

	if (end_stage<=3)
	{
		write_model();
		exit(0);
	}



	cout << endl << "STAGE 4: Train SQS models" << endl;
	cout <<         "*************************" << endl << endl;
	if (start_stage>4)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only specified charge " <<  specific_charge << endl << endl;

		vector< vector<float> > weights;
		int max_c = 4;
		if (fm.get_max_charge()+1>max_c)
			max_c = fm.get_max_charge()+1;
		weights.resize(max_c);
		
		int i;
		for (i=1; i<max_c; i++)
			weights[i].resize(3,0);
	
		weights[1][0] = 0.1; weights[1][1] = 0.1;  weights[1][2] = 0.4;
		weights[2][0] = 0.6; weights[2][1] = 0.75; weights[2][2] = 0.5;
		weights[3][0] = 0.3; weights[3][1] = 0.15; weights[3][2] = 0.1;
		for (i=4; i<max_c; i++)
			weights[i]=weights[3];

		train_sqs(fm,neg_sqs_list,specific_charge,&weights);
	}

	if (end_stage<=4)
	{
		write_model();
		exit(0);
	}

	cout << endl << "STAGE 5: Train PMCR models" << endl;
	cout <<         "**************************" << endl << endl;
	if (start_stage>5)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only specified charge " <<  specific_charge << endl << endl;


		train_pmc_rank_models(fm,specific_charge);
	} 


	if (end_stage<=5)
	{
		write_model();
		exit(0);
	}

	cout << endl << "STAGE 6: Train edge models" << endl;
	cout <<         "**************************" << endl << endl;
	if (start_stage>6)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only specified charge " <<  specific_charge << endl << endl;


		edge_model.train_all_edge_models(fm,this,specific_charge);
	}

	cout << endl << "STAGE 7: Train Amino Acid models" << endl;
	cout <<         "********************************" << endl << endl;
	if (start_stage>7)
	{
		cout << endl << "Already done." << endl;
	}
	else
	{
		if (specific_charge>0)
			cout << "+++ Only specified charge " <<  specific_charge << endl << endl;

		amino_acid_probs.train_amino_acid_prob_models(fm,this,specific_charge,specific_size);
	} 


	if (end_stage<=7)
	{
		write_model();
		exit(0);
	}


	exit(0);
}





/******************************************************************************
	This model selects the fragment types that will take part in the models.
	The fragments are selected according to the offset frequency function.
	If there isn't a suffcient number of spectra from the desired charge,
	the fragment selection is skipped.
*******************************************************************************/
bool Model::select_fragments(const char *name, 
							 const FileManager& fm, 
							 int   max_num_frags, 
							 float min_prob)
{
	FragmentTypeSet fragment_types;

	int c;
	cout << "Training set consists of:" << endl;
	for (c=fm.get_min_charge(); c<=fm.get_max_charge(); c++)
		cout << "Charge " << c <<"  " << fm.get_num_spectra(c) << " spectra."<< endl;
	cout<<endl;

	// select potential fragment type using the fragment offset test
	select_frags_using_frag_offset_counts(fm,&config, fragment_types, min_prob);

	// add these fragments to the existing set
	config.add_fragment_types(fragment_types);

	for (c=1; c<=fm.get_max_charge(); c++)
	{
	
		// check that there is a minimal number of files...
		int num_charge_spectra = fm.get_num_spectra(c);

		if (num_charge_spectra<MINIMAL_NUMBER_SPECTRA_FOR_FRAGMENT_SELECTION)
			continue;

		config.init_regional_fragment_set_defaults(0,c);

		select_regional_fragments(fm,&config,c,true);

		config.select_fragments_in_sets(1.5,max_num_frags);

		// select strong, combos...
		int max_num_combos = max_num_frags > 0 ? max_num_frags : 2;
		if (max_num_combos>3)
			max_num_combos = 3;
		
		config.select_strong_fragments(c,0.5,3);

		select_frag_combos(fm,&config,c,max_num_combos);
	}

	string fragments_file = config.get_resource_dir() + "/" + string(name) + "_fragments.txt";
	ofstream os(fragments_file.c_str(),ios::out);
	config.print_fragments(os);
	config.set_fragments_file(fragments_file);
	os.close();

	string regional_fragment_sets_file = config.get_resource_dir() + "/" + string(name) + "_fragment_sets.txt";
	os.open(regional_fragment_sets_file.c_str(),ios::out);
	config.print_regional_fragment_sets(os);
	config.set_regional_fragment_sets_file(regional_fragment_sets_file);
	os.close();

	return true;
}




// determines the tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t calc_tolerance_distribution(Model *model, 
								   const FileManager& fm, 
								   mass_t max_tolerance,
								   float cutoff_prob)
{
	FileSet fs;
	Config *config = model->get_config();
	FragmentTypeSet frags;
	vector<string> file_list;

	fs.select_all_files(fm);

	vector<int> test_frag_idxs;
	test_frag_idxs.clear();
	
	if (config->get_strong_type1_idx()>=0)
		test_frag_idxs.push_back(config->get_strong_type1_idx());

	if (config->get_strong_type2_idx()>=0)
		test_frag_idxs.push_back(config->get_strong_type2_idx());

	if (test_frag_idxs.size()==0)
	{
		cout << endl <<"Warning: no strong fragments selected, using maximal tolerance!!!" << endl << endl;
		return max_tolerance;
	}

	mass_t tol_increment = max_tolerance * 0.05;
	vector<float> offset_bins;
	vector<float> offsets;

	offset_bins.resize(41,0);

	int total_frag_count=0;
	while(1)
	{
		Spectrum s;
		vector<mass_t> break_masses;
		mass_t true_mass_with_19,true_mass;
	
		if (! fs.get_next_spectrum(fm,config,&s))
			break;
		
		s.init_spectrum();
		s.get_peptide().calc_expected_breakage_masses(config,break_masses);
		true_mass=s.get_peptide().get_mass();
		true_mass_with_19 =  true_mass + MASS_OHHH;

		if (break_masses.size()<3)
			continue;		
	
		// loop on fragments first, so high count fragments get precedence over
		// low count fragments that are actually due to b/y ions of previous or
		// next amino acids
		int f;
		for (f=0; f<test_frag_idxs.size(); f++)
		{
			const FragmentType& frag = config->get_fragment(test_frag_idxs[f]);
			int b;

			
			for (b=1; b<break_masses.size()-1; b++)
			{
				mass_t break_mass = break_masses[b];

				const mass_t exp_mass = frag.calc_expected_mass(break_mass,true_mass_with_19);
				const int p_idx = s.get_max_inten_peak(exp_mass,max_tolerance);

				if (p_idx>=0)
				{
					
					total_frag_count++;
					mass_t offset =  s.get_peak_mass(p_idx) - exp_mass;

					int bin_idx = 20 + (int)((offset / max_tolerance)*20);
					if (bin_idx<0)
						bin_idx=0;
					if (bin_idx>40)
						bin_idx=40;

					offset_bins[bin_idx]++;
					offsets.push_back(offset);
				}
			}
		}
	}

	int i;
	cout << "bin histogram: " << endl;
	for (i=0; i<=40; i++)
		cout << setprecision(4) << (20-i)*tol_increment << " " << 
			    offset_bins[i]/total_frag_count << endl;

	// find the offset that keeps the desired proportion of fragments
	sort(offsets.begin(),offsets.end());
	int count=0;
	int target_count = (int)((1.0 - cutoff_prob)*total_frag_count);
	int left_idx=0;
	int right_idx=offsets.size()-1;
	mass_t cutoff_offset=-1;
	while (count<target_count)
	{
		if (fabs(offsets[left_idx])>offsets[right_idx])
		{
			left_idx++;
		}
		else
			right_idx--;

		if (++count == target_count)
		{
			if (fabs(offsets[left_idx])>fabs(offsets[right_idx]))
			{
				cutoff_offset = fabs(offsets[left_idx]); 
			}
			else
				cutoff_offset = fabs(offsets[right_idx]);

			break;
		}
	}

	cout << "offset for " << cutoff_prob << " is " << cutoff_offset << endl;
	return cutoff_offset;
}


// determines the parent mass tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t calc_parent_mass_tolerance_distribution(Model *model,  const FileManager& fm, 
											   float cutoff_prob)
{
	FileSet fs;
	Config *config = model->get_config();
	FragmentTypeSet frags;
	
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	vector<float> offsets;

	int total_frag_count=0;
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		SingleSpectrumFile *ssf = all_ssf[i];
		if (ssf->peptide.get_num_aas()<3)
			continue;

		total_frag_count++;

		ssf->pm_with_19 = ssf->charge * ssf->m_over_z - MASS_PROTON * (ssf->charge -1 );
		mass_t offset =  ssf->pm_with_19  - ssf->peptide.get_mass_with_19();
	//	cout << setprecision(3) << offset << " " << ssf->charge << " " << ssf->peptide.as_string(config) << endl;
		offsets.push_back(offset);
	}



	// find the offset that keeps the desired proportion of fragments
	sort(offsets.begin(),offsets.end());
	int count=0;
	int target_count = (int)((1.0 - cutoff_prob)*total_frag_count);
	int left_idx=0;
	int right_idx=offsets.size()-1;
	mass_t cutoff_offset=-1;
	while (count<target_count)
	{
		if (fabs(offsets[left_idx])>offsets[right_idx])
		{
			left_idx++;
		}
		else
			right_idx--;

		if (++count == target_count)
		{
			if (fabs(offsets[left_idx])>fabs(offsets[right_idx]))
			{
				cutoff_offset = fabs(offsets[left_idx]); 
			}
			else
				cutoff_offset = fabs(offsets[right_idx]);

			break;
		}
	}

	cout << "Parent mass offset for " << setprecision(4) << cutoff_prob << " is " << cutoff_offset << endl;
	return cutoff_offset;
}







