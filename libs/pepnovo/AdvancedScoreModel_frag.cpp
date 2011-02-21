#include "AdvancedScoreModel.h"


bool StrongFragModel::write_model(ostream &os) const
{
	if (! ind_has_models)
		return false;

	os << model_frag_idx << " " << mirror1_idx << " " << mirror2_idx << " " <<
		parent1_idx << " " << parent2_idx << endl;

	os << model_frag_charge << " " << mirror1_charge << " " <<  mirror2_charge << " " << 
		parent1_charge << " " << parent2_charge << endl;

	os << setprecision(4) << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		cout << "Model for frag " << config->get_fragment(model_frag_idx).label << " had scaling problems: " << endl;
		cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		cout << "NOT writing the model, fix it!" << endl;
		exit(1);
	}

	inten_model.write_regression_model(os);
	no_inten_model.write_regression_model(os);

	return true;
}


bool StrongFragModel::read_model(istream& is, bool silent_ind)
{
	char buff[256];
	is.getline(buff,256);
	istringstream iss(buff);
	iss >> model_frag_idx >> mirror1_idx >> mirror2_idx >> parent1_idx >> parent2_idx;

	is.getline(buff,256);
	iss.str(buff);
	iss >> model_frag_charge >> mirror1_charge >> mirror2_charge >> parent1_charge >> parent2_charge;

	is.getline(buff,256);
	sscanf(buff,"%f %f",&inten_log_scaling_factor,&no_inten_log_scaling_factor);

	inten_model.read_regression_model(is);
	no_inten_model.read_regression_model(is);

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		cout << "Model for frag " << config->get_fragment(model_frag_idx).label << " had scaling problems: " << endl;
		cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		exit(1);
	}

	if (! inten_model.get_has_weights() || ! no_inten_model.get_has_weights())
	{
		
		cout << "Model for frag " << config->get_fragment(model_frag_idx).label << " has no weights!" << endl;
		exit(1);
	}

	ind_has_models = true;

	return true;
}


bool RegularFragModel::write_model(ostream& os) const
{
	if (! ind_has_models)
		return false;

	os << model_frag_idx << " " << model_frag_charge << endl;

	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		os <<  NEG_INF << " " << NEG_INF << endl;
	}
	else
		os << setprecision(4) << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;

	os << parent_idxs.size() << " " << parent_idx_with_same_charge_ori << endl;
	int i;
	for (i=0; i<parent_idxs.size() && i<num_parents; i++)
		os << parent_idxs[i] << " ";
	os << endl;

	inten_model.write_regression_model(os);
	no_inten_model.write_regression_model(os);

	return true;
}


bool RegularFragModel::read_model(istream& is, bool silent_ind)
{
	char buff[256];

	is.getline(buff,256);
	istringstream iss(buff);
	iss >> model_frag_idx >> model_frag_charge;
	
	is.getline(buff,256);
	sscanf(buff,"%f %f",&inten_log_scaling_factor,&no_inten_log_scaling_factor);

	is.getline(buff,256);
	iss.str(buff);
	iss >> num_parents >> parent_idx_with_same_charge_ori;

	is.getline(buff,256);
	iss.str(buff);
	int i;
	parent_idxs.resize(num_parents);
	for (i=0; i<num_parents; i++)
		iss >> parent_idxs[i];

	inten_model.read_regression_model(is);
	no_inten_model.read_regression_model(is);

	ind_has_models = true;

	// check if something was wrong with this fragment
	if (fabs(inten_log_scaling_factor)>5 || fabs(no_inten_log_scaling_factor)>5)
	{
		if (! silent_ind)
		{
			cout << "Model for frag " << model_frag_idx << " (" << 
				config->get_fragment(model_frag_idx).label << ") had scaling problems: " << endl;
			cout << inten_log_scaling_factor << " " << no_inten_log_scaling_factor << endl;
		}
		ind_has_models=false;
	}

	if (! inten_model.get_has_weights() || ! no_inten_model.get_has_weights())
		ind_has_models=false;

	return true;
}


string RegionalScoreModel::make_model_file_name(const char *name) const
{
	char dir_path[256];
	char model_name[64];
	strcpy(dir_path,config->get_resource_dir().c_str());
	strcat(dir_path,"/");
	strcat(dir_path,name);
	strcat(dir_path,"_SCORE/");
	sprintf(model_name,"%s_%d_%d_%d.txt",name,charge,size_idx,region_idx);
	string file_path = dir_path;
	file_path += model_name;

	return file_path;
}






bool RegionalScoreModel::write_regional_score_model(const char *name) const
{
	if (! was_initialized)
		return false;

	string file_path = make_model_file_name(name);
	ofstream os(file_path.c_str());

	if (! os.good() || ! os.is_open())
	{
		cout << "Error: couldn't write model file." << endl;
		cout << "Make sure the following path can be written: " << file_path << endl;
		exit(1);
	}

	int i;
	for (i=0; i<strong_models.size(); i++)
		if (! strong_models[i].write_model(os))
		{
			cout << "Error: no data exists for " << config->get_fragment(strong_models[i].model_frag_idx).label << endl;
			exit(1);
		}

	for (i=0; i<regular_models.size(); i++)
		if (! regular_models[i].write_model(os))
		{
			cout << "Warning: no data exists for " << config->get_fragment(regular_models[i].model_frag_idx).label << endl;
			exit(1);
		}

	os.close();

	return true;
}


bool RegionalScoreModel::read_regional_score_model(const char *name, bool silent_ind)
{
	string file_path = make_model_file_name(name);
	ifstream is(file_path.c_str());

	if (! is.good() || ! is.is_open())
	{
		is.close();
		return false;
	}

	int i;
	for (i=0; i<strong_models.size(); i++)
	{
		strong_models[i].read_model(is, silent_ind);
		if (! strong_models[i].ind_has_models)
			return false;
	}
	
	for (i=0; i<regular_models.size(); i++)
	{
		regular_models[i].read_model(is, silent_ind);
	}

	is.close();
	
	was_initialized = true;

	return true;
}

