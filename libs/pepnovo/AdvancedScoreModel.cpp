#include "AdvancedScoreModel.h"
#include "DeNovoRankScore.h"


void AdvancedScoreModel::read_score_model(const char *name, bool silent_ind)
{
	rank_models[0]=NULL;
	rank_models[1]=NULL;
	rank_models[2]=NULL;

	int i;
	for (i=0; i<10; i++)
		rank_tag_models[i]=NULL;

	// resize regional breakage score models according to regional fragment sets
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config.get_regional_fragment_sets();

	int c;
	regional_breakage_score_models.resize(all_rfs.size());
	for (c=0; c<all_rfs.size(); c++)
	{
		regional_breakage_score_models[c].resize(all_rfs[c].size());
		int s;
		for (s=0; s<all_rfs[c].size(); s++)
		{
			regional_breakage_score_models[c][s].resize(all_rfs[c][s].size());
			int r;
			for (r=0; r<regional_breakage_score_models[c][s].size(); r++)
				if (! regional_breakage_score_models[c][s][r].was_initialized)
					regional_breakage_score_models[c][s][r].init(&config,c,s,r);
		}
	}
	
	int num_models_read=0;
	for (c=0; c<regional_breakage_score_models.size(); c++)
	{
		int s;
		for (s=0; s<regional_breakage_score_models[c].size(); s++)
		{
			int r;
			for (r=0; r<regional_breakage_score_models[c][s].size(); r++)
			{
				if (regional_breakage_score_models[c][s][r].read_regional_score_model(name,silent_ind))
				{
					if (! silent_ind)
						cout << "Read breakage model: " << c << " " << s << " " << r << endl;
					num_models_read++;
				}
			}
		}
	}

	if (! silent_ind)
		cout << "Read " << num_models_read << " regional breakage models.." << endl;

	// reaf prm_normalize model file if it exists
	read_prm_normalizer_values();

}
	
void AdvancedScoreModel::write_score_model(const char *name) const
{
	int num_models_read=0;
	int c;
	for (c=0; c<regional_breakage_score_models.size(); c++)
	{
		int s;
		for (s=0; s<regional_breakage_score_models[c].size(); s++)
		{
			int r;
			for (r=0; r<regional_breakage_score_models[c][s].size(); r++)
			{
				if (regional_breakage_score_models[c][s][r].was_initialized)
					regional_breakage_score_models[c][s][r].write_regional_score_model(name);
				
			}
		}
	}
}


void AdvancedScoreModel::read_rank_models(const char *name, bool silent_ind)
{
	const string suffixes[]={"DB","DNVPART","DNVCOMP"};

	const int num_models = sizeof(rank_models)/sizeof(void *);

	int i;
	for (i=0; i<num_models; i++)
	{
		const string suffix = suffixes[i];
		DeNovoRankScorer *rank_model = (DeNovoRankScorer *)rank_models[i];
	
		string rank_name = string(name) + "_" + suffix;
		string path = config.get_resource_dir() + "/" + rank_name + "/" + suffix + "_rank_model.txt";

		ifstream fs(path.c_str());
		if (! fs.is_open() || ! fs.good())
		{
			if (! silent_ind)
				cout << "No " << path << endl;
			continue;
		}

		fs.close();
		rank_model = new DeNovoRankScorer;
		rank_model->set_type(i);
		rank_model->set_model(this);
		rank_model->read_denovo_rank_scorer_model(path.c_str(),suffix,silent_ind);
		rank_models[i] = (void *)rank_model;
	
		if (! silent_ind)
			cout << "Read " << path << endl;

	}

	// read tag models
	for (i=3; i<10; i++)
	{
		DeNovoRankScorer *rank_model = (DeNovoRankScorer *)rank_tag_models[i];
	
		
		ostringstream oss;
		oss << "TAG" << i;
		string suffix = oss.str();
		string rank_name = string(name) + "_" + suffix;
		string path = config.get_resource_dir() + "/" + rank_name + "/" + suffix + "_rank_model.txt";

		ifstream fs(path.c_str());
		if (! fs.is_open() || ! fs.good())
		{
			if (! silent_ind)
				cout << "No " << path << endl;
			continue;
		}

		fs.close();
		rank_model = new DeNovoRankScorer;
		rank_model->set_type(3);
		rank_model->set_model(this);
		rank_model->read_denovo_rank_scorer_model(path.c_str(),suffix,silent_ind);
		rank_tag_models[i] = (void *)rank_model;
	
		if (! silent_ind)
			cout << "Read " << path << endl;

	}

//	cout << endl;
}

// simple Dancik-style score
void AdvancedScoreModel::score_breakage(Spectrum *spec, Breakage *breakage, bool verbose) const
{
	const RegionalScoreModel& breakage_model = regional_breakage_score_models[breakage->parent_charge]
														  [breakage->parent_size_idx][breakage->region_idx];

	if (0 || breakage_model.has_all_breakage_models)
	{
	}
	else
	{
		score_t breakage_score=0;
		int i;
		for (i=0; i<breakage_model.frag_type_idxs.size(); i++)
		{
			const int frag_idx = breakage_model.frag_type_idxs[i];
			if (breakage->is_frag_type_visible(frag_idx))
			{
				if (breakage->get_position_of_frag_idx(frag_idx)>=0)
				{
					breakage_score += breakage_model.frag_inten_scores[frag_idx];
				}
				else
					breakage_score += breakage_model.frag_no_inten_scores[frag_idx];
			}
		}
		breakage->score = breakage_score;
	}
}
	


void AdvancedScoreModel::score_graph_edges(PrmGraph& prm) const
{
	edge_model.score_graph_edges(prm);
}



int AdvancedScoreModel::get_max_score_model_charge() const
{
	return regional_breakage_score_models.size();
}


void AdvancedScoreModel::train_score_model(const char *name, const FileManager& fm, 
						int charge, int size_idx, int region_idx)
{
	// resize regional breakage score models according to regional fragment sets
	const vector< vector< vector< RegionalFragments > > >& all_rfs = config.get_regional_fragment_sets();

	int c;
	regional_breakage_score_models.resize(all_rfs.size());
	for (c=0; c<all_rfs.size(); c++)
	{
		regional_breakage_score_models[c].resize(all_rfs[c].size());
		int s;
		for (s=0; s<all_rfs[c].size(); s++)
		{
			regional_breakage_score_models[c][s].resize(all_rfs[c][s].size());
			int r;
			for (r=0; r<regional_breakage_score_models[c][s].size(); r++)
				if (! regional_breakage_score_models[c][s][r].was_initialized)
					regional_breakage_score_models[c][s][r].init(&config,c,s,r);
		}
	}


	// train models
	for (c=0; c<regional_breakage_score_models.size(); c++)
	{
		if (regional_breakage_score_models.size() == 0 || (charge>0 && charge != c))
			continue;

		int s;
		for (s=0; s<regional_breakage_score_models[c].size(); s++)
		{
			if (size_idx>=0 && s != size_idx)
				continue;

			int r;
			for (r=0; r<regional_breakage_score_models[c][s].size(); r++)
			{
				if (region_idx>=0 && r != region_idx)
					continue;
				
				regional_breakage_score_models[c][s][r].train_regional_score_model(this,name,fm);
			}
		}
	}

	// train PRM normalizer values
	cout << endl << "Training PRM normalizer vlaues..." << endl;

	learn_prm_normalizer_values(fm);
}



/***********************************************************************
Creates a score table for each node that has an enry for each possible
entrance and exit combination of amino aicds.
************************************************************************/
void AdvancedScoreModel::score_all_node_combos(PrmGraph *prm) const
{
	const vector<int>& org_aas = config.get_org_aa();
	const vector<MultiEdge>& multi_edges = prm->get_multi_edges();
	const int num_nodes = prm->get_num_nodes();
	const mass_t mid_mass = prm->get_pm_with_19() * 0.5;
	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = prm->get_non_const_node(i);
		const RegionalScoreModel& score_model = 
			regional_breakage_score_models[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];

		vector<BreakageInfo> infos;

		BreakageInfo double_gap_info;
		prm->fill_breakage_info(this,&double_gap_info,i,NEG_INF,NEG_INF,NEG_INF,NEG_INF);
		infos.push_back(double_gap_info);

		int j;
		for (j=0; j<node.in_edge_idxs.size(); j++)
		{
			const int in_edge_idx = node.in_edge_idxs[j];
			const MultiEdge& in_edge = multi_edges[in_edge_idx];
			const int num_in_variants = in_edge.get_num_variants();

			int k;
			for (k=0; k<num_in_variants; k++)
			{
				BreakageInfo c_gap_info;
				prm->fill_breakage_info(this,&c_gap_info,i,in_edge_idx,k,NEG_INF,NEG_INF);
				infos.push_back(c_gap_info);
			}
		}

		for (j=0; j<node.out_edge_idxs.size(); j++)
		{
			const int out_edge_idx = node.out_edge_idxs[j];
			const MultiEdge& out_edge = multi_edges[out_edge_idx];
			const int num_out_variants = out_edge.get_num_variants();

			int k;
			for (k=0; k<num_out_variants; k++)
			{
				BreakageInfo n_gap_info;
				prm->fill_breakage_info(this,&n_gap_info,i,NEG_INF,NEG_INF,out_edge_idx,k);
				infos.push_back(n_gap_info);
			}
		}
		
		for (j=0; j<node.in_edge_idxs.size(); j++)
		{
			const int in_edge_idx = node.in_edge_idxs[j];
			const MultiEdge& in_edge = multi_edges[in_edge_idx];
			const int num_in_variants = in_edge.get_num_variants();

			int k;
			for (k=0; k<num_in_variants; k++)
			{
				BreakageInfo c_gap_info;
				prm->fill_breakage_info(this,&c_gap_info,i,in_edge_idx,k,NEG_INF,NEG_INF);
				infos.push_back(c_gap_info);

				int a;
				for (a=0; a<node.out_edge_idxs.size(); a++)
				{
					const int out_edge_idx = node.out_edge_idxs[a];
					const MultiEdge& out_edge = multi_edges[out_edge_idx];
					const int num_out_variants = out_edge.get_num_variants();

					int b;
					for (b=0; b<num_out_variants; b++)
					{
						BreakageInfo info;
						prm->fill_breakage_info(this,&info,i,in_edge_idx,k,out_edge_idx,b);
						infos.push_back(info);
					}
				}
			}
		}

	
		node.score_combos.clear();

		for (j=0; j<infos.size(); j++)
			infos[j].score = score_model.score_a_single_breakage_combo(prm, node, &node.breakage,infos[j]);

		for (j=0; j<infos.size(); j++)
			node.score_combos[ScoreComboLoc(infos[j])]=infos[j].score;

		score_t max_score = NEG_INF;
	
		if (i>0 && i<num_nodes-1)
		{
			if (node.mass>mid_mass)
			{
				int j;
				for (j=0; j<infos.size(); j++)
					if (infos[j].n_aa >= Ala && infos[j].c_aa <= Gap && infos[j].score>max_score)
						max_score=infos[j].score;
			}
			else
			{
				int j;
				for (j=0; j<infos.size(); j++)
					if (infos[j].n_aa <= Gap && infos[j].c_aa >= Ala && infos[j].score>max_score)
						max_score=infos[j].score;
			}
		}
		else
			max_score = infos[0].score;

		node.score = max_score; 
		node.breakage.score = node.score;
	}
	prm->set_has_node_combo_scores(true);
}



void AdvancedScoreModel::score_peptide_node_combos(PrmGraph *prm, const Peptide& peptide ) const
{
	const vector<int>& org_aas		= config.get_org_aa();
	const vector<mass_t>& aa2mass	= config.get_aa2mass();
	const vector<MultiEdge>& multi_edges = prm->get_multi_edges();
	const int num_nodes		   = prm->get_num_nodes();
	const vector<int>& pep_aas = peptide.get_amino_acids();
	const int num_pep_aas = pep_aas.size();

	mass_t p_mass=0;
	int aa_idx=0;

	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = prm->get_non_const_node(i);
		const RegionalScoreModel& score_model = 
			regional_breakage_score_models[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];

		int in_edge_idx=NEG_INF,  in_edge_variant=NEG_INF;
		int out_edge_idx=NEG_INF, out_edge_variant=NEG_INF;

	//	cout << "N: " << node.mass << endl;
		while (aa_idx<pep_aas.size() && fabs(p_mass-node.mass)>0.1)
		{
			p_mass += aa2mass[pep_aas[aa_idx]];
			aa_idx++;
	//		cout << aa_idx << "\t" << p_mass << endl;
		}
		
		if (aa_idx == pep_aas.size() && i != num_nodes-1)
		{
			int j;
			for (j=0; j<num_nodes; j++)
				cout << j << "\t" << prm->get_node(j).mass << endl;

			cout << endl << "PEP:" << endl;
			vector<mass_t> exp_masses;
			peptide.calc_expected_breakage_masses((Config *)&config,exp_masses);
			for (j=0; j<exp_masses.size(); j++)
				cout << j << "\t" << exp_masses[j] << endl;

			cout << "Error: mismatch between nodes and peptide!" << endl;
			exit(1);
		}

		if (node.in_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[edge_idx];
				const int num_aa = in_edge.num_aa;

				if (num_aa>aa_idx)
					continue;

				const int var_idx = in_edge.get_variant_idx(num_aa,&pep_aas[aa_idx-num_aa]);
				if (var_idx<0)
					continue;

				in_edge_idx = edge_idx;
				in_edge_variant = var_idx;
				break;
			}
		}

		if (node.out_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				const MultiEdge& out_edge = multi_edges[edge_idx];
				const int num_aa = out_edge.num_aa;

				if (num_aa + aa_idx >num_pep_aas)
					continue;

				const int var_idx = out_edge.get_variant_idx(num_aa,&pep_aas[aa_idx]);
				if (var_idx<0)
					continue;

				out_edge_idx = edge_idx;
				out_edge_variant = var_idx;
				break;
			}
		}


		BreakageInfo info;
		prm->fill_breakage_info(this,&info,i,in_edge_idx,in_edge_variant,out_edge_idx,out_edge_variant);
	
		node.score_combos.clear();

	//	cout << in_edge_idx << " " << in_edge_variant << " " << out_edge_idx << " " << out_edge_variant << "\t";

		info.score = score_model.score_a_single_breakage_combo(prm, node, &node.breakage, info);
	
		node.score_combos[ScoreComboLoc(info)]=info.score;
		node.score = info.score; 
		node.breakage.score = node.score;

	//	cout << node.score << endl;
	}
	prm->set_has_node_combo_scores(true);
}



/***************************************************************************
Scores the Gap-Gap combination for all nodes and sets it as the node score
****************************************************************************/
void AdvancedScoreModel::initial_combos_score(PrmGraph *prm) const
{
	const vector<int>& org_aas = config.get_org_aa();
	const vector<MultiEdge>& multi_edges = prm->get_multi_edges();
	const int num_nodes = prm->get_num_nodes();
	const mass_t mid_mass = prm->get_pm_with_19()*0.5;
	
	int i;
	for (i=0; i<num_nodes; i++)
	{
		Node& node = prm->get_non_const_node(i);
		const RegionalScoreModel& score_model = 
			regional_breakage_score_models[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];

		vector<BreakageInfo> infos;

		BreakageInfo double_gap_info;
		prm->fill_breakage_info(this,&double_gap_info,i,NEG_INF,NEG_INF,NEG_INF,NEG_INF);
		infos.push_back(double_gap_info);

		// fill all combos for a digest node
		if (node.type == NODE_DIGEST)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int in_edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[in_edge_idx];
				const int num_in_variants = in_edge.get_num_variants();

				int k;
				for (k=0; k<num_in_variants; k++)
				{
					int a;
					for (a=0; a<node.out_edge_idxs.size(); a++)
					{
						const int out_edge_idx = node.out_edge_idxs[a];
						const MultiEdge& out_edge = multi_edges[out_edge_idx];
						const int num_out_variants = out_edge.get_num_variants();

						int b;
						for (b=0; b<num_out_variants; b++)
						{
							BreakageInfo info;
							prm->fill_breakage_info(this,&info,i,in_edge_idx,k,out_edge_idx,b);
							if (info.connects_to_C_term || info.connects_to_N_term)
								infos.push_back(info);
						}
					}
				}
			}
		}
		else
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int in_edge_idx = node.in_edge_idxs[j];
				const MultiEdge& in_edge = multi_edges[in_edge_idx];
				const int num_in_variants = in_edge.get_num_variants();

				int k;
				for (k=0; k<num_in_variants; k++)
				{
					BreakageInfo c_gap_info;
					prm->fill_breakage_info(this,&c_gap_info,i,in_edge_idx,k,NEG_INF,NEG_INF);
					infos.push_back(c_gap_info);
				}
			}
		
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int out_edge_idx = node.out_edge_idxs[j];
				const MultiEdge& out_edge = multi_edges[out_edge_idx];
				const int num_out_variants = out_edge.get_num_variants();

				int k;
				for (k=0; k<num_out_variants; k++)
				{
					BreakageInfo n_gap_info;
					prm->fill_breakage_info(this,&n_gap_info,i,NEG_INF,NEG_INF,out_edge_idx,k);
					infos.push_back(n_gap_info);
				}
			}
		}

		node.score_combos.clear();
		int j;
		for (j=0; j<infos.size(); j++)
			infos[j].score = score_model.score_a_single_breakage_combo(prm, node, &node.breakage, infos[j]);

		for (j=0; j<infos.size(); j++)
			node.score_combos[ScoreComboLoc(infos[j])]=infos[j].score;
		
		score_t max_score=NEG_INF;
		for (j=1; j<infos.size(); j++)
			if (infos[j].score>max_score)
				max_score=infos[j].score;

		node.score = max_score; 
		node.breakage.score = node.score;
	}

	prm->set_has_node_combo_scores(true);
}

// performs scoring on demand (if the combo was not previously scored, calculates
// values, otherwise returns hashed value
score_t AdvancedScoreModel::get_node_combo_score(PrmGraph *prm, int node_idx, 
										 int in_edge_idx,  int in_var_idx, 
										 int out_edge_idx, int out_var_idx) const
{
	Config *config = prm->get_config();
	const Node& node = prm->get_node(node_idx);
	ScoreComboLoc loc(in_edge_idx,out_edge_idx,in_var_idx,out_var_idx);

	const map<ScoreComboLoc,score_t>::const_iterator it = node.score_combos.find(loc);
	if (it != node.score_combos.end())
		return it->second;

	if (node.type == NODE_N_TERM || node.type == NODE_C_TERM)
	{
		const score_t terminal_score=config->get_terminal_score();
		Node& non_const_node = prm->get_non_const_node(node_idx);
		non_const_node.score_combos[loc]=terminal_score;
		return terminal_score;
	}

	const RegionalScoreModel& score_model = 
			regional_breakage_score_models[prm->get_charge()][prm->get_size_idx()][node.breakage.region_idx];
	Node& non_const_node = prm->get_non_const_node(node_idx);

	BreakageInfo info;
	prm->fill_breakage_info(this,&info,node_idx,in_edge_idx,in_var_idx,out_edge_idx,out_var_idx);

	// score not here, calculate the combo score, store it, and return it
	const int charge = prm->get_charge();
	const int size_idx = prm->get_size_idx();
	Spectrum *spec = prm->get_source_spectrum();
	mass_t pm_with_19 = prm->get_pm_with_19();
	const RegionalScoreModel& rsm = regional_breakage_score_models[charge][size_idx][info.breakage->region_idx];

	score_t combo_score = rsm.score_a_single_breakage_combo(prm, non_const_node, info.breakage,info);

	non_const_node.score_combos[loc]=combo_score;
	non_const_node.breakage.score = combo_score;

	return combo_score;

}







/****************************************************************************
Returns the score of the constant element in the breakage:
strong features tha don't depend on aa, and all weak features
*****************************************************************************/
void RegionalScoreModel::calc_constant_element(
							   Node& node,
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage) const
{
	node.const_strong_exps.clear();
	node.const_strong_exps.resize(strong_models.size(),0);

	int f;
	for (f=0; f<this->strong_models.size(); f++)
	{
		const StrongFragModel& strong_model = strong_models[f];
		const int frag_idx = strong_model.model_frag_idx;

		if (! strong_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;	
		strong_model.fill_constant_vals(spec,pm_with_19,breakage,sam.f_vals);

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			node.const_strong_exps[f]=strong_model.inten_model.get_sum_exp(sam);
		}
		else
			node.const_strong_exps[f]=strong_model.no_inten_model.get_sum_exp(sam);
	}

	node.const_regular_exps.clear();
	node.const_regular_exps.resize(regular_models.size(),0);

	for (f=0; f<regular_models.size(); f++)
	{
		const RegularFragModel& regular_model = regular_models[f];
		const int frag_idx = regular_model.model_frag_idx;

		if (! regular_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		regular_model.fill_constant_vals(spec,pm_with_19,breakage,sam.f_vals);

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			node.const_regular_exps[f] = regular_model.inten_model.get_sum_exp(sam);
		}
		else
			node.const_regular_exps[f] = regular_model.no_inten_model.get_sum_exp(sam);
	}
}


/****************************************************************************
Returns the score of the variable element in the breakages:
the strong features that are aa dependant.
*****************************************************************************/
score_t RegionalScoreModel::score_a_single_breakage_combo(
							   PrmGraph *prm,  
							   Node& node,
							   const Breakage *breakage,
							   BreakageInfo& info,
							   bool verbose) const
{
	if (node.type == NODE_N_TERM || node.type == NODE_C_TERM)
	{
		const score_t terminal_score=config->get_terminal_score();
		info.score = terminal_score;
		return terminal_score;
	}
	
	Spectrum *spec = prm->get_source_spectrum();
	const mass_t pm_with_19 = prm->get_pm_with_19();

	if (node.const_strong_exps.size()==0)
		calc_constant_element(node,spec,pm_with_19,breakage);

	const bool print_all = false;

	score_t score=0;

	int f;
	for (f=0; f<this->strong_models.size(); f++)
	{
		const StrongFragModel& strong_model = strong_models[f];
		const int frag_idx = strong_model.model_frag_idx;

		if (! strong_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		sam.f_vals.clear();
		strong_model.fill_aa_variable_vals(spec,pm_with_19,breakage,&info,sam.f_vals);
	
		score_t prev = score;
		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			const float var_exp = strong_model.inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_strong_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.001)
				prob=0.001;

			const float log_random_peak = spec->get_peak(breakage->fragments[pos].peak_idx).log_rand_prob;
			score += (strong_model.inten_log_scaling_factor + log(prob) - log_random_peak);

			if (print_all)
			{
				cout << "viz >> SCORE " << score << "\t(lrp= " << log_random_peak << " , ilsf=" << strong_model.inten_log_scaling_factor <<
					" , lp=" << log(prob) << " , lrandp=" << log_random_peak << endl;
			}
		}
		else
		{
			const float var_exp = strong_model.no_inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_strong_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.01)
				prob=0.01;

			score += (strong_model.no_inten_log_scaling_factor + log(prob) - log_one_minus_random);

			if (print_all)
			{
				cout << "no viz >> SCORE " << score << "\t(nilsf=" << strong_model.inten_log_scaling_factor <<
					" , lp=" << log(prob) << " , loneminusrandp=" << log_one_minus_random << endl;
			}
		}

		if (verbose)
		{
			cout << setprecision(4) << fixed << f << "\t" << score-prev << endl;
		}
	}

	for (f=0; f<regular_models.size(); f++)
	{
		const RegularFragModel& regular_model = regular_models[f];
		const int frag_idx = regular_model.model_frag_idx;

		if (! regular_model.ind_has_models)
			continue;

		if (! breakage->is_frag_type_visible(frag_idx))
			continue;

		ME_Regression_Sample sam;
		sam.f_vals.clear();
		regular_model.fill_aa_variable_vals(spec,pm_with_19,breakage,&info,sam.f_vals);
		
		score_t prev= score;

		const int pos = breakage->get_position_of_frag_idx(frag_idx);
		if (pos>=0)
		{
			const float var_exp = regular_model.inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_regular_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.01)
				prob=0.01;

			const float log_random_peak = spec->get_peak(breakage->fragments[pos].peak_idx).log_rand_prob;
			score += (regular_model.inten_log_scaling_factor + log(prob) - log_random_peak);
		}
		else
		{
			const float var_exp = regular_model.no_inten_model.get_sum_exp(sam);
			const float e = exp(var_exp + node.const_regular_exps[f]);

			float prob = e/(1.0 + e);
			if (prob>0.99)
				prob=0.99;
			if (prob<0.02)
				prob=0.02;

			score += (regular_model.no_inten_log_scaling_factor + log(prob) - log_one_minus_random);
		}

		if (verbose)
		{
			cout << f << "\t" << score-prev << endl;
		}
	}

	// correct for digest node scores
	if (node.type == NODE_DIGEST)
	{
		const score_t digest_score=config->get_digest_score();

		if ((info.connects_to_N_term && info.n_edge_idx<0) || 
			(info.connects_to_C_term && info.c_edge_idx<0) )
		{
			score -= digest_score; // penalty for ending here!
		}
			
		if (info.connects_to_N_term && info.preferred_digest_aa_N_term && score<0)
			score = 0;
		if (info.connects_to_C_term && info.preferred_digest_aa_C_term && score<0)
			score = 0;
	}

	return score;
}



struct NodeType {
	int   type;
	int   region;
	float org_score;
	float mod_score;
};




void AdvancedScoreModel::learn_prm_normalizer_values(const FileManager& fm)
{
	const float step = 0.5;
	const float min_delta = -1.0;
	const float max_delta = 7.0;
	const float target_mid_ratio = 0.96;
	const float target_side_ratio = 0.94;


	config.set_use_spectrum_charge(1);

	regional_prm_normalizers.resize(regional_breakage_score_models.size());
	int c;
	for (c=0; c<regional_breakage_score_models.size(); c++)
	{
		regional_prm_normalizers[c].resize(regional_breakage_score_models[c].size());
		int s;
		for (s=0; s<regional_breakage_score_models[c].size(); s++)
			regional_prm_normalizers[c][s].resize(regional_breakage_score_models[c][s].size(),0);
	}
	

	const vector< vector<mass_t> >& mass_threshes = config.get_size_thresholds();
	for (c=1; c<regional_prm_normalizers.size(); c++)
	{
		int s;
		for (s=0; s<regional_prm_normalizers[c].size(); s++)
		{
			const mass_t min_mass = (s == 0 ? 0 : mass_threshes[c][s-1]);
			const mass_t max_mass =  mass_threshes[c][s];
			const int num_regions = regional_prm_normalizers[c][s].size();
			
			cout << "Finding normalizers for charge " << c << " size " << s << "  (masses " << min_mass << " - " <<
				max_mass << ")" << endl;

			FileSet fs;
			BasicSpecReader bsr;

			fs.select_files_in_mz_range(fm,min_mass/c,max_mass/c,c);
			fs.randomly_reduce_ssfs(2000);

			vector< vector< NodeType > > all_prms;
			const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

			if (fs.get_total_spectra()<50)
			{
				cout << "Insufficient number of spectra... skipping" << endl;
				continue;
			}

			int sc;
			for (sc=0; sc<all_ssf.size(); sc++)
			{
				PrmGraph prm;
				static vector<QCPeak> peaks;
				SingleSpectrumFile *ssf = all_ssf[sc];
				if (peaks.size()<ssf->num_peaks)
				{
					int new_size = ssf->num_peaks*2;
					if (new_size<2500)
						new_size=2500;
					peaks.resize(new_size); 
				}

				const int num_peaks = bsr.read_basic_spec(&config,fm,ssf,&peaks[0]);	
				if (num_peaks<5)
					continue;

				// convert peak list ot a spectrum with charge (if original charge ==0)
				// the spectrum gets charge 2, but the true charge is computed from the data
			
				Spectrum s;
				s.init_from_QCPeaks(&config,&peaks[0],num_peaks,ssf);

				vector<mass_t> pms_with_19;
				vector<int>    charges;

				pms_with_19.clear();
				charges.clear();		
				
				BasicSpectrum bs;
				bs.ssf = ssf;
				bs.peaks = &peaks[0];
				bs.num_peaks = num_peaks;

				// output m/z and prob values for the different charge states
				
				select_pms_and_charges(&config,bs,pms_with_19,charges);
				if (pms_with_19.size()<=0)
					continue;
			
				s.set_charge(charges[0]);
				init_model_for_scoring_spectrum(&s);
				prm.create_graph_from_spectrum(this,&s,pms_with_19[0]);

				vector<NodeType> spec_prms;
				vector<mass_t>   exp_masses;
				const mass_t true_mass_with_19 = s.get_true_mass_with_19();
				s.get_peptide().calc_expected_breakage_masses(&config,exp_masses);

				int i;
				for (i=1; i<prm.get_num_nodes()-1; i++)
				{
					const Node& node = prm.get_node(i);
					if (node.score == 0)
						continue;
					
					NodeType nt;

					nt.type = 0;
					int j;
					for (j=0; j<exp_masses.size(); j++)
						if (fabs(exp_masses[j]-node.mass)<config.get_tolerance())
						{
							nt.type=1;
							break;
						}
					
					if (nt.type<=0)
					{
						int j;
						for (j=0; j<exp_masses.size(); j++)
							if (fabs(true_mass_with_19 - exp_masses[j] -node.mass-MASS_PROTON)<config.get_tolerance())
							{
								nt.type=2;
								break;
							}
					}

					nt.org_score = node.score;
					nt.mod_score = node.score;
					nt.region = node.breakage.region_idx;
					spec_prms.push_back(nt);
				}
				all_prms.push_back(spec_prms);
			}
		
	
			vector< vector< double > > per_pre, per_suf, per_covered;
			vector<float> deltas;

			per_pre.resize(num_regions);
			per_suf.resize(num_regions);
			per_covered.resize(num_regions);

			float delta;
			for (delta = min_delta; delta<=max_delta; delta+=step )
			{
				// perform mods
				int a;
				for (a=0; a<all_prms.size(); a++)
				{
					int b;
					for (b=0; b<all_prms[a].size(); b++)
					{
						NodeType& nt = all_prms[a][b];
						if (nt.org_score< -delta)
						{
							nt.mod_score = NEG_INF;
							continue;
						}
						
					/*	if (nt.org_score>delta)
						{
							nt.mod_score = nt.org_score ;
						}
						else
							nt.mod_score = nt.org_score + (delta-nt.org_score)*0.5;*/
						nt.mod_score = nt.org_score + delta;
					}
				}

				// compute stats (if score is negative treat as 0)
				vector<double> num_pre,num_suf;
				vector<double> num_pre_wpos, num_suf_wpos;
				vector<double> score_pre, score_suf, total_score;
			

				num_pre.resize(num_regions,0);
				num_suf.resize(num_regions,0);
				num_pre_wpos.resize(num_regions,0);
				num_suf_wpos.resize(num_regions,0);
				score_pre.resize(num_regions,0);
				score_suf.resize(num_regions,0);
				total_score.resize(num_regions,0);
				
				for (a=0; a<all_prms.size(); a++)
				{
					int b;
					for (b=0; b<all_prms[a].size(); b++)
					{
						const int   type =    all_prms[a][b].type;
						const float score =   all_prms[a][b].mod_score;
						const int   region =  all_prms[a][b].region;

						if (type == 1)
						{
							num_pre[region]++;
							if (score>0)
							{
								num_pre_wpos[region]++;
								score_pre[region]+= score;
							}
						}

						if (type == 2)
						{
							num_suf[region]++;
							if (score>0)
							{
								num_suf_wpos[region]++;
								score_suf[region]+=score;
							}
						}

						if (score>0)
							total_score[region]+=score;
					}
				}

				
				deltas.push_back(delta);
				int r;
				for (r=0; r<num_regions; r++)
				{
					per_pre[r].push_back(num_pre_wpos[r]/num_pre[r]);
					per_suf[r].push_back(num_suf_wpos[r]/num_suf[r]);
					per_covered[r].push_back((score_pre[r]+score_suf[r])/total_score[r]);
				}
			}

			// report
			int r;
			for (r=0; r<num_regions; r++)
			{
				cout << endl << "Region " << r << endl;
				int d;
				for (d=0; d<deltas.size(); d++)
					cout << "\t" << deltas[d];
				cout << endl << "% Pre";
				for (d=0; d<per_pre[r].size(); d++)
					cout << "\t" << per_pre[r][d];
				cout << endl << "% Suf";
				for (d=0; d<per_suf[r].size(); d++)
					cout << "\t" << per_suf[r][d];
				cout << endl << "% Cov";
				for (d=0; d<per_covered[r].size(); d++)
					cout << "\t" << per_covered[r][d];
				cout << endl;

				// select
				float target_val = target_mid_ratio;
				if (r==0 || r == num_regions-1)
					target_val = target_side_ratio;

				float best_val=POS_INF;
				float best_delta=0;

				for (d=0; d<deltas.size(); d++)
					if (fabs(per_pre[r][d]-target_val)<best_val)
					{
						best_val = fabs(per_pre[r][d]-target_val);
						best_delta = deltas[d];
					}
				
				cout << "Chose delta = " << best_delta << endl << endl;
				regional_prm_normalizers[c][s][r]=best_delta;
			}	
		}
	}
}



bool AdvancedScoreModel::read_prm_normalizer_values()
{

	string file_path = config.get_resource_dir() + "/" + this->model_name + "_prm_norm.txt";

	regional_prm_normalizers.resize(regional_breakage_score_models.size());
	int c;
	for (c=0; c<regional_breakage_score_models.size(); c++)
	{
		regional_prm_normalizers[c].resize(regional_breakage_score_models[c].size());
		int s;
		for (s=0; s<regional_breakage_score_models[c].size(); s++)
			regional_prm_normalizers[c][s].resize(regional_breakage_score_models[c][s].size(),0);
	}


	ifstream in_stream(file_path.c_str());
	if (! in_stream.good())
	{
	//	cout << "Error: couldn't open prm_norm file for reading: " << file_path << endl;
	//	exit(1);
		return false;
	}


	while (! in_stream.eof())
	{
		char line_buff[256];
		in_stream.getline(line_buff,256);

		istringstream iss(line_buff);
		int c=0,s=-1,n=0;
		iss >> c >> s >> n;

		if (c<=0 || s<0 || n<=0)
			continue;

		if (c > regional_prm_normalizers.size() ||
			s > regional_prm_normalizers[c].size() ||
			n > regional_prm_normalizers[c][s].size() )
		{
			cout << "Error: mismatch between regional models and prm normalizers: " << c << " " << s << " " << n << endl;
			exit(1);
		}

		int r;
		for (r=0; r<n; r++)
		{
			iss >> regional_prm_normalizers[c][s][r];
		}
	}

	in_stream.close();
	return true;
}

void AdvancedScoreModel::write_prm_normalizer_values() const
{
	string file_path = config.get_resource_dir() + "/" + this->model_name + "_prm_norm.txt";

	ofstream out_stream(file_path.c_str());
	if (! out_stream.good())
	{
		cout << "Error: couldn't open prm_norm file for writing: " << file_path << endl;
		exit(1);
	}

	out_stream << setprecision(3);
	int c;
	for (c=1; c<regional_prm_normalizers.size(); c++)
	{
		if (regional_prm_normalizers[c].size()==0)
			continue;
		int s;
		for (s=0; s<regional_prm_normalizers[c].size(); s++)
		{
			out_stream << c << " " << s << " " << regional_prm_normalizers[c][s].size();
			int r;
			for (r=0; r<regional_prm_normalizers[c][s].size(); r++)
				out_stream << " " << regional_prm_normalizers[c][s][r];
			out_stream << endl;
		}
	}
	out_stream.close();
}


/****************************************************************************
This function should be used for adjusting the scores of the PrmGraph
when printing of the node scores only
*****************************************************************************/
void AdvancedScoreModel::normalize_prm_scores(PrmGraph &prm) const
{
	const int charge   = prm.get_charge();
	const int size_idx = prm.get_size_idx();
	const vector<score_t>& region_norms = this->regional_prm_normalizers[charge][size_idx];

	const int num_nodes = prm.get_num_nodes();
	int i;
	for (i=1; i<num_nodes-1; i++)
	{
		Node& node = prm.get_non_const_node(i);

		if (node.type == NODE_DIGEST)
			continue;

		const int region = node.breakage.region_idx;
		const score_t norm_val = region_norms[region];
		if (node.score<-norm_val)
		{
			node.score=NEG_INF;
			continue;
		}

		if  (node.score>norm_val)
			continue;

		node.score += (norm_val - node.score)*0.5;
	}
}


struct PeakTuple {
	bool operator< (const PeakTuple& rhs) const
	{
		return (score>rhs.score);
	}
	size_t frag_idx;
	size_t pos;
	int	   rank;
	float  score;
};


void predict_fragmentation(AdvancedScoreModel* model, const char* input_file, size_t num_peaks)
{
	FILE* stream = fopen(input_file,"r");
	if (! stream)
	{
		cout << "Error: couldn't open file for reading: " << input_file << endl;
		exit(1);
	}

	Config *config = model->get_config();
	DeNovoRankScorer *dnv_rank = (DeNovoRankScorer *)model->get_rank_model_ptr(1);
	PeakRankModel *prm = dnv_rank->get_peak_prediction_model(3);


	char buffer[128];
	char pep_str[128];
	while (fgets(buffer,128,stream))
	{
		int charge;

		if (sscanf(buffer,"%s %d",pep_str,&charge) != 2)
			continue;

		cout << ">> " << pep_str << "\t" << charge << endl;
		if (charge<1 || charge>=prm->get_size_thresholds().size())
		{
			cout << "Invalid charge!" << endl;
			continue;
		}

		Peptide pep;
		pep.parse_from_string(config,static_cast<string>(pep_str));


		PeptideSolution sol;
		sol.pep = pep;
		sol.reaches_n_terminal=true;
		sol.reaches_c_terminal=true;
		sol.charge = charge;
		sol.pm_with_19 = pep.get_mass_with_19();

		PeptidePeakPrediction ppp;
		prm->calc_peptide_predicted_scores(sol, ppp);

		const size_t num_frags = ppp.frag_idxs.size();
		vector< vector<int> > predicted_ranks;
		calc_combined_peak_ranks(ppp.rank_scores, predicted_ranks);

		vector<PeakTuple> tuples;
		for (size_t f=0; f<num_frags; f++)
			for (size_t i=0; i<ppp.rank_scores[f].size(); i++)
				if (predicted_ranks[f][i]<999)
				{
					PeakTuple pt;
					pt.frag_idx = f;
					pt.pos =i;
					pt.rank = predicted_ranks[f][i];
					pt.score = ppp.rank_scores[f][i];
					tuples.push_back(pt);
				}
		
		sort(tuples.begin(),tuples.end());

		if (tuples.size()<1)
			continue;

		const size_t num_aas = pep.get_num_aas();
		vector<mass_t> breakage_masses;
		pep.calc_expected_breakage_masses(config, breakage_masses);

		cout << fixed << "Rank\tIon\tm/z\tScore" << endl;
		for (size_t i=0; i<num_peaks && i<tuples.size(); i++)
		{
			PeakTuple pt = tuples[i];
			cout << i+1 << "\t";
			const FragmentType& ft = config->get_fragment(ppp.frag_idxs[pt.frag_idx]);
			cout << ft.label << ":" << (ft.orientation == PREFIX ? pt.pos : num_aas - pt.pos) << "\t";

			mass_t mz =  ft.calc_expected_mass(breakage_masses[pt.pos],pep.get_mass_with_19());
			cout << setprecision(2);
			if (mz<100)
				cout << " ";
			if (mz<1000)
				cout << " ";

			cout << mz << "\t";
			cout << setprecision(3) << pt.score << endl;
		}
		cout << endl;
	}

	fclose(stream);
}

