#ifndef __DISCRETEPEAKMDODEL_H__
#define __DISCRETEPEAKMDODEL_H__

#include "Model.h"
#include "AnnotatedSpectrum.h"
#include "FileManagement.h"
#include "RegionalPepNovoModel.h"
#include "EdgeModel.h"
#include "PrmGraph.h"
#include "includes.h"

#define NUM_PREVIOUS_FRAG_VALS 4

struct Edge; // forward declr
// struct SeqPath;
class  PrmGraph;

class DiscretePeakModel : public Model {
public:

	virtual void write_model_specifics(ostream& os) const =0;

	virtual void read_model_specifics(istream& is) =0;

	void write_tables(ostream& os) const;

	void read_tables(ostream& os) const;

	void clone_charge_model(int source_charge, int target_charge);

	void init_score_model() { init_default_thresholds(); }

	virtual void init_default_thresholds() =0;

	// this function should always be run once before each new spectrum is 
	void init_model_for_scoring_spectrum(Spectrum *spec) =0;
	

	// prints joint scores of two top  fragments
	void print_joint_scores(ostream& os = cout) const;

	void print_report(int charge, int size, int region, ostream& os = cout) const;
	
	void print_frag_probs(const vector<int>& frag_type_idxs, 
		int charge, int size, int region, ostream& os) const;

	void print_all_table_names(ostream& os = cout) const;
	virtual void print_level_legend( ostream& os) const =0;

	// This is the interface of the model with the rest of the program
	// all peaks get their breakages set by the model
	virtual void set_breakage_peak_levels(Breakage *breakage) const =0;

	void score_breakage(Spectrum *spec, Breakage *breakage, bool verbose=false) const
	{
		set_breakage_peak_levels(breakage);

		score_t b_score = regional_models[breakage->parent_charge]
			[breakage->parent_size_idx][breakage->region_idx].calc_breakage_score(breakage,
			verbose, spec->get_config());

		breakage->score = b_score;

		add_breakage_iso_score(spec, breakage);
	}


	void score_all_node_combos(PrmGraph *prm) const {};
	void initial_combos_score(PrmGraph *prm) const {};
	void score_peptide_node_combos(PrmGraph *prm, const Peptide& peptide ) const {};

	score_t get_node_combo_score(PrmGraph *prm, int node_idx,
								 int in_edge_idx, int in_var_idx, 
								int out_edge_idx, int out_var_idx) const
	{
		return prm->get_node(node_idx).score;
	}


	score_t get_missing_breakage_score(int charge, int size_idx, int region_idx) const
	{
		return  regional_models[charge][size_idx][region_idx].get_missed_cleavage_score();
	}

	
	void score_graph_edges(PrmGraph& prm) const
	{
		edge_model.score_graph_edges(prm);
	}


	int get_max_score_model_charge() const
	{
		return max_score_model_charge;
	}

	void normalize_prm_scores(PrmGraph &prm) const
	{
		return;
	}

protected:
	
	int max_score_model_charge;

	int num_peak_levels;          // the number of levels + 1 (for level 0)

	vector<double> q_rand;           // rank (rank 0 = no peak);
	vector< vector<double> > q_frag; // rank (rank 0 = no peak), frag_type_idx

	Spectrum *current_spectrum; // the current spectrum being ranked

	vector<int> level_thresholds; // rank levels for rank models


	vector< vector< vector<RegionalPepNovoModel> > > regional_models; // charge, size, region

	EdgeModel edge_model;


	// calcs the score of participating fragments that have istopic peak scores
	void add_breakage_iso_score(Spectrum *spec, Breakage *breakage) const;


	// creates models allowing each fragment to have upto two parents
	// (chooses the ones that give the best difference in probability
	// in terms of DKL). Parent fragments must appear before the current
	// fragment in the order of the regional fragment set
	void learn_parent_dependencies(const FileManager& fm, int charge = 2);

	// learns the probabilities of randomly observing peaks
	void learn_q_rand(const FileManager& fm, int charge);


	// this function is the interface for querying on peak's intensity level,
	// rank etc. The peak's level might depend on the type of fragment that
	// is inteded to be used, so there is an option to add this to the query
	// (though systems like regular rank won't use the second parameter).
	virtual int get_peak_level(int p_idx, int f_idx = -1) const =0;


	virtual int get_lowest_level_in_spectrum(Spectrum *spec) const =0;



	void train_score_model(const char *name, const FileManager& fm, int charge=0,
						   int size_idx=-1, int region_idx=-1);


	// reads all relevant info: 
	// thresholds and tables
	void read_score_model(const char *name, bool silent_ind = false);

	// writes all relevant info:
	// thresholds and tables
	void write_score_model(const char *name) const;


	// converts all the probabilities to socres in the tables
	void convert_probs_to_scores();


};


#endif



