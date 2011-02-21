#ifndef __MODEL_H__
#define __MODEL_H__

/*************************************************************************

  This is the top contianer class for PepNovo's models.
  All models must inherit this class. The acutal model class that is used
  is RegularRankModel which is inherits from DiscretePeakModel, which inherits from Model.

  Training/reading/writing of all models is done through the functions of this class.

  The model class contains virtual functions for scoring breakages and edges.

  Before scoring a spectrum the model should be initialized for that spectrum.
**************************************************************************/
#include "Config.h"
#include "BasicDataStructs.h"
#include "Spectrum.h"
#include "FileManagement.h"
#include "PMCSQS.h"
#include "PeptideComp.h"
#include "EdgeModel.h"
#include "AminoAcidProbs.h"
#include "includes.h"



#define MINIMAL_NUMBER_SPECTRA_FOR_FRAGMENT_SELECTION 100
 
struct Edge; // forward declr
struct SeqPath;
class PrmGraph;

class Model {
public:
	
	string get_model_name() const { return model_name; }

	void set_model_name(string _name) { model_name = _name; config.set_model_name(_name); }

	Config * get_config() { return &config; }

	void set_config(Config *_config) { config= *_config; }

	// this function performs the entire training process of the model
	void train_model_in_stages(const char *name, 
		const FileManager& fm, 
		mass_t tolerance, 
		int start_stage = 0,
		int end_stage = POS_INF,
		int specific_charge=-1, 
		int specific_size=-1, 
		int specific_region=-1,
		char *neg_sqs_list = NULL);



	void train_pmc_rank_models(const FileManager& fm, int charge = 0)
	{
		pmcsqs.train_pmc_rank_models(&config,fm,charge);
	}

	void test_pmc(char *specs_file, int charge, mass_t min_mass=0,
		mass_t max_mass = POS_INF)
	{
		pmcsqs.test_pmc(&config,specs_file,charge,min_mass,max_mass);
	}

	void compute_sqs_cum_stats_for_ided(char *list)
	{
		pmcsqs.compute_sqs_cum_stats_for_ided(&config, list);
	}

	void compute_sqs_cum_stats_for_crap(char *list)
	{
		pmcsqs.compute_sqs_cum_stats_for_crap(&config, list);
	}

	void benchmark_sqs(char *list, char *anns)
	{
		pmcsqs.benchmark_sqs(&config, list, anns);
	}


	void train_sqs(const FileManager& pos_fm, char *neg_list, int specific_charge, vector< vector<float> > *weights=NULL)
	{
		pmcsqs.train_sqs_models(&config,pos_fm,neg_list,specific_charge, weights);
	}

	float get_best_mz_charge(Config *config, const BasicSpectrum& bs, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res = NULL)
	{
		return pmcsqs.get_best_mz_charge(config,bs,mz1,charge1,prob1,mz2,charge2,prob2,all_res);
	}

	void select_pms_and_charges(Config *config, 
								const BasicSpectrum& bs,
								vector<mass_t>& pms_with_19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* all_res = NULL)
	{
		pmcsqs.select_pms_and_charges(config,bs,pms_with_19,charges,all_res);
	}



	float get_sqs_only(Config *config, const BasicSpectrum& bs, int *max_charge=NULL)
	{
		return pmcsqs.get_sqs_for_spectrum(config,bs,max_charge);
	}

	bool get_ind_pmcsqs_was_intialized() const { return (pmcsqs.get_ind_initialized_pmcr() &&
		pmcsqs.get_ind_initialized_sqs()); }

	const PMCSQS_Scorer* get_pmcsqs_ptr() const { return &pmcsqs; }

	const AminoAcidProbs* get_amino_acid_probs_ptr() const { return &amino_acid_probs; }

	int get_aa_category(int num_aa, int *aas, bool n_term, bool c_term) const
	{
		const vector<int>& org_aa = config.get_org_aa();
		int conv_aas[6];
		int i;
		for (i=0; i<num_aa; i++)
			conv_aas[i]=org_aa[aas[i]];

		return comp_assigner.get_aa_category(num_aa,conv_aas,n_term,c_term);
	}
	

	
	// uses the current model fragments to test what tolerance can be used that 
	// catches at least *percen_frags_caught* of the dominant fragments
	mass_t calculate_tolerance(const FileManager& fm, mass_t max_tolerance,
							   float percent_frags_caught = 0.95);
	
	// 
	void perform_offset_frequency_function(const FileManager& fm, mass_t tolerance,
										   int max_charge);

	void read_model(const char* model_name, bool silent = false);

	void clone_charge_model(int source_charge, int target_charge);


	void write_model();


	virtual void init_score_model() = 0;

	// in case there are som preliminary actions the model needs to 
	// do before it can score a spectrum (most models might not need this)
	virtual void init_model_for_scoring_spectrum(Spectrum *spec) = 0;

	virtual void score_breakage(Spectrum *spec, Breakage *breakage, bool verbose=false) const =0;

	virtual void score_all_node_combos(PrmGraph *prm) const =0;

	virtual void initial_combos_score(PrmGraph *prm) const =0;

	virtual void score_peptide_node_combos(PrmGraph *prm, const Peptide& peptide) const =0;

	virtual score_t get_node_combo_score(PrmGraph *prm, int node_idx, 
										 int in_edge_idx, int in_var_idx, 
										 int out_edge_idx, int out_var_idx) const =0;
	
	virtual score_t get_missing_breakage_score(int charge, int size_idx, int region_idx) const =0;

	virtual void score_graph_edges(PrmGraph& prm) const =0;

	virtual void normalize_prm_scores(PrmGraph &prm) const =0;

	virtual int get_max_score_model_charge() const =0;


protected:

	string model_name;
	
	Config config;

	PMCSQS_Scorer pmcsqs;

	PeptideCompAssigner comp_assigner;

	EdgeModel edge_model;

	AminoAcidProbs amino_acid_probs;

	vector<bool> stages_intialized;


	bool select_fragments(const char *name, const FileManager& fm,
		int max_num_frags =0, float min_prob = 0.5);


	// The actual training and writing of score models depends on the type 
	// of model that is used

	virtual void train_score_model(const char *name, const FileManager& fm, 
		int charge=0, int size_idx =-1, int region_idx=-1) =0;

	virtual void read_score_model(const char *name, bool silent_ind = false) =0;

	virtual void write_score_model(const char *name) const =0;
};




// determines the tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t calc_tolerance_distribution(Model *model,  const FileManager& fm, mass_t max_tolerance,
								   float cutoff_prob=0.96);

// determines the parent mass tolerance for which *cuttoff_prob* of the abundant fragments
// are caught
mass_t calc_parent_mass_tolerance_distribution(Model *model,  const FileManager& fm, 
											   float cutoff_prob=0.98);


#endif

