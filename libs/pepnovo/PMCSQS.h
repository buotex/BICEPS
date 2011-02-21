#ifndef __PMSQS_H__
#define __PMSQS_H__

#include "ME_REG.h"
#include "RankBoost.h"
#include "QuickClustering.h"
#include "includes.h"

#define MIN_SPECTRA_FOR_PMCSQS_MODEL 100

// The 3 SQS model makes binary decisions: good charge c spectrum or not good for charges 
// c=1,2,3.
// The PMC model assumes that we have a good spectrum and in the correct charge.
// 




typedef enum SQS_Fields 
{ 
	//0
	SQS_CONST, SQS_PEAK_DENSITY,

	SQS_PROP_UPTO2G,  SQS_PROP_UPTO5G, 
	SQS_PROP_UPTO10G, SQS_PROP_MORE10G,

	//14
	SQS_PROP_INTEN_UPTO2G, SQS_PROP_INTEN_UPTO5G, SQS_PROP_INTEN_MORE5G,

	SQS_PROP_ISO_PEAKS, SQS_PROP_STRONG_WITH_ISO_PEAKS,

	//19
	SQS_PROP_ALL_WITH_H2O_LOSS,    SQS_PROP_ALL_WITH_NH3_LOSS,    SQS_PROP_ALL_WITH_CO_LOSS,
	SQS_PROP_STRONG_WITH_H2O_LOSS, SQS_PROP_STRONG_WITH_NH3_LOSS, SQS_PROP_STRONG_WITH_CO_LOSS,

	SQS_C2_PROP_ALL_WITH_H2O_LOSS,    SQS_C2_PROP_ALL_WITH_NH3_LOSS,    SQS_C2_PROP_ALL_WITH_CO_LOSS,
	SQS_C2_PROP_STRONG_WITH_H2O_LOSS, SQS_C2_PROP_STRONG_WITH_NH3_LOSS, SQS_C2_PROP_STRONG_WITH_CO_LOSS,

	SQS_DIFF_ALL_WITH_H2O_LOSS,    SQS_DIFF_ALL_WITH_NH3_LOSS,    SQS_DIFF_ALL_WITH_CO_LOSS,
	SQS_DIFF_STRONG_WITH_H2O_LOSS, SQS_DIFF_STRONG_WITH_NH3_LOSS, SQS_DIFF_STRONG_WITH_CO_LOSS,
	
	SQS_PROP_PEAKS_WITH_C1C2, SQS_PROP_STRONG_PEAKS_WITH_C1C2, SQS_PROP_INTEN_WITH_C1C2,

	SQS_IND_MAX_TAG_LENGTH_ABOVE_4, SQS_IND_MAX_TAG_LENGTH_BELOW_4,
	SQS_MAX_TAG_LENGTH_ABOVE_4,     SQS_MAX_TAG_LENGTH_BELOW_4,

	SQS_PROP_INTEN_IN_TAGS,
	SQS_PROP_TAGS1, SQS_PROP_STRONG_PEAKS_IN_TAG1, SQS_PROP_INTEN_TAG1, 
	SQS_IND_PROP_STRONG_BELOW30_TAG1,
	SQS_PROP_TAGS2, SQS_PROP_STRONG_PEAKS_IN_TAG2, SQS_PROP_INTEN_TAG2,	
	SQS_IND_PROP_STRONG_BELOW20_TAG2,
	SQS_PROP_TAGS3, SQS_PROP_STRONG_PEAKS_IN_TAG3, SQS_PROP_INTEN_TAG3, 
	SQS_IND_PROP_STRONG_BELOW10_TAG3,

	SQS_C2_IND_MAX_TAG_LENGTH_ABOVE_4, SQS_C2IND_MAX_TAG_LENGTH_BELOW_4,
	SQS_C2MAX_TAG_LENGTH_ABOVE_4,      SQS_C2MAX_TAG_LENGTH_BELOW_4,

	
	SQS_C2PROP_INTEN_IN_TAGS,
	SQS_C2PROP_TAGS1, SQS_C2PROP_STRONG_PEAKS_IN_TAG1, SQS_C2PROP_INTEN_TAG1, 
	SQS_IND_C2PROP_STRONG_BELOW30_TAG1,
	SQS_C2PROP_TAGS2, SQS_C2PROP_STRONG_PEAKS_IN_TAG2, SQS_C2PROP_INTEN_TAG2,	
	SQS_IND_C2PROP_STRONG_BELOW20_TAG2,
	SQS_C2PROP_TAGS3, SQS_C2PROP_STRONG_PEAKS_IN_TAG3, SQS_C2PROP_INTEN_TAG3, 
	SQS_IND_C2PROP_STRONG_BELOW10_TAG3,

	SQS_DIFF_MAX_TAG_LENGTH, SQS_DIFF_PROP_INTEN_IN_TAGS,
	SQS_DIFF_PROP_TAGS1, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG1, SQS_DIFF_PROP_INTEN_TAG1,
	SQS_DIFF_PROP_TAGS2, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG2, SQS_DIFF_PROP_INTEN_TAG2,
	SQS_DIFF_PROP_TAGS3, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG3, SQS_DIFF_PROP_INTEN_TAG3,
	
	SQS_PEAK_DENSE_T1,  SQS_PEAK_DENSE_T2,  SQS_PEAK_DENSE_T3,
	SQS_INTEN_DENSE_T1, SQS_INTEN_DENSE_T2, SQS_INTEN_DENSE_T3,

	SQS_PEAK_DENSE_H1,  SQS_PEAK_DENSE_H2,
	SQS_INTEN_DENSE_H1, SQS_INTEN_DENSE_H2,

	SQS_PROP_MZ_RANGE_WITH_33_INTEN,
	SQS_PROP_MZ_RANGE_WITH_50_INTEN,
	SQS_PROP_MZ_RANGE_WITH_75_INTEN,
	SQS_PROP_MZ_RANGE_WITH_90_INTEN,

	// get these from the PM features (max values found)
	// all values are after subtracting the maximum background values
	// obtained by using erroneous parent masses

	SQS_NUM_FRAG_PAIRS_1,			SQS_NUM_STRONG_FRAG_PAIRS_1,
	SQS_NUM_C2_FRAG_PAIRS_1,		SQS_NUM_STRONG_C2_FRAG_PAIRS_1,

	SQS_NUM_FRAG_PAIRS_2,			SQS_NUM_STRONG_FRAG_PAIRS_2,
	SQS_NUM_C2_FRAG_PAIRS_2,		SQS_NUM_STRONG_C2_FRAG_PAIRS_2,

	SQS_NUM_FRAG_PAIRS_3,			SQS_NUM_STRONG_FRAG_PAIRS_3,
	SQS_NUM_C2_FRAG_PAIRS_3,		SQS_NUM_STRONG_C2_FRAG_PAIRS_3,

	SQS_PROP_OF_MAX_FRAG_PAIRS_1,			SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_1,
	SQS_PROP_OF_MAX_C2_FRAG_PAIRS_1,		SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_1,

	SQS_PROP_OF_MAX_FRAG_PAIRS_2,			SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_2,
	SQS_PROP_OF_MAX_C2_FRAG_PAIRS_2,		SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_2,

	SQS_PROP_OF_MAX_FRAG_PAIRS_3,			SQS_PROP_OF_MAX_STRONG_FRAG_PAIRS_3,
	SQS_PROP_OF_MAX_C2_FRAG_PAIRS_3,		SQS_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_3,

	SQS_PROP_FRAG_PAIRS_1,			SQS_PROP_STRONG_FRAG_PAIRS_1,
	SQS_PROP_C2_FRAG_PAIRS_1,		SQS_PROP_STRONG_C2_FRAG_PAIRS_1,

	SQS_PROP_FRAG_PAIRS_2,			SQS_PROP_STRONG_FRAG_PAIRS_2,
	SQS_PROP_C2_FRAG_PAIRS_2,		SQS_PROP_STRONG_C2_FRAG_PAIRS_2,

	SQS_PROP_FRAG_PAIRS_3,			SQS_PROP_STRONG_FRAG_PAIRS_3,
	SQS_PROP_C2_FRAG_PAIRS_3,		SQS_PROP_STRONG_C2_FRAG_PAIRS_3,

	SQS_DIFF_NUM_FRAG_PAIRS_23,					SQS_DIFF_NUM_STRONG_FRAG_PAIRS_23,
	SQS_DIFF_NUM_C2_FRAG_PAIRS_23,				SQS_DIFF_NUM_STRONG_C2_FRAG_PAIRS_23,
	SQS_DIFF_PROP_OF_MAX_FRAG_PAIRS_23,			SQS_DIFF_PROP_OF_MAX_STRONG_FRAG_PAIRS_23,
	SQS_DIFF_PROP_OF_MAX_C2_FRAG_PAIRS_23,		SQS_DIFF_PROP_OF_MAX_STRONG_C2_FRAG_PAIRS_23,
	SQS_DIFF_PROP_FRAG_PAIRS_23,				SQS_DIFF_PROP_STRONG_FRAG_PAIRS_23,
	SQS_DIFF_PROP_C2_FRAG_PAIRS_23,				SQS_DIFF_PROP_STRONG_C2_FRAG_PAIRS_23,


	SQS_NUM_FIELDS 
} SQS_Fields;







struct PMCRankStats {

	void clear();

	PMCRankStats() { clear(); }
	PMCRankStats(const PMCRankStats& other);

	mass_t m_over_z;

	score_t rank_score;

	float num_frag_pairs;
	float num_strong_frag_pairs;
	float num_c2_frag_pairs;
	float num_strong_c2_frag_pairs;
	float num_h2o_loss_frag_pairs;
	float num_h2o_loss_c2_frag_pairs;

	float inten_frag_pairs;
	float inten_strong_pairs;
	float inten_c2_pairs;
	float inten_c2_strong_pairs;
	float inten_h2o_loss_frag_pairs;
	float itnen_h2o_loss_c2_frag_pairs;

	float mean_offset_pairs;
	float mean_offset_strong_pairs;
	float mean_offset_c2_pairs;
	float mean_offset_c2_strong_pairs;
	float mean_offset_h2o_pairs;
	float mean_offset_c2_h2o_pairs;

	bool ind_pairs_with_min_tol;
	bool ind_strong_pairs_with_min_tol;
	bool ind_c2_pairs_with_min_tol;
	bool ind_c2_strong_pairs_with_min_tol;

	float log_dis_from_pairs_min_tol;
	float log_dis_from_strong_pairs_min_tol;
	float log_dis_from_c2_pairs_min_tol;
	float log_dis_from_c2_strong_pairs_min_tol;

	vector<float> offset_pairs_ordered_by_inten;
	vector<float> strong_offset_pairs_ordered_by_inten;
	vector<float> c2_offset_pairs_ordered_by_inten;

	float num_strict_pairs0, inten_strict_pairs0;
	float num_strict_pairs1, inten_strict_pairs1;
	float num_strict_pairs2, inten_strict_pairs2;

	float c2_num_strict_pairs0, c2_inten_strict_pairs0;
	float c2_num_strict_pairs1, c2_inten_strict_pairs1;
	float c2_num_strict_pairs2, c2_inten_strict_pairs2;
};


// per charge
struct PmcSqsChargeRes {
	
	PmcSqsChargeRes() : mz1(-1), score1(NEG_INF), mz2(-1), score2(NEG_INF), min_comp_prob(-1), sqs_prob(0) {};

	float mz1;
	float score1;
	float mz2;
	float score2;

	float min_comp_prob; // the minimal probaiblity when comparing to other charges
	float sqs_prob;
};

struct  DPColumn {
	int pointers[Val+1];
};


class PMCSQS_Scorer {
public:

	PMCSQS_Scorer() : max_model_charge(0),
				frag_pair_sum_offset(NEG_INF), bin_increment(NEG_INF), 
				ind_initialized_pmcr(false), ind_initialized_sqs(false), config(NULL), 
				curr_spec_total_intensity(0), curr_spec_strong_intensity(0),
				curr_spec_num_strong(0) {};
	~PMCSQS_Scorer();
	bool init_for_current_spec(Config *config, const BasicSpectrum& bs);

	void   set_frag_pair_sum_offset(mass_t offset) { frag_pair_sum_offset = offset; }
	mass_t get_frag_pair_sum_offset() const { return frag_pair_sum_offset; }

	void   set_bin_increment(mass_t inc) { bin_increment = inc; }
	mass_t get_bin_increment() const { return bin_increment; }


	const vector< vector<PMCRankStats> >& get_curr_spec_rank_pmc_tables() const { return curr_spec_rank_pmc_tables; }

	void find_best_mz_values_from_rank_model(const BasicSpectrum& bs, 
											 int charge,
											 mass_t pm_tolerance,
											 PmcSqsChargeRes& res);


	float get_pmcsqs_results_for_spectrum(Config *config, const BasicSpectrum& bs,
						vector<PmcSqsChargeRes>& res);

	float get_best_mz_charge(Config *config, const BasicSpectrum& bs, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res = NULL);

	void select_pms_and_charges(Config *config, 
								const BasicSpectrum& bs,
								vector<mass_t>& pms_with_19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* all_res = NULL);

	void benchmark_pm_selection(Config *config, FileManager& fm, mass_t pm_val_tol);


	float get_sqs_for_spectrum(Config *config, const BasicSpectrum& bs, 
							   int *charge = NULL, bool verbose = false);

	mass_t get_charge_mz_bias(int charge) const { return charge_mz_biases[charge]; }

	void print_spec(const BasicSpectrum& bs) const;

	void test_pmc(Config *config, char *specs_file, int charge, mass_t min_mass=0,
		mass_t max_mass = POS_INF);

	void output_pmc_rank_results(const FileManager& fm, int charge,
				const vector<SingleSpectrumFile *>& test_ssfs); 


	void train_sqs_models(Config *config, const FileManager& pos_fm, char *neg_list, 
		int specific_charge=-1, vector<vector<float> > *inp_weights = NULL);

	void train_pmc_rank_models(Config *config, const FileManager& fm, int sel_charge = 0, bool overwrite=true);

	bool read_pmc_rank_models(Config *config, char *file);

	void write_pmc_rank_models(const char *path) const;

	bool read_sqs_models(Config *config, char *file);

	void write_sqs_models(const char *path) const;

	void compute_sqs_cum_stats_for_ided(Config *config, char *list);
	void compute_sqs_cum_stats_for_crap(Config *config, char *list);

	void benchmark_sqs(Config *config, char *list, char *anns);

	bool get_ind_initialized_pmcr() const { return ind_initialized_pmcr; }
	bool get_ind_initialized_sqs() const { return ind_initialized_sqs; }

	int get_max_model_charge() const { return max_model_charge; }

	void output_filtered_spectra_to_mgfs(Config *config,
									 const vector<string>& files,
									 char *out_dir,
									 float filter_prob, 
									 int& num_written, 
									 int& num_read);

private:

	vector< vector<RankBoostModel * > >  pmc_rank_models; // charge / size thresholds
	vector< vector<mass_t> >             pmc_rank_mass_thresholds;
	vector< vector<mass_t> >			 pmc_charge_mz_biases;


	vector<mass_t> sqs_mass_thresholds;
	vector< vector<float> > sqs_correction_factors; // sqs prob += this[charge][size_idx]
	vector< vector<float> > sqs_mult_factors;      // sqs prob *= this[charge][size_idx]

	vector< vector< vector<ME_Regression_Model *> > > sqs_models; // charge1, charge2 (charge2=0 is pos 
																  // vs negative classes) / size_idx



	vector<mass_t> charge_mz_biases; // add this to the predicted m/z

	int max_model_charge; // the maximal charge for which we have a model

	mass_t frag_pair_sum_offset; // b+y or c+z - (PM+19)

	mass_t bin_increment;

	bool ind_initialized_pmcr;

	bool ind_initialized_sqs;


	Config *config;

	float curr_spec_total_intensity;
	float curr_spec_strong_intensity;
	int   curr_spec_num_strong;

	vector<bool>  curr_spec_strong_inds;
	vector<float> curr_spec_iso_levels;
	vector<bool>  curr_spec_strict_iso_ind; // holds for each spectrum ind if there is any peak at -1

//	vector< vector<PmcStats> > curr_spec_pmc_tables;
//	vector<PmcStats>		   curr_spec_background_stats;
//	vector<PmcStats>		   curr_spec_maximal_values;

	vector< vector<PMCRankStats> > curr_spec_rank_pmc_tables;
	vector<PMCRankStats>		curr_spec_rank_background_stats;
	vector<PMCRankStats>		curr_spec_rank_maximal_values;
	
	void fill_SQS_DP(const BasicSpectrum& bs, vector<DPColumn>& dp, int farge_charge ) const;


	void fill_fval_vector_with_SQS(const BasicSpectrum& bs,
								   ME_Regression_Sample& sam) const;

	
	void fill_RankBoost_smaples_with_PMC(const BasicSpectrum& bs,
										 int charge,
										 vector<RankBoostSample>& samples) const;

	
	int get_optimal_bin(int true_mz_bin, int charge) const;

	void select_training_sample_idxs(int charge,
		const vector<ME_Regression_Sample>& spec_samples,
		const BasicSpectrum& bs,
		int& correct_idx,
		vector<int>& bad_pmc_idxs) const;

	void select_training_sample_idxs(int charge,
		const vector<RankBoostSample>& spec_samples,
		const BasicSpectrum& bs,
		int& correct_idx,
		vector<int>& bad_pmc_idxs) const;

	int get_rank_model_size_idx(int charge, mass_t pm_with_19) const
	{
		if (pmc_rank_mass_thresholds.size()<=charge)
		{
			cout << "Error: PMC does not support charge " << charge << endl;
			exit(1);
		}

		int i;
		for (i=0; i<pmc_rank_mass_thresholds[charge].size(); i++)
			if (pm_with_19<=pmc_rank_mass_thresholds[charge][i])
				break;
		return i;
	}



	void calculate_curr_spec_pmc_values( const BasicSpectrum& bs, mass_t bin_increment);

	void get_sqs_features_from_pmc_tables(const BasicSpectrum& bs,
										  vector< vector<float> >& sqs_featrues) const;



	void set_pmc_mass_thresholds(int option = 0); // option 1 , fixed for it data

	void set_sqs_mass_thresholds();

	void set_default_sqs_correct_factors();

	void init_sqs_correct_factors(int max_charge, int num_sizes);

	int  get_sqs_size_idx(mass_t mass) const
	{
		int i;
		for (i=0; i<sqs_mass_thresholds.size(); i++)
			if (mass<=sqs_mass_thresholds[i])
				break;
		return i;
	}

	void create_filtered_peak_list_for_sqs(
									  QCPeak *org_peaks, int num_org_peaks,
									  QCPeak *new_peaks, int& num_new_peaks) const;

};







void fill_rank_PMC_stats(int charge,
						  const mass_t single_charge_pair_offset, // the sum of b+y or c+z
						  mass_t minus_range, 
						  mass_t plus_range,
						  mass_t increment,
						  Config *config,
						  const BasicSpectrum& bs,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  vector<PMCRankStats>& pmc_stats_vec);

void calc_pmc_rank_stats_for_mass(const QCPeak *peaks, 
										  int num_peaks, 
										  mass_t single_charge_pair_sum,
										  mass_t tolerance, 
										  const vector<float>& iso_levels,
										  const vector<bool>& strong_inds,
										  const vector<bool>& strict_iso_inds,
										  PMCRankStats& stats);



void create_training_files(Config *config);

void init_PMC_feature_names(vector<string>& names);




#endif

