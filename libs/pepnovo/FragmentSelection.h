#ifndef __FRAGMENTSELECTION_H__
#define __FRAGMENTSELECTION_H__

#include "FileManagement.h"
#include "AnnotatedSpectrum.h"
#include "Config.h"
#include "includes.h"


// 
void select_frags_using_frag_offset_counts(const FileManager& fm, 
										   Config *config,
										   FragmentTypeSet& fts, 
										   float min_frag_prob = 0.05);


// returns the probablility of observing a the different types of
// fragments for different charges/sizes/regions
void collect_probs_of_known_fragments(const FileManager& fm, Config *config,
				vector< vector< vector< vector<double> > > >& frag_probs, // charge , size, region, frag_idx
				vector< vector< vector< vector<double> > > >& in_range_counts,
				double &avg_rand, int &num_files_used, int charge = 0);


// Selects for a charge state's regional models the set of fragments
// that should be used.
void select_regional_fragments(const FileManager& fm, Config *config, int charge,
						bool verbose = true);



// chooses for each regional model the pair of fragments that identify
// the largest number of cuts not already found using strong fragments
void select_frag_combos(const FileManager& fm, Config *config, int charge,
						int max_num_combos=4);


/*********************************************************************
Adds the counts for peaks around a breakage to their respective bins

**********************************************************************/
void add_offset_counts_arround_mass(vector<int>& counts, Spectrum *spec,
									mass_t min_mass, mass_t max_mass, 
									mass_t bin_coef, mass_t break_mass,
									int charge);




void explore_fragment_set(FileManager& fm, Config *config);

void show_occurences(FileManager& fm, Config *config, string label);

void make_frag_rank_histogram(FileManager& fm, Config *config);

// calculates the average random probability according to the offset frequency function
double calc_avg_rand(FileManager& fm, Config *config);

void make_bin_histograms(FileManager& fm, Config *config);

void select_fragments_from_bins(vector<int>& counts, FragmentTypeSet& fts, int max_num_frags, 
								int charge, int orientation, mass_t min_offset_mass, 
								mass_t bin_coef, mass_t tolerance);

void calculate_true_fragment_probabilities(const FileManager& fm, Config *config,
							 FragmentTypeSet& fts,  float min_prob);

#endif

