#include "QuickClustering.h"
#include "RegularRankModel.h"
#include "PMCSQS.h"
#include "auxfun.h"

// the sim matrix stores the similarity distances computed between clusters
int     num_sim_matrix_spectra = 0;
unsigned char  * sim_matrix = NULL;
unsigned char  * max_sim_addr = NULL;



void print_byte(unsigned char byte)
{
	int i;
	unsigned char mask=0X1;

	for (i=0; i<8; i++)
	{
		cout << ( (byte & mask) ? '1' : '0');
		mask <<= 1;
	}
}


void mark_bit_zero(unsigned char *addr, int position)
{
	const unsigned char ANDmasks[]={254,253,251,247,239,223,191,127};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);
	
	*(addr+byte_offset) &= ANDmasks[bit_offset];
}


void mark_bit_one(unsigned char *addr, int position)
{
	const unsigned char ORmasks[] ={1,2,4,8,16,32,64,128};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);

	*(addr+byte_offset) |= ORmasks[bit_offset];
}


// returns full int assumes we are at position 31 in the 4 bytes
int get_matrix_32_bits(unsigned char *row_start, int position)
{
	int cell_off = (position >> 5);
	return (*((int *)row_start+cell_off)); 
//	return 1;
}

int get_matrix_val(unsigned char *row_start, int position)
{
	const unsigned char masks[] ={1,2,4,8,16,32,64,128};
	int bit_offset = (position & 0X7);
	int byte_offset = (position >> 3);
	
	return ((*(row_start+byte_offset) & masks[bit_offset]));
//	return 1;
}


// global debugging variable
int wrongly_filtered_spectra;

/***************************************************************************
	This function creates clusters from a list of files containing spectra
	(possibly different file types).
	The cluster spectra are outputted as mgf files in the output dir (x spectra
	per file). In addition, for each cluster file there is a map file that holds
	the indices (position in list, and idx in file) of the original spectra
	that are part of the cluster.
****************************************************************************/
void cluster_full_dataset(Config *config,
							  char *list_file,
							  const string& out_dir,
							  const string& clust_name,
							  int batch_idx,
							  int specs_per_slice,
							  mass_t min_m_over_z,
							  mass_t max_m_over_z,
							  float  min_similarity, 
							  int min_cluster_size,
							  int max_cluster_size,
							  bool verbose,
							  int  max_small_cluster_size,
							  int  k_value,
							  void *pmcsqs,
							  bool ind_pkl_mode,
							  float filter_prob,
							  bool ind_assign_charges,
							  int  max_mzxml_idx,
							  char *good_anns_file)
{
	const mass_t pm_tolerance = config->get_pm_tolerance();
	const mass_t double_pm_tolerance = pm_tolerance * 2.0;
	const mass_t additional_pm_tolerance = pm_tolerance * 2.5;
	const mass_t tolerance = (config->get_tolerance() <= 0.01) ? config->get_tolerance() :
							  config->get_tolerance() * 0.75;
	FileManager fm;
	FileSet all_spec_fs;
	vector<QCPeak> basic_peaks;      // bulk allocation

	// for cluster histograms
	const int clust_vals[]={1,2,5,10,20,50,100,200,500};
	const int  num_clust_vals = sizeof(clust_vals)/sizeof(int);
	vector<int> clust_counts;

	clust_counts.resize(num_clust_vals+1,0);

	double avg_cluster_size=0;
	int num_clusters=0;
	int total_spectra_in_clusters = 0;
	int spec_idx=0;
	int next_cluster_idx = 0;

	wrongly_filtered_spectra = 0;

	QCOutputter qco;
	qco.init(clust_name,out_dir, batch_idx, min_m_over_z, max_m_over_z,
			 min_similarity, min_cluster_size);

	ostringstream oss;
	oss << batch_idx;
	string batch_str = oss.str();
	string last_good_mass_name = out_dir + "/" + clust_name + "_" + batch_str + "_lastmass.txt";

	ClusterSpectrum::init_statics(config);

//	ClusterSpectrum::set_num_top_peaks_per_1000_da(k_value);

	fm.init_from_list_file(config,list_file,min_m_over_z,max_m_over_z);

	all_spec_fs.select_all_files(fm);

	if (max_mzxml_idx>=0)
	{
		const int org_num = all_spec_fs.get_total_spectra();
		cout << "Filtering for max mzxml idx " << max_mzxml_idx << endl;
		all_spec_fs.filter_dat_spectra_by_mzxml_idx(max_mzxml_idx);
		cout << "#spectra went from " << org_num << " to " << all_spec_fs.get_total_spectra() << endl;
	}

	all_spec_fs.sort_according_to_m_over_z();
	const int total_spectra = all_spec_fs.get_total_spectra();
	const vector<SingleSpectrumFile *>& all_ssf = all_spec_fs.get_ssf_pointers();

	if (all_ssf.size()==0)
	{
		cout <<"Warning: no files were selected for clusering in mass range: " <<
			min_m_over_z << " - " << max_m_over_z << endl;
		return;
	}




	// used for debugging purposes
	map<mzXML_annotation,int> ann_map;
	map<mzXML_annotation,int> *ann_map_ptr=NULL;
	if (good_anns_file)
	{
//		read_mzXML_annotations_to_map(good_anns_file,ann_map);
		ann_map_ptr = &ann_map;
	}
	
	int total_spectra_read=0;
	int total_mismatched = 0;

	// set the sizes for the static arrays of peaks, sim matrix and determine
	// the maximal slice size for the clustering
	// vlaues are set according to defaults of:
	// 4M spectra
	// 4M peaks
	// Slice of 20000 spectra per iteration
	// Slice width 3 Da.
	

	float multiplier = (float)specs_per_slice / 25000.0;

	int total_spectra_to_cluster = all_ssf.size();
	int num_total_peaks = (int)(4000000 * multiplier); 
	int slice_size      =   specs_per_slice;
	int sim_matrix_size =   ((specs_per_slice+7) / 8);
	int sim_matrix_bytes = sim_matrix_size * sim_matrix_size * 4 + specs_per_slice;



	cout << "Going to cluster " << all_ssf.size() << " spectra in this run." << endl;
	cout << "Need to allocate following memory: " << endl;
	cout << "peaks            " << right << num_total_peaks << endl;
	cout << "spectrum headers " << right << slice_size << endl;
	cout << "sim matrix       " << right << sim_matrix_bytes << " (bytes) " << endl;
	cout << "Using " << k_value << " peaks per 1000 Da for similarity computations." << endl;

	if (pmcsqs)
		cout << "Filtering with quality threshold: " << filter_prob << endl;

	cout << endl;

	// adjusting the quality filter probability (to improve the default behaviour)
	// singletons should generally have a higher prob for accepting them 
	// this setting is used only for large clustering jobs

	float adjusted_prob_for_min_size_clusters = filter_prob;

	const mass_t mz_range = all_ssf[all_ssf.size()-1]->m_over_z - all_ssf[0]->m_over_z;
	if (mz_range>0)
	{
		float spectra_density =  (all_ssf.size()/mz_range);
		if (spectra_density>1200.0)
		{
			float mult_val = (spectra_density/1600.0)*3.0;
			if (mult_val>3)
				mult_val = 3;

			adjusted_prob_for_min_size_clusters = (filter_prob * mult_val);
					
			if (adjusted_prob_for_min_size_clusters > 0.3)
				adjusted_prob_for_min_size_clusters = 0.3;

			if (adjusted_prob_for_min_size_clusters<filter_prob)
				adjusted_prob_for_min_size_clusters=filter_prob;
		}
	}
	

	const float filter_prob_for_min_size_clusters = adjusted_prob_for_min_size_clusters;
	const float filter_prob_for_larger_than_min_clusters = filter_prob;	


	basic_peaks.resize(num_total_peaks+20000);

	sim_matrix = new unsigned char [sim_matrix_bytes];
	max_sim_addr = sim_matrix + sim_matrix_bytes;

	if (! sim_matrix)
	{
		cout << "Error: couldn't allocate memory for sim matrix!" << endl;
		exit(1);
	}

	double total_sims =0;

	// collect spectra into two ssf vectors
	while (spec_idx<total_spectra)
	{
		const mass_t max_m_over_z = all_ssf[spec_idx]->m_over_z + double_pm_tolerance;
		const mass_t additional_m_over_z = all_ssf[spec_idx]->m_over_z + additional_pm_tolerance;
		int num_used_peaks=0;
		int num_peaks_first_stage = -1;
		int start_idx = spec_idx;
		int end_stage_one_idx = -1;
		int end_additional_idx = -1;

		int start_spec_idx = spec_idx;
		bool add_additional_spectra = false; // flag if spectra from the next mass range
											 // should be added at the later stage
		int i;

		// add spectra wto the 
		int max_spec_idx = spec_idx + slice_size;
		if (max_spec_idx> total_spectra)
			max_spec_idx = total_spectra;

		while (spec_idx < max_spec_idx && 
			   num_used_peaks<num_total_peaks && 
			   all_ssf[spec_idx]->m_over_z<max_m_over_z)
					num_used_peaks+=all_ssf[spec_idx++]->num_peaks;

		end_stage_one_idx = spec_idx-1;
		num_peaks_first_stage = num_used_peaks;

		if (spec_idx < max_spec_idx)
		{
			add_additional_spectra = true;

			while (spec_idx < max_spec_idx && 
			   num_used_peaks<num_total_peaks && 
			   all_ssf[spec_idx]->m_over_z<additional_m_over_z)
					num_used_peaks+=all_ssf[spec_idx++]->num_peaks;

			end_additional_idx = spec_idx-1;
		}

		FileSet cluster_fs, additional_fs;
		cluster_fs.init_from_another_fs(all_spec_fs,start_idx,end_stage_one_idx);

		if (add_additional_spectra)
			additional_fs.init_from_another_fs(all_spec_fs,end_stage_one_idx+1,end_additional_idx);

		vector<ClusterSpectrum> clusters;
		clusters.clear();

		cout << fixed << setprecision(3) << "Clustering: " << all_ssf[start_idx]->m_over_z << " - " << 
			(add_additional_spectra ? all_ssf[end_additional_idx]->m_over_z : 
					all_ssf[end_stage_one_idx]->m_over_z )  << "  (" <<
			spec_idx - start_spec_idx  << "  spectra,  " << num_used_peaks << " peaks)" << endl;

		int num_in_clusters=0;
		total_spectra_read+= cluster_spec_in_file_set(
										config, 
										fm, 
										cluster_fs, 
										tolerance,
										&basic_peaks[0], 
										clusters, 
										min_similarity, 
										max_small_cluster_size,
										max_cluster_size,
										k_value,
										false,
										pmcsqs,
										filter_prob_for_larger_than_min_clusters,
										ind_assign_charges,
										ann_map_ptr); 


		// join singletons from the next half of the clustering window

		if (add_additional_spectra && additional_fs.get_total_spectra()>0) 
		{
			int num_added =  add_additional_spectra_to_existing_clusters(config,fm,additional_fs,tolerance,
				&basic_peaks[num_peaks_first_stage], clusters, min_similarity,max_cluster_size,
				pmcsqs,filter_prob_for_larger_than_min_clusters, ind_assign_charges, ann_map_ptr, false);

			total_spectra_read += num_added;
		} 

		// update cluster info
		for (i=0; i<clusters.size(); i++)
		{
			ClusterSpectrum& cluster = clusters[i];
			if (cluster.get_tmp_cluster_idx()<0)
				continue;

			if	(cluster.get_num_basic_spectra()<min_cluster_size)
				continue;

			// check if sqs is high enough
			if (pmcsqs && cluster.get_num_basic_spectra() == 1 &&
				cluster.get_basic_spectrum(0).ssf->sqs < filter_prob_for_min_size_clusters)
				continue;

			total_spectra_in_clusters += cluster.get_num_basic_spectra();

//			cluster.set_charge();
//			cluster.set_cluster_m_over_z();
	
		/*	if (ind_pkl_mode)
			{
				qco.output_cluster_spectrum_as_single_pkl(cluster);
			}
			else
				qco.output_cluster_spectrum(cluster,ind_assign_charges);*/

			if (verbose)
			{
				cout << num_clusters<< " " << cluster.get_num_basic_spectra() << endl;
//				cluster.print_cluster_similarities();
				cout << endl;
			}

			const int num_spec_in_cluster = cluster.get_num_basic_spectra();
			if (num_spec_in_cluster>=min_cluster_size)
			{
				num_clusters++;
				avg_cluster_size += num_spec_in_cluster;
			}

			// add counts to histogram
			int j;
			for (j=0; j<num_clust_vals; j++)
				if (num_spec_in_cluster<= clust_vals[j])
					break;
			clust_counts[j]++;

			// check for mixed clusters
//			int n_mis = clusters[i].get_num_misassigned_spectra();
//			total_mismatched += n_mis;
//			if (n_mis>0)
//				clusters[i].print_cluster_peptides();
		}
		spec_idx = end_stage_one_idx + 1; // go back to the end of stage one

		// update the file which holds the last mass clustered
		ofstream last_mass_stream(last_good_mass_name.c_str(),ios::out);
		last_mass_stream << fixed << setprecision(3) << max_m_over_z << endl;
		last_mass_stream.close();
	}

	
	if (sim_matrix)
		delete [] sim_matrix;

	cout << endl << endl << "Total spectra read and clustered: " << total_spectra_read << " (" <<
		total_spectra << ")" << endl;

	cout << "# spectra in clusters: " << total_spectra_in_clusters << endl;
	cout << "% spectra in clusters: " << setprecision(3) << 
		(double)total_spectra_in_clusters / (double)total_spectra << endl;

	cout << "# clusters: " << num_clusters << "  , " << "Avg cluster size: " << 
		avg_cluster_size/(double)num_clusters << endl;

	if (total_mismatched>0)
		cout << "Total mismatched spectra: " << total_mismatched << "  (" <<
			(double)total_mismatched/total_spectra_read << ")" << endl;

	// cluster histogram
	cout << "Histogram of clusters: " << endl;
	cout << "max size     count" << endl;
	int i;
	for (i=0; i<num_clust_vals; i++)
	{
		cout << setw(8) << left << clust_vals[i] << clust_counts[i] << endl;
	}
	cout << ">" << setw(7) << left << clust_vals[num_clust_vals-1] <<
			clust_counts[num_clust_vals] << endl;

	double ratio = (double)total_spectra_in_clusters / (double)total_spectra;
	if ((filter_prob_for_min_size_clusters<0.2 && ratio<0.45) ||
		 (filter_prob_for_min_size_clusters>0.2 &&  ratio < 0.4))
	{
		cout << endl << "WARNING: only " << fixed << setprecision(2) << ratio*100 << "%" <<
			" of the spectra found there way into clusters. This might lead to loss of idnetifications." << endl;
		cout << "If you are filtering, you might consider reducing the quality threshold to a lower value,";
		cout << "e.g., by using the flag \"-min_filter_prob " << filter_prob*0.5 << endl;
	}

}





/**********************************************************************
	Creates cluster spec for the set of basic spectra.
	First reads the spectra and copies the peaks into the bulk peak allocation
	returns number of spectra actually read (does not read spectra that
	were already assigned to a cluster).

	The clustering is done in two phases. First a tight distance threshold
	is implemented, and in the second phase it is relaxed (this way the clusters
	should be more homogeneous).
***********************************************************************/
int cluster_spec_in_file_set(Config *config, 
							 const FileManager& fm, 
							 FileSet& cluster_fs,
							 mass_t tolerance,
							 QCPeak *basic_peaks, 
							 vector<ClusterSpectrum>& clusters, 
							 float	min_similarity,
							 int	max_small_cluster_size,
							 int	max_cluster_size,
							 int	num_top_peaks_per_1000_da,
							 bool	verbose,
							 void	*pmcsqs_ptr,
							 float	filter_prob,
							 bool	ind_assign_charges,
							 map<mzXML_annotation,int> *ann_map_ptr)
{
	// set clustering similarity thresholds
	vector<float> similarity_vals;
	int   num_rounds;
	if (min_similarity >= 0.9)
	{
		similarity_vals.push_back(min_similarity);
		num_rounds=1;
	}
	else
	{
		similarity_vals.push_back(0.9);
		if (min_similarity>=0.8)
		{
			similarity_vals.push_back(min_similarity);
			num_rounds=2;
		}
		else
		{
			similarity_vals.push_back((min_similarity+0.9)/2.0);
			similarity_vals.push_back(min_similarity);
			num_rounds=3;
		}
	}
	const float min_similarity_thresh    = (min_similarity <0.2 ? min_similarity : 0.2); // don't test similarity if a previously recored																				 
																						// similarity between clusters is less than this value

	PMCSQS_Scorer * pmcsqs = (PMCSQS_Scorer *)pmcsqs_ptr;

	BasicSpecReader bsr;
	const int num_spectra = cluster_fs.get_total_spectra();
	vector<SingleSpectrumFile *>& all_ssf = cluster_fs.get_non_const_ssf_pointers();
	static vector<BasicSpectrum> basic_spectra;
	int i;

	// set max_small_cluster_size
	if (max_small_cluster_size<0)
		max_small_cluster_size = 10 + (int)(0.8*log((float)all_ssf.size()));

	if (max_small_cluster_size > max_cluster_size)
		max_small_cluster_size = max_cluster_size;

	// read all the basic spectra into a central spectra repository
	int total_peaks_read=0;
	mass_t min_m_over_z = 1E7;
	mass_t max_m_over_z = 0;

	int num_filtered=0;

	basic_spectra.clear();
	basic_spectra.reserve(num_spectra);
	for (i=0; i<num_spectra; i++)
	{
		if (! all_ssf[i]) // invalidated
			continue;

		if (! config->get_use_spectrum_charge())
		{
			if (i>0 && all_ssf[i-1] &&
				all_ssf[i]->get_scan()  == all_ssf[i-1]->get_scan() && 
				all_ssf[i-1]->num_peaks == all_ssf[i]->num_peaks)
				continue;
		}

		const int num_spec_peaks = bsr.read_basic_spec(config,fm,all_ssf[i], basic_peaks + total_peaks_read);
		if (num_spec_peaks<5)
		{
			all_ssf[i]=NULL;
			continue;
		}

		BasicSpectrum bs;
		bs.num_peaks = num_spec_peaks;
		bs.peaks = basic_peaks + total_peaks_read;
		bs.ssf = all_ssf[i];

		if (pmcsqs && bs.ssf->sqs<0)
		{
			int max_charge=0;
			const float prob = pmcsqs->get_sqs_for_spectrum(config,bs,&max_charge);
			if (prob<filter_prob)
			{
				if (ann_map_ptr)
				{
					DAT_single *dat_ssf = (DAT_single *)bs.ssf;

					map<mzXML_annotation,int>::const_iterator it;
					mzXML_annotation ann_pos;
					ann_pos.mzXML_file_idx = dat_ssf->mzxml_file_idx;
					ann_pos.scan = dat_ssf->scan_number;

					it = ann_map_ptr->find(ann_pos);
					if (it != ann_map_ptr->end())
					{
						
						cout << ">> " << ++wrongly_filtered_spectra << "\t" <<
						dat_ssf->mzxml_file_idx << " " << dat_ssf->scan_number << 
						" --> " << prob << endl;
					}
				}
				all_ssf[i]=NULL; // this specturm was filtered!
				num_filtered++;
				continue;
			}

			// update m/z and charge state (yes it is supposed to be const...)	
			if (ind_assign_charges)
				bs.ssf->charge = max_charge;
			bs.ssf->sqs = prob;
		}
		

		basic_spectra.push_back(bs);

		total_peaks_read += num_spec_peaks;

		mass_t& m_over_z = bs.ssf->m_over_z;
		if (m_over_z<min_m_over_z)
			min_m_over_z = m_over_z;
		if (m_over_z>max_m_over_z)
			max_m_over_z = m_over_z;
	}


	// First stage, compare the basic spectra with clusters
	// Use high similarity threshold
	// If no cluster is found, create a new clusters for the spectrum
	// the calculated simlarities are stored is the sim matrix and can be used
	// in later stages to detect the need to re test the similarity

	const float first_stage_sim = similarity_vals[0];
	
	unsigned char * start_pos = sim_matrix;

	static vector<int> idx_permutations;
	idx_permutations.resize(basic_spectra.size());
	for (i=0; i<basic_spectra.size(); i++)
		idx_permutations[i]=i;

	permute_vector(idx_permutations);

	static vector<int> spec_top_idxs;
	static float top_x_masses[NUM_TOP_CLUSTER_PEAKS];
	for (i=0; i<basic_spectra.size(); i++)
	{
		const int spec_idx = idx_permutations[i];
		BasicSpectrum& spec = basic_spectra[spec_idx];
		const int spec_charge = spec.ssf->charge;
	
		set_adjusted_inten(spec.peaks,spec.num_peaks);
	//	select_top_peak_idxs(spec.peaks,spec.num_peaks,spec.ssf->m_over_z,
	//		tolerance, spec_top_idxs, top_x_masses, 
	//		ClusterSpectrum::get_num_top_peaks_per_1000_da(), config);

		// compare to previous clusters
		
		int j;
		for (j=0; j<clusters.size(); j++)
		{ 
			ClusterSpectrum& curr_clust = clusters[j];
	
			if (! curr_clust.find_match_in_top_masses(top_x_masses))
			{
				mark_bit_zero(start_pos,j);
				continue;
			}

			if (curr_clust.get_num_basic_spectra()>= max_cluster_size)
				continue;

			if (config->get_use_spectrum_charge() && 
				curr_clust.get_charge()>0 && 
				spec_charge>0 && 
				curr_clust.get_charge() !=spec_charge)
			{
				mark_bit_zero(start_pos,j);
				continue;
			}

			const float sim = calc_selected_dot_prod(tolerance,
				spec.peaks,spec.num_peaks, spec_top_idxs,
				curr_clust.get_peaks_pointer(),curr_clust.get_num_peaks(), 
				curr_clust.get_top_ranked_idxs(),verbose);

			if (sim >= min_similarity_thresh)
			{
				mark_bit_one(start_pos,j);
			}
			else
				mark_bit_zero(start_pos,j);
			

			// add this spectrum to an existing cluster
			if (sim > first_stage_sim)
			{
//				curr_clust.add_spectrum_to_cluster(spec, spec_top_idxs,top_x_masses);
				break;
			}
		}

		if (j<clusters.size())  // we added the spectrum to an existing cluster
			continue;

	
		// create new cluster from spectrum
		
		clusters.resize(clusters.size()+1);
		ClusterSpectrum& cs = clusters[clusters.size()-1];

//		cs.create_new_cluster(config, spec, clusters.size()-1);
//		cs.set_charge(spec.ssf->charge);
		cs.set_top_ranked_idxs(spec_top_idxs);
		cs.set_top_masses(top_x_masses);
		cs.set_sim_matrix_row_start(start_pos);
		
		// round off the start position to the next byte
		unsigned char *old = start_pos;
		start_pos += ((j+7) >> 3);
	}



	// second stage try joining clusters
	// first start with the joining larger clusters
	// use lower threshold

	int round;
	for (round=1; round<num_rounds; round++)
	{
		const float round_similarity = similarity_vals[round];
		const float tighter_similarity = 1.0 - (1.0 - similarity_vals[round])/2.0;

		// join larger clusters, use tighter similarity
		for (i=clusters.size()-1; i>0; i--)
		{
			if (clusters[i].get_tmp_cluster_idx()<0)
				continue;

			if (clusters[i].get_num_basic_spectra()>=max_cluster_size)
				continue;

			unsigned char *sim_row_start = clusters[i].get_sim_matrix_row_start();
			const int num_spec_i = clusters[i].get_num_basic_spectra();
			const int clust_charge = clusters[i].get_charge();

			int j;
			for (j=i-1; j>=0; j--)
			{
				const int size_sum = clusters[j].get_num_basic_spectra() + num_spec_i;
				if (clusters[j].get_tmp_cluster_idx()<0 ||
					size_sum <= max_small_cluster_size ||
					size_sum >= max_cluster_size) 
					continue;

				// skip 32 places if the matirx is all zeros in that area
				if ((j % 32 == 31) &&  ! get_matrix_32_bits(sim_row_start,j))
				{
					j-=31;
					continue;
				}

				if (! get_matrix_val(sim_row_start,j))
					continue;

				if (config->get_use_spectrum_charge() && clust_charge>0)
				{
					const int other_clust_charge = clusters[j].get_charge();
					if (other_clust_charge>0 && clust_charge != other_clust_charge)
					{
						mark_bit_zero(sim_row_start,j);
						continue;
					}
				}

				const float sim = calc_selected_dot_prod(tolerance,
					clusters[j].get_peaks_pointer(),clusters[j].get_num_peaks(), 
					clusters[j].get_top_ranked_idxs(),
					clusters[i].get_peaks_pointer(),clusters[i].get_num_peaks(), 
					clusters[i].get_top_ranked_idxs());


				if (sim >= min_similarity_thresh)
				{
					mark_bit_one(sim_row_start,j);
				}
				else
					mark_bit_zero(sim_row_start,j);
				

				if (sim > tighter_similarity)
				{
			//		if (clusters[j].add_cluster(clusters[i], tighter_similarity))
			//		{
			//			clusters[i].set_tmp_cluster_idx(-1);
			//			break;
			//		}
				}
			}
		} 


		// join smaller clusters, use the round similarity
		for (i=clusters.size()-1; i>0; i--)
		{
			if (clusters[i].get_tmp_cluster_idx()<0 || 
				clusters[i].get_num_basic_spectra() >= max_small_cluster_size)
				continue;

			unsigned char * sim_row_start = clusters[i].get_sim_matrix_row_start();
			const int clust_charge = clusters[i].get_charge();
			const int num_spec_i = clusters[i].get_num_basic_spectra();
			int j;
			for (j=i-1; j>=0; j--)
			{
				if (clusters[j].get_tmp_cluster_idx()<0 ||
					num_spec_i + clusters[j].get_num_basic_spectra() >= max_small_cluster_size)
					continue;

				// skip 32 places if the matirx is all zeros in that area
				if ((j % 32 == 31) &&  ! get_matrix_32_bits(sim_row_start,j))
				{
					j-=31;
					continue;
				}

				if ( ! get_matrix_val(sim_row_start,j))
					continue;

				if (clusters[j].get_num_basic_spectra()>max_cluster_size)
					continue;

				if (config->get_use_spectrum_charge() && clust_charge>0)
				{
					const int other_clust_charge = clusters[j].get_charge();
					if (other_clust_charge>0 && clust_charge != other_clust_charge)
					{
						mark_bit_zero(sim_row_start,j);
						continue;
					}
				}

				const float sim = calc_selected_dot_prod(tolerance,
					clusters[j].get_peaks_pointer(),clusters[j].get_num_peaks(), 
					clusters[j].get_top_ranked_idxs(),
					clusters[i].get_peaks_pointer(),clusters[i].get_num_peaks(), 
					clusters[i].get_top_ranked_idxs());


				if (sim >= min_similarity_thresh)
				{
					mark_bit_one(sim_row_start,j);
				}
				else
					mark_bit_zero(sim_row_start,j);
				
				if (sim > round_similarity)
				{
//					if (clusters[j].add_cluster(clusters[i],round_similarity))
//					{
//						clusters[i].set_tmp_cluster_idx(-1);
//						break;
//					}
				}
			}
		}
	} 



//	cout << "Filtered " << num_filtered << "/" << num_spectra << endl;

	if (verbose)
	{
		for (i=0; i<clusters.size(); i++)
			if (clusters[i].get_tmp_cluster_idx() >=0 &&
				clusters[i].get_basic_spectra().size() == 1 && 
				clusters[i].get_basic_spectra()[0].ssf->peptide.get_num_aas()>0)
				cout << clusters[i].get_basic_spectra()[0].ssf->peptide.as_string(config) << endl;
	}


	
	return num_spectra;
}



/**********************************************************************
	Adds spectra from the additional set to the existing clusters.
	If they are added they are invalidated from further clusetering.
***********************************************************************/
int add_additional_spectra_to_existing_clusters(
							Config *config, 
							const FileManager& fm, 
							FileSet& additional_fs, 
							mass_t tolerance, 
							QCPeak *basic_peaks, 
							vector<ClusterSpectrum>& clusters, 
							float min_similarity, 
							int max_cluster_size,
							void *pmcsqs_ptr,
							float filter_prob,
							bool  ind_assign_charges,
							map<mzXML_annotation,int> *ann_map_ptr,
							bool verbose)
{
	float spectrum_join_similarity = 0.8;
	if (min_similarity>spectrum_join_similarity)
		spectrum_join_similarity = min_similarity;

	// read spectra
	BasicSpecReader bsr;
	const int num_spectra = additional_fs.get_total_spectra();
	vector<SingleSpectrumFile *>& all_ssf = additional_fs.get_non_const_ssf_pointers();
	
	static vector<BasicSpectrum> basic_spectra;

	basic_spectra.clear();
	basic_spectra.reserve(num_spectra);

	PMCSQS_Scorer * pmcsqs = (PMCSQS_Scorer *)pmcsqs_ptr;

	int i,total_peaks_read=0;

	mass_t min_m_over_z = 1E7;
	mass_t max_m_over_z = 0;

	for (i=0; i<num_spectra; i++)
	{
		if (! all_ssf[i] || all_ssf[i]->assigned_cluster>=0)
			continue;

		const int num_spec_peaks = bsr.read_basic_spec(config,fm,all_ssf[i], basic_peaks + total_peaks_read);
		
		if (num_spec_peaks<5)
		{
			all_ssf[i]=NULL;
			continue;
		}

		BasicSpectrum bs;
		bs.num_peaks = num_spec_peaks;
		bs.peaks = basic_peaks + total_peaks_read;
		bs.ssf = all_ssf[i];

		if (pmcsqs && bs.ssf->sqs<0)
		{
			int max_charge=0;
			const float prob = pmcsqs->get_sqs_for_spectrum(config,bs,&max_charge);
			if (prob<filter_prob)
			{
				if (ann_map_ptr)
				{
					DAT_single *dat_ssf = (DAT_single *)bs.ssf;

					map<mzXML_annotation,int>::const_iterator it;
					mzXML_annotation ann_pos;
					ann_pos.mzXML_file_idx = dat_ssf->mzxml_file_idx;
					ann_pos.scan = dat_ssf->scan_number;

					it = ann_map_ptr->find(ann_pos);
					if (it != ann_map_ptr->end())
					{
						
						cout << ">> " << ++wrongly_filtered_spectra << "\t" <<
						dat_ssf->mzxml_file_idx << " " << dat_ssf->scan_number << 
						" --> " << prob << endl;
					}
				}
				all_ssf[i]=NULL; // this specturm was filtered!
				continue; 
			}

			// update m/z and charge state (yes it is supposed to be const...)	
			if (ind_assign_charges)
				bs.ssf->charge = max_charge;
			//ssf->m_over_z = res[max_charge].mz1;
			bs.ssf->sqs = prob;
		}

		basic_spectra.push_back(bs);

		total_peaks_read += num_spec_peaks;

		mass_t& m_over_z = bs.ssf->m_over_z;
		if (m_over_z<min_m_over_z)
			min_m_over_z = m_over_z;
		if (m_over_z>max_m_over_z)
			max_m_over_z = m_over_z;
	}


	static vector<int> idx_permutations;
	idx_permutations.resize(basic_spectra.size());
	for (i=0; i<basic_spectra.size(); i++)
		idx_permutations[i]=i;

	permute_vector(idx_permutations);

	// cluster the spectra
	static float top_x_masses[NUM_TOP_CLUSTER_PEAKS];
	static vector<int> spec_top_idxs;
	int num_added=0;
	for (i=0; i<basic_spectra.size(); i++)
	{
		const int spec_idx = idx_permutations[i];
		BasicSpectrum& spec = basic_spectra[spec_idx];
		const float spec_retention_time = spec.ssf->retention_time;
	
		set_adjusted_inten(spec.peaks,spec.num_peaks);
	//	select_top_peak_idxs(spec.peaks,spec.num_peaks,spec.ssf->m_over_z,
	//		tolerance,spec_top_idxs, top_x_masses, 
	//		ClusterSpectrum::get_num_top_peaks_per_1000_da(),
	//		config);

		// compare to previous clusters
		int j;
		for (j=0; j<clusters.size(); j++)
		{
			if (clusters[j].get_tmp_cluster_idx() < 0)
				continue;

			if (! clusters[j].find_match_in_top_masses(top_x_masses))
				continue;

			if (clusters[j].get_num_basic_spectra()>=max_cluster_size)
				continue;

			const float sim = calc_selected_dot_prod(tolerance,
				spec.peaks,spec.num_peaks, spec_top_idxs,
				clusters[j].get_peaks_pointer(),clusters[j].get_num_peaks(), 
				clusters[j].get_top_ranked_idxs());

			if (sim > spectrum_join_similarity)
			{
//				clusters[j].add_spectrum_to_cluster(spec, spec_top_idxs, top_x_masses);
				num_added++;
				break;
			}
		}
	}

	return num_added;
}



void check_fs_for_missing_anns_and_test_sqs(Config *config,
											const vector<SingleSpectrumFile *>& ssfs,
											const FileManager& fm,
											void *pmcsqs_ptr,
											char *anns_file)
{
	const mass_t tolerance = config->get_tolerance();

	map<mzXML_annotation,int> ann_map, ssf_map;
//	read_mzXML_annotations_to_map(anns_file,ann_map);
	
	int i;
	for (i=0; i<ssfs.size(); i++)
	{
		DAT_single *dat_ssf = (DAT_single *)ssfs[i];
		mzXML_annotation ann;
		ann.charge = 0;
		ann.scan = dat_ssf->scan_number;
		ann.mzXML_file_idx = dat_ssf->mzxml_file_idx;


		ssf_map.insert(make_pair(ann,i));
	}

	// check if anns are in map
	int num_read=0;
	int num_not_read=0;
	vector<SingleSpectrumFile *> miss_ssfs;

	map<mzXML_annotation,int>::const_iterator it;
	for (it = ann_map.begin(); it != ann_map.end(); it++)
	{
		map<mzXML_annotation,int>::const_iterator ssf_it;
		mzXML_annotation ann_pos;
		ann_pos.mzXML_file_idx = it->first.mzXML_file_idx;
		ann_pos.scan = it->first.scan;

		ssf_it = ssf_map.find(ann_pos);
		if (ssf_it != ssf_map.end())
		{
			num_read++;
			miss_ssfs.push_back(ssfs[ssf_it->second]);
			ssfs[ssf_it->second]->peptide.parse_from_string(config,it->first.pep);
		}
		else
			num_not_read++;
	}
				

	cout << "NUM anns read in SSFs:    " << num_read << endl;
	cout << "NUM in the twilight zone: " << num_not_read << endl;

	PMCSQS_Scorer * pmcsqs = (PMCSQS_Scorer *)pmcsqs_ptr;

	BasicSpecReader bsr;
	vector<QCPeak> peaks;
	peaks.resize(6000);

	QCOutputter qco;
	qco.init("missed_xxx","out",0);

	int num_with_little_peaks=0;
	for (i=0; i<miss_ssfs.size(); i++)
	{
		const int num_spec_peaks = bsr.read_basic_spec(config,fm,miss_ssfs[i],
											 &peaks[0]);
		if (num_spec_peaks<3)
		{
			num_with_little_peaks++;
			continue;
		}

		BasicSpectrum bs;
		
		bs.num_peaks = num_spec_peaks;
		bs.peaks = &peaks[0];
		bs.ssf = miss_ssfs[i];


		if (pmcsqs && bs.ssf->sqs<0)
		{
			static QCPeak *tmp_peaks = NULL;
			static int num_tmp_peaks = 0;
			if (! tmp_peaks || bs.num_peaks>num_tmp_peaks)
			{
				if (tmp_peaks)
					delete [] tmp_peaks;
				num_tmp_peaks = (int)(bs.num_peaks * 1.5 + 1);
				tmp_peaks = new QCPeak[num_tmp_peaks];
			}

			// need to create a special peak list that passes filtering
			QCPeak *org_peaks = bs.peaks;
			int    num_org_peaks = bs.num_peaks;
			BasicSpectrum sqs_bs = bs; 
			sqs_bs.peaks = tmp_peaks;
			sqs_bs.ssf   = bs.ssf;

			const mass_t join_tolerance = (tolerance < 0.05 ? tolerance : 0.6 * tolerance);
			int p_idx=0;
			int j=1;
			while (j<num_org_peaks)
			{
				if (org_peaks[j].mass - org_peaks[p_idx].mass<=join_tolerance)
				{
					intensity_t inten_sum = org_peaks[i].intensity + org_peaks[p_idx].intensity;
					mass_t new_mass = (org_peaks[j].intensity * org_peaks[j].mass + 
									   org_peaks[p_idx].intensity * org_peaks[p_idx].mass ) / inten_sum;

					org_peaks[p_idx].mass = new_mass;
					org_peaks[p_idx].intensity = inten_sum;	
				}
				else
				{
					org_peaks[++p_idx]=org_peaks[j];
				}
				j++;
			}
			num_org_peaks = p_idx+1;

			vector<bool> indicators;
			mark_top_peaks_with_sliding_window(bs.peaks, 
											   bs.num_peaks,
											   config->get_local_window_size(),
											   config->get_max_number_peaks_per_local_window(),
											   indicators);
		
			if (bs.ssf->precursor_intensity <=0)
			{
				float inten=0;
				int j;
				for (j=0; j<num_org_peaks; i++)
					inten+=org_peaks[j].intensity;
				bs.ssf->precursor_intensity=inten;
			}

			p_idx=0;
			for (j=0; j<num_org_peaks; j++)
				if (indicators[j])
					tmp_peaks[p_idx++]=org_peaks[j];
			sqs_bs.num_peaks = p_idx;
		//	cout << num_org_peaks << "->" << sqs_bs.num_peaks << endl;

			int max_charge=0;
			const float prob = pmcsqs->get_sqs_for_spectrum(config,sqs_bs,&max_charge);
			cout << i << "\t" << setprecision(3) << prob << "\t";
			bs.ssf->type = DAT;
			bs.ssf->print_ssf_stats(config);
		}

		ClusterSpectrum cs;
//		cs.create_new_cluster(config,bs,i);
//		qco.output_cluster_spectrum(cs);
	}
	cout << "NUM WITH LITTLE PEAKS: " << num_with_little_peaks << endl;
}




