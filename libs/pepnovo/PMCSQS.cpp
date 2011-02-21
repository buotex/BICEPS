#include "PMCSQS.h"
#include "auxfun.h"


PMCRankStats::PMCRankStats(const PMCRankStats & Source)
{
	*this = Source;
//	offset_pairs_ordered_by_inten		 = Source.offset_pairs_ordered_by_inten;
//	strong_offset_pairs_ordered_by_inten = Source.strong_offset_pairs_ordered_by_inten;
//  c2_offset_pairs_ordered_by_inten	 =  Source.c2_offset_pairs_ordered_by_inten;
}


void PMCSQS_Scorer::set_sqs_mass_thresholds()
{
	sqs_mass_thresholds.clear();
	sqs_mass_thresholds.push_back(800);
	sqs_mass_thresholds.push_back(1200);
}


void PMCSQS_Scorer::set_default_sqs_correct_factors()
{
	sqs_correction_factors.resize(4);
	sqs_mult_factors.resize(4);

	sqs_correction_factors[1].push_back(0);
	sqs_correction_factors[1].push_back(-0.1);
	sqs_correction_factors[1].push_back(-0.15);

	sqs_correction_factors[2].push_back(0);
	sqs_correction_factors[2].push_back(0.05);
	sqs_correction_factors[2].push_back(0.05);

	sqs_correction_factors[2].push_back(0);
	sqs_correction_factors[2].push_back(0.05);
	sqs_correction_factors[2].push_back(0.05);

	int i;
	for (i=0; i<4; i++)
	{
		int j;
		for (j=0; j<sqs_correction_factors.size(); j++)
			sqs_mult_factors[i].push_back(1.0/(1.0+sqs_correction_factors[i][j]));
	}
}


void PMCSQS_Scorer::init_sqs_correct_factors(int max_charge, int num_sizes)
{
	sqs_correction_factors.resize(max_charge+max_charge);
	sqs_mult_factors.resize(max_charge+1);

	int i;
	for (i=1; i<=max_charge; i++)
	{
		sqs_correction_factors[i].clear();
		sqs_correction_factors[i].resize(num_sizes,0);
		sqs_mult_factors[i].clear();
		sqs_mult_factors[i].resize(num_sizes,1);
	}
}





/***************************************************************************
Gives detailed m/z and prob values for the different charges.
Returns the highest sqs probability found for a spectrum.
****************************************************************************/
float PMCSQS_Scorer::get_pmcsqs_results_for_spectrum(Config *config, 
													 const BasicSpectrum& bs,
													 vector<PmcSqsChargeRes>& res)
{
	static ME_Regression_Sample sqs_sam;

	if (! init_for_current_spec(config,bs))
		return 0; // corrupt file

	calculate_curr_spec_pmc_values(bs,bin_increment);

	const int num_charges = pmc_rank_models.size();
	const int sqs_size_idx = this->get_sqs_size_idx(bs.ssf->m_over_z);

	if (res.size()<num_charges)
		res.resize(num_charges);

	
	vector<float> min_comp_probs;		// the minimal probability when used in a comparison model
	min_comp_probs.resize(num_charges,1.0);
	min_comp_probs[0]=0;

	if (ind_initialized_sqs)
	{
		fill_fval_vector_with_SQS(bs, sqs_sam);
	
		int charge;
		for (charge=1; charge<num_charges; charge++)
		{
			
			float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

			// correct prob
			prob += sqs_correction_factors[charge][sqs_size_idx];
			prob *= sqs_mult_factors[charge][sqs_size_idx];
			if (prob<0.0)
				prob=0.0;

			res[charge].sqs_prob=prob;
		}
		
		int c;
		for (c=1; c<num_charges-1; c++)
		{
			float comp_prob = 2.0;
			int d;
			for (d=c+1; d<num_charges; d++)
			{
				if (! sqs_models[d][c][sqs_size_idx])
					continue;

				comp_prob=sqs_models[d][c][sqs_size_idx]->p_y_given_x(0,sqs_sam);

				if (comp_prob<min_comp_probs[d])
					min_comp_probs[d]=comp_prob;

				const float one_minus_prob = 1.0-comp_prob;
				if (one_minus_prob<min_comp_probs[c])
					min_comp_probs[c]=one_minus_prob;
			}
		}	
	}
	else // give the max prob charge to the input charge, the rest get 0
	{
		if (bs.ssf->charge<=0)
		{
			cout << "Error: no SQS model was read, so charge must be supplied in the spectrum!" << endl;
			exit(1);
		}

		int c;
		for (c=1; c<num_charges; c++)
			if (c != bs.ssf->charge)
				min_comp_probs[c]=0;
	}

	int charge;
	for (charge=0; charge<min_comp_probs.size(); charge++)
	{
	//	if (min_comp_probs[charge]>1.5)
	//		min_comp_probs[charge]=0;
//		cout << charge << "\t" << min_comp_probs[c] << endl;
	}


	float max_prob=0;
	for (charge=1; charge<num_charges; charge++)
	{
		

		const float prob = min_comp_probs[charge];
		if (prob>max_prob)
			max_prob=prob;

		res[charge].min_comp_prob = prob;

		if (prob>0.02)
		{
			find_best_mz_values_from_rank_model(bs,charge,config->get_pm_tolerance(),res[charge]);
		}
		else // only give one guess, the second one will come from another charge
		{
			res[charge].mz1=bs.ssf->m_over_z;
			res[charge].score1=0;
			res[charge].mz2=NEG_INF;
			res[charge].score2=NEG_INF;
		}
	}
	

	return max_prob;
}




/********************************************************************
Computes the best sqs value for the spectrum. It is faster than the charge
m/z function since it only looks at a few mass positions. 
*********************************************************************/
float PMCSQS_Scorer::get_sqs_for_spectrum(Config *config, 
										  const BasicSpectrum& bs, 
										  int *max_charge,
										  bool verbose)
{
	static ME_Regression_Sample sqs_sam;
	static vector<QCPeak> sqs_peaks;
	static int num_sqs_peaks = 0;
	BasicSpectrum sqs_bs;

	if (verbose)
	{
		cout << "ORG spectrum: " << bs.num_peaks << endl;
//		bs.print_peaks();
		cout << endl;
	}

	if (num_sqs_peaks< 2000 || num_sqs_peaks<bs.num_peaks)
	{
		num_sqs_peaks = (int)(bs.num_peaks * 1.5);
		if (num_sqs_peaks<2000)
			num_sqs_peaks = 2000;
		
		if (num_sqs_peaks>100000)
		{
			cout << "Error: too many peak in spectrum: " << bs.num_peaks << endl;
			exit(1);
		}
		sqs_peaks.resize(num_sqs_peaks);
	}

	sqs_bs.peaks = &sqs_peaks[0];
	sqs_bs.ssf = bs.ssf;

	int best_charge=0;

	// filter peaks to acheive the required peak density for the models

	int new_num_peaks =0;
	create_filtered_peak_list_for_sqs(bs.peaks,bs.num_peaks,sqs_bs.peaks,new_num_peaks);
	sqs_bs.num_peaks= new_num_peaks;

	if (verbose)
	{
		cout << "ATER FILTERING: " << new_num_peaks << endl;
//		sqs_bs.print_peaks();
		cout << endl;
	}

	if (! init_for_current_spec(config,sqs_bs))
		return 0; // corrupt spectrum

	calculate_curr_spec_pmc_values(sqs_bs,bin_increment*3);

	const int num_charges = pmc_rank_models.size();
	const int sqs_size_idx = this->get_sqs_size_idx(sqs_bs.ssf->m_over_z);

	float max_prob=-1;
	
	fill_fval_vector_with_SQS(sqs_bs, sqs_sam);
	
	int charge;
	for (charge=1; charge<num_charges; charge++)
	{
			
		float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

		// correct prob
		prob += sqs_correction_factors[charge][sqs_size_idx];
		prob *= sqs_mult_factors[charge][sqs_size_idx];
		if (prob<0.0)
			prob=0.0;

		if (prob>max_prob)
		{
			max_prob=prob;
			best_charge=charge;
		}
	}

	if (max_charge)
		*max_charge = best_charge;

	return max_prob;
}



/******************************************************************************
Selects the the two best values of charges 1,2,3
returns the max prob found
*******************************************************************************/
float PMCSQS_Scorer::get_best_mz_charge(Config *config, const BasicSpectrum& bs, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res)
{
	static vector<PmcSqsChargeRes> res;

	if (! this->ind_initialized_pmcr)
	{
		cout << "Error: no PMC model was read!" << endl;
		if (bs.ssf->charge <=0)
			cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
		exit(1);
	}

	get_pmcsqs_results_for_spectrum(config,bs,res);

	float best_prob=-1;
	int best_charge=0;

	int charge;
	for (charge =1; charge<=max_model_charge; charge++)
		if (res[charge].min_comp_prob>best_prob)
		{
			best_charge = charge;
			best_prob = res[charge].min_comp_prob;
		}

	const PmcSqsChargeRes& cr = res[best_charge];
	float second_best_prob=1.0 - (cr.score1-cr.score2)/(fabs(cr.score1)+fabs(cr.score2)+2.0);
	int   second_best_charge = best_charge;
	
	for (charge =1; charge<=max_model_charge; charge++)
		if (charge != best_charge && res[charge].min_comp_prob>second_best_prob)
		{
			second_best_charge = charge;
			second_best_prob = res[charge].min_comp_prob;
		}

	*mz1 = (mass_t)res[best_charge].mz1;
	*charge1 = best_charge;
	*prob1 = best_prob;

	if (mz2 && charge2)
	{
		*charge2 = second_best_charge;
		*prob2   = second_best_prob;
		if (second_best_charge == best_charge)
		{
			*mz2 = (mass_t)res[best_charge].mz2;
		}
		else
			*mz2 = (mass_t)res[second_best_charge].mz1;
	}

	if (all_res)
		*all_res = res;

	return best_prob;
}



/********************************************************************
Select a set of mzs nd charges that have high probabilities.
If unless the use_spectrum_mz is set to 1 in the config, it adds -1,+1  
(if there are too many mzs, not all these offfsets will be added)
Uses emprically set probability threhsolds to chose corrected pms.
*********************************************************************/
void PMCSQS_Scorer::select_pms_and_charges(Config *config, 
								const BasicSpectrum& bs,
								vector<mass_t>& pms_with_19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* all_res)
{
	static const mass_t offsets[]={-MASS_ISO,MASS_ISO};
	static const int num_all_offsets = sizeof(offsets)/sizeof(mass_t);
	static const float min_comp_prob_for_adding_second = 0.04;
	static const float min_sqs_prob_for_adding_second = 0.05;
	const mass_t half_tol = config->get_tolerance() * 0.5;

	int spec_charge = bs.ssf->charge;
	mass_t spec_mz = bs.ssf->m_over_z;
	int org_spec_charge = spec_charge;
	mass_t org_spec_mz = spec_mz;
	
	pms_with_19.clear();
	charges.clear();
	int    specific_charge=0;
	mass_t specific_mz=0;
	

	if (config->get_use_spectrum_charge())
		specific_charge=spec_charge;

	static vector<PmcSqsChargeRes> res;
	int max_charge=0;
	float max_prob = 0;
	if (spec_charge>0 && 
		config->get_use_spectrum_charge() && 
		config->get_use_spectrum_mz())
	{
		pms_with_19.push_back(spec_mz * spec_charge - (spec_charge-1)*MASS_PROTON);
		charges.push_back(spec_charge);
		max_charge=spec_charge;
	}
	else
	{
		if (! this->ind_initialized_pmcr)
		{
			cout << "Error: no PMC model was read!" << endl;
			if (bs.ssf->charge <=0)
				cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
			exit(1);
		}

		if (config->get_use_spectrum_charge() && spec_charge>= this->pmc_rank_models.size())
		{
			// adjust charge and m/z to lowest charge modeled
			spec_charge = pmc_rank_models.size()-1;
			mass_t mass_with_19 = org_spec_charge * spec_mz - (org_spec_charge-1)*MASS_PROTON;
			spec_mz = (mass_with_19 + (spec_charge -1 ) * MASS_PROTON) / spec_charge;

			bs.ssf->charge = spec_charge;
			bs.ssf->m_over_z = spec_mz;
		}
		

		get_pmcsqs_results_for_spectrum(config,bs,res);

		if (specific_charge>0)
		{
			max_charge = spec_charge;
			max_prob = res[spec_charge].min_comp_prob;
		}
		else
		{
			int c;
			for (c=1; c<res.size(); c++)
			{
				if (res[c].min_comp_prob>max_prob)
				{
					max_prob = res[c].min_comp_prob;
					max_charge = c;
				}
			}
		}

		pms_with_19.push_back(res[max_charge].mz1 * max_charge - (max_charge-1)*MASS_PROTON);
		charges.push_back(max_charge);
		if (res[max_charge].mz2>0)
		{
			pms_with_19.push_back(res[max_charge].mz2 * max_charge - (max_charge-1)*MASS_PROTON);
			charges.push_back(max_charge);
		}

	}

	// only added to the first pm_with_19, check that it doesn't overlap with others
	if (! config->get_use_spectrum_mz())
	{
		int num_offsets = num_all_offsets;
		int i;
		for (i=0; i<num_offsets; i++)
		{
			const mass_t pm_with_19 = pms_with_19[0] + offsets[i];
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge && fabs(pm_with_19 - pms_with_19[j])<half_tol)
					break;
	
			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}

		// add tol to -1 of 2nd
		if (charges[0]>=2)
		{
			const mass_t pm_with_19 = pms_with_19[1] - MASS_PROTON;
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge &&
					fabs(pm_with_19 - pms_with_19[j])<half_tol)
						break;
			
			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}
		

		if (charges[0]>=3)
		{
			const mass_t pm_with_19 = pms_with_19[1] + MASS_PROTON;
			int j;
			for (j=0; j<pms_with_19.size(); j++)
				if (charges[j]==max_charge &&
					fabs(pm_with_19 - pms_with_19[j])<half_tol)
						break;

			if (j==pms_with_19.size())
			{
				pms_with_19.push_back(pm_with_19);
				charges.push_back(max_charge);
			}
		}
	}

	// add other charges if their comp probability and sqs probs are high enough
	if (specific_charge==0)
	{
		// find best charge
		int c; 
		float max_prob=-1.0;

		for (c=1; c<res.size(); c++)
		{
			if (c==max_charge)
				continue;

			if (res[c].min_comp_prob > min_comp_prob_for_adding_second &&
				res[c].sqs_prob > min_sqs_prob_for_adding_second)
			{
				pms_with_19.push_back(res[c].mz1 * c - (c-1)*MASS_PROTON);
				charges.push_back(c);
			}
		}
	}

	if (all_res)
		*all_res = res;

	if (spec_charge != org_spec_charge)
	{
		bs.ssf->charge = org_spec_charge;
		bs.ssf->m_over_z = org_spec_mz;
	}
}



void PMCSQS_Scorer::benchmark_pm_selection(Config *config, FileManager& fm, mass_t pm_val_tol)
{
	const vector< vector< mass_t > >& threshes = config->get_size_thresholds();
	int c=1;

	for (c=1; c<threshes.size(); c++)
	{
		int s;
		for (s=0; s<threshes[c].size(); s++)
		{
			mass_t min_mz = (s>0 ? threshes[c][s-1]/c : 0);
			mass_t max_mz = threshes[c][s]/c;

			FileSet fs;
			fs.select_files_in_mz_range(fm,min_mz,max_mz,c);

			if (fs.get_total_spectra()<200)
				continue;

			cout << "CHARGE " << c <<" size " << s << "  (" << fs.get_total_spectra() << " spectra)" << endl;

			const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
			BasicSpecReader bsr;
			vector<QCPeak> peaks;
			peaks.resize(10000);

			vector<int> correct_counts;
			correct_counts.resize(8,0);

			int num_correct=0;
			int num_wrong_charge=0;
			int num_diff_charge=0;

			int i;
			for (i=0; i<all_ssfs.size(); i++)
			{
				SingleSpectrumFile *ssf=all_ssfs[i];
				BasicSpectrum bs;

				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
				bs.peaks = &peaks[0];
				bs.ssf = ssf;

				const mass_t true_mass = ssf->peptide.get_mass_with_19();

				vector<mass_t> pms_with_19;
				vector<int>    charges;
				select_pms_and_charges(config,bs,pms_with_19,charges);

				if (charges[0] != c)
					num_wrong_charge++;

				bool got_diff_charge=false;
				int j;
				for (j=1; j<charges.size(); j++)
					if (charges[j] != charges[0])
						got_diff_charge=true;

				if (got_diff_charge)
					num_diff_charge++;

				for (j=0; j<pms_with_19.size(); j++)
					if (fabs(pms_with_19[j]-true_mass)<pm_val_tol)
						break;
				
				if (j==pms_with_19.size())
					continue;
				

				num_correct++;

				if (got_diff_charge && j == pms_with_19.size()-1)
				{
					correct_counts[7]++;
				}
				else
					correct_counts[j]++;
			}

			double num_total = (double)all_ssfs.size();
			cout << "Had correct       " << fixed << setprecision(4) << num_correct/num_total << endl;
			cout << "First correct     " << correct_counts[0]/num_total << endl;
			cout << "Second correct    " << correct_counts[1]/num_total << endl;
			cout << "Off-1  correct    " << correct_counts[2]/num_total << endl;
			cout << "Off+1  correct    " << correct_counts[3]/num_total << endl;
			cout << "Off-1  2nd        " << correct_counts[4]/num_total << endl;
			cout << "Off+1  2nd        " << correct_counts[5]/num_total << endl;
			cout << "Diff Ch correct   " << correct_counts[7]/num_total << endl;
			cout << "With wrong charge " << num_wrong_charge/num_total << endl;
			cout << "With diff charge  " << num_diff_charge/num_total << endl << endl;
		}
	}
}


const int DPColumnBytes = sizeof(int)*(Val+1);

// for each peak,aa holds the idx of the previous aa if they have a mass
// diff of that aa (within tolerance)
// entry 0 in each column holds an indicator if peak is in aa diff
// entry 1 in each column holds an indicator if peak has a 
void PMCSQS_Scorer::fill_SQS_DP(const BasicSpectrum& bs, vector<DPColumn>& dp, int frag_charge ) const
{
	const QCPeak *peaks = bs.peaks;
	const int num_peaks = bs.num_peaks;
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const mass_t tag_tolerance = (config->get_tolerance()<0.15 ? 
									config->get_tolerance() : config->get_tolerance()*0.33);

	
	const mass_t mult_val = (1.0 / frag_charge);
	dp.resize(num_peaks);
	int i;
	for (i=0; i<num_peaks; i++)
	{
		dp[i].pointers[0]=0;
		int j;
		for (j=1; j<=Val; j++)
			dp[i].pointers[j]=-1;
	}

	int aa;
	for (aa=Ala; aa<=Val; aa++)
	{
		if (aa==Ile || aa==Xle)
			continue;

		const mass_t aamass = mult_val * aa2mass[aa];
		const mass_t min_offset = aamass-tag_tolerance;
		const mass_t max_offset = aamass+tag_tolerance;

		int trail_idx=0;
		int lead_idx=1;

		while (lead_idx<num_peaks)
		{
			if (curr_spec_iso_levels[lead_idx]>0)
			{
				lead_idx++;
				continue;
			}

			while (peaks[lead_idx].mass-peaks[trail_idx].mass>max_offset)
				trail_idx++;

			if (curr_spec_iso_levels[trail_idx]==0 && 
				peaks[lead_idx].mass-peaks[trail_idx].mass>min_offset)
			{
				dp[lead_idx].pointers[aa]=trail_idx;
				dp[lead_idx].pointers[0]=1;
				dp[trail_idx].pointers[0]=1;

			//	cout << "Off: " << peaks[lead_idx].mass-peaks[trail_idx].mass - min_offset<< endl;
			}
			else 
				dp[lead_idx].pointers[aa]=-1;

			lead_idx++;
		}

	}
}


/*****************************************************************************
Does the raw calculations required for processing the spectrum.
******************************************************************************/
bool PMCSQS_Scorer::init_for_current_spec(Config *_config, 
										  const BasicSpectrum& bs)
{
	config = _config;

	bs.calc_peak_isotope_levels(config->get_tolerance(),this->curr_spec_iso_levels);
	if (! bs.select_strong_peak_idxs(this->curr_spec_iso_levels,this->curr_spec_strong_inds))
		return false;

	bs.mark_all_possible_isotope_peaks(config->get_tolerance(),curr_spec_strict_iso_ind);

	curr_spec_total_intensity=0;
	curr_spec_strong_intensity=0;
	curr_spec_num_strong=0;
	int i;
	for (i=0; i<bs.num_peaks; i++)
	{
		if (curr_spec_iso_levels[i]>0)
			continue;

		curr_spec_total_intensity+=bs.peaks[i].intensity;
		if (curr_spec_strong_inds[i])
		{
			curr_spec_strong_intensity+=bs.peaks[i].intensity;
			curr_spec_num_strong++;
		}
	}


	if (curr_spec_rank_pmc_tables.size()<=max_model_charge+1)
	{
		curr_spec_rank_pmc_tables.clear();
		curr_spec_rank_pmc_tables.resize(max_model_charge+1);

		curr_spec_rank_background_stats.clear();
		curr_spec_rank_background_stats.resize(max_model_charge+1);

		curr_spec_rank_maximal_values.clear();
		curr_spec_rank_maximal_values.resize(max_model_charge+1);
	}

	return true;
}


/****************************************************************************
*****************************************************************************/
void PMCSQS_Scorer::fill_fval_vector_with_SQS(const BasicSpectrum& bs,
								   ME_Regression_Sample& sam) const
{
	const mass_t tolerance = config->get_tolerance();
	const mass_t frag_tolerance = (tolerance<0.15 ? tolerance : tolerance * 0.5);
	const QCPeak *peaks  = bs.peaks;
	const int  num_peaks = bs.num_peaks;
	const mass_t m_over_z = bs.ssf->m_over_z;
	const int num_strong=curr_spec_num_strong;
	const float total_intensity=curr_spec_total_intensity;
	const float one_over_np = 10.0 / (num_peaks+1);
	const float one_over_ns = 5.0/ (num_strong + 1);
	const float one_over_total_intensity = 1.0 / (total_intensity + 1.0);
	const mass_t max_peak_mass = (num_peaks <5 ? POS_INF : peaks[num_peaks-1].mass);

	sam.f_vals.clear();
	sam.f_vals.push_back(fval(SQS_CONST,1.0));
	sam.f_vals.push_back(fval(SQS_PEAK_DENSITY, (float)num_peaks/max_peak_mass));

	float grass_level_inten=NEG_INF;
	
	// calculate grass level peaks
	if (1)
	{
		int i;
		vector<float> peak_intens;
		peak_intens.resize(num_peaks);
		for (i=0; i<num_peaks; i++)
			peak_intens[i]=peaks[i].intensity;

		sort(peak_intens.begin(),peak_intens.end());

		int idx_G = num_peaks/3;
		grass_level_inten = peak_intens[idx_G];

		float cum_intensity2G=0;
		float inten_2G=2.0 * grass_level_inten;
		int idx_2G = idx_G;
		while (idx_2G<num_peaks && peak_intens[idx_2G]<inten_2G)
		{
			cum_intensity2G+=peak_intens[idx_2G];
			idx_2G++;
		}

		float cum_intensity5G=0;
		float inten_5G = 5.0 * grass_level_inten;
		int idx_5G = idx_2G;
		while (idx_5G<num_peaks && peak_intens[idx_5G]<inten_5G)
		{
			cum_intensity5G+=peak_intens[idx_5G];
			idx_5G++;
		}

		float inten_10G = 10.0 * grass_level_inten;
		int idx_10G = idx_5G;
		while (idx_10G<num_peaks && peak_intens[idx_10G]<inten_10G)
			idx_10G++;

		sam.f_vals.push_back(fval(SQS_PROP_UPTO2G,(float)idx_2G*one_over_np));
		sam.f_vals.push_back(fval(SQS_PROP_UPTO5G,(float)(idx_5G-idx_2G)*one_over_np));
		sam.f_vals.push_back(fval(SQS_PROP_UPTO10G,(float)(idx_10G-idx_5G)*one_over_np)); 
		sam.f_vals.push_back(fval(SQS_PROP_MORE10G,(float)(num_peaks-idx_10G)*one_over_np));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO2G,cum_intensity2G*one_over_total_intensity));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_UPTO5G,cum_intensity5G*one_over_total_intensity));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_MORE5G, (total_intensity-cum_intensity2G-cum_intensity5G)*one_over_total_intensity));
	}

	// isotope features
	if (1)
	{
		int i;
		int num_with_iso=0;
		int strong_with_iso=0;
		for (i=1; i<num_peaks; i++)
			if (curr_spec_iso_levels[i]>0)
			{
				num_with_iso++;
				if (curr_spec_iso_levels[i-1]==0 && curr_spec_strong_inds[i-1])
					strong_with_iso++;
			}

		sam.f_vals.push_back(fval(SQS_PROP_ISO_PEAKS,num_with_iso*one_over_np)); 
		sam.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_ISO_PEAKS,strong_with_iso*one_over_ns)); 

	//	cout << "PROP WITH ISO  : " << num_with_iso/(float)num_peaks << endl;
	//	cout << "STRONG WITH IOS: " << strong_with_iso/(float)num_strong << endl;
	}


	// neutral loss features
	if (1)
	{
		vector<float> tmp_vals;
		tmp_vals.clear();
		int frag_charge;
		for (frag_charge=1; frag_charge<=2; frag_charge++)
		{
			const mass_t offsets[3]={MASS_H2O/frag_charge, MASS_NH3/frag_charge, MASS_CO/frag_charge};
			const int num_offsets = 3;
			
			const int fc_off = (frag_charge-1)*num_offsets * 2;
			int i;
			for (i=0; i<num_offsets; i++)
			{
				const mass_t min_offset = offsets[i]-frag_tolerance;
				const mass_t max_offset = offsets[i]+frag_tolerance;

				int num_pairs=0;
				int num_strong_pairs=0;

				int trail_idx=0;
				int lead_idx=1;

				while (lead_idx<num_peaks)
				{
					while (peaks[lead_idx].mass-peaks[trail_idx].mass>max_offset)
						trail_idx++;

					if (peaks[lead_idx].mass-peaks[trail_idx].mass>min_offset)
					{
						num_pairs++;
						if (curr_spec_strong_inds[lead_idx])
							num_strong_pairs++;
					}
					lead_idx++;
				}

				const float prop_peaks = num_pairs*one_over_np;
				const float prop_strong = num_strong_pairs*one_over_ns;
				sam.f_vals.push_back(fval(SQS_PROP_ALL_WITH_H2O_LOSS+fc_off+i,prop_peaks));
				sam.f_vals.push_back(fval(SQS_PROP_STRONG_WITH_H2O_LOSS+fc_off+i,prop_strong));

				tmp_vals.push_back(prop_peaks);
				tmp_vals.push_back(prop_strong);

			//	cout << "OFF REG " << i << " " << num_pairs/(float)num_peaks << endl;
			//	cout << "OFF STR " << i << " " << num_strong_pairs/(float)num_strong << endl;
			}
		}

		int i;
		const int half_size = tmp_vals.size()/2; 
		for (i=0; i<half_size; i++)
			sam.f_vals.push_back(fval(SQS_DIFF_ALL_WITH_H2O_LOSS+i,tmp_vals[i]-tmp_vals[half_size+i]));

	//	curr_spec_strong_inds
		
		

	//	
		const mass_t half_max = max_peak_mass*0.5;
		const mass_t half_tol = tolerance * 0.66;
		
		int num_pairs=0;
		int num_strong_pairs=0;
		float inten_both=0;
		int didx=0;
		for (i=0; i<num_peaks; i++)
		{
			const mass_t peak_mass = peaks[i].mass;
			if (peak_mass>half_max)
				break;

			const mass_t doub_mass = peak_mass * 2 - MASS_PROTON;
			const mass_t min_mass = doub_mass - half_tol;
			const mass_t max_mass = doub_mass + half_tol;
			while (didx<num_peaks && peaks[didx].mass<min_mass)
				didx++;

			if (didx==num_peaks)
				break;
			
			if (peaks[didx].mass<max_mass)
			{
				num_pairs++;
				if (curr_spec_strong_inds[i] || curr_spec_strong_inds[didx])
					num_strong_pairs++;
				inten_both+= peaks[i].intensity + peaks[didx].intensity;
			}
		}

		sam.f_vals.push_back(fval(SQS_PROP_PEAKS_WITH_C1C2,num_pairs*one_over_np));
		sam.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_WITH_C1C2,num_strong_pairs*one_over_ns));
		sam.f_vals.push_back(fval(SQS_PROP_INTEN_WITH_C1C2,inten_both*one_over_total_intensity));
	}

	// tag features
	if (1)
	{
		vector<DPColumn> dp;

		vector<float> tmp_vals;
		int frag_charge;
		for (frag_charge=1; frag_charge<=2; frag_charge++)
		{
			fill_SQS_DP(bs,dp,frag_charge);
			float inten_in_tags[3]={0,0,0};
			int num_peaks_in_tags[3]={0,0,0};
			int num_strong_peaks_in_tags[3]={0,0,0};

			vector<int> tag_lengths;
			tag_lengths.resize(num_peaks,0);

			int longest=0;
			int i;
			for (i=1; i<num_peaks; i++)
			{
				int max_tl=-1;
				int aa;
				for (aa=Ala; aa<=Val; aa++)
				{
					const int prev = dp[i].pointers[aa];
					if (prev>=0 && tag_lengths[prev]>max_tl)
						max_tl = tag_lengths[prev];
				}

				tag_lengths[i]=max_tl+1;
				if (tag_lengths[i]>longest)
					longest = tag_lengths[i];

				int k;
				for (k=0; k<3; k++)
				{
					if (tag_lengths[i]==k)
						break;
					
					inten_in_tags[k]+=peaks[i].intensity;
					num_peaks_in_tags[k]++;
					if (curr_spec_strong_inds[i])
						num_strong_peaks_in_tags[k]++;
				}		
			}

			
			const int idx_off = (SQS_C2_IND_MAX_TAG_LENGTH_ABOVE_4 - SQS_IND_MAX_TAG_LENGTH_ABOVE_4)*(frag_charge-1);
			if (longest>=4)
			{
				sam.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_ABOVE_4+idx_off,1.0));
				sam.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_ABOVE_4+idx_off,(float)longest-4.0));
			}
			else
			{
				sam.f_vals.push_back(fval(SQS_IND_MAX_TAG_LENGTH_BELOW_4+idx_off,1.0));
				sam.f_vals.push_back(fval(SQS_MAX_TAG_LENGTH_BELOW_4+idx_off,(float)longest));
			}

	

			float inten_in_tags_both_sides=0;
			for (i=0; i<num_peaks; i++)
				if (dp[i].pointers[0])
					inten_in_tags_both_sides+=peaks[i].intensity;

			sam.f_vals.push_back(fval(SQS_PROP_INTEN_IN_TAGS+idx_off,inten_in_tags_both_sides*one_over_total_intensity));

			tmp_vals.push_back(longest);
			tmp_vals.push_back(inten_in_tags_both_sides*one_over_total_intensity);
		//	cout << "MAX TAG: " << longest << endl;
		//	cout << "Inten in TAG: " << inten_in_tags_both_sides/total_intensity << endl;

			const float strong_threshes[]={0.3,0.2,0.1};
			int k;
			for (k=0; k<3; k++)
			{
				const int pos_off = 4*k;
				const float prop_tags = (float)num_peaks_in_tags[k]*one_over_np;
				const float prop_strong = (float)num_strong_peaks_in_tags[k]*one_over_ns;
				const float prop_inten  = (float)inten_in_tags[k]*one_over_total_intensity;
				sam.f_vals.push_back(fval(SQS_PROP_TAGS1+idx_off+pos_off,prop_tags));
				sam.f_vals.push_back(fval(SQS_PROP_STRONG_PEAKS_IN_TAG1+idx_off+pos_off,prop_strong));
				sam.f_vals.push_back(fval(SQS_PROP_INTEN_TAG1+idx_off+pos_off,prop_inten));

				if (prop_strong<strong_threshes[k])
					sam.f_vals.push_back(fval(SQS_IND_PROP_STRONG_BELOW30_TAG1+idx_off+pos_off,1.0));

				// save vals for diff features
				tmp_vals.push_back(prop_tags);
				tmp_vals.push_back(prop_strong);
				tmp_vals.push_back(prop_inten);

			}
		}

		int pos_off = tmp_vals.size()/2;
		int i;
		for (i=0; i<pos_off; i++)
			sam.f_vals.push_back(fval(SQS_DIFF_MAX_TAG_LENGTH+i,tmp_vals[i]-tmp_vals[pos_off+i]));
		/*
		SQS_DIFF_MAX_TAG_LENGTH, SQS_DIFF_PROP_INTEN_IN_TAGS,
		SQS_DIFF_PROP_TAGS1, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG1, SQS_DIFF_PROP_INTEN_TAG1,
		SQS_DIFF_PROP_TAGS2, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG2, SQS_DIFF_PROP_INTEN_TAG2,
		SQS_DIFF_PROP_TAGS3, SQS_DIFF_PROP_STRONG_PEAKS_IN_TAG3, SQS_DIFF_PROP_INTEN_TAG3,
		*/
	}


	// density features
	if (1)
	{
		mass_t t1=m_over_z*0.6666;
		mass_t t2=m_over_z*1.3333;
		mass_t h1=m_over_z;

		float inten_t1=0,inten_t2=0, inten_h1=0;
		int	  num_peaks_t1=0, num_peaks_t2=0, num_peaks_h1=0;

		int i;
		for (i=0; i<num_peaks && peaks[i].mass<t1; i++)
			inten_t1+=peaks[i].intensity;
		
		num_peaks_t1=i;

		for ( ; i<num_peaks && peaks[i].mass<h1; i++)
			inten_h1+=peaks[i].intensity;

		num_peaks_h1=i;

		inten_t2=inten_h1;

		for ( ; i<num_peaks && peaks[i].mass<t2; i++)
			inten_t2+=peaks[i].intensity;

		num_peaks_t2=i-num_peaks_t1;

		int num_peaks_t3= num_peaks - num_peaks_t2 - num_peaks_t1;
		float inten_t3 = total_intensity - inten_t2 - inten_t1;

	
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T1,(float)(num_peaks_t1)*one_over_np));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T2,(float)(num_peaks_t2)*one_over_np));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_T3,(float)(num_peaks_t3)*one_over_np));

	//	cout << "T1 PEAKS: " << (float)(num_peaks_t1)/(float)num_peaks << endl;
	//	cout << "T2 PEAKS: " << (float)(num_peaks_t2)/(float)num_peaks << endl;
	//	cout << "T3 PEAKS: " << (float)(num_peaks_t3)/(float)num_peaks << endl;


		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T1,inten_t1*one_over_total_intensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T2,inten_t2*one_over_total_intensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_T3,inten_t3*one_over_total_intensity));

	//	cout << "T1 INTEN: " << inten_t1/total_intensity << endl;
	//	cout << "T2 INTEN: " << inten_t2/total_intensity << endl;
	//	cout << "T3 INTEN: " << inten_t3/total_intensity << endl;


		int num_peaks_h2 = num_peaks - num_peaks_h1;
		float inten_h2 = total_intensity - inten_h1;

		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_H1,(num_peaks_h1)*one_over_np));
		sam.f_vals.push_back(fval(SQS_PEAK_DENSE_H2,(num_peaks_h2)*one_over_np));

		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_H1,inten_h1*one_over_total_intensity));
		sam.f_vals.push_back(fval(SQS_INTEN_DENSE_H2,inten_h2*one_over_total_intensity));
	
		const float inten33=0.333* total_intensity;
		const float inten50=0.5  * total_intensity;
		const float inten75=0.75 * total_intensity;
		const float inten90=0.90 * total_intensity;

		float cum_inten=0;
		const mass_t one_over_mz = 1.0 / (m_over_z + 0.1);

		for (i=0; i<num_peaks && cum_inten<inten33; i++)
			cum_inten+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_33_INTEN,(peaks[i].mass*one_over_mz)-0.333));

		for ( ; i<num_peaks && cum_inten<inten50; i++)
			cum_inten+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_50_INTEN,(peaks[i].mass*one_over_mz)-0.5));

		for ( ; i<num_peaks && cum_inten<inten75; i++)
			cum_inten+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_75_INTEN,(peaks[i].mass*one_over_mz)-0.75));

		for ( ; i<num_peaks && cum_inten<inten90; i++)
			cum_inten+=peaks[i].intensity;

		if (i==num_peaks)
			i--;

		sam.f_vals.push_back(fval(SQS_PROP_MZ_RANGE_WITH_90_INTEN,(peaks[i].mass*one_over_mz)-0.9));		
	}

	// pmc features
	if (1)
	{
		static vector< vector<float> > sqs_features;
		
		get_sqs_features_from_pmc_tables(bs,sqs_features);

		vector<float> max_vals;
		max_vals.resize(4,0.000001);

		vector< vector<float> > tmp_vals;
		tmp_vals.resize(4);

		// first give absolute counts
		int charge;
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqs_features[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				sam.f_vals.push_back(fval(SQS_NUM_FRAG_PAIRS_1+idx_off+i,vals[i]));	
				tmp_vals[charge].push_back(vals[i]);
				if (vals[i]>max_vals[i])
					max_vals[i]=vals[i];
			}	
		}

		// add features for prop of max
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqs_features[charge];
			const int idx_off = 4 *( charge-1);

			int i;
			for (i=0; i<4; i++)
			{
				const float ratio = (max_vals[i]>0 ? vals[i]/max_vals[i] : 0);
				sam.f_vals.push_back(fval(SQS_PROP_OF_MAX_FRAG_PAIRS_1+idx_off+i,ratio));			
				tmp_vals[charge].push_back(ratio);
			}
		}

		// conver to proportions by dividing by the number of peaks/strong_peaks
		// and subtract the background levels (first)
		const float one_over_peaks = 1.0/((float)num_peaks+0.1);
		const float one_over_strong = 1.0/((float)num_strong+0.1);

		
	
		for (charge=1; charge<=3; charge++)
		{
			vector<float>& vals = sqs_features[charge];
			const int idx_off = 4 *( charge-1);

			// normalize to get proportions of total number of peaks/strong
			vals[0]-= this->curr_spec_rank_background_stats[charge].num_frag_pairs;
			vals[1]-= this->curr_spec_rank_background_stats[charge].num_strong_frag_pairs;
			vals[2]-= this->curr_spec_rank_background_stats[charge].num_c2_frag_pairs;
			vals[3]-= this->curr_spec_rank_background_stats[charge].num_strong_c2_frag_pairs;

			vals[0]*=one_over_peaks;
			vals[1]*=one_over_strong;
			vals[2]*=one_over_peaks;
			vals[3]*=one_over_strong;

			int i;
			for (i=0; i<4; i++)
			{
				sam.f_vals.push_back(fval(SQS_PROP_FRAG_PAIRS_1+idx_off+i,vals[i]));
				tmp_vals[charge].push_back(vals[i]);
			}
		}

		// diff is between values for c2 and c3
		int i;
		for (i=0; i<tmp_vals[2].size(); i++)
			sam.f_vals.push_back(fval(SQS_DIFF_NUM_FRAG_PAIRS_23+i,tmp_vals[2][i]-tmp_vals[3][i]));
	}

	
	sort(sam.f_vals.begin(),sam.f_vals.end());
}






/*******************************************************************************

  Calculates the average statistics observed for the background using different
  mass offsets

********************************************************************************/
void calc_background_stats(const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  Config *config,
						  const QCPeak *peaks, 
						  const int num_peaks,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  PMCRankStats& pmc_stats_total)
{
	const mass_t tolerance = config->get_tolerance();
	const mass_t offsets[]={-22.0,-10.0,-8.5,8.5,12.0,22.5};
	const int num_offsets = sizeof(offsets)/sizeof(mass_t);
	const float one_over_num_offsets = 1.0 / (float)num_offsets;

	int i;
	for (i=0; i<num_offsets; i++)
	{
		static PMCRankStats pmc_stats; 
		calc_pmc_rank_stats_for_mass(peaks,num_peaks, single_charge_pair_sum+offsets[i], 
			1.5*tolerance,iso_levels,strong_inds, strict_iso_inds,pmc_stats);

		if (i==0)
		{
			pmc_stats_total = pmc_stats;
		}
		else
		{
			pmc_stats_total.num_c2_frag_pairs        += pmc_stats.num_c2_frag_pairs;
			pmc_stats_total.num_frag_pairs           += pmc_stats.num_frag_pairs;
			pmc_stats_total.num_strong_c2_frag_pairs += pmc_stats.num_strong_c2_frag_pairs;
			pmc_stats_total.num_strong_frag_pairs	 += pmc_stats.num_strong_frag_pairs;

			pmc_stats_total.inten_frag_pairs		 += pmc_stats.inten_frag_pairs;
			pmc_stats_total.inten_strong_pairs		 += pmc_stats.inten_strong_pairs;
			pmc_stats_total.inten_c2_pairs			 += pmc_stats.inten_c2_pairs;
			pmc_stats_total.inten_c2_strong_pairs	 += pmc_stats.inten_c2_strong_pairs;

			
		}
	}

	pmc_stats_total.num_c2_frag_pairs		 *= one_over_num_offsets;
	pmc_stats_total.num_frag_pairs			 *= one_over_num_offsets;
	pmc_stats_total.num_strong_c2_frag_pairs *= one_over_num_offsets;
	pmc_stats_total.num_strong_frag_pairs	 *= one_over_num_offsets;

	pmc_stats_total.inten_frag_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_strong_pairs		 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_pairs			 *= one_over_num_offsets;
	pmc_stats_total.inten_c2_strong_pairs	 *= one_over_num_offsets;


}





/*********************************************************************************
Fills the PMC sample features
Assumes the frag_pair_sum_offset is set.
Sets the PMC features and also creates vetors of features that are to be added
to the SQS sampels.
**********************************************************************************/
void PMCSQS_Scorer::calculate_curr_spec_pmc_values(const BasicSpectrum& bs,
												   mass_t increment)
{
	const mass_t m_over_z = bs.ssf->m_over_z;


	if (frag_pair_sum_offset<-999)
	{
		cout << "Error: must first set the expected frag pair offset!" << endl;
		exit(1);
	}

	int charge;
	for (charge=1; charge<=max_model_charge; charge++)
	{
		const mass_t org_pm_with_19 = bs.ssf->m_over_z * charge - MASS_PROTON*(charge - 1);
		const mass_t frag_pair_sum = org_pm_with_19 + frag_pair_sum_offset;
		const int size_idx = this->get_rank_model_size_idx(charge,org_pm_with_19);

		float bias = 0;
		if (pmc_charge_mz_biases[charge].size()>size_idx)
			bias = this->pmc_charge_mz_biases[charge][size_idx];

		if (config->get_pm_tolerance()<0.075) // don't use bias if dealing with accurate pm
			bias=0;

		int bin_range =3*charge;
		if (charge>=2)
			bin_range = int(2*charge);
		if (charge>=3)
			bin_range = int(1.5*charge);

		fill_rank_PMC_stats(  charge,
							  frag_pair_sum + bias,
							  -bin_range-1, 
							  bin_range-1,
							  increment,
							  config,
							  bs,
							  curr_spec_strong_inds,
							  curr_spec_iso_levels,
							  curr_spec_strict_iso_ind,
							  curr_spec_rank_pmc_tables[charge]);


		calc_background_stats(frag_pair_sum,config,bs.peaks,bs.num_peaks, 
			curr_spec_strong_inds, curr_spec_iso_levels, 
			this->curr_spec_strict_iso_ind, curr_spec_rank_background_stats[charge]);


		// find maximal values
		int i;
		PMCRankStats& maximal = curr_spec_rank_maximal_values[charge];
		maximal = curr_spec_rank_pmc_tables[charge][0];

		for (i=1; i<curr_spec_rank_pmc_tables[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = curr_spec_rank_pmc_tables[charge][i];

			// mathced pairs look for maximum number

			if (curr_stats.num_c2_frag_pairs>maximal.num_c2_frag_pairs)
				maximal.num_c2_frag_pairs = curr_stats.num_c2_frag_pairs;

			if (curr_stats.num_frag_pairs>maximal.num_frag_pairs)
				maximal.num_frag_pairs = curr_stats.num_frag_pairs;

			if (curr_stats.num_strong_frag_pairs>maximal.num_strong_frag_pairs)
				maximal.num_strong_frag_pairs = curr_stats.num_strong_frag_pairs;

			if (curr_stats.num_strong_c2_frag_pairs > maximal.num_strong_c2_frag_pairs)
				maximal.num_strong_c2_frag_pairs = curr_stats.num_strong_c2_frag_pairs;

			// intensity look for maximal numbers

			if (curr_stats.inten_frag_pairs > maximal.inten_frag_pairs)
				maximal.inten_frag_pairs = curr_stats.inten_frag_pairs;
			
			if (curr_stats.inten_strong_pairs> maximal.inten_strong_pairs)
				maximal.inten_strong_pairs = curr_stats.inten_strong_pairs;

			if (curr_stats.inten_c2_pairs> maximal.inten_c2_pairs)
				maximal.inten_c2_pairs = curr_stats.inten_c2_pairs;

			if (curr_stats.inten_c2_strong_pairs > maximal.inten_c2_strong_pairs)
				maximal.inten_c2_strong_pairs = curr_stats.inten_c2_strong_pairs;


	
		}

		// find indicators for min tolerance and max pairs
		float tol_pairs =           POS_INF;
		float tol_strong_pairs =    POS_INF;
		float tol_c2_pairs =        POS_INF;
		float tol_c2_strong_pairs = POS_INF;

		int   idx_pairs=0;
		int	  idx_strong_pairs=0;
		int   idx_c2_pairs=0;
		int   idx_c2_strong_pairs=0;

		for (i=0; i<curr_spec_rank_pmc_tables[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = curr_spec_rank_pmc_tables[charge][i];

			if (curr_stats.num_frag_pairs == maximal.num_frag_pairs &&
				curr_stats.mean_offset_pairs < tol_pairs)
			{
				idx_pairs=i;
				tol_pairs=curr_stats.mean_offset_pairs;
			}

			if (curr_stats.num_strong_frag_pairs == maximal.num_strong_frag_pairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_strong_pairs)
			{
				idx_strong_pairs=i;
				tol_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

			if (curr_stats.num_c2_frag_pairs == maximal.num_c2_frag_pairs &&
				curr_stats.mean_offset_c2_pairs < tol_c2_pairs)
			{
				idx_c2_pairs=i;
				tol_c2_pairs=curr_stats.mean_offset_c2_pairs;
			}

			if (curr_stats.num_strong_c2_frag_pairs == maximal.num_strong_c2_frag_pairs &&
				curr_stats.mean_offset_c2_strong_pairs < tol_c2_strong_pairs)
			{
				idx_c2_strong_pairs=i;
				tol_c2_strong_pairs=curr_stats.mean_offset_c2_strong_pairs;
			}

		}

		curr_spec_rank_pmc_tables[charge][idx_pairs].ind_pairs_with_min_tol=true;
		curr_spec_rank_pmc_tables[charge][idx_strong_pairs].ind_strong_pairs_with_min_tol=true;
		curr_spec_rank_pmc_tables[charge][idx_c2_pairs].ind_c2_pairs_with_min_tol=true;
		curr_spec_rank_pmc_tables[charge][idx_c2_strong_pairs].ind_c2_strong_pairs_with_min_tol=true;


		static vector<float> log_distances;
		if (log_distances.size()<curr_spec_rank_pmc_tables[charge].size())
		{
			log_distances.resize(curr_spec_rank_pmc_tables[charge].size(),0);
			int i;
			for (i=1; i<log_distances.size(); i++)
				log_distances[i]=log(1.0+(float)i);
		}

		for (i=0; i<curr_spec_rank_pmc_tables[charge].size(); i++)
		{
			PMCRankStats& curr_stats = curr_spec_rank_pmc_tables[charge][i];
			curr_stats.log_dis_from_pairs_min_tol = log_distances[abs(i-idx_pairs)];
			curr_stats.log_dis_from_strong_pairs_min_tol = log_distances[abs(i-idx_strong_pairs)];
			curr_stats.log_dis_from_c2_pairs_min_tol = log_distances[abs(i-idx_c2_pairs)];
			curr_stats.log_dis_from_c2_strong_pairs_min_tol = log_distances[abs(i-idx_c2_strong_pairs)];
		}

	}
}
	


/*****************************************************************************************

******************************************************************************************/
void PMCSQS_Scorer::get_sqs_features_from_pmc_tables(const BasicSpectrum& bs,
						vector< vector<float> >& sqs_features) const
{
	const mass_t org_m_over_z = bs.ssf->m_over_z;
	float max_num_strong_pairs=0;
	int   best_table_idx=0;
	mass_t mz_diff = 99999;

	if (sqs_features.size() != max_model_charge +1)
		sqs_features.resize(max_model_charge+1);

	int charge;
	for (charge=1; charge<=max_model_charge; charge++)
	{
		// find entry which has the maximal number of strong pairs, while maining the minimial
		// m/z distance shift from the original
		int i;
		for (i=1; i<this->curr_spec_rank_pmc_tables[charge].size(); i++)
		{
			const float curr_num_strong = curr_spec_rank_pmc_tables[charge][i].num_strong_frag_pairs;
			
			if (curr_num_strong>=max_num_strong_pairs)
			{
				const mass_t distance = fabs(curr_spec_rank_pmc_tables[charge][i].m_over_z - org_m_over_z);
				if (curr_num_strong == max_num_strong_pairs && distance>= mz_diff)
					continue;
				
				max_num_strong_pairs = curr_num_strong;
				mz_diff = distance;
				best_table_idx = i;
			}
		}

		const PMCRankStats& best_stats = curr_spec_rank_pmc_tables[charge][best_table_idx];

		sqs_features[charge].clear();

		sqs_features[charge].push_back(best_stats.num_frag_pairs);

		sqs_features[charge].push_back(best_stats.num_strong_frag_pairs);

		sqs_features[charge].push_back(best_stats.num_c2_frag_pairs);

		sqs_features[charge].push_back(best_stats.num_strong_c2_frag_pairs);
	}
}



void PMCSQS_Scorer::compute_sqs_cum_stats_for_ided(Config *config, char *list)
{
	const float max_prob=0.5;
	const int max_bin_idx = int(max_prob*100)+1;
	vector< vector< vector<int> > > bin_counts; //charge / size_idx / bin
	vector< vector< int > > totals;

	bin_counts.resize(4);
	totals.resize(4);
	int i;
	for (i=0; i<4; i++)
	{
		int n = this->sqs_mass_thresholds.size() + 1;
		totals[i].resize(n,0);
		bin_counts[i].resize(n);
		int j;
		for (j=0; j<n; j++)
			bin_counts[i][j].resize(max_bin_idx+1,0);
	}

	// prefrom sqs on all

		

	FileManager fm;
	FileSet fs;

	fm.init_from_list_file(config,list);
	fs.select_all_files(fm);

	const int max_to_read_per_file = 100000;

	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	int low_count = 0;
	for (i=0; i<all_ssf.size(); i++)
	{
		SingleSpectrumFile* ssf = all_ssf[i];
		BasicSpectrum bs;
	
		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

		const int org_charge=ssf->peptide.calc_charge(ssf->m_over_z);
		const int size_idx = get_sqs_size_idx(ssf->m_over_z);
		int charge=0;
		float prob=get_sqs_for_spectrum(config,bs,&charge);

		int bin_idx=int(100*prob);
		if (bin_idx>max_bin_idx)
			bin_idx = max_bin_idx;

		totals[org_charge][size_idx]++;
		int j;
		for (j=0; j<=bin_idx; j++)
			bin_counts[org_charge][size_idx][j]++;

		if (prob<0.1)
		{
			cout << ++low_count << "\t" << ssf->single_name << "\t" << org_charge << " : " << 
				charge << " --> " << prob << endl;
		}
	}

	int c;
	for (c=0; c<4; c++)
	{
		int specs=0;
		int i;
		for (i=0; i<totals[c].size(); i++)
			specs+=totals[c][i];
		if (specs<100)
			continue;

		cout << endl << "CHARGE " << c << endl;
		for (i=0; i<totals[c].size(); i++)
		{
			if (totals[c][i]<50)
				continue;
			cout << "SIZE " << i << " ( charge " << c << " )" << endl;
			int j;
			cout << setprecision(4);
			for (j=0; j<=max_bin_idx; j++)
			{
				cout << j*0.01 << "\t" << bin_counts[c][i][j] << "\t" << 
					bin_counts[c][i][j]/(float)totals[c][i] << endl;
			}
			cout << endl;
		}
	}
}






void PMCSQS_Scorer::create_filtered_peak_list_for_sqs(
									  QCPeak *org_peaks, int num_org_peaks,
									  QCPeak *new_peaks, int& num_new_peaks) const
{
	static vector<MassInten> peak_list;
	static int peak_list_size =0;
	int i;

	if (num_org_peaks>peak_list_size)
	{
		peak_list_size = (int)(num_org_peaks * 1.5);
		if (peak_list_size<2000)
			peak_list_size = 2000;

		peak_list.resize(peak_list_size);
	}

	// copy org_peaks to the temporary peak_list
	int f_idx=0;
	for (i=0; i<num_org_peaks; i++)
	{
		peak_list[i].mass=org_peaks[i].mass;
		peak_list[i].intensity=org_peaks[i].intensity;
	}

	const mass_t tolerance = config->get_tolerance();
	const mass_t join_tolerance = (tolerance < 0.05 ? tolerance : 0.5 * tolerance);
	int p_idx=0;
	i=1;
	while (i<num_org_peaks)
	{
		if (peak_list[i].mass - peak_list[p_idx].mass<=join_tolerance)
		{
			intensity_t inten_sum = peak_list[i].intensity + peak_list[p_idx].intensity;
			mass_t new_mass = (peak_list[i].intensity * peak_list[i].mass + 
							   peak_list[p_idx].intensity * peak_list[p_idx].mass ) / inten_sum;

			peak_list[p_idx].mass = new_mass;
			peak_list[p_idx].intensity = inten_sum;	
		}
		else
		{
			peak_list[++p_idx]=peak_list[i];
		}
		i++;
	}
	int num_peaks = p_idx+1;


	// filter low intensity noise
	const mass_t half_window_size = 0.5 * config->get_local_window_size();
	const int num_peaks_in_window = config->get_max_number_peaks_per_local_window();
	const int max_peak_idx = num_peaks -1;
	int min_idx=1;
	int max_idx=1;
	p_idx =1;

	
	new_peaks[0].mass      = peak_list[0].mass;
	new_peaks[0].intensity = peak_list[0].intensity;
	f_idx=1;

	// check the rest of the peaks
	for (i=1; i<max_peak_idx; i++)
	{
		const mass_t& peak_mass=peak_list[i].mass;
		mass_t min_mass = peak_list[min_idx].mass;
		mass_t max_mass = peak_list[max_idx].mass;

		// advance min/max pointers
		while (peak_mass-min_mass > half_window_size)
			min_mass=peak_list[++min_idx].mass;

		while (max_idx < max_peak_idx && max_mass - peak_mass <= half_window_size)
			max_mass=peak_list[++max_idx].mass;

		if (max_mass - peak_mass > half_window_size)
			max_idx--;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < num_peaks_in_window)
		{
			new_peaks[f_idx].mass = peak_list[i].mass;
			new_peaks[f_idx].intensity = peak_list[i].intensity;
			f_idx++;
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peak_list[j].intensity > peak_list[i].intensity)
				higher_count++;

		if (higher_count < num_peaks_in_window)
		{
			new_peaks[f_idx].mass = peak_list[i].mass;
			new_peaks[f_idx].intensity = peak_list[i].intensity;
			f_idx++;
		}
	}
	new_peaks[f_idx].mass = peak_list[i].mass;
	new_peaks[f_idx].intensity = peak_list[i].intensity;
	f_idx++;


	num_new_peaks = f_idx;

	// normalize intensities

	if (1)
	{
		intensity_t total_inten=0;

		for (i=1; i<num_new_peaks; i++)
			total_inten+=new_peaks[i].intensity;

		const mass_t one_over_total_inten = (1000.0 / total_inten);

		for (i=1; i<num_new_peaks; i++)
			new_peaks[i].intensity *= one_over_total_inten; 
	}
}



/*****************************************************************
Creates for each input file an mgf file that holds the spectra
that passed quality filtering does not correct PM and charge in the 
mgf files.
******************************************************************/
void PMCSQS_Scorer::output_filtered_spectra_to_mgfs(
									 Config *config,
									 const vector<string>& files,
									 char *out_dir,
									 float filter_prob, 
									 int& total_num_written, 
									 int& total_num_read)
{
	total_num_read = 0;
	total_num_read = 0;
	int f;
	for (f=0; f<files.size(); f++)
	{
		const char *spectra_file = files[f].c_str();
		FileManager fm;
		FileSet fs;
		BasicSpecReader bsr;

		string fname, mgf_name, map_name;
		get_file_name_without_extension(files[f],fname);

		mgf_name = string(out_dir) + "/" + fname + "_fil.mgf";
		map_name = string(out_dir) + "/" + fname + "_map.txt";

		///////////////////////////////////////////////
		// Quick read, get all pointers to begining of spectra
		if (get_file_extension_type(files[f]) != MZXML)
		{
			fm.init_from_file(config,spectra_file);
		}
		else
			fm.init_and_read_single_mzXML(config,spectra_file,f);

		fs.select_all_files(fm);

		const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
		int sc;
		int  num_spec_written=0;
		bool first=true;
		ofstream out_stream, map_stream;
		for (sc=0; sc<all_ssf.size(); sc++)
		{
			static vector<QCPeak> peaks;
			SingleSpectrumFile *ssf = all_ssf[sc];
			
			if (peaks.size()<ssf->num_peaks)
			{
				int new_size = ssf->num_peaks*2;
				if (new_size<2500)
					new_size=2500;
				peaks.resize(new_size);
			}

			// read without processing peaks
			const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0],false,true);
			ssf->file_idx = f;

			BasicSpectrum bs;
			bs.peaks = &peaks[0];
			bs.num_peaks = num_peaks;
			bs.ssf = all_ssf[sc];

			int max_charge;
			float prob = get_sqs_for_spectrum(config,bs,&max_charge);
			if (prob<filter_prob)
			{
				continue;
			}

			if (first)
			{
				out_stream.open(mgf_name.c_str(),ios::out);
				map_stream.open(map_name.c_str(),ios::out);

				cout << "Filtering spectra to minumum quality score: " << filter_prob << endl;
				cout << "Writing spectra info to:" << endl;
				cout << mgf_name << endl << map_name << endl;

				if (! out_stream.is_open() || ! out_stream.good())
				{
					cout << "Error: couldn\'t open for out mgf stream for writing: " <<
						endl << mgf_name << endl;
					exit(1);
				}
				first = false;
			}

			char single_name[64];
			sprintf(single_name,"%d:%d",f,bs.ssf->get_scan());
			bs.ssf->single_name = single_name;
			bs.output_to_mgf(out_stream,config);
			if (prob>1.0)
				prob=1.0;
			map_stream << num_spec_written++ << "\t" << all_ssf[sc]->get_scan() << "\t" << fixed << prob << endl;

			
		}
		out_stream.close();
		map_stream.close();

		total_num_read+= all_ssf.size();
		total_num_written += num_spec_written;
	}
}
	
