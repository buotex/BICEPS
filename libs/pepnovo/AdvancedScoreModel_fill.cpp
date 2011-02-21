#include "AdvancedScoreModel.h"


struct FragStats {
	FragStats() : frag_idx(NEG_INF), is_viz(false), has_intensity(false), peak_idx(NEG_INF),
					   mass(NEG_INF), log_intensity(NEG_INF), log_local_rank(NEG_INF), log_global_rank(NEG_INF) {};

	void fill_from_breakage(const Breakage *breakage, Spectrum *spec, int f)
	{
		frag_idx = f;
		if (breakage->is_frag_type_visible(f))
		{
			is_viz=true;
			const int pos = breakage->get_position_of_frag_idx(f);
			if (pos>=0)
			{
				has_intensity=true;
				peak_idx = breakage->fragments[pos].peak_idx;

				const Peak& peak = spec->get_peak(peak_idx);
				mass		    = peak.mass;
				iso_level	    = peak.iso_level;
				log_intensity   = peak.log_intensity;
				log_local_rank  = peak.log_local_rank;
				log_global_rank = log(1.0+(float)peak.rank);
			}
		}
	}

	int frag_idx;
	bool is_viz;
	bool has_intensity;
	int    peak_idx;
	mass_t mass;
	float  iso_level;
	float  log_intensity;
	float  log_local_rank;
	float  log_global_rank;
};

void StrongFragModel::fill_constant_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage, 
							   vector<fval>& f_vals) const
{
	const mass_t iso_tolerance = spec->get_config()->get_tolerance()*0.5;
	FragStats frag_stats;
	FragStats parent1_stats, parent2_stats;
	FragStats mirror1_stats, mirror2_stats;

	frag_stats.fill_from_breakage(breakage,spec,model_frag_idx);
	if (parent1_idx>=0)
		parent1_stats.fill_from_breakage(breakage,spec,parent1_idx);
	if (parent2_idx>=0)
		parent2_stats.fill_from_breakage(breakage,spec,parent2_idx);
	if (mirror1_idx>=0)
		mirror1_stats.fill_from_breakage(breakage,spec,mirror1_idx);
	if (mirror2_idx>=0)
		mirror2_stats.fill_from_breakage(breakage,spec,mirror2_idx);

	f_vals.clear();
	if (frag_stats.has_intensity) // fill features for visible fragment
	{
		f_vals.push_back(fval(SI_CONST,1.0));

		const float log_inten = frag_stats.log_intensity;

		// Mirror1 features
		if (mirror1_idx>=0)
		{
			if (mirror1_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_MIRROR1_VIZ,1.0));
				if (mirror1_stats.has_intensity)
				{
					f_vals.push_back(fval(SI_IND_HAS_MIRROR1_INTEN,1.0));
					if (mirror1_stats.iso_level!=0)
						f_vals.push_back(fval(SI_MIRROR1_ISO_LEVEL,mirror1_stats.iso_level));
				
					vector<float> iso_intens;
					spec->get_iso_intens(mirror1_stats.peak_idx, iso_intens, iso_tolerance, mirror1_charge);	
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR1_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR1_MINUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR1_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR1_MINUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR1_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR1_PLUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR1_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR1_PLUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[3]));	
						}
					}
						

					const mass_t sum_masses = frag_stats.mass * model_frag_charge + 
										mirror1_stats.mass * mirror1_charge;

					const mass_t offset = fabs(pm_with_19 - sum_masses + (model_frag_charge + mirror1_charge -1)*MASS_PROTON);
					const float offset_level = offset * one_over_tolerance;

					if (offset_level<0.25)
					{
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF25,1.0));
					} 
					else if (offset_level<0.75)
					{
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF75,1.0));
					}
					else
						f_vals.push_back(fval(SI_MIRROR1_MASS_DIFF_LARGE,1.0));
				}
			}
			else
				f_vals.push_back(fval(SI_IND_MIRROR1_NOT_VIZ,1.0));
		}

		// Mirror2 features
		if (mirror2_idx>=0)
		{
			if (mirror2_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_MIRROR2_VIZ,1.0));
				if (mirror2_stats.has_intensity)
				{
					f_vals.push_back(fval(SI_IND_HAS_MIRROR2_INTEN,1.0));
					if (mirror2_stats.iso_level != 0)
						f_vals.push_back(fval(SI_MIRROR2_ISO_LEVEL,mirror2_stats.iso_level));

					vector<float> iso_intens;
					spec->get_iso_intens(mirror2_stats.peak_idx, iso_intens, iso_tolerance, mirror2_charge);
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR2_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR2_MINUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR2_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR2_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SI_IND_MIRROR2_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SI_MIRROR2_PLUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SI_IND_MIRROR2_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SI_MIRROR2_PLUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[3]));	
						}
					}



					const mass_t sum_masses = frag_stats.mass * model_frag_charge + 
										mirror2_stats.mass * mirror2_charge;

					const mass_t offset = fabs(pm_with_19 - sum_masses + (model_frag_charge + mirror2_charge -1)*MASS_PROTON);
					const float offset_level = offset * one_over_tolerance;

					if (offset_level<0.25)
					{
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF25,1.0));
					} 
					else if (offset_level<0.75)
					{
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF75,1.0));
					}
					else
						f_vals.push_back(fval(SI_MIRROR2_MASS_DIFF_LARGE,1.0));
				}
			}
			else
				f_vals.push_back(fval(SI_IND_MIRROR2_NOT_VIZ,1.0));
		}

		// Parent 1 features
		if (parent1_idx>=0)
		{
			if (parent1_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_PARENT1_VIZ,1.0));
				if (parent1_stats.has_intensity)
				{
					const float inten_diff = parent1_stats.log_intensity - log_inten;
					const mass_t dis_min = parent1_stats.mass - spec->get_min_peak_mass();
					const mass_t dis_max = spec->get_max_peak_mass() - parent1_stats.mass;
					const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);
					
					if (parent1_stats.iso_level != 0)
						f_vals.push_back(fval(SI_PARENT1_ISO_LEVEL,parent1_stats.iso_level));
					if (dis<100)
					{
						f_vals.push_back(fval(SI_IND_PARENT1_LESS_THAN_100_MIN_MAX,1.0));
					} 
					else if (dis<200)
						f_vals.push_back(fval(SI_IND_PARENT1_LESS_THAN_200_MIN_MAX,1.0));

					if (inten_diff>0)
					{
						f_vals.push_back(fval(SI_IND_PARENT1_INTEN_MORE,1.0));
						f_vals.push_back(fval(SI_PARENT1_INTEN_DIFF_MORE,inten_diff));
					}
					else
					{
						f_vals.push_back(fval(SI_IND_PARENT1_INTEN_LESS,1.0));
						f_vals.push_back(fval(SI_PARENT1_INTEN_DIFF_LESS,inten_diff));
					}
				}
				else
					f_vals.push_back(fval(SI_IND_PARENT1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SI_IND_PARENT1_NOT_VIZ,1.0));
		}

		// Parent 2 features
		if (parent2_idx>=0)
		{
			if (parent2_stats.is_viz)
			{
				f_vals.push_back(fval(SI_IND_PARENT2_VIZ,1.0));
				if (parent2_stats.has_intensity)
				{
					const float inten_diff = parent2_stats.log_intensity - log_inten;
					const mass_t dis_min = parent2_stats.mass - spec->get_min_peak_mass();
					const mass_t dis_max = spec->get_max_peak_mass() - parent2_stats.mass;
					const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);
					
					if (parent2_stats.iso_level != 0)
						f_vals.push_back(fval(SI_PARENT2_ISO_LEVEL,parent2_stats.iso_level));
					if (dis<100)
					{
						f_vals.push_back(fval(SI_IND_PARENT2_LESS_THAN_100_MIN_MAX,1.0));
					} 
					else if (dis<200)
						f_vals.push_back(fval(SI_IND_PARENT2_LESS_THAN_200_MIN_MAX,1.0));

					if (inten_diff>0)
					{
						f_vals.push_back(fval(SI_IND_PARENT2_INTEN_MORE,1.0));
						f_vals.push_back(fval(SI_PARENT2_INTEN_DIFF_MORE,inten_diff));
					}
					else
					{
						f_vals.push_back(fval(SI_IND_PARENT2_INTEN_LESS,1.0));
						f_vals.push_back(fval(SI_PARENT2_INTEN_DIFF_LESS,inten_diff));
					}
				}
				else
					f_vals.push_back(fval(SI_IND_PARENT2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SI_IND_PARENT2_NOT_VIZ,1.0));
		}

		// self intensity
		f_vals.push_back(fval(SI_LOG_LOCAL_RANK,frag_stats.log_local_rank));
		f_vals.push_back(fval(SI_LOG_GLOBAL_RANK,frag_stats.log_global_rank));
		if (frag_stats.iso_level != 0)
			f_vals.push_back(fval(SI_ISO_LEVEL,frag_stats.iso_level));
		
		if (log_inten<1.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS1,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS1,log_inten));
		} 
		else if (log_inten<2.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS2,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS2,log_inten-1.0));
		}
		else if (log_inten<3.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS3,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS3,log_inten-2.0));
		}
		else if (log_inten<4.0)
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_LESS4,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_LESS4,log_inten-3.0));
		}
		else
		{
			f_vals.push_back(fval(SI_IND_LOG_INTEN_MORE,1.0));
			f_vals.push_back(fval(SI_LOG_INTEN_MORE,log_inten-4.0));
		}

		// self distance
		const mass_t dis_min = frag_stats.mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - frag_stats.mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX0,dis));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_50,log_inten));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX50,dis-50.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_150,log_inten));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX150,dis-150.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_250,log_inten));
		}
		else
		{
			f_vals.push_back(fval(SI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(SI_DIS_FROM_MINMAX250,dis-250.0));
			f_vals.push_back(fval(SI_LOG_INTEN_DIS_MORE,log_inten));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(SI_REL_POS0+rel_pos,1.0));

		const float one_over_model_frag_charge = 1.0 / model_frag_charge;
		const mass_t forward_offsets[]={MASS_NH3,MASS_H2O,MASS_CO,MASS_H2ONH3,MASS_H2OH2O};
		const int num_forward_offsets = sizeof(forward_offsets)/sizeof(mass_t);
		const mass_t peak_mass = frag_stats.mass;

		int t;
		for (t=0; t<num_forward_offsets; t++)
		{
			const int forward_idx = spec->get_max_inten_peak(peak_mass + forward_offsets[t]*one_over_model_frag_charge,frag_tolerance);
			if (forward_idx>0)
			{
				f_vals.push_back(fval(SI_IND_HAS_PLUS_NH3+2*t,1.0));
				f_vals.push_back(fval(SI_IND_HAS_PLUS_NH3+2*t+1,spec->get_peak(forward_idx).log_intensity-log_inten));
			}
		}

		
		const int plus_idx = spec->get_max_inten_peak((peak_mass * model_frag_charge + MASS_PROTON)/(model_frag_charge+1.0),frag_tolerance);
		if (plus_idx>=0)
		{
			f_vals.push_back(fval(SI_IND_HAS_CHARGE_PLUS1,1.0));
			f_vals.push_back(fval(SI_CHARGE_PLUS1_INTEN_DIFF,spec->get_peak(plus_idx).log_intensity-log_inten));
		}

		if (model_frag_charge>1)
		{
			const int minus_idx = spec->get_max_inten_peak((peak_mass * model_frag_charge - MASS_PROTON)/(model_frag_charge-1.0),frag_tolerance);
			if (minus_idx>=0)
			{
				f_vals.push_back(fval(SI_IND_HAS_CHARGE_MINUS1,1.0));
				f_vals.push_back(fval(SI_CHARGE_MINUS1_INTEN_DIFF,spec->get_peak(minus_idx).log_intensity-log_inten));
			}
		}
	}
	else // Fill features for non-visible
	{
		f_vals.push_back(fval(SNI_CONST,1.0));

		// Mirror1 features
		if (mirror1_idx>=0)
		{
			if (mirror1_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_MIRROR1_VIZ,1.0));
				if (mirror1_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_HAS_MIRROR1_INTEN,1.0));
					if (mirror1_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_MIRROR1_ISO_LEVEL,mirror1_stats.iso_level));

					vector<float> iso_intens;
					spec->get_iso_intens(mirror1_stats.peak_idx, iso_intens, iso_tolerance, mirror1_charge);
					
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_MINUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR1_MINUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_PLUS_1_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR1_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR1_PLUS_2_INTEN_DIFF,mirror1_stats.log_intensity-iso_intens[3]));	
						}
					}
				}
				else
					f_vals.push_back(fval(SNI_IND_MIRROR1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_MIRROR1_NOT_VIZ,1.0));
		}

		// Mirror2 features
		if (mirror2_idx>=0)
		{
			if (mirror2_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_MIRROR2_VIZ,1.0));
				if (mirror2_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_HAS_MIRROR2_INTEN,1.0));
					if (mirror2_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_MIRROR2_ISO_LEVEL,mirror2_stats.iso_level));

					vector<float> iso_intens;
					spec->get_iso_intens(mirror2_stats.peak_idx, iso_intens, iso_tolerance, mirror2_charge);
					
					if (iso_intens[1]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_MINUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR1_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[1]));
						if (iso_intens[0]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_MINUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR2_MINUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[0]));	
						}
					}
					if (iso_intens[2]>=0)
					{
						f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_PLUS_1,1.0));
						f_vals.push_back(fval(SNI_MIRROR2_PLUS_1_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[2]));
						if (iso_intens[3]>=0)
						{
							f_vals.push_back(fval(SNI_IND_MIRROR2_HAS_PLUS_2,1.0));
							f_vals.push_back(fval(SNI_MIRROR2_PLUS_2_INTEN_DIFF,mirror2_stats.log_intensity-iso_intens[3]));	
						}
					}
				}
				else
					f_vals.push_back(fval(SNI_IND_MIRROR2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_MIRROR2_NOT_VIZ,1.0));
		}

		// Parent 1 features
		if (parent1_idx>=0)
		{
			if (parent1_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_PARENT1_VIZ,1.0));
				if (parent1_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_PARENT1_INTEN,1.0));
					if (parent1_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_PARENT1_ISO_LEVEL,parent1_stats.iso_level));
					f_vals.push_back(fval(SNI_PARENT1_LOG_INTEN,parent1_stats.log_intensity));
					f_vals.push_back(fval(SNI_PARENT1_LOG_GLOBAL_RANK,parent1_stats.log_global_rank));
				}
				else
					f_vals.push_back(fval(SNI_IND_PARENT1_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_PARENT1_NOT_VIZ,1.0));
		}

		// Parent 2 features
		if (parent2_idx>=0)
		{
			if (parent2_stats.is_viz)
			{
				f_vals.push_back(fval(SNI_IND_PARENT2_VIZ,1.0));
				if (parent2_stats.has_intensity)
				{
					f_vals.push_back(fval(SNI_IND_PARENT2_INTEN,1.0));
					if (SNI_PARENT2_ISO_LEVEL,parent2_stats.iso_level != 0)
						f_vals.push_back(fval(SNI_PARENT2_ISO_LEVEL,parent2_stats.iso_level));
					f_vals.push_back(fval(SNI_PARENT2_LOG_INTEN,parent2_stats.log_intensity));
					f_vals.push_back(fval(SNI_PARENT2_LOG_GLOBAL_RANK,parent2_stats.log_global_rank));
				}
				else
					f_vals.push_back(fval(SNI_IND_PARENT2_NO_INTEN,1.0));
			}
			else
				f_vals.push_back(fval(SNI_IND_PARENT2_NOT_VIZ,1.0));
		}

		// self distance
		const mass_t expected_mass = config->get_fragment(model_frag_idx).calc_expected_mass(breakage->mass,pm_with_19);
		const mass_t dis_min = expected_mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - expected_mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX0,dis));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX50,dis-50.0));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX150,dis-150.0));
		}
		else
		{
			f_vals.push_back(fval(SNI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(SNI_DIS_FROM_MINMAX250,dis-250.0));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(SNI_REL_POS0+rel_pos,1.0));
	}
}



/**********************************************************************************
***********************************************************************************/
void StrongFragModel::fill_aa_variable_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const
{
	const vector<int>& org_aas = config->get_org_aa();
	const int n_aa = (info->n_aa>=0 ? org_aas[info->n_aa] : Gap);
	const int c_aa = (info->c_aa>=0 ? org_aas[info->c_aa] : Gap);
	const int pos = breakage->get_position_of_frag_idx(model_frag_idx);
	
	const bool do_n_features = (n_aa != Gap);
	const bool do_c_features = (c_aa != Gap);

	const FragmentType& fragment = config->get_fragment(model_frag_idx);
	const mass_t exp_peak_mass = (info->breakage ? 
		fragment.calc_expected_mass(info->breakage->mass,pm_with_19) : NEG_INF);
	const mass_t exp_n_peak_mass = (do_n_features && info->n_break ? 
		fragment.calc_expected_mass(info->n_break->mass,pm_with_19) : NEG_INF);
	const mass_t exp_c_peak_mass = (do_c_features && info->c_break ? 
		fragment.calc_expected_mass(info->c_break->mass,pm_with_19) : NEG_INF);	
	
	static vector<int> threshes;
	if (threshes.size()==0)
	{
		threshes.push_back(18);
		threshes.push_back(12);
		threshes.push_back(8);
		threshes.push_back(4);
		threshes.push_back(2);
		threshes.push_back(NEG_INF);
	}

	// fill intensity
	if (pos>=0)
	{
		if (! do_n_features)
			f_vals.push_back(fval(SI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(SI_IND_C_IS_GAP,1.0));

		// aa category feature
		if (do_n_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->n_side_cat>threshes[k])
					break;
			if (info->connects_to_N_term)
			{
				f_vals.push_back(fval(SI_N_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SI_N_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;
		}

		if (do_c_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;

			if (info->connects_to_C_term)
			{
				f_vals.push_back(fval(SI_C_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SI_C_EDGE_CAT20+k,1.0));

			if (do_n_features)
			{
				int k;
				for (k=0; k<threshes.size(); k++)
					if (info->span_cat>threshes[k])
						break;
				f_vals.push_back(fval(SI_SPAN_CAT20+k,1.0));

				if (info->n_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->n_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SI_ND_SPAN_CAT20+k,1.0));
				}

				if (info->c_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->c_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SI_CD_SPAN_CAT20+k,1.0));
				}
			}
		}

		const Peak& peak = spec->get_peak(breakage->fragments[pos].peak_idx);

		if (do_n_features)
		{
			if (info->connects_to_N_term)
				f_vals.push_back(fval(SI_IND_CONNECTS_TO_N_TERM,1.0));

			if (info->preferred_digest_aa_N_term)
				f_vals.push_back(fval(SI_IND_PREFERRED_DIGEST_AA_N_TERM,1.0));
		}

		if (do_c_features)
		{
			if (info->connects_to_C_term)
				f_vals.push_back(fval(SI_IND_CONNECTS_TO_C_TERM,1.0));

			if (info->preferred_digest_aa_C_term)
				f_vals.push_back(fval(SI_IND_PREFERRED_DIGEST_AA_C_TERM,1.0));
		}

		if (!info->connects_to_N_term && ! info->connects_to_C_term)
			f_vals.push_back(fval(SI_IND_NOT_CONNECTED_TO_TERMS,1.0));
		
		if (info->missed_cleavage)
			f_vals.push_back(fval(SI_IND_MISSED_CLEAVAGE,1.0));
	

		if (do_n_features)
		{
			if (! info->n_break || ! info->n_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SI_IND_N_FRAG_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SI_IND_N_N_TERM,1.0));
				f_vals.push_back(fval(SI_IND_N_N_TERM+n_aa,1.0));
				
				const int n_pos=info->n_break->get_position_of_frag_idx(model_frag_idx);
				if (n_pos>=0)
				{
					f_vals.push_back(fval(SI_IND_N_INTEN,1.0));
					f_vals.push_back(fval(SI_IND_N_N_TERM_INTEN,peak.log_intensity));
					f_vals.push_back(fval(SI_IND_N_N_TERM_INTEN+n_aa,peak.log_intensity));
				}
				else
					f_vals.push_back(fval(SI_IND_N_NO_INTEN,1.0));

				const mass_t dis_min =exp_n_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_n_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->n_edge_is_single)
				{
					f_vals.push_back(fval(SI_IND_N_SE,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->get_peak(info->n_break->fragments[n_pos].peak_idx);
						const mass_t n_diff = fabs(fabs((peak.mass - n_peak.mass)* model_frag_charge) - info->exp_n_edge_mass);

						f_vals.push_back(fval(SI_SE_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (n_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_01,1.0));
						}
						else if (n_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_SE_IND_N_FRAG_DIFF_LARGE,1.0));

						const float diff_inten = peak.log_intensity - n_peak.log_intensity;
						f_vals.push_back(fval(SI_SE_IND_N_N_TERM_DIFF_INTEN,diff_inten));
						f_vals.push_back(fval(SI_SE_IND_N_N_TERM_DIFF_INTEN+n_aa,diff_inten));
					}
					else
					{
					
						f_vals.push_back(fval(SI_SE_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SI_IND_N_ME,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->get_peak(info->n_break->fragments[n_pos].peak_idx);
						const mass_t n_diff = fabs(fabs((peak.mass - n_peak.mass)* model_frag_charge) - info->exp_n_edge_mass);

						f_vals.push_back(fval(SI_ME_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (n_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_01,1.0));
						}
						else if (n_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_ME_IND_N_FRAG_DIFF_LARGE,1.0));
					}
					else
					{
						f_vals.push_back(fval(SI_ME_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}

		if (do_c_features)
		{
			if (! info->c_break || ! info->c_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SI_IND_C_FRAG_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SI_IND_C_N_TERM,1.0));
				f_vals.push_back(fval(SI_IND_C_N_TERM+c_aa,1.0));
				
				const int c_pos=info->c_break->get_position_of_frag_idx(model_frag_idx);
				if (c_pos>=0)
				{
					f_vals.push_back(fval(SI_IND_C_INTEN,1.0));
					f_vals.push_back(fval(SI_IND_C_N_TERM_INTEN,peak.log_intensity));
					f_vals.push_back(fval(SI_IND_C_N_TERM_INTEN+c_aa,peak.log_intensity));
				}
				else
					f_vals.push_back(fval(SI_IND_C_NO_INTEN,1.0));

				const mass_t dis_min =exp_c_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_c_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->c_edge_is_single)
				{
					f_vals.push_back(fval(SI_IND_C_SE,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->get_peak(info->c_break->fragments[c_pos].peak_idx);
						const mass_t c_diff = fabs(fabs((peak.mass - c_peak.mass)* model_frag_charge) - info->exp_c_edge_mass);

						f_vals.push_back(fval(SI_SE_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (c_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_01,1.0));
						}
						else if (c_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_SE_IND_C_FRAG_DIFF_LARGE,1.0));

						const float diff_inten = peak.log_intensity - c_peak.log_intensity;
						f_vals.push_back(fval(SI_SE_IND_C_N_TERM_DIFF_INTEN,diff_inten));
						f_vals.push_back(fval(SI_SE_IND_C_N_TERM_DIFF_INTEN+c_aa,diff_inten));
					}
					else
					{
						f_vals.push_back(fval(SI_SE_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SI_IND_C_ME,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->get_peak(info->c_break->fragments[c_pos].peak_idx);
						const mass_t c_diff = fabs(fabs((peak.mass - c_peak.mass)* model_frag_charge) - info->exp_c_edge_mass);

						f_vals.push_back(fval(SI_ME_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
				
						if (c_diff<exact_peak_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_01,1.0));
						}
						else if (c_diff<frag_tolerance)
						{
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_05,1.0));
						}
						else
							f_vals.push_back(fval(SI_ME_IND_C_FRAG_DIFF_LARGE,1.0));
					}
					else
					{
						f_vals.push_back(fval(SI_ME_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}
	}
	else  // fill features for no intensity
	{
		if (! do_n_features)
			f_vals.push_back(fval(SNI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(SNI_IND_C_IS_GAP,1.0));

		// aa category feature
		if (do_n_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->n_side_cat>threshes[k])
					break;
			if (info->connects_to_N_term)
			{
				f_vals.push_back(fval(SNI_N_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SNI_N_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;
		}

		if (do_c_features)
		{
			int k;
			for (k=0; k<threshes.size(); k++)
				if (info->c_side_cat>threshes[k])
					break;

			if (info->connects_to_C_term)
			{
				f_vals.push_back(fval(SNI_C_TERM_CAT20+k,1.0));
			}
			else
				f_vals.push_back(fval(SNI_C_EDGE_CAT20+k,1.0));

			for (k=0; k<threshes.size(); k++)
				if (info->span_cat>threshes[k])
					break;
		
			if (do_n_features)
			{
				f_vals.push_back(fval(SNI_SPAN_CAT20+k,1.0));

				if (info->n_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->n_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SNI_ND_SPAN_CAT20+k,1.0));
				}

				if (info->c_double_span_cat>NEG_INF)
				{
					int k;
					for (k=0; k<threshes.size(); k++)
						if (info->c_double_span_cat>threshes[k])
							break;
					f_vals.push_back(fval(SNI_CD_SPAN_CAT20+k,1.0));
				}
			}
		}

		if (do_n_features)
		{
			if (info->connects_to_N_term)
				f_vals.push_back(fval(SNI_IND_CONNECTS_TO_N_TERM,1.0));

			if (info->preferred_digest_aa_N_term)
				f_vals.push_back(fval(SNI_IND_PREFERRED_DIGEST_AA_N_TERM,1.0));
		}

		if (do_c_features)
		{
			if (info->connects_to_C_term)
				f_vals.push_back(fval(SNI_IND_CONNECTS_TO_C_TERM,1.0));

			if (info->preferred_digest_aa_C_term)
				f_vals.push_back(fval(SNI_IND_PREFERRED_DIGEST_AA_C_TERM,1.0));
		}


		if (! info->connects_to_N_term && ! info->connects_to_C_term)
			f_vals.push_back(fval(SNI_IND_NOT_CONNECTED_TO_TERMS,1.0));

		if (info->missed_cleavage)
			f_vals.push_back(fval(SNI_IND_MISSED_CLEAVAGE,1.0));
	

		if (do_n_features)
		{
			if (! info->n_break || ! info->n_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SNI_IND_N_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SNI_IND_N_N_TERM,1.0));
				f_vals.push_back(fval(SNI_IND_N_N_TERM+n_aa,1.0));
			
				const int n_pos=info->n_break->get_position_of_frag_idx(model_frag_idx);
				if (n_pos>=0)
				{
					f_vals.push_back(fval(SNI_IND_N_INTEN,1.0));
				}
				else
					f_vals.push_back(fval(SNI_IND_N_NO_INTEN,1.0));

				const mass_t dis_min = exp_n_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_n_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->n_edge_is_single)
				{
					f_vals.push_back(fval(SNI_IND_N_SE,1.0));
					if (n_pos>=0)
					{
						const Peak& n_peak = spec->get_peak(info->n_break->fragments[n_pos].peak_idx);

						f_vals.push_back(fval(SNI_SE_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));

						f_vals.push_back(fval(SNI_SE_IND_N_N_TERM_DIFF_INTEN,n_peak.log_intensity));
						f_vals.push_back(fval(SNI_SE_IND_N_N_TERM_DIFF_INTEN+n_aa,n_peak.log_intensity));
					}
					else
					{
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SNI_IND_N_ME,1.0));
					if (n_pos>=0)
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
					}
					else
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_NO_N_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}


		if (do_c_features)
		{
			if (! info->c_break || ! info->c_break->is_frag_type_visible(model_frag_idx))
			{
				f_vals.push_back(fval(SNI_IND_C_NOT_VIZ,1.0));
			}
			else
			{
				f_vals.push_back(fval(SNI_IND_C_N_TERM,1.0));
				f_vals.push_back(fval(SNI_IND_C_N_TERM+c_aa,1.0));
			
				const int c_pos=info->c_break->get_position_of_frag_idx(model_frag_idx);
				if (c_pos>=0)
				{
					f_vals.push_back(fval(SNI_IND_C_INTEN,1.0));
				}
				else
					f_vals.push_back(fval(SNI_IND_C_NO_INTEN,1.0));

				const mass_t dis_min = exp_c_peak_mass - spec->get_min_peak_mass();
				const mass_t dis_max = spec->get_max_peak_mass() - exp_c_peak_mass;
				const mass_t min_dis = (dis_min<dis_max ? dis_min : dis_max);

				if (info->c_edge_is_single)
				{
					f_vals.push_back(fval(SNI_IND_C_SE,1.0));
					if (c_pos>=0)
					{
						const Peak& c_peak = spec->get_peak(info->c_break->fragments[c_pos].peak_idx);
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));

						f_vals.push_back(fval(SNI_SE_IND_C_N_TERM_DIFF_INTEN,c_peak.log_intensity));
						f_vals.push_back(fval(SNI_SE_IND_C_N_TERM_DIFF_INTEN+c_aa,c_peak.log_intensity));
					}
					else
					{
					
						f_vals.push_back(fval(SNI_SE_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
				else // multiple edge
				{
					f_vals.push_back(fval(SNI_IND_C_ME,1.0));
					if (c_pos>=0)
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,1.0));
					}
					else
					{
						f_vals.push_back(fval(SNI_ME_IND_HAS_NO_C_FRAG_INTEN,1.0));
						if (min_dis<50)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN,1.0));
						} 
						else if (min_dis<150)
						{
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN,1.0));
						} 
						else if (min_dis<250)
							f_vals.push_back(fval(SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,1.0));
					}
				}
			}
		}
	}
}



void RegularFragModel::fill_constant_vals(
							Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage, 
							vector<fval>& f_vals) const
{
	const int num_parents = parent_idxs.size();
	FragStats frag_stats;
	vector<FragStats> parent_frag_stasts;

	parent_frag_stasts.resize(num_parents);
	frag_stats.fill_from_breakage(breakage,spec,model_frag_idx);
	int i;
	for (i=0; i<num_parents; i++)
		parent_frag_stasts[i].fill_from_breakage(breakage,spec,parent_idxs[i]);
	
	f_vals.clear();

	if (frag_stats.has_intensity) // fill features for visible fragment
	{
		f_vals.push_back(fval(RI_CONST,1.0));

		const float log_inten = frag_stats.log_intensity;
		f_vals.push_back(fval(RI_LOG_LOCAL_RANK,frag_stats.log_local_rank));
		f_vals.push_back(fval(RI_LOG_GLOBAL_RANK,frag_stats.log_global_rank));
		if (frag_stats.iso_level != 0)
			f_vals.push_back(fval(RI_ISO_LEVEL,frag_stats.iso_level));

		if (log_inten<1.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS1,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS1,log_inten));
		} 
		else if (log_inten<2.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS2,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS2,log_inten-1.0));
		}
		else if (log_inten<3.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS3,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS3,log_inten-2.0));
		}
		else if (log_inten<4.0)
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_LESS4,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_LESS4,log_inten-3.0));
		}
		else
		{
			f_vals.push_back(fval(RI_IND_LOG_INTEN_MORE,1.0));
			f_vals.push_back(fval(RI_LOG_INTEN_MORE,log_inten-4.0));
		}


			// self distance
		const mass_t dis_min = frag_stats.mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - frag_stats.mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		if (dis<50)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX0,dis));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS50,log_inten));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX50,dis-50.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS150,log_inten));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX150,dis-150.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DIS250,log_inten));
		}
		else
		{
			f_vals.push_back(fval(RI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(RI_DIS_FROM_MINMAX250,dis-250.0));
			f_vals.push_back(fval(RI_LOG_INTEN_DISMORE,log_inten));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(RI_REL_POS0+rel_pos,1.0));

		bool got_prefix = false;
		bool got_suffix = false;
		int combo_idx=0;
		int pow=1;
		int num_parents_with_intensity=0;
		int i;
		for (i=0; i<num_parents; i++)
		{
			if (parent_frag_stasts[i].has_intensity)
			{
				num_parents_with_intensity++;

				if (config->get_fragment(parent_idxs[i]).orientation == PREFIX)
				{
					got_prefix = true;
				}
				else
					got_suffix = true;

				if (i<4)
					combo_idx += pow;
			}
			pow*=2;
		}

		if (num_parents_with_intensity>6)
			num_parents_with_intensity=7;

		f_vals.push_back(fval(RI_IND_NUM_PARENTS_WITH_INTEN_IS_0+num_parents_with_intensity,1.0));

		f_vals.push_back(fval(RI_IND_PARENT_COMBO_0+combo_idx,1.0));

		if (got_prefix && got_suffix)
		{
			f_vals.push_back(fval(RI_IND_GOT_BOTH_ORIS,1.0));	   
		}
		else if (got_prefix)
		{
			f_vals.push_back(fval(RI_IND_GOT_PREFIX,1.0)); 
		}
		else if (got_suffix)
			f_vals.push_back(fval(RI_IND_GOT_SUFFIX,1.0));

		for (i=0; i<num_parents; i++)
		{
			const int offset_idx = 7*i;
			if (! parent_frag_stasts[i].is_viz)
			{
				f_vals.push_back(fval(RI_IND_PARENT1_NOT_VIZ+offset_idx,1.0));
			}
			else
			{
				if (! parent_frag_stasts[i].has_intensity)
				{
					f_vals.push_back(fval(RI_IND_PARENT1_NO_INTEN+offset_idx,1.0));
				}
				else
				{
					if (parent_frag_stasts[i].iso_level != 0)
						f_vals.push_back(fval(RI_PARENT1_ISO_LEVEL+offset_idx,parent_frag_stasts[i].iso_level));

					const float parent_log = parent_frag_stasts[i].log_intensity;
					if (parent_log>log_inten)
					{
						f_vals.push_back(fval(RI_IND_PARENT1_INTEN_MORE+offset_idx,1.0));
						f_vals.push_back(fval(RI_PARENT1_INTEN_DIFF_MORE+offset_idx,parent_log-log_inten));
					}
					else
					{
						f_vals.push_back(fval(RI_IND_PARENT1_INTEN_LESS+offset_idx,1.0));
						f_vals.push_back(fval(RI_PARENT1_INTEN_DIFF_LESS+offset_idx,parent_log-log_inten));
					}
				}
			}
		}

		for (i=0; i<f_vals.size(); i++)
		if (f_vals[i].f_idx>RI_NUM_FEATURES)
		{
			int qq=1;
		}
	
	}
	else // Fill features for non-visible
	{
		const mass_t expected_mass = config->get_fragment(model_frag_idx).calc_expected_mass(breakage->mass,pm_with_19);
		const mass_t dis_min = expected_mass - spec->get_min_peak_mass();
		const mass_t dis_max = spec->get_max_peak_mass() - expected_mass;
		const mass_t dis = (dis_min<dis_max ? dis_min : dis_max);

		f_vals.push_back(fval(RNI_CONST,1.0));

		if (dis<50)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_50,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX0,dis));
		}
		else if (dis<150)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_150,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX50,dis-50.0));
		}
		else if (dis<250)
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_LESS_250,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX150,dis-150.0));
		}
		else
		{
			f_vals.push_back(fval(RNI_IND_DIS_FROM_MINMAX_MORE,1.0));
			f_vals.push_back(fval(RNI_DIS_FROM_MINMAX250,dis-250.0));
		}

		const int rel_pos = int(10*breakage->mass/pm_with_19);
		f_vals.push_back(fval(RNI_REL_POS0+rel_pos,1.0));


		bool got_prefix=false;
		bool got_suffix=false;
		int combo_idx=0;
		int pow=1;
		int num_parents_with_intensity=0;
		int i;
		for (i=0; i<num_parents; i++)
		{
			if (parent_frag_stasts[i].has_intensity)
			{
				num_parents_with_intensity++;

				if (config->get_fragment(parent_idxs[i]).orientation == PREFIX)
				{
					got_prefix = true;
				}
				else
					got_suffix = true;

				if (i<4)
					combo_idx += pow;
			}
			pow*=2;
		}

		if (num_parents_with_intensity>6)
			num_parents_with_intensity=7;

		f_vals.push_back(fval(RNI_IND_NUM_PARENTS_WITH_INTEN_IS_0+num_parents_with_intensity,1.0));

		f_vals.push_back(fval(RNI_IND_PARENT_COMBO_0+combo_idx,1.0));

		if (got_prefix && got_suffix)
		{
			f_vals.push_back(fval(RI_IND_GOT_BOTH_ORIS,1.0));	   
		}
		else if (got_prefix)
		{
			f_vals.push_back(fval(RI_IND_GOT_PREFIX,1.0)); 
		}
		else if (got_suffix)
			f_vals.push_back(fval(RI_IND_GOT_SUFFIX,1.0));

		for (i=0; i<num_parents; i++)
		{
			const int offset_idx = 5*i;
			if (! parent_frag_stasts[i].is_viz)
			{
				f_vals.push_back(fval(RNI_IND_PARENT1_NOT_VIZ+offset_idx,1.0));
			}
			else
			{
				if (! parent_frag_stasts[i].has_intensity)
				{
					f_vals.push_back(fval(RNI_IND_PARENT1_NO_INTEN+offset_idx,1.0));
				}
				else
				{
					f_vals.push_back(fval(RNI_PARENT1_LOG_INTEN+offset_idx,parent_frag_stasts[i].log_intensity));
					f_vals.push_back(fval(RNI_PARENT1_LOG_GLOBAL_RANK+offset_idx,parent_frag_stasts[i].log_global_rank));
					if (parent_frag_stasts[i].iso_level != 0)
						f_vals.push_back(fval(RNI_PARENT1_ISO_LEVEL+offset_idx,parent_frag_stasts[i].iso_level));
				}
			}
		}	
	}
}


/**********************************************************************************
***********************************************************************************/
void RegularFragModel::fill_aa_variable_vals(
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const
{
	const vector<int>& org_aas = config->get_org_aa();
	const int n_aa = (info->n_aa>=0 ? org_aas[info->n_aa] : Gap);
	const int c_aa = (info->c_aa>=0 ? org_aas[info->c_aa] : Gap);
	const int pos = breakage->get_position_of_frag_idx(model_frag_idx);
	
	const bool do_n_features = (n_aa != Gap);
	const bool do_c_features = (c_aa != Gap);


	// fill intensity
	if (pos>=0)
	{
		if (! do_n_features)
			f_vals.push_back(fval(RI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(RI_IND_C_IS_GAP,1.0));

		const float log_inten = spec->get_peak(breakage->fragments[pos].peak_idx).log_intensity;

		// add aa presence indicators
		if (do_n_features)
		{
			int *var_ptr = info->n_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int n_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RI_IND_N_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RI_IND_N_HAS_N_TERM+n_aa,1.0));
				f_vals.push_back(fval(RI_N_N_TERM_SELF_INTEN,log_inten));
				f_vals.push_back(fval(RI_N_N_TERM_SELF_INTEN+n_aa,log_inten));
			}
		}

		if (do_c_features)
		{
			int *var_ptr = info->c_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int c_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RI_IND_C_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RI_IND_C_HAS_N_TERM+c_aa,1.0));
				f_vals.push_back(fval(RI_C_N_TERM_SELF_INTEN,log_inten));
				f_vals.push_back(fval(RI_C_N_TERM_SELF_INTEN+c_aa,log_inten));
			}
		}
	}
	else // Fill no frag features
	{
		if (! do_n_features)
			f_vals.push_back(fval(RNI_IND_N_IS_GAP,1.0));
		
		if (! do_c_features)
			f_vals.push_back(fval(RNI_IND_C_IS_GAP,1.0));


		// add aa presence indicators
		if (do_n_features)
		{
			int *var_ptr = info->n_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int n_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RNI_IND_N_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RNI_IND_N_HAS_N_TERM+n_aa,1.0));
			}
		}

		if (do_c_features)
		{
			int *var_ptr = info->c_var_ptr;
			const int num_aa = *var_ptr++;
			int k;
			for (k=0; k<num_aa; k++)
			{
				const int c_aa = org_aas[var_ptr[k]];
				f_vals.push_back(fval(RNI_IND_C_HAS_N_TERM,1.0));
				f_vals.push_back(fval(RNI_IND_C_HAS_N_TERM+c_aa,1.0));
			}
		}
	}
}





void StrongFragModel::fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const
{
	vector<fval> const_vals;

	fill_constant_vals(spec,pm_with_19,breakage,const_vals);
	samples.resize(infos.size());
	int i;
	for (i=0; i<infos.size(); i++)
	{
		vector< fval > var_vals;
		fill_aa_variable_vals(spec,pm_with_19,breakage,&infos[i],var_vals);

		samples[i].f_vals = const_vals;
		int j;
		for (j=0; j<var_vals.size(); j++)
			samples[i].f_vals.push_back(var_vals[j]);
	}
}


void StrongFragModel::fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const
{

	fill_constant_vals(spec,pm_with_19,breakage,f_vals);
	fill_aa_variable_vals(spec,pm_with_19,breakage,&info,f_vals);
	sort(f_vals.begin(),f_vals.end());
}

void RegularFragModel::fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const
{
	vector<fval> const_vals;
	fill_constant_vals(spec,pm_with_19,breakage,const_vals);
	samples.resize(infos.size());
	int i;
	for (i=0; i<infos.size(); i++)
	{
		vector< fval > var_vals;
		fill_aa_variable_vals(spec,pm_with_19,breakage,&infos[i],var_vals);

		samples[i].f_vals = const_vals;
		int j;
		for (j=0; j<var_vals.size(); j++)
			samples[i].f_vals.push_back(var_vals[j]);
	}
}


void RegularFragModel::fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const
{

	fill_constant_vals(spec,pm_with_19,breakage,f_vals);
	fill_aa_variable_vals(spec,pm_with_19,breakage,&info,f_vals);
	sort(f_vals.begin(),f_vals.end());
}
