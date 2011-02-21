#include "Spectrum.h"
#include "Isotopes.h"
#include "QuickClustering.h"
#include "FileManagement.h"
#include "auxfun.h"



void Spectrum::init_spectrum(bool perform_filtering)
{
	join_adjacent_peaks();
	if (perform_filtering)
		filter_peaks();

	init_index_array();
	calc_ranks();
	calc_log_local_ranks();
	normalize_intensities();
	calc_isotope_levels();
	select_strong_peaks();
	set_log_random_probs();
}





/*********************************************************************
Initializes the index array.
For each rounded off Dalton m, it gives the index of the closest peak i
with mass m_i>= m.
**********************************************************************/
void Spectrum::init_index_array()
{
	int i,c,size=(int)(max_mass+27.0);
	const int max_peak_idx = peaks.size()-1;
	
	index_array.clear();
	index_array.resize(size+1,max_peak_idx);
	
	i=0;
	int m=(int)peaks[0].mass;
	while (i<m)
		index_array[i++]=0;

	c=0;
	while (c< max_peak_idx)
	{
		int curr_m=(int)peaks[c].mass;
		int next_m = curr_m;
		int next_c = c;

		while (next_m == curr_m && next_c<max_peak_idx)
			next_m=(int)peaks[++next_c].mass;

		if (next_m>=size)
			next_m = size-1;

		while (i<next_m)
			index_array[i++]=c;
		
		c=next_c;
	}
}




/**********************************************************************
It then joins all pairs that are less than max_proximity away from each
other.
**********************************************************************/
void Spectrum::join_adjacent_peaks()
{
	vector<Peak> new_peaks;
	mass_t max_proximity = config->get_tolerance() * 0.5;
	int i;

	if (peaks.size() == 0)
		return;

	new_peaks.clear();
	new_peaks.push_back(peaks[0]);
	int prev_idx=0;
	for (i=1; i<peaks.size(); i++)
	{
		if 	(peaks[i].mass - new_peaks[prev_idx].mass < max_proximity)
		{
			// join peaks with proportion to their intensities

			int tt = new_peaks.size();

			intensity_t inten_sum=(new_peaks[prev_idx].intensity + peaks[i].intensity);
			mass_t ratio = new_peaks[prev_idx].intensity/inten_sum;
			mass_t new_mass = ratio *new_peaks[prev_idx].mass + (1-ratio)*peaks[i].mass;
			
			new_peaks[prev_idx].intensity = inten_sum;
			new_peaks[prev_idx].mass = new_mass;
		}
		else
		{
			new_peaks.push_back(peaks[i]);
			prev_idx++;
		}
	}

	// realocate peak memory
	if (new_peaks.size() != peaks.size())
		peaks = new_peaks;

}



/*******************************************************************
Filters the number of peaks in the spectra. Keeps only the highest
intensity peaks in the window. Uses values for number of peaks per
windows of 100 Da. from the Config.
Also keeps pairs of peaks that add up to mass+20 (b,y pairs)
Also keeps neutral losses for kept peaks (-H2O and -NH3)
********************************************************************/
void Spectrum::filter_peaks()
{
	const mass_t max_allowed_peak_mass = max_mass;
	const mass_t window_size = 0.5 * config->get_local_window_size();
	const int num_peaks_in_window = config->get_max_number_peaks_per_local_window();
	const mass_t tolerance = config->get_tolerance();
	int i,j,min_idx,max_idx;
	vector<bool> keep_peaks;	
	int new_num_peaks=0;
	vector<Peak> new_peaks;

	if (peaks.size()<5)
		return;

	int max_peak_idx = peaks.size() -1;

	keep_peaks.resize(peaks.size(),false);

	// keep first peak and last peak
	keep_peaks[0]=true;
	if (peaks[max_peak_idx].mass<max_allowed_peak_mass)
		keep_peaks[max_peak_idx]=true;
	
	min_idx=1;
	max_idx=1;

	// check the rest of the peaks
	for (i=1; i<max_peak_idx; i++)
	{
		mass_t peak_mass=peaks[i].mass;
		mass_t min_mass=peaks[min_idx].mass;
		mass_t max_mass=peaks[max_idx].mass;

		if (peaks[i].mass > max_allowed_peak_mass)
			break;

		// advance min/max pointers
		while (peak_mass-min_mass > window_size)
			min_mass=peaks[++min_idx].mass;

		while (max_idx < max_peak_idx && max_mass - peak_mass <= window_size)
			max_mass=peaks[++max_idx].mass;

		if (max_mass - peak_mass > window_size)
			max_idx--;

		// this peak might already be marked for keeping (isotpoic peak)
		if (keep_peaks[i])
			continue;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < num_peaks_in_window)
		{
			keep_peaks[i]=true;
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peaks[j].intensity > peaks[i].intensity)
				higher_count++;

		if (higher_count < num_peaks_in_window)
		{
			keep_peaks[i]=true;
		}
	}



	// look for b/y pairs
	mass_t pm_with_20 = (corrected_pm_with_19>0 ? corrected_pm_with_19 : 
						 org_pm_with_19)+MASS_PROTON ;

	mass_t pm_with_20_upper  = pm_with_20 + tolerance;
	mass_t pm_with_20_lower = pm_with_20 - tolerance;

	int f_idx =0;
	int b_idx = peaks.size()-1;
	while (f_idx<peaks.size() && b_idx>=0)
	{
		if (! keep_peaks[f_idx])
		{
			f_idx++;
			continue;
		}

		while (b_idx>=0 && peaks[f_idx].mass + peaks[b_idx].mass > pm_with_20_upper )
			b_idx--;

		if (b_idx<0)
			break;

		mass_t mass_sum = peaks[f_idx].mass + peaks[b_idx].mass;
		if (mass_sum > pm_with_20_lower && mass_sum < pm_with_20_upper)
		{
			keep_peaks[f_idx]=true;
			keep_peaks[b_idx]=true;
		}
		f_idx++;
	}


	// look for -H2O -NH3 peaks
	// 17.0265, 18.0105
	const mass_t frag_tolerance = tolerance * 0.6;
	int p_idx = 0;
	while (p_idx<peaks.size())
	{
		if (! keep_peaks[p_idx])
		{
			p_idx++;
			continue;
		}

		const mass_t upper_H2O = peaks[p_idx].mass + MASS_H2O + frag_tolerance;
		const mass_t lower_H2O = peaks[p_idx].mass + MASS_H2O - frag_tolerance;
		const mass_t upper_NH3 = peaks[p_idx].mass + MASS_NH3 + frag_tolerance;
		const mass_t lower_NH3 = peaks[p_idx].mass + MASS_NH3 - frag_tolerance;

		int f_idx;
		for (f_idx=p_idx+1; f_idx<peaks.size(); f_idx++)
		{
			if (peaks[f_idx].mass>upper_H2O)
				break;

			if (peaks[f_idx].mass>lower_H2O)
			{
				keep_peaks[f_idx]=true;
				break;
			}

			if (peaks[f_idx].mass>=lower_NH3 && peaks[f_idx].mass<=upper_NH3)
				keep_peaks[f_idx]=true;
		}
		p_idx++;
	}




	new_num_peaks=0;
	for (i=0; i<peaks.size(); i++)
		if (keep_peaks[i])
			new_num_peaks++;

	new_peaks.resize(new_num_peaks);
	
	j=0;
	for (i=0; i<peaks.size(); i++)
		if (keep_peaks[i])
			new_peaks[j++]=peaks[i];
		
	peaks=new_peaks;

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; 
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
}



struct peak_pair {
	bool operator< (const peak_pair& other) const
	{
		return inten>other.inten;
	}
	int idx;
	intensity_t inten;
};

void Spectrum::calc_ranks()
{
	int i;
	vector<peak_pair> pairs;
	pairs.resize(peaks.size());

	for (i=0; i<peaks.size(); i++)
	{
		pairs[i].idx=i;
		pairs[i].inten=peaks[i].intensity;
	}

	sort(pairs.begin(), pairs.end());

	for (i=0; i<pairs.size(); i++)
		peaks[pairs[i].idx].rank=i+1;
}


/*********************************************************************
// gives each peak it log local rank
// good be done more effciently...
**********************************************************************/
void Spectrum::calc_log_local_ranks()
{
	const mass_t half_window_size = config->get_local_window_size() * 0.5;
	const int num_peaks = peaks.size();
	int i;
	
	for (i=0; i<num_peaks; i++)
	{
		
		const mass_t peak_mass = peaks[i].mass;

		const PeakRange pr= get_peaks_in_range(peak_mass - half_window_size,
											   peak_mass + half_window_size);


		if (pr.num_peaks==0)
		{

		}

		int j, above=0;

	/*	if (pr.high_idx>=num_peaks)
		{
			cout << i << " MASS " << peak_mass << endl;
			this->print_spectrum();
			cout << endl;
			for (j=0; j<index_array.size(); j++)
				cout << j << "\t" << index_array[j] << endl;
			exit(0);
		}*/
		for (j=pr.low_idx; j<=pr.high_idx && j<num_peaks; j++)
			if (peaks[j].intensity>peaks[i].intensity)
				above++;

		peaks[i].log_local_rank = log(1.0 + (float)above);
	}	
}



/**************************************************************
Calculates the normalized intensity for each peaks. 
The normalization is done so that the new sum of all intensities 
equals 1000
***************************************************************/
void Spectrum::normalize_intensities()
{
	int i;

	if (! config->get_need_to_normalize())
		return;

	// remove the intensity of pm+20 / 2 and pm+2/2 from the total intensity
	// if this is an issue, normalization might have to be done twice,
	// before and after parent mass correction
	intensity_t total_intensity = this->total_original_intensity;
	
	const float norm_val = 1000.0 / total_intensity;
	for (i=0; i<peaks.size(); i++)
	{
		peaks[i].intensity *= norm_val;
		peaks[i].log_intensity = log(1.0 + peaks[i].intensity);
	}
}


struct rank_pair {
	bool operator< (const rank_pair& other) const
	{
		return inten > other.inten;
	}

	int idx;
	intensity_t inten;
};


/**************************************************************
Sets the iso_level for each peak. Iso level 0 means that there 
is no evidence that this peaks is an isotopic peak. The higher
the level, the more this looks like an isotopic peak
***************************************************************/
void Spectrum::calc_isotope_levels()
{
	int i;
	const mass_t iso_tolerance = (config->get_tolerance()<0.2) ?
								  config->get_tolerance() : 0.2;

	const int last_peak_idx = peaks.size()-1;

	for (i=0; i<last_peak_idx; i++)
	{	
		if (peaks[i].iso_level>0)
			continue;

		// look for +1 peak
		int idx1 = get_max_inten_peak(peaks[i].mass + 1.0033, iso_tolerance);
		if (idx1<0)  
			continue;

		float one_over_intensity = 1.0 / peaks[i].intensity;
		float ratio1 = peaks[idx1].intensity * one_over_intensity;

		// ignore strong +1
		if ( ratio1 > 3)
			continue;

		// examine ratios
		vector<float> expected_ratios, observed_ratios, relative_ratios;
		vector<int> iso_idxs;
		observed_ratios.resize(6);
		observed_ratios[0]=1.0;
		observed_ratios[1]= ratio1;

		iso_idxs.resize(6);
		iso_idxs[0]=i;
		iso_idxs[1]=idx1;

		// find additional peaks
		int j;
		for (j=2; j<=5; j++)
		{
			int idx = get_max_inten_peak(peaks[i].mass + j, iso_tolerance);
			if (idx<0)
				break;
			observed_ratios[j] = peaks[idx].mass * one_over_intensity;
			iso_idxs[j]=idx;
		}
		int last_iso = j-1;

		// get expected iso ratios
		calc_expected_iso_ratios(peaks[i].mass,expected_ratios,j);

		// calc ratios between observed and expected		
		relative_ratios.resize(j);
		relative_ratios[0]=1;
		for (j=1; j<=last_iso; j++)
			relative_ratios[j]=observed_ratios[j] / expected_ratios[j];

		peaks[i].iso_level=0;

		float level_mul=1.0;
		for (j=1; j<= last_iso; j++)
		{
			float iso_level;

			if (relative_ratios[j]>= 0.75 && relative_ratios[j]<=1.333)
			{
				iso_level=2.0;
			}
			else if (relative_ratios[j] >= 0.5 && relative_ratios[j] <=2)
			{
				iso_level=1.3333;
			}
			else if (relative_ratios[j] >= 0.3333 && relative_ratios[j] <=3)
			{
				iso_level=0.6666;
			}
			else if (relative_ratios[j] >= 0.25 && relative_ratios[j] <= 4)
			{
				iso_level=0.3333;
			}
			else
				break;

		//	if (relative_ratios[j] / relative_ratios[j-1] > 3)
		//		break;
			
			peaks[iso_idxs[j]].iso_level = peaks[iso_idxs[j-1]].iso_level + level_mul * iso_level;
			level_mul *= 0.5;
		}
	}

	vector<rank_pair> mono_pairs, iso_pairs;
	for (i=0; i<peaks.size(); i++)
	{
		rank_pair p;
		p.idx=i;
		p.inten = peaks[i].intensity;
		if (peaks[i].iso_level == 0)
		{
			mono_pairs.push_back(p);
		}
		else
			iso_pairs.push_back(p);

	}

	sort(mono_pairs.begin(),mono_pairs.end());
	sort(iso_pairs.begin(),iso_pairs.end());

	ranked_mono_peaks.resize(peaks.size());
	
	for (i=0; i<mono_pairs.size(); i++)
		ranked_mono_peaks[i]=mono_pairs[i].idx;

	for (i=0; i<iso_pairs.size(); i++)
		ranked_mono_peaks[mono_pairs.size()+i]=iso_pairs[i].idx;

}



/*****************************************************************
Find the strongest peaks by compairng their log local rank
to the number from the config file
******************************************************************/
void Spectrum::select_strong_peaks()
{
	int i;
	float thresh_log_level = log((float)(config->get_number_of_strong_peaks_per_local_window()));

	strong_peak_idxs.clear();
	for (i=0 ;i<peaks.size(); i++)
		if (peaks[i].log_local_rank<=thresh_log_level && peaks[i].iso_level == 0)
			strong_peak_idxs.push_back(i);
}





/***********************************************************************
calcs for each peak the probability of observing it at random (based on
 the neighbor's distribution. Assumes the log_intens are distributed
 according to a normal disribution.
************************************************************************/
void Spectrum::set_log_random_probs()
{
	if (peaks.size()<2)
		return;
	const float log_add = log(1.2);
	const float one_over_sqr_2pi = 1.0 / sqrt(2*3.1415927);
	const mass_t peak_window_size = 0.6; // this is fixed and independent of the tolerance!
	const mass_t margin      = 25.0;
	const mass_t window_size = 100.0;
	const mass_t min_mass = peaks[0].mass;
	const mass_t max_mass = peaks[peaks.size()-1].mass;
	const mass_t viz_range = (max_mass - min_mass);

	
	int i;
	for (i=0; i<peaks.size(); i++)
	{
		const mass_t peak_mass    = peaks[i].mass;
		const mass_t rel_position = (peak_mass-min_mass)/viz_range;
		const mass_t left_window  = margin + rel_position * window_size;
		const mass_t right_window = margin + window_size - left_window;
		const PeakRange pr = get_peaks_in_range(peak_mass-left_window,peak_mass+right_window);
		const float peak_window_prob = peak_window_size /(left_window + right_window);

		// some freak cases have 0 peak counts (only in unix)
		const int num_peaks_in_range = (pr.num_peaks>0 ? pr.num_peaks : 1);
		const float zero_prob = pow((1.0 - peak_window_prob),num_peaks_in_range); 

		if (pr.num_peaks<5)
		{
			peaks[i].log_rand_prob = log(1.0-zero_prob) + log_add;
		}
		else
		{
			vector<float> log_intens;
			int j;
			for (j=pr.low_idx; j<=pr.high_idx; j++)
				log_intens.push_back(peaks[j].log_intensity);
		
			float mean=0,sd=1;
			calc_mean_sd(log_intens,&mean,&sd);
			const float e = (peaks[i].log_intensity - mean)/sd;
			if (e<0)
			{
				peaks[i].log_rand_prob = log(1 - zero_prob) + log_add;
			}
			else
			{
				const float norm = (one_over_sqr_2pi/ sd) * exp(-0.5*e*e);
				const float norm_const = (1 - zero_prob) / (one_over_sqr_2pi/ sd); //
				peaks[i].log_rand_prob = log(norm*norm_const);
			}
		} 
	}
}
								 


/**************************************************************
Sets all peaks that have iso level=0 as strong
***************************************************************/
void Spectrum::mark_all_non_iso_peaks_as_strong()
{
	int i;
	strong_peak_idxs.clear();
	for (i=0 ;i<peaks.size(); i++)
		if (peaks[i].iso_level == 0)
			strong_peak_idxs.push_back(i);
}



/****************************************************************
Calculates different pm_with_19 values, for different charges
(calculates assumed pm_with_19 from the m/z value)
*****************************************************************/
void Spectrum::calc_original_pm_with_19_for_different_charges(
								   vector<mass_t>& pms_with_19) const
{
	int i;
	mass_t m_over_z = (org_pm_with_19 - MASS_PROTON) / charge;

	pms_with_19.resize(5);
	pms_with_19[0]=-1;
	for (i=1; i<=4; i++)
		pms_with_19[i] = m_over_z * i + MASS_PROTON;
}


/*********************************************************
Reads the spectrum info from a dta file.
returns false if there was a problem
**********************************************************/
bool Spectrum::read_from_dta( Config* _config, const char *dta_name, 
							 int _charge , bool verbose)
{
	ifstream fs(dta_name,ios::in);
	if (! fs)
		return false;

	peaks.clear();
	config = _config;
	file_name = dta_name;

	if (verbose)
		cout << "Reading: " << file_name << endl;

	char buff[1024];

	fs.getline(buff,1024);
	if (fs.bad())
		return false;

	while (buff[0] =='#') 
	{
		istringstream is(buff);
		string arg1,arg2;

		is >> arg1 >> arg2;

		if (! arg1.compare(0,5,"#FILE"))
		{
			file_name = arg2;
			if (verbose)
				cout << "#FILE " << file_name << endl;
		}
		else 
		if (! arg1.compare(0,4,"#SEQ"))
		{
			peptide.parse_from_string(config,arg2);
			if (verbose)
				cout << "#SEQ " << peptide.as_string(config) << endl;

		}
		fs.getline(buff,1024);
	}
	
	istringstream is(buff);
	is >> org_pm_with_19 >> charge;

	if (_charge>=0)
		charge = _charge;

//	charge =0;

	if (charge == 0)
	{
		m_over_z = org_pm_with_19;
	}
	else
		m_over_z = (org_pm_with_19 + (charge - 1) * 1.0078)  / charge;

	if (verbose)
		cout << org_pm_with_19 << " " << charge << endl;


	while (fs.getline(buff,1024))
	{
		istringstream is(buff);
		Peak p;
		is >> p.mass >> p.intensity;
	
		if (p.mass <0 || p.intensity==0)   // the peak probably got rounded down
			continue;

		peaks.push_back(p);
		total_original_intensity+=p.intensity;

		if (verbose)
			cout << p.mass << " " << p.intensity << endl;
	}

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; // margin 1 Dalton
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
	else
		return false;


	if (charge>=2)
	{
		max_mass =  org_pm_with_19 > max_peak_mass ? 
					org_pm_with_19 : max_peak_mass;
	}
	else
		max_mass = max_peak_mass;

	max_mass += 11.0; // margin of error

	size_idx = config->calc_size_idx(charge,org_pm_with_19);

	return true;
}





bool Spectrum::read_from_MGF_stream(Config *_config, FILE *stream, int _charge, bool verbose)
{
	config = _config;
	total_original_intensity=0;
	char buff[1024];

	peaks.clear();
	
/*	while (1)
	{
		if( ! fgets(buff, 1024, stream) )
			return false;

		if (strlen(buff)<3)
			continue;

		if (strncmp(buff,"BEGIN IONS",10) )
			continue;

		break;
	}*/
	scan_number=-1;
	retention_time=-1;
	cluster_size=-1;
	charge=-1;
	m_over_z=-1;
	org_pm_with_19=-1;

	peaks.clear();

	// read header info and first peak
	while (1)
	{
		if( ! fgets(buff, 1024, stream))
			return false;

		if (! strncmp(buff,"END IONS",7))
			return false;

		if (! strncmp(buff,"TITLE=",6) )
		{
			buff[strlen(buff)-1]='\0';
			file_name = buff+6;
			continue;
		}
		else
		if (! strncmp(buff,"SEQ=",4) )
		{
			string seq_string(buff+4);
			peptide.parse_from_string(config,seq_string);
			peptide.calc_mass(config);
			continue;		
		}
		else
		if (! strncmp(buff,"SCAN=",5) )
		{
			if (sscanf(buff+5,"%d",&scan_number) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		if (! strncmp(buff,"SCANS=",6) )
		{
			if (sscanf(buff+6,"%d",&scan_number) != 1)
			{
				cout << "Error: couldn't read scan number!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"RT=",3) )
		{
			if (sscanf(buff+3,"%f",&retention_time) != 1)
			{
				cout << "Error: couldn't read retention_time!" << endl;
				exit(1);
			}
			continue;
		}
		else
		if (! strncmp(buff,"CLUSTER_SIZE=",13) )
		{
			if (sscanf(buff+13,"%d",&cluster_size) != 1)
			{
				cout << "Error: couldn't read cluster size!" << endl;
				exit(1);
			}
			continue;
		}
		else	
		if ( ! strncmp(buff,"CHARGE=",6))
		{
			if (sscanf(buff,"CHARGE=%d",&charge) != 1)
			{
				cout <<  "Error: couldn't read charge!" << endl;
				return false;
			}
		}
		else
		if (! strncmp(buff,"PEPMASS=",8))
		{
			istringstream is(buff+8);
			is >> m_over_z;
			
			if (m_over_z < 0)
			{
				cout << "Error: reading pepmass:" << m_over_z << endl;
				return false;
			}		
		}
		else
		{
			istringstream is(buff);
			Peak p;
	
			is >> p.mass >> p.intensity;

			if (p.mass >0 && p.intensity>0)
			{
				peaks.push_back(p);
				total_original_intensity+=p.intensity;
				break;
			}
		}
	}



	if (_charge > 0)
			charge = _charge;

	org_pm_with_19 = m_over_z * charge + MASS_PROTON * (1 - charge);
	if (org_pm_with_19<0)
		org_pm_with_19=-1;



	// read rest of peaks
	
	while (fgets(buff, 1024, stream))
	{
		istringstream is(buff);
		Peak p;

		if (! strncmp("END IONS",buff,7))
			break;

		is >> p.mass >> p.intensity;
	
		if (p.mass <0 || p.intensity==0)   // the peak probably got rounded down
		{
			continue;
		}

		peaks.push_back(p);
		total_original_intensity+=p.intensity;

		if (verbose)
			cout << p.mass << " " << p.intensity << endl;
	}

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; // margin 1 Dalton
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
	else
		return false;


	max_mass =      org_pm_with_19 > max_peak_mass ? 
					org_pm_with_19 : max_peak_mass;

	max_mass += 11.0; // margin of error

	size_idx = config->calc_size_idx(charge,org_pm_with_19);

	if (peptide.get_num_aas()>0 && charge>0)
	{
		mass_t org_diff = org_pm_with_19 - peptide.get_mass() - 19.0183;

		if (fabs(org_diff)>7.0)
		{
			// try and correct charge!
			for (charge=1; charge<=4; charge++)
			{

				mass_t new_org_pm_with_19 =  m_over_z * charge + MASS_PROTON * (1 - charge);
				mass_t diff = fabs(new_org_pm_with_19 - peptide.get_mass() - MASS_OHHH);
				if (diff<6)
				{
					org_pm_with_19 = new_org_pm_with_19;
					return true;
				}
			}
		
			cout << "Error: sequence mass doesn't add up: " << file_name << " (diff: "
				 << setprecision(3) << fixed << org_diff << ")" << endl;
			cout << "Pepitde: " << peptide.as_string(config) << endl;
			cout << "Mass Cys = " << config->get_session_tables().get_aa2mass(Cys) << endl;
			exit(1);
		}
	}

	return true;
}


bool Spectrum::read_from_peak_arrays(Config* _config, Peptide& _pep, mass_t _pm_with_19, 
		int _charge, int _num_peaks, mass_t *_masses, intensity_t *_intensities)
{
	peaks.clear();
	config = _config;
	file_name = "xxx";

	peptide =_pep;
	org_pm_with_19 = _pm_with_19;
	charge = _charge;
	m_over_z = (org_pm_with_19 + (charge - 1) * MASS_PROTON ) / charge;

	peaks.resize(_num_peaks);
	int i;
	for (i=0; i<_num_peaks; i++)
	{
		peaks[i].mass = _masses[i];
		peaks[i].intensity = _intensities[i];
		total_original_intensity+=_intensities[i];
	}

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; // margin 1 Dalton
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
	else
		return false;


	if (charge>=2)
	{
		max_mass =  org_pm_with_19 > max_peak_mass ? 
					org_pm_with_19 : max_peak_mass;
	}
	else
		max_mass = max_peak_mass;

	max_mass += 11.0; // margin of error

	size_idx = config->calc_size_idx(charge,org_pm_with_19);

	return true;
}

bool Spectrum::init_from_QCPeaks(Config* _config, 
								 void *QCPeaks_ptr, 
								 int num_peaks, 
								 void *ssf_ptr, 
								 bool ind_perfrom_filtering)
{
	SingleSpectrumFile *ssf = (SingleSpectrumFile *)ssf_ptr;
	QCPeak     *QCpeaks = (QCPeak *)QCPeaks_ptr;

	peaks.clear();
	config = _config;
	file_name = "xxx";

	charge = (ssf->charge>0 ? ssf->charge : 2);
	m_over_z = ssf->m_over_z;
	org_pm_with_19 = m_over_z * charge - ( charge - 1)* MASS_PROTON;
	peptide=ssf->peptide;

	peaks.resize(num_peaks);
	int i;
	for (i=0; i<num_peaks; i++)
	{
		peaks[i].mass = QCpeaks[i].mass;
		peaks[i].intensity = QCpeaks[i].intensity;
		total_original_intensity+=peaks[i].intensity;
	}

	if (peaks.size() >0)
	{
		min_peak_mass = peaks[0].mass -1; // margin 1 Dalton
		max_peak_mass = peaks[peaks.size()-1].mass +1;
	}
	else
		return false;


//	if (charge>=2)
//	{
		max_mass =  org_pm_with_19 > max_peak_mass ? 
					org_pm_with_19 : max_peak_mass;
//	}
//	else
//		max_mass = max_peak_mass;

	max_mass += 2.0; // margin of error

	size_idx = 0;

	init_spectrum(ind_perfrom_filtering);

	return true;
}


void Spectrum::print_spectrum(ostream& os) const
{
	int i;
	if (file_name.length() > 0)
		os << "#FILE " << file_name << endl;
	if (peptide.get_num_aas()>0)
		os << "#SEQ " << peptide.as_string(config) << endl;

	os << fixed << setprecision(2) << org_pm_with_19 << " " << charge << endl;

	for (i=0; i<peaks.size(); i++)
	{
		os << left << setw(5) << i << setw(8) << setprecision(NUM_SIG_DIGITS) 
			<< fixed << right  << peaks[i].mass;
		os << setw(12) << right  << setprecision(1) << peaks[i].intensity << " " << setw(4) 
			<< setprecision(1) << peaks[i].iso_level <<  setw(4) 
			<< setprecision(1) << peaks[i].log_local_rank << "\t" << setprecision(3) << exp(peaks[i].log_rand_prob) << endl;
	}
}


void Spectrum::output_as_MGF(ostream& os) const
{
	os << "BEGIN IONS" << endl;
	os << "TITLE=" << file_name << endl;
	
	if (peptide.get_num_aas()>0)
		os << "SEQ=" << peptide.as_string(config) << endl;
	
	if (scan_number>=0 && cluster_size == 1)
		os << "SCANS=" << scan_number << endl;
	
//	if (retention_time >=0)
//		os << "RT=" << retention_time << endl;

	if (cluster_size>1)
		os << "CLUSTER_SIZE=" << cluster_size << endl;

	os << "CHARGE=+" << charge << endl;

	os << "PEPMASS=" << fixed << setprecision(NUM_SIG_DIGITS) << m_over_z << endl;

	int i;
	for (i=0; i<peaks.size(); i++)
		os << fixed << setprecision(NUM_SIG_DIGITS) << peaks[i].mass << " " 
		   << fixed << setprecision(NUM_SIG_DIGITS) << peaks[i].intensity << endl;
	
	os << "END IONS" << endl << endl;
}



void Spectrum::print_expected_by(ostream& os) const
{
	vector<string> labels;
	labels.push_back("b");
	labels.push_back("y");
	print_expected_fragment_peaks(labels,os);
}

void Spectrum::print_expected_fragment_peaks(vector<string>& frag_labels, ostream& os) const
{
	int i;
	const vector<FragmentType>& all_fragments = config->get_all_fragments();
	vector<mass_t> break_masses;
	
	const int frag_label_width =20;

	if (peptide.get_num_aas() == 0)
		return;

	vector<string> pre_frags, suf_frags;
	vector<int> pre_frag_idxs, suf_frag_idxs;
	
	for (i=0; i<frag_labels.size(); i++)
	{
		if (frag_labels[i].length()>0 &&
			(frag_labels[i][0] == 'a' || frag_labels[i][0]== 'b' || frag_labels[i][0] == 'c') )
		{
			pre_frags.push_back(frag_labels[i]);
			pre_frag_idxs.push_back(config->get_frag_idx_from_label(frag_labels[i]));
		}
		else
		{
			suf_frags.push_back(frag_labels[i]);
			suf_frag_idxs.push_back(config->get_frag_idx_from_label(frag_labels[i]));
		}
	}


	peptide.calc_expected_breakage_masses(config,break_masses);
	mass_t true_mass_with_19 = peptide.get_mass() + MASS_OHHH;
	os << peptide.as_string(config) << " (" << fixed << setprecision(4) << true_mass_with_19 << ")" << endl;

	if ( suf_frags.size() == 0)
	{
		os << setw(frag_label_width*frag_labels.size()+3)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<frag_labels.size(); i++)
		{
			int w = (frag_label_width+6- frag_labels[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-w-1) << left << frag_labels[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+3)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		
		for (i=0; i<break_masses.size(); i++)
		{
			int region = config->calc_region_idx(break_masses[i],
				true_mass_with_19, charge, min_peak_mass, max_peak_mass);

			const RegionalFragments& rf = config->get_regional_fragments(charge,size_idx,region);

			if (rf.get_frag_type_idxs().size() == 0)
			{
				cout << "Error: no fragments selected for region " << charge << " " <<
						size_idx << " " << region << endl;
				exit(1);
			}

			cout << "|" << setw(2) << left << i <<"|";
			
			int j;
			for (j=0; j<frag_labels.size(); j++)
			{
				int idx=pre_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;
				
				if (idx<0)
				{
					os << setw(frag_label_width-2) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=get_max_inten_peak(exp_mass,config->get_tolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}
			os << endl;
		}
		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}
	else
	if ( pre_frags.size() == 0)
	{
		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<frag_labels.size(); i++)
		{
			int w = (frag_label_width - 1- frag_labels[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-1-w) << left << frag_labels[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		for (i=0; i<break_masses.size(); i++)
		{
			int region = config->calc_region_idx(break_masses[i],
				true_mass_with_19, charge, min_peak_mass, max_peak_mass);

			const RegionalFragments& rf = config->get_regional_fragments(charge,size_idx,region);

			cout << "|" << setw(2) << left << break_masses.size() - i - 1<<"|";
			
			int j;
			for (j=0; j<frag_labels.size(); j++)
			{
				int idx=suf_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=get_max_inten_peak(exp_mass,config->get_tolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}
			os << endl;
		}

		os << setw(frag_label_width*frag_labels.size()+5)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}
	else // have both prefix and suffix fragments, separate between them
	{
		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		os << "|  |";
		for (i=0; i<pre_frags.size(); i++)
		{
			int w = (frag_label_width - 2 - pre_frags[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-2-w) << left << pre_frags[i] << "|";
		}
		os << "   |";
		for (i=0; i<suf_frags.size(); i++)
		{
			int w = (frag_label_width - 2 - suf_frags[i].length())/2;
			os << setw(w) << " " << setw(frag_label_width-2-w) << left << suf_frags[i] << "|";
		}
		os<<endl;

		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');

		for (i=0; i<break_masses.size(); i++)
		{
			int region = config->calc_region_idx(break_masses[i],
				true_mass_with_19, charge, min_peak_mass, max_peak_mass);

			const RegionalFragments& rf = config->get_regional_fragments(charge,size_idx,region);


			cout << "|" << setw(2) << left << i <<"|";
			
			int j;
			for (j=0; j<pre_frags.size(); j++)
			{
				int idx=pre_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=get_max_inten_peak(exp_mass,config->get_tolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}

			// output suffix fragments

			cout << " " << setw(2) << left << break_masses.size() - i - 1<<"|";
			for (j=0; j<suf_frags.size(); j++)
			{
				int idx=suf_frag_idxs[j];
				if (rf.get_position_of_frag_type_idx(idx)<0)
					idx=-1;

				if (idx<0)
				{
					os << setw(13) <<  left << "      - ";
				}
				else
				{
					mass_t exp_mass = all_fragments[idx].calc_expected_mass(break_masses[i],
										true_mass_with_19);
					int p_idx=get_max_inten_peak(exp_mass,config->get_tolerance());
					os << " " << setw(6) << right << fixed << setw(9) << setprecision(3) << exp_mass;
					if (p_idx>=0)
					{
						os << " " << setw(6) << fixed << setprecision(3) << peaks[p_idx].mass - exp_mass;
					}
					else
						os << "       ";
				}
				os << " |";
			}


			os << endl;
		}

		os << setw(frag_label_width*frag_labels.size()+7)<< setfill('-') << right << " " << endl;
		os << setfill(' ');
	}

	// print pm stats

	

	os << "true - org_pm:    " << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		org_pm_with_19 << "  = " << true_mass_with_19 - org_pm_with_19 << endl;
	
	if (this->corrected_pm_with_19>0)
	{
		os << "true - cor1_pm:   "  << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		corrected_pm_with_19 << "  = " << true_mass_with_19 - corrected_pm_with_19 << endl;
	}
	if (this->secondary_pm_with_19>0)
	{
		os << "true - cor2_pm: "  << fixed << setprecision(4) << true_mass_with_19 << " - " <<
		secondary_pm_with_19 << "  = " << true_mass_with_19 - secondary_pm_with_19 << endl;
	}
	os << endl;
	
}


// checks several charges to see which one has a good m over z match
// writes the expected b/y ions if finds a match
bool Spectrum::check_m_over_z_and_sequence(ostream& os)
{
	int c;
	mass_t pep_mass = peptide.get_mass() + MASS_OHHH;

	for (c=1; c<=4; c++)
	{
		mass_t pm_19 = m_over_z * c - (c-1)*MASS_PROTON;

	//	cout <<"c: " <<c << "   pm19: " << pm_19 << endl;

		if (fabs(pm_19-pep_mass)<7)
		{
			set_org_pm_with_19(pm_19);
			set_charge(c);
		//	print_expected_by(os);
		//	os << endl;
			return true;
		}
	}

	for (c=1; c<=4; c++)
	{
		mass_t pm_19 = m_over_z * c - (c-1)*MASS_PROTON;

		cout <<"c: " <<c << "   pm19: " << pm_19 << endl;
	}

	os << "No charge matches: " << peptide.as_string(config) << "  (" <<
		peptide.get_mass() + 19.083 << ") !!!! " << endl << endl;

	return false;
}


ostream& operator << (ostream& os, const Spectrum& spec)
{
	const vector<Peak>& peaks = spec.get_peaks();
	const string& file_name = spec.get_file_name();
	const Peptide& peptide = spec.get_peptide();
	const Config *config = spec.get_config();

	if (file_name.length() > 0)
		os << "#FILE " << file_name << endl;
	if (peptide.get_num_aas()>0)
		os << "#SEQ " << peptide.as_string(config) << endl;

	os << spec.get_org_pm_with_19() << " " << spec.get_charge() << endl;
	int i;
	for (i=0; i<peaks.size(); i++)
		os << setw(8) << left << peaks[i].mass << " \t" << peaks[i].intensity << 
		" \t" << exp(peaks[i].log_rand_prob)<< endl;

	return os;
}


/*****************************************************************
 Returns the indices of all peaks that are within the mass range
 ******************************************************************/
PeakRange  Spectrum::get_peaks_in_range(mass_t min_r, mass_t max_r) const
{
	PeakRange pr;
	const int max_peak_idx = peaks.size()-1;
    
	if (max_r>= max_mass+0.1)
		max_r = max_mass+0.1;
    
	if (min_r<0)
		min_r=0;
    
	if (max_r<=min_r)
		return pr;
    
	int i_min=index_array[(int)(min_r)];
	int i_max=index_array[(int)(max_r)];
    
	if (i_max<max_peak_idx)
		i_max++;
    
	while (peaks[i_min].mass<min_r && i_min<i_max)
		i_min++;
    
	while (peaks[i_max].mass>max_r && i_max>0)
		i_max--;
    
	if (peaks[i_min].mass > max_r || peaks[i_max].mass<min_r)
		return pr;
    
	pr.num_peaks = i_max - i_min+1;
	pr.low_idx = i_min;
	pr.high_idx = i_max;
    
	return pr;
}


/***********************************************************
 returns the idx of the peak that is closes to the mass
 (-1 is returned if no peak is found)
 ************************************************************/
int Spectrum::get_min_dis_peak(mass_t mass, mass_t tolerance) const
{
	PeakRange pr = get_peaks_in_range(mass-tolerance,mass+tolerance);
	if (pr.num_peaks==0)
		return -1;
	mass_t min_dis=fabs(mass - peaks[pr.low_idx].mass);
	int min_idx=pr.low_idx;
	int j;
    
	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t dis = fabs(mass - peaks[pr.low_idx+j].mass);
		if (dis<min_dis)
		{
			min_dis=dis;
			min_idx = pr.low_idx+j;
		}
	}
	return min_idx;
}

/***********************************************************
 returns the idx of the peak that has the highest intensity
 in the range (-1 is returned if no peak is found)
 Balances between intensity and proximity to the expected position
 A peak at the edge needs to be at least 2 times stronger than
 a peak exactly at the middle to be chosen.
 ************************************************************/
int Spectrum::get_max_inten_peak(mass_t exp_mass, mass_t tolerance) const
{
	PeakRange pr = get_peaks_in_range(exp_mass-tolerance,exp_mass+tolerance);
	if (pr.num_peaks==0)
		return -1;
    
	if (pr.num_peaks ==1)
		return pr.low_idx;
    
	mass_t max_val = (2*tolerance - fabs(peaks[pr.low_idx].mass - exp_mass)) *
    peaks[pr.low_idx].intensity;
    
	int max_idx=pr.low_idx;
	int j;
	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t peak_val = (2*tolerance - fabs(peaks[pr.low_idx].mass - exp_mass)) *
        peaks[pr.low_idx].intensity;
		if (peak_val>max_val)
		{
			max_val=peak_val;
			max_idx = pr.low_idx+j;
		}
		else
			break;
	}
	return max_idx;
}



// returns intensities at positions -2,-1,+1,+2
void Spectrum::get_iso_intens(int p_idx, vector<float>& iso_intens, mass_t iso_tolerance, int charge) const
{
	const mass_t one_over = 1.0/charge;
	const mass_t one_proton = MASS_PROTON*one_over;
	const mass_t two_proton = 2*one_proton;
	const mass_t p_mass=peaks[p_idx].mass;
	const mass_t min_minus2_mass =  p_mass - two_proton - iso_tolerance;
	const mass_t min_minus1_mass =  p_mass - one_proton - iso_tolerance;
	const mass_t min_plus1_mass  =  p_mass + one_proton - iso_tolerance;
	const mass_t min_plus2_mass  =  p_mass + two_proton - iso_tolerance;
	const int num_peaks = peaks.size();
    
	iso_intens.clear();
	iso_intens.resize(4,NEG_INF);
    
	int p=p_idx;
	while (p>=0 && peaks[p].mass>min_minus2_mass)
		p--;
    
	p++;
    
	if (fabs(peaks[p].mass- p_mass + two_proton)<=iso_tolerance)
		iso_intens[0]=peaks[p].log_intensity;
    
	while (p<p_idx && peaks[p].mass<min_minus1_mass)
		p++;
    
	if (fabs(peaks[p].mass - p_mass + one_proton)<=iso_tolerance)
		iso_intens[1]=peaks[p].log_intensity;
    
	p=p_idx+1;
    
	while (p<num_peaks && peaks[p].mass<min_plus1_mass)
		p++;
    
	if (p<num_peaks && fabs(peaks[p].mass - p_mass - one_proton)<=iso_tolerance)
		iso_intens[2]=peaks[p].log_intensity;
    
	while (p<num_peaks && peaks[p].mass<min_plus2_mass)
		p++;
    
	if (p<num_peaks && fabs(peaks[p].mass - p_mass - two_proton)<=iso_tolerance)
		iso_intens[3]=peaks[p].log_intensity;
}





