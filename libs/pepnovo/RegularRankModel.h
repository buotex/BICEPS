#ifndef __REGULARRANKMODEL_H__
#define __REGULARRANKMODEL_H__

#include "DiscretePeakModel.h"
#include "AnnotatedSpectrum.h"
#include "FileManagement.h"

#include "includes.h"

class RegularRankModel : public DiscretePeakModel {
public:

	void read_model_specifics(istream& is);

	void write_model_specifics(ostream& os) const;

	void init_model_for_scoring_spectrum(class Spectrum *spec) { current_spectrum = spec; };

	void init_default_thresholds(); 

	void set_breakage_peak_levels(Breakage *breakage) const
	{
		int i;
		for (i=0; i<breakage->fragments.size(); i++)
			breakage->fragments[i].peak_level = get_peak_level(breakage->fragments[i].peak_idx);
	}

	int get_peak_level(int p_idx, int f_idx=-1) const
	{
		const int rank = current_spectrum->get_peak_rank(p_idx);
		int i;
		for (i=0; i<level_thresholds.size(); i++)
			if (rank<=level_thresholds[i])
				return i;
		return level_thresholds.size();
	}


	void print_level_legend( ostream& os) const;

private:

	int get_lowest_level_in_spectrum(Spectrum *spec) const
	{
		const int num_peaks = spec->get_num_peaks();
		int i;
		for (i=0; i<level_thresholds.size(); i++)
			if (num_peaks<=level_thresholds[i])
				return i;
		return level_thresholds.size();
	}


};


#endif

