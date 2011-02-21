#include "Peak.h"
namespace Pepsplice{
Peak::Peak(int bin0, float mz0, float intensity0)
{
	bin = bin0;
	mz = mz0;
	intensity = intensity0;

	isotopesim = 0;
	intensityNormLocAvg = 0;
	intensityNormInvert = 0;
	above_isotope_cutoff = true;
	above_intensity_cutoff = true;
	above_inv_intensity_cutoff = true;
	above_glob_intensity_cutoff = true;
	
}

Peak::~Peak()
{
}
}