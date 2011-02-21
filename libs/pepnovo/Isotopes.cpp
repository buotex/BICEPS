#include "Isotopes.h"

// These are the ratios for isotpic peaks
// each line represents 100 Da (starting from 0 (or 12 to be precise)

const int MAX_ISO_ROW = 35; 
const float AVG_ratios[MAX_ISO_ROW+1][6]={
{1.0000,0.0112,0.0000,0.0000,0.0000,0.0000},
{1.0000,0.0521,0.0032,0.0001,0.0000,0.0000},
{1.0000,0.1116,0.0118,0.0009,0.0000,0.0000},
{1.0000,0.1653,0.0210,0.0020,0.0001,0.0000},
{1.0000,0.2269,0.0348,0.0041,0.0004,0.0000},
{1.0000,0.2770,0.0512,0.0071,0.0008,0.0001},
{1.0000,0.3385,0.0719,0.0114,0.0014,0.0001},
{1.0000,0.3923,0.0933,0.0166,0.0024,0.0003},
{1.0000,0.4517,0.1221,0.0245,0.0040,0.0005},
{1.0000,0.5038,0.1489,0.0323,0.0058,0.0008},
{1.0000,0.5557,0.1782,0.0418,0.0079,0.0012},
{1.0000,0.6153,0.2168,0.0557,0.0116,0.0021},
{1.0000,0.6689,0.2531,0.0695,0.0153,0.0028},
{1.0000,0.7304,0.2979,0.0876,0.0206,0.0041},
{1.0000,0.7837,0.3860,0.1418,0.0421,0.0104},
{1.0000,0.8452,0.4380,0.1687,0.0522,0.0136},
{1.0000,0.8990,0.4867,0.1948,0.0629,0.0171},
{1.0000,0.9583,0.5454,0.2291,0.0775,0.0221},
{1.0000,1.0105,0.5987,0.2607,0.0911,0.0270},
{1.0000,1.0626,0.6543,0.2951,0.1069,0.0326},
{1.0000,1.1223,0.7234,0.3403,0.1284,0.0408},
{1.0000,1.1759,0.7867,0.3828,0.1489,0.0488},
{1.0000,1.2371,0.8627,0.4355,0.1756,0.0596},
{1.0000,1.2875,0.9299,0.4854,0.2022,0.0705},
{1.0000,1.3494,1.0128,0.5474,0.2355,0.0849},
{1.0000,1.4026,1.0884,0.6060,0.2683,0.0996},
{1.0000,1.4621,1.1772,0.6792,0.3108,0.1191},
{1.0000,1.5139,1.2565,0.7451,0.3497,0.1374},
{1.0000,1.5756,1.3531,0.8278,0.4006,0.1620},
{1.0000,1.6256,1.4370,0.9037,0.4490,0.1863},
{1.0000,1.6795,1.5277,0.9861,0.5024,0.2139},
{1.0000,1.7413,1.6348,1.0866,0.5693,0.2484},
{1.0000,1.7908,1.7265,1.1774,0.6317,0.2823},
{1.0000,1.8528,1.8405,1.2908,0.7114,0.3261},
{1.0000,1.9069,1.9436,1.3959,0.7864,0.3684},
{1.0000,1.9659,2.0626,1.5215,0.8806,0.4240} };


// returns the ratios of the expected isotopic peaks
// bases numbers on the averagine statistics
void calc_expected_iso_ratios(mass_t peak_mass, vector<float>& ratios,
						int num_ratios)
{
	int i;
	int i1;
	float w1,w2;

	i1 = (int)(peak_mass * 0.01);

	if (i1>MAX_ISO_ROW)
	{
		i1=MAX_ISO_ROW-1;
		w1=0;
		w2=1;
	}
	else
	{
		w1  = 1.0 - (peak_mass - i1*100.0)*0.01;
		w2  = 1.0 - w1;
	}

	ratios.resize(num_ratios+1);
	ratios[0]=1.0;
	for (i=1; i<=num_ratios; i++)
		ratios[i]=w1*AVG_ratios[i1][i] + w2*AVG_ratios[i1+1][i];	
}


