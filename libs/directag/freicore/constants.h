#ifndef _CONSTANTS_H
#define _CONSTANTS_H

// constants derived from http://www.lbl.gov/abc/wallchart/chapters/appendix/appendixc.html

#define PROTON				1.007276f
#define NEUTRON				1.008665f
#define ELECTRON			0.000549f

#define NUM_ISOTOPES		5
extern const float CARBON_ISOTOPES[NUM_ISOTOPES];
extern const float HYDROGEN_ISOTOPES[NUM_ISOTOPES];
extern const float OXYGEN_ISOTOPES[NUM_ISOTOPES];
extern const float NITROGEN_ISOTOPES[NUM_ISOTOPES];
extern const float SULFUR_ISOTOPES[NUM_ISOTOPES];

//extern const float C12_MASS;
//extern const float SCALE_FACTOR;
extern const float CARBON_MONO;
extern const float CARBON_AVG;
extern const float HYDROGEN_MONO;
extern const float HYDROGEN_AVG;
extern const float OXYGEN_MONO;
extern const float OXYGEN_AVG;
extern const float NITROGEN_MONO;
extern const float NITROGEN_AVG;
extern const float SULFUR_MONO;
extern const float SULFUR_AVG;

extern const float WATER_MONO;
extern const float WATER_AVG;
#define WATER(useAvg) (useAvg?WATER_AVG:WATER_MONO)

extern const float AMMONIA_MONO;
extern const float AMMONIA_AVG;
#define AMMONIA(useAvg) (useAvg?AMMONIA_AVG:AMMONIA_MONO)

#define PEPTIDE_N_TERMINUS_SYMBOL	'('
#define PEPTIDE_C_TERMINUS_SYMBOL	')'
extern const char PEPTIDE_N_TERMINUS_STRING[2];
extern const char PEPTIDE_C_TERMINUS_STRING[2];

#define PROTEIN_N_TERMINUS_SYMBOL	'['
#define PROTEIN_C_TERMINUS_SYMBOL	']'
extern const char PROTEIN_N_TERMINUS_STRING[2];
extern const char PROTEIN_C_TERMINUS_STRING[2];

#endif
