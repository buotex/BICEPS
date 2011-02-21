#ifndef __INCLUDES_H__
#define __INCLUDES_H__


#pragma warning (disable:4786)
#pragma warning (disable:4305)
#pragma warning (disable:4503)

#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <time.h>



using namespace std;
using std::vector;



#define BYTEORDER_LITTLE_ENDIAN

#define NEG_INF -999999999
#define POS_INF  999999999

#define MASS_PROTON 1.00728
#define MASS_ISO 1.0033
#define MASS_ELECTRON 0.00055
#define MASS_H2O 18.01056 
#define MASS_NH3 17.02655
#define MASS_CO	 27.99492
#define MASS_H2ONH3 35.03711
#define MASS_H2OH2O 36.02112
#define MASS_NH3NH3 34.0531
#define MASS_OHHH   19.0184

#define NUM_SIG_DIGITS 3

#define ESI_MASS_SPEC 1

typedef enum AminoAcids {N_TERM, C_TERM, Gap,Xle,Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,
						 Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val} AminoAcids;


// Data types for common variables

typedef float  mass_t;
typedef float  intensity_t;
typedef float  score_t;
typedef double weight_t;

typedef map< string , int , less<string> > STRING2INT_MAP;
typedef map< mass_t , mass_t, less<mass_t> > MASS_T_MAP;
typedef map< int , int , less<int> > INT_MAP;
typedef map< mass_t , int   , less<mass_t> > MASS_T2INT_MAP;


#endif


