#ifndef __BICEPS_INCLUDE_BICEPSDEFINITIONS_H__
#define __BICEPS_INCLUDE_BICEPSDEFINITIONS_H__
#include <string>
#define PEN_PTM 1.2f
#define PEN_SNP 2.0f
#define PEN_METHOX 0.3f
#define PEN_TRYPTIC 1.2f
#define PEN_GENOMIC 1.2f
#define PEN_SPLICED 2.0f
#define PEN_MISSCLEAV 0.3f
#define TUPLEBUFFER 50000 //Amount of tuples till they get scored
#define POOLSIZE 51000 //Size of the memory pool for the tuples, has to be > tuplebuffer as extra space has to be reserved for the Tuple-Matches:w
#define MAXSEQUENCESIZE 142 //how long a resulting sequence can be.
namespace biceps{
static std::string bicepsconfigpath("${BICEPS_CONFIG_PATH}");
}
#endif
