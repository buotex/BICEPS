#ifndef __BREAKAGE_SCORE_MODEL_H__
#define __BREAKAGE_SCORE_MODEL_H__

#include "ME_REG.h"
#include "BasicDataStructs.h"
#include "Spectrum.h"
#include "FileManagement.h"
#include "Model.h"
#include "PrmGraph.h"
#include "includes.h"


// model for strong fragment with intensity
typedef enum ScoreModelFields_SI {

	SI_CONST,                          
	
	SI_IND_MIRROR1_VIZ,			SI_IND_MIRROR1_NOT_VIZ,
	SI_IND_HAS_MIRROR1_INTEN,   SI_MIRROR1_ISO_LEVEL,
	SI_IND_MIRROR1_HAS_MINUS_1,	SI_MIRROR1_MINUS_1_INTEN_DIFF, // with mirror
	SI_IND_MIRROR1_HAS_MINUS_2, SI_MIRROR1_MINUS_2_INTEN_DIFF,
	SI_IND_MIRROR1_HAS_PLUS_1,	SI_MIRROR1_PLUS_1_INTEN_DIFF,
	SI_IND_MIRROR1_HAS_PLUS_2,	SI_MIRROR1_PLUS_2_INTEN_DIFF,
	SI_MIRROR1_MASS_DIFF25,     SI_MIRROR1_MASS_DIFF75,		SI_MIRROR1_MASS_DIFF_LARGE,

	SI_IND_MIRROR2_VIZ,			SI_IND_MIRROR2_NOT_VIZ,
	SI_IND_HAS_MIRROR2_INTEN,	SI_MIRROR2_ISO_LEVEL, 
	SI_IND_MIRROR2_HAS_MINUS_1,	SI_MIRROR2_MINUS_1_INTEN_DIFF,
	SI_IND_MIRROR2_HAS_MINUS_2, SI_MIRROR2_MINUS_2_INTEN_DIFF,
	SI_IND_MIRROR2_HAS_PLUS_1,	SI_MIRROR2_PLUS_1_INTEN_DIFF,
	SI_IND_MIRROR2_HAS_PLUS_2,	SI_MIRROR2_PLUS_2_INTEN_DIFF,
	SI_MIRROR2_MASS_DIFF25,		SI_MIRROR2_MASS_DIFF75,		SI_MIRROR2_MASS_DIFF_LARGE,

	SI_IND_PARENT1_VIZ, SI_IND_PARENT1_NOT_VIZ,	SI_IND_PARENT1_NO_INTEN,	
	SI_PARENT1_ISO_LEVEL,
	SI_IND_PARENT1_LESS_THAN_100_MIN_MAX, SI_IND_PARENT1_LESS_THAN_200_MIN_MAX,   
	SI_IND_PARENT1_INTEN_MORE, SI_PARENT1_INTEN_DIFF_MORE,
	SI_IND_PARENT1_INTEN_LESS, SI_PARENT1_INTEN_DIFF_LESS,
	
	SI_IND_PARENT2_VIZ, SI_IND_PARENT2_NOT_VIZ,	SI_IND_PARENT2_NO_INTEN,	
	SI_PARENT2_ISO_LEVEL,
	SI_IND_PARENT2_LESS_THAN_100_MIN_MAX, SI_IND_PARENT2_LESS_THAN_200_MIN_MAX,   
	SI_IND_PARENT2_INTEN_MORE, SI_PARENT2_INTEN_DIFF_MORE,
	SI_IND_PARENT2_INTEN_LESS, SI_PARENT2_INTEN_DIFF_LESS,
	
	SI_LOG_LOCAL_RANK,		 SI_LOG_GLOBAL_RANK,  SI_ISO_LEVEL,
	SI_IND_HAS_MINUS_1,		 SI_IND_HAS_MINUS_2,
	SI_IND_HAS_PLUS_1,		 SI_IND_HAS_PLUS_2,
	
	SI_IND_LOG_INTEN_LESS1,  SI_LOG_INTEN_LESS1,
	SI_IND_LOG_INTEN_LESS2,  SI_LOG_INTEN_LESS2,
	SI_IND_LOG_INTEN_LESS3,  SI_LOG_INTEN_LESS3,
	SI_IND_LOG_INTEN_LESS4,  SI_LOG_INTEN_LESS4,
	SI_IND_LOG_INTEN_MORE,   SI_LOG_INTEN_MORE,

	SI_IND_DIS_FROM_MINMAX_LESS_50,		SI_DIS_FROM_MINMAX0,	SI_LOG_INTEN_DIS_50,
	SI_IND_DIS_FROM_MINMAX_LESS_150,	SI_DIS_FROM_MINMAX50,	SI_LOG_INTEN_DIS_150,
	SI_IND_DIS_FROM_MINMAX_LESS_250,	SI_DIS_FROM_MINMAX150,	SI_LOG_INTEN_DIS_250,
	SI_IND_DIS_FROM_MINMAX_MORE,		SI_DIS_FROM_MINMAX250,	SI_LOG_INTEN_DIS_MORE,

	SI_REL_POS0, SI_REL_POS1, SI_REL_POS2, SI_REL_POS3, SI_REL_POS4,
	SI_REL_POS5, SI_REL_POS6, SI_REL_POS7, SI_REL_POS8, SI_REL_POS9,

	// check for forward offsets from this frag 
	SI_IND_HAS_PLUS_NH3,	SI_PLUS_NH3_INTEN_DIFF,
	SI_IND_HAS_PLUS_H2O,	SI_PLUS_H2O_INTEN_DIFF,
	SI_IND_HAS_PLUS_CO,		SI_PLUS_CO_INTEN_DIFF,
	SI_IND_HAS_PLUS_H2ONH3,	SI_PLUS_H2ONH3_INTEN_DIFF,
	SI_IND_HAS_PLUS_H2OH2O,	SI_PLUS_H2OH2O_INTEN_DIFF,

	SI_IND_HAS_CHARGE_PLUS1,	SI_CHARGE_PLUS1_INTEN_DIFF,
	SI_IND_HAS_CHARGE_MINUS1,	SI_CHARGE_MINUS1_INTEN_DIFF,

	// variable aa features disntiguish between single aa edge and multiple aa edge

	SI_IND_N_IS_GAP,	SI_IND_C_IS_GAP,

	SI_N_TERM_CAT20, SI_N_TERM_CAT18, SI_N_TERM_CAT12, SI_N_TERM_CAT8, SI_N_TERM_CAT4, SI_N_TERM_CAT2,
	SI_N_EDGE_CAT20, SI_N_EDGE_CAT18, SI_N_EDGE_CAT12, SI_N_EDGE_CAT8, SI_N_EDGE_CAT4, SI_N_EDGE_CAT2,
	SI_C_TERM_CAT20, SI_C_TERM_CAT18, SI_C_TERM_CAT12, SI_C_TERM_CAT8, SI_C_TERM_CAT4, SI_C_TERM_CAT2,
	SI_C_EDGE_CAT20, SI_C_EDGE_CAT18, SI_C_EDGE_CAT12, SI_C_EDGE_CAT8, SI_C_EDGE_CAT4, SI_C_EDGE_CAT2,
	SI_SPAN_CAT20,	 SI_SPAN_CAT18,   SI_SPAN_CAT12,   SI_SPAN_CAT8,   SI_SPAN_CAT4,   SI_SPAN_CAT2,
	SI_ND_SPAN_CAT20, SI_ND_SPAN_CAT18, SI_ND_SPAN_CAT12, SI_ND_SPAN_CAT8, SI_ND_SPAN_CAT4, SI_ND_SPAN_CAT2,
	SI_CD_SPAN_CAT20, SI_CD_SPAN_CAT18, SI_CD_SPAN_CAT12, SI_CD_SPAN_CAT8, SI_CD_SPAN_CAT4, SI_CD_SPAN_CAT2,

	SI_IND_CONNECTS_TO_N_TERM, SI_IND_PREFERRED_DIGEST_AA_N_TERM, SI_PREFERRED_DIGEST_AA_N_TERM_INTEN, SI_NON_PREFERRED_DIGEST_AA_N_TERM_INTEN,
	SI_IND_CONNECTS_TO_C_TERM, SI_IND_PREFERRED_DIGEST_AA_C_TERM, SI_PREFERRED_DIGEST_AA_C_TERM_INTEN, SI_NON_PREFERRED_DIGEST_AA_C_TERM_INTEN,

	SI_IND_NOT_CONNECTED_TO_TERMS,   SI_IND_MISSED_CLEAVAGE,

	SI_IND_N_FRAG_NOT_VIZ,			 SI_IND_C_FRAG_NOT_VIZ,
	SI_IND_N_INTEN,					 SI_IND_C_INTEN,
	SI_IND_N_NO_INTEN,				 SI_IND_C_NO_INTEN,

	// N-term features
	SI_IND_N_N_TERM, SI_IND_N_C_TERM, SI_IND_N_Gap, SI_IND_N_Xle, SI_IND_N_Ala, SI_IND_N_Arg, SI_IND_N_Asn,
	SI_IND_N_Asp,    SI_IND_N_Cys,    SI_IND_N_Gln, SI_IND_N_Glu, SI_IND_N_Gly, SI_IND_N_His, SI_IND_N_Ile,
	SI_IND_N_Leu,    SI_IND_N_Lys,    SI_IND_N_Met, SI_IND_N_Phe, SI_IND_N_Pro, SI_IND_N_Ser, SI_IND_N_Thr,
	SI_IND_N_Trp,    SI_IND_N_Tyr,    SI_IND_N_Val,

	// self inten
	SI_IND_N_N_TERM_INTEN, SI_IND_N_C_TERM_INTEN, SI_IND_N_Gap_INTEN, SI_IND_N_Xle_INTEN, SI_IND_N_Ala_INTEN, SI_IND_N_Arg_INTEN, SI_IND_N_Asn_INTEN,
	SI_IND_N_Asp_INTEN,    SI_IND_N_Cys_INTEN,    SI_IND_N_Gln_INTEN, SI_IND_N_Glu_INTEN, SI_IND_N_Gly_INTEN, SI_IND_N_His_INTEN, SI_IND_N_Ile_INTEN,
	SI_IND_N_Leu_INTEN,    SI_IND_N_Lys_INTEN,    SI_IND_N_Met_INTEN, SI_IND_N_Phe_INTEN, SI_IND_N_Pro_INTEN, SI_IND_N_Ser_INTEN, SI_IND_N_Thr_INTEN,
	SI_IND_N_Trp_INTEN,    SI_IND_N_Tyr_INTEN,    SI_IND_N_Val_INTEN,

	// The following are only filled if the N are in the visible range
	SI_IND_N_SE, 

	SI_SE_IND_HAS_N_FRAG_INTEN, 
	SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN, SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN, SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,
	SI_SE_IND_N_FRAG_DIFF_01,  SI_SE_IND_N_FRAG_DIFF_05,  SI_SE_IND_N_FRAG_DIFF_LARGE,

	SI_SE_IND_N_N_TERM_DIFF_INTEN, SI_SE_IND_N_C_TERM_DIFF_INTEN, SI_SE_IND_N_Gap_DIFF_INTEN, SI_SE_IND_N_Xle_DIFF_INTEN, SI_SE_IND_N_Ala_DIFF_INTEN, SI_SE_IND_N_Arg_DIFF_INTEN, SI_SE_IND_N_Asn_DIFF_INTEN,
	SI_SE_IND_N_Asp_DIFF_INTEN,    SI_SE_IND_N_Cys_DIFF_INTEN,    SI_SE_IND_N_Gln_DIFF_INTEN, SI_SE_IND_N_Glu_DIFF_INTEN, SI_SE_IND_N_Gly_DIFF_INTEN, SI_SE_IND_N_His_DIFF_INTEN, SI_SE_IND_N_Ile_DIFF_INTEN,
	SI_SE_IND_N_Leu_DIFF_INTEN,    SI_SE_IND_N_Lys_DIFF_INTEN,    SI_SE_IND_N_Met_DIFF_INTEN, SI_SE_IND_N_Phe_DIFF_INTEN, SI_SE_IND_N_Pro_DIFF_INTEN, SI_SE_IND_N_Ser_DIFF_INTEN, SI_SE_IND_N_Thr_DIFF_INTEN,
	SI_SE_IND_N_Trp_DIFF_INTEN,    SI_SE_IND_N_Tyr_DIFF_INTEN,    SI_SE_IND_N_Val_DIFF_INTEN,

	SI_SE_IND_HAS_NO_N_FRAG_INTEN,
	SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN, SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN, SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	SI_IND_N_ME, 
	SI_ME_IND_HAS_N_FRAG_INTEN, 
	SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN, SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN, SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,
	SI_ME_IND_N_FRAG_DIFF_01,  SI_ME_IND_N_FRAG_DIFF_05,  SI_ME_IND_N_FRAG_DIFF_LARGE,

	SI_ME_IND_HAS_NO_N_FRAG_INTEN,
	SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN, SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN, SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	// C-term features

	SI_IND_C_N_TERM, SI_IND_C_C_TERM, SI_IND_C_Gap, SI_IND_C_Xle, SI_IND_C_Ala, SI_IND_C_Arg, SI_IND_C_Asn,
	SI_IND_C_Asp,    SI_IND_C_Cys,    SI_IND_C_Gln, SI_IND_C_Glu, SI_IND_C_Gly, SI_IND_C_His, SI_IND_C_Ile,
	SI_IND_C_Leu,    SI_IND_C_Lys,    SI_IND_C_Met, SI_IND_C_Phe, SI_IND_C_Pro, SI_IND_C_Ser, SI_IND_C_Thr,
	SI_IND_C_Trp,    SI_IND_C_Tyr,    SI_IND_C_Val,

	// self inten
	SI_IND_C_N_TERM_INTEN, SI_IND_C_C_TERM_INTEN, SI_IND_C_Gap_INTEN, SI_IND_C_Xle_INTEN, SI_IND_C_Ala_INTEN, SI_IND_C_Arg_INTEN, SI_IND_C_Asn_INTEN,
	SI_IND_C_Asp_INTEN,    SI_IND_C_Cys_INTEN,    SI_IND_C_Gln_INTEN, SI_IND_C_Glu_INTEN, SI_IND_C_Gly_INTEN, SI_IND_C_His_INTEN, SI_IND_C_Ile_INTEN,
	SI_IND_C_Leu_INTEN,    SI_IND_C_Lys_INTEN,    SI_IND_C_Met_INTEN, SI_IND_C_Phe_INTEN, SI_IND_C_Pro_INTEN, SI_IND_C_Ser_INTEN, SI_IND_C_Thr_INTEN,
	SI_IND_C_Trp_INTEN,    SI_IND_C_Tyr_INTEN,    SI_IND_C_Val_INTEN,

	// The following are only filled if the C are in the visible range
	SI_IND_C_SE, 

	SI_SE_IND_HAS_C_FRAG_INTEN, 
	SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN, SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN, SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,
	SI_SE_IND_C_FRAG_DIFF_01,  SI_SE_IND_C_FRAG_DIFF_05,  SI_SE_IND_C_FRAG_DIFF_LARGE,

	SI_SE_IND_C_N_TERM_DIFF_INTEN, SI_SE_IND_C_C_TERM_DIFF_INTEN, SI_SE_IND_C_Gap_DIFF_INTEN, SI_SE_IND_C_Xle_DIFF_INTEN, SI_SE_IND_C_Ala_DIFF_INTEN, SI_SE_IND_C_Arg_DIFF_INTEN, SI_SE_IND_C_Asn_DIFF_INTEN,
	SI_SE_IND_C_Asp_DIFF_INTEN,    SI_SE_IND_C_Cys_DIFF_INTEN,    SI_SE_IND_C_Gln_DIFF_INTEN, SI_SE_IND_C_Glu_DIFF_INTEN, SI_SE_IND_C_Gly_DIFF_INTEN, SI_SE_IND_C_His_DIFF_INTEN, SI_SE_IND_C_Ile_DIFF_INTEN,
	SI_SE_IND_C_Leu_DIFF_INTEN,    SI_SE_IND_C_Lys_DIFF_INTEN,    SI_SE_IND_C_Met_DIFF_INTEN, SI_SE_IND_C_Phe_DIFF_INTEN, SI_SE_IND_C_Pro_DIFF_INTEN, SI_SE_IND_C_Ser_DIFF_INTEN, SI_SE_IND_C_Thr_DIFF_INTEN,
	SI_SE_IND_C_Trp_DIFF_INTEN,    SI_SE_IND_C_Tyr_DIFF_INTEN,    SI_SE_IND_C_Val_DIFF_INTEN,

	SI_SE_IND_HAS_NO_C_FRAG_INTEN,
	SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN, SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN, SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	SI_IND_C_ME, 
	SI_ME_IND_HAS_C_FRAG_INTEN, 
	SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN, SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN, SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,
	SI_ME_IND_C_FRAG_DIFF_01,  SI_ME_IND_C_FRAG_DIFF_05,  SI_ME_IND_C_FRAG_DIFF_LARGE,

	SI_ME_IND_HAS_NO_C_FRAG_INTEN,
	SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN, SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN, SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,


	SI_NUM_FEATURES

} ScoreModelFields_SI;


// Strong no intensity
typedef enum ScoreModelFields_SNI {
	SNI_CONST,

	SNI_IND_MIRROR1_VIZ,	SNI_IND_MIRROR1_NOT_VIZ, 
	SNI_IND_HAS_MIRROR1_INTEN,  SNI_IND_MIRROR1_NO_INTEN, SNI_MIRROR1_ISO_LEVEL,
	SNI_IND_MIRROR1_HAS_MINUS_1, SNI_MIRROR1_MINUS_1_INTEN_DIFF, // with mirror
	SNI_IND_MIRROR1_HAS_MINUS_2, SNI_MIRROR1_MINUS_2_INTEN_DIFF,
	SNI_IND_MIRROR1_HAS_PLUS_1,	SNI_MIRROR1_PLUS_1_INTEN_DIFF,
	SNI_IND_MIRROR1_HAS_PLUS_2,	SNI_MIRROR1_PLUS_2_INTEN_DIFF,

	SNI_IND_MIRROR2_VIZ,	SNI_IND_MIRROR2_NOT_VIZ, 
	SNI_IND_HAS_MIRROR2_INTEN,  SNI_IND_MIRROR2_NO_INTEN, SNI_MIRROR2_ISO_LEVEL,
	SNI_IND_MIRROR2_HAS_MINUS_1, SNI_MIRROR2_MINUS_1_INTEN_DIFF, // with mirror
	SNI_IND_MIRROR2_HAS_MINUS_2, SNI_MIRROR2_MINUS_2_INTEN_DIFF,
	SNI_IND_MIRROR2_HAS_PLUS_1,	SNI_MIRROR2_PLUS_1_INTEN_DIFF,
	SNI_IND_MIRROR2_HAS_PLUS_2,	SNI_MIRROR2_PLUS_2_INTEN_DIFF,

	SNI_IND_PARENT1_VIZ,	SNI_IND_PARENT1_NOT_VIZ, 
	SNI_IND_PARENT1_INTEN,  SNI_IND_PARENT1_NO_INTEN, SNI_PARENT1_ISO_LEVEL,
	SNI_PARENT1_LOG_INTEN,	SNI_PARENT1_LOG_GLOBAL_RANK,

	SNI_IND_PARENT2_VIZ,	SNI_IND_PARENT2_NOT_VIZ, 
	SNI_IND_PARENT2_INTEN,  SNI_IND_PARENT2_NO_INTEN, SNI_PARENT2_ISO_LEVEL,
	SNI_PARENT2_LOG_INTEN,	SNI_PARENT2_LOG_GLOBAL_RANK,

	SNI_IND_DIS_FROM_MINMAX_LESS_50,	SNI_DIS_FROM_MINMAX0,
	SNI_IND_DIS_FROM_MINMAX_LESS_150,	SNI_DIS_FROM_MINMAX50,
	SNI_IND_DIS_FROM_MINMAX_LESS_250,	SNI_DIS_FROM_MINMAX150,
	SNI_IND_DIS_FROM_MINMAX_MORE,		SNI_DIS_FROM_MINMAX250,

	SNI_REL_POS0, SNI_REL_POS1, SNI_REL_POS2, SNI_REL_POS3, SNI_REL_POS4,
	SNI_REL_POS5, SNI_REL_POS6, SNI_REL_POS7, SNI_REL_POS8, SNI_REL_POS9,

	SNI_IND_N_IS_GAP,	SNI_IND_C_IS_GAP,

	SNI_N_TERM_CAT20, SNI_N_TERM_CAT18,	SNI_N_TERM_CAT12, SNI_N_TERM_CAT8, SNI_N_TERM_CAT4, SNI_N_TERM_CAT2,
	SNI_N_EDGE_CAT20, SNI_N_EDGE_CAT18, SNI_N_EDGE_CAT12, SNI_N_EDGE_CAT8, SNI_N_EDGE_CAT4, SNI_N_EDGE_CAT2,
	SNI_C_TERM_CAT20, SNI_C_TERM_CAT18, SNI_C_TERM_CAT12, SNI_C_TERM_CAT8, SNI_C_TERM_CAT4, SNI_C_TERM_CAT2,
	SNI_C_EDGE_CAT20, SNI_C_EDGE_CAT18, SNI_C_EDGE_CAT12, SNI_C_EDGE_CAT8, SNI_C_EDGE_CAT4, SNI_C_EDGE_CAT2,
	SNI_SPAN_CAT20,	  SNI_SPAN_CAT18,   SNI_SPAN_CAT12,   SNI_SPAN_CAT8,   SNI_SPAN_CAT4,   SNI_SPAN_CAT2,
	SNI_ND_SPAN_CAT20, SNI_ND_SPAN_CAT18, SNI_ND_SPAN_CAT12, SNI_ND_SPAN_CAT8, SNI_ND_SPAN_CAT4, SNI_ND_SPAN_CAT2,
	SNI_CD_SPAN_CAT20, SNI_CD_SPAN_CAT18, SNI_CD_SPAN_CAT12, SNI_CD_SPAN_CAT8, SNI_CD_SPAN_CAT4, SNI_CD_SPAN_CAT2,


	SNI_IND_CONNECTS_TO_N_TERM,			SNI_IND_CONNECTS_TO_C_TERM,
	SNI_IND_PREFERRED_DIGEST_AA_C_TERM, SNI_IND_PREFERRED_DIGEST_AA_N_TERM,
	SNI_IND_NOT_CONNECTED_TO_TERMS,     SNI_IND_MISSED_CLEAVAGE, 


	SNI_IND_N_NOT_VIZ,	SNI_IND_C_NOT_VIZ,
	SNI_IND_N_INTEN,	SNI_IND_N_NO_INTEN,
	SNI_IND_C_INTEN,	SNI_IND_C_NO_INTEN,

	// For N-terminal

	SNI_IND_N_N_TERM, SNI_IND_N_C_TERM, SNI_IND_N_Gap, SNI_IND_N_Xle, SNI_IND_N_Ala, SNI_IND_N_Arg, SNI_IND_N_Asn,
	SNI_IND_N_Asp,    SNI_IND_N_Cys,    SNI_IND_N_Gln, SNI_IND_N_Glu, SNI_IND_N_Gly, SNI_IND_N_His, SNI_IND_N_Ile,
	SNI_IND_N_Leu,    SNI_IND_N_Lys,    SNI_IND_N_Met, SNI_IND_N_Phe, SNI_IND_N_Pro, SNI_IND_N_Ser, SNI_IND_N_Thr,
	SNI_IND_N_Trp,    SNI_IND_N_Tyr,    SNI_IND_N_Val,

	// The following are only filled if the N are in the visible range
	SNI_IND_N_SE, 

	SNI_SE_IND_HAS_N_FRAG_INTEN, 
	SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN, SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN, SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,

	SNI_SE_IND_N_N_TERM_DIFF_INTEN, SNI_SE_IND_N_C_TERM_DIFF_INTEN, SNI_SE_IND_N_Gap_DIFF_INTEN, SNI_SE_IND_N_Xle_DIFF_INTEN, SNI_SE_IND_N_Ala_DIFF_INTEN, SNI_SE_IND_N_Arg_DIFF_INTEN, SNI_SE_IND_N_Asn_DIFF_INTEN,
	SNI_SE_IND_N_Asp_DIFF_INTEN,    SNI_SE_IND_N_Cys_DIFF_INTEN,    SNI_SE_IND_N_Gln_DIFF_INTEN, SNI_SE_IND_N_Glu_DIFF_INTEN, SNI_SE_IND_N_Gly_DIFF_INTEN, SNI_SE_IND_N_His_DIFF_INTEN, SNI_SE_IND_N_Ile_DIFF_INTEN,
	SNI_SE_IND_N_Leu_DIFF_INTEN,    SNI_SE_IND_N_Lys_DIFF_INTEN,    SNI_SE_IND_N_Met_DIFF_INTEN, SNI_SE_IND_N_Phe_DIFF_INTEN, SNI_SE_IND_N_Pro_DIFF_INTEN, SNI_SE_IND_N_Ser_DIFF_INTEN, SNI_SE_IND_N_Thr_DIFF_INTEN,
	SNI_SE_IND_N_Trp_DIFF_INTEN,    SNI_SE_IND_N_Tyr_DIFF_INTEN,    SNI_SE_IND_N_Val_DIFF_INTEN,

	SNI_SE_IND_HAS_NO_N_FRAG_INTEN,
	SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN, SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN, SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	SNI_IND_N_ME, 
	SNI_ME_IND_HAS_N_FRAG_INTEN, 
	SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN, SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN, SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN,

	SNI_ME_IND_HAS_NO_N_FRAG_INTEN,
	SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN, SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN, SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	// For C-terminal

	SNI_IND_C_N_TERM, SNI_IND_C_C_TERM, SNI_IND_C_Gap, SNI_IND_C_Xle, SNI_IND_C_Ala, SNI_IND_C_Arg, SNI_IND_C_Asn,
	SNI_IND_C_Asp,    SNI_IND_C_Cys,    SNI_IND_C_Gln, SNI_IND_C_Glu, SNI_IND_C_Gly, SNI_IND_C_His, SNI_IND_C_Ile,
	SNI_IND_C_Leu,    SNI_IND_C_Lys,    SNI_IND_C_Met, SNI_IND_C_Phe, SNI_IND_C_Pro, SNI_IND_C_Ser, SNI_IND_C_Thr,
	SNI_IND_C_Trp,    SNI_IND_C_Tyr,    SNI_IND_C_Val,

	// The following are only filled if the N are in the visible range
	SNI_IND_C_SE, 

	SNI_SE_IND_HAS_C_FRAG_INTEN, 
	SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN, SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN, SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,

	SNI_SE_IND_C_N_TERM_DIFF_INTEN, SNI_SE_IND_C_C_TERM_DIFF_INTEN, SNI_SE_IND_C_Gap_DIFF_INTEN, SNI_SE_IND_C_Xle_DIFF_INTEN, SNI_SE_IND_C_Ala_DIFF_INTEN, SNI_SE_IND_C_Arg_DIFF_INTEN, SNI_SE_IND_C_Asn_DIFF_INTEN,
	SNI_SE_IND_C_Asp_DIFF_INTEN,    SNI_SE_IND_C_Cys_DIFF_INTEN,    SNI_SE_IND_C_Gln_DIFF_INTEN, SNI_SE_IND_C_Glu_DIFF_INTEN, SNI_SE_IND_C_Gly_DIFF_INTEN, SNI_SE_IND_C_His_DIFF_INTEN, SNI_SE_IND_C_Ile_DIFF_INTEN,
	SNI_SE_IND_C_Leu_DIFF_INTEN,    SNI_SE_IND_C_Lys_DIFF_INTEN,    SNI_SE_IND_C_Met_DIFF_INTEN, SNI_SE_IND_C_Phe_DIFF_INTEN, SNI_SE_IND_C_Pro_DIFF_INTEN, SNI_SE_IND_C_Ser_DIFF_INTEN, SNI_SE_IND_C_Thr_DIFF_INTEN,
	SNI_SE_IND_C_Trp_DIFF_INTEN,    SNI_SE_IND_C_Tyr_DIFF_INTEN,    SNI_SE_IND_C_Val_DIFF_INTEN,

	SNI_SE_IND_HAS_NO_C_FRAG_INTEN,
	SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN, SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN, SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,

	SNI_IND_C_ME, 
	SNI_ME_IND_HAS_C_FRAG_INTEN, 
	SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN, SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN, SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN,

	SNI_ME_IND_HAS_NO_C_FRAG_INTEN,
	SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN, SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN, SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN,


	SNI_NUM_FEATURES

} ScoreModelFields_SNI;


// a regular fragment with intensity
typedef enum ScoreModelFields_RI {
	RI_CONST,                          

	RI_LOG_LOCAL_RANK,		RI_LOG_GLOBAL_RANK,		RI_ISO_LEVEL,
	
	RI_IND_LOG_INTEN_LESS1,  RI_LOG_INTEN_LESS1,
	RI_IND_LOG_INTEN_LESS2,  RI_LOG_INTEN_LESS2,
	RI_IND_LOG_INTEN_LESS3,  RI_LOG_INTEN_LESS3,
	RI_IND_LOG_INTEN_LESS4,  RI_LOG_INTEN_LESS4,
	RI_IND_LOG_INTEN_MORE,   RI_LOG_INTEN_MORE,

	RI_IND_DIS_FROM_MINMAX_LESS_50,		RI_DIS_FROM_MINMAX0,	RI_LOG_INTEN_DIS50,
	RI_IND_DIS_FROM_MINMAX_LESS_150,	RI_DIS_FROM_MINMAX50,	RI_LOG_INTEN_DIS150,
	RI_IND_DIS_FROM_MINMAX_LESS_250,	RI_DIS_FROM_MINMAX150,	RI_LOG_INTEN_DIS250,
	RI_IND_DIS_FROM_MINMAX_MORE,		RI_DIS_FROM_MINMAX250,	RI_LOG_INTEN_DISMORE,

	RI_REL_POS0, RI_REL_POS1, RI_REL_POS2, RI_REL_POS3, RI_REL_POS4,
	RI_REL_POS5, RI_REL_POS6, RI_REL_POS7, RI_REL_POS8, RI_REL_POS9,

	RI_IND_NUM_PARENTS_WITH_INTEN_IS_0,	RI_IND_NUM_PARENTS_WITH_INTEN_IS_1,
	RI_IND_NUM_PARENTS_WITH_INTEN_IS_2,	RI_IND_NUM_PARENTS_WITH_INTEN_IS_3,
	RI_IND_NUM_PARENTS_WITH_INTEN_IS_4,	RI_IND_NUM_PARENTS_WITH_INTEN_IS_5,
	RI_IND_NUM_PARENTS_WITH_INTEN_IS_6,	RI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6,

	// exact combo for top 4 parents (if not viz, assume frag is present)
	RI_IND_PARENT_COMBO_0,	RI_IND_PARENT_COMBO_1, RI_IND_PARENT_COMBO_2, RI_IND_PARENT_COMBO_3,
	RI_IND_PARENT_COMBO_4,	RI_IND_PARENT_COMBO_5, RI_IND_PARENT_COMBO_6, RI_IND_PARENT_COMBO_7,
	RI_IND_PARENT_COMBO_8,	RI_IND_PARENT_COMBO_9, RI_IND_PARENT_COMBO_10, RI_IND_PARENT_COMBO_11,
	RI_IND_PARENT_COMBO_12,	RI_IND_PARENT_COMBO_13, RI_IND_PARENT_COMBO_14, RI_IND_PARENT_COMBO_15,

	RI_IND_GOT_BOTH_ORIS,	   RI_IND_GOT_PREFIX,	   RI_IND_GOT_SUFFIX,

	RI_IND_PARENT1_NOT_VIZ,	   RI_IND_PARENT1_NO_INTEN,	RI_PARENT1_ISO_LEVEL,
	RI_IND_PARENT1_INTEN_MORE, RI_PARENT1_INTEN_DIFF_MORE,
	RI_IND_PARENT1_INTEN_LESS, RI_PARENT1_INTEN_DIFF_LESS,

	RI_IND_PARENT2_NOT_VIZ,	   RI_IND_PARENT2_NO_INTEN,	RI_PARENT2_ISO_LEVEL,
	RI_IND_PARENT2_INTEN_MORE, RI_PARENT2_INTEN_DIFF_MORE,
	RI_IND_PARENT2_INTEN_LESS, RI_PARENT2_INTEN_DIFF_LESS,

	RI_IND_PARENT3_NOT_VIZ,	   RI_IND_PARENT3_NO_INTEN,	RI_PARENT3_ISO_LEVEL,
	RI_IND_PARENT3_INTEN_MORE, RI_PARENT3_INTEN_DIFF_MORE,
	RI_IND_PARENT3_INTEN_LESS, RI_PARENT3_INTEN_DIFF_LESS,

	RI_IND_PARENT4_NOT_VIZ,	   RI_IND_PARENT4_NO_INTEN, RI_PARENT4_ISO_LEVEL,	
	RI_IND_PARENT4_INTEN_MORE, RI_PARENT4_INTEN_DIFF_MORE,
	RI_IND_PARENT4_INTEN_LESS, RI_PARENT4_INTEN_DIFF_LESS,

	RI_IND_PARENT5_NOT_VIZ,	   RI_IND_PARENT5_NO_INTEN,	RI_PARENT5_ISO_LEVEL,
	RI_IND_PARENT5_INTEN_MORE, RI_PARENT5_INTEN_DIFF_MORE,
	RI_IND_PARENT5_INTEN_LESS, RI_PARENT5_INTEN_DIFF_LESS,

	RI_IND_PARENT6_NOT_VIZ,	   RI_IND_PARENT6_NO_INTEN,	RI_PARENT6_ISO_LEVEL,
	RI_IND_PARENT6_INTEN_MORE, RI_PARENT6_INTEN_DIFF_MORE,
	RI_IND_PARENT6_INTEN_LESS, RI_PARENT6_INTEN_DIFF_LESS,

	RI_IND_PARENT7_NOT_VIZ,	   RI_IND_PARENT7_NO_INTEN,	RI_PARENT7_ISO_LEVEL,
	RI_IND_PARENT7_INTEN_MORE, RI_PARENT7_INTEN_DIFF_MORE,
	RI_IND_PARENT7_INTEN_LESS, RI_PARENT7_INTEN_DIFF_LESS,

	RI_IND_PARENT8_NOT_VIZ,	   RI_IND_PARENT8_NO_INTEN,	RI_PARENT8_ISO_LEVEL,
	RI_IND_PARENT8_INTEN_MORE, RI_PARENT8_INTEN_DIFF_MORE,
	RI_IND_PARENT8_INTEN_LESS, RI_PARENT8_INTEN_DIFF_LESS,

	// variable aa features
	RI_IND_N_IS_GAP,	RI_IND_C_IS_GAP,

	RI_IND_N_HAS_N_TERM, RI_IND_N_HAS_C_TERM, RI_IND_N_HAS_Gap, RI_IND_N_HAS_Xle, RI_IND_N_HAS_Ala, RI_IND_N_HAS_Arg, RI_IND_N_HAS_Asn,
	RI_IND_N_HAS_Asp,    RI_IND_N_HAS_Cys,    RI_IND_N_HAS_Gln, RI_IND_N_HAS_Glu, RI_IND_N_HAS_Gly, RI_IND_N_HAS_His, RI_IND_N_HAS_Ile,
	RI_IND_N_HAS_Leu,    RI_IND_N_HAS_Lys,    RI_IND_N_HAS_Met, RI_IND_N_HAS_Phe, RI_IND_N_HAS_Pro, RI_IND_N_HAS_Ser, RI_IND_N_HAS_Thr,
	RI_IND_N_HAS_Trp,    RI_IND_N_HAS_Tyr,    RI_IND_N_HAS_Val,

	RI_N_N_TERM_SELF_INTEN, RI_N_C_TERM_SELF_INTEN, RI_N_Gap_SELF_INTEN, RI_N_Xle_SELF_INTEN, RI_N_Ala_SELF_INTEN, RI_N_Arg_SELF_INTEN, RI_N_Asn_SELF_INTEN,
	RI_N_Asp_SELF_INTEN,    RI_N_Cys_SELF_INTEN,    RI_N_Gln_SELF_INTEN, RI_N_Glu_SELF_INTEN, RI_N_Gly_SELF_INTEN, RI_N_His_SELF_INTEN, RI_N_Ile_SELF_INTEN,
	RI_N_Leu_SELF_INTEN,    RI_N_Lys_SELF_INTEN,    RI_N_Met_SELF_INTEN, RI_N_Phe_SELF_INTEN, RI_N_Pro_SELF_INTEN, RI_N_Ser_SELF_INTEN, RI_N_Thr_SELF_INTEN,
	RI_N_Trp_SELF_INTEN,    RI_N_Tyr_SELF_INTEN,    RI_N_Val_SELF_INTEN,

	RI_IND_C_HAS_N_TERM, RI_IND_C_HAS_C_TERM, RI_IND_C_HAS_Gap, RI_IND_C_HAS_Xle, RI_IND_C_HAS_Ala, RI_IND_C_HAS_Arg, RI_IND_C_HAS_Asn,
	RI_IND_C_HAS_Asp,    RI_IND_C_HAS_Cys,    RI_IND_C_HAS_Gln, RI_IND_C_HAS_Glu, RI_IND_C_HAS_Gly, RI_IND_C_HAS_His, RI_IND_C_HAS_Ile,
	RI_IND_C_HAS_Leu,    RI_IND_C_HAS_Lys,    RI_IND_C_HAS_Met, RI_IND_C_HAS_Phe, RI_IND_C_HAS_Pro, RI_IND_C_HAS_Ser, RI_IND_C_HAS_Thr,
	RI_IND_C_HAS_Trp,    RI_IND_C_HAS_Tyr,    RI_IND_C_HAS_Val,

	RI_C_N_TERM_SELF_INTEN, RI_C_C_TERM_SELF_INTEN, RI_C_Gap_SELF_INTEN, RI_C_Xle_SELF_INTEN, RI_C_Ala_SELF_INTEN, RI_C_Arg_SELF_INTEN, RI_C_Asn_SELF_INTEN,
	RI_C_Asp_SELF_INTEN,    RI_C_Cys_SELF_INTEN,    RI_C_Gln_SELF_INTEN, RI_C_Glu_SELF_INTEN, RI_C_Gly_SELF_INTEN, RI_C_His_SELF_INTEN, RI_C_Ile_SELF_INTEN,
	RI_C_Leu_SELF_INTEN,    RI_C_Lys_SELF_INTEN,    RI_C_Met_SELF_INTEN, RI_C_Phe_SELF_INTEN, RI_C_Pro_SELF_INTEN, RI_C_Ser_SELF_INTEN, RI_C_Thr_SELF_INTEN,
	RI_C_Trp_SELF_INTEN,    RI_C_Tyr_SELF_INTEN,    RI_C_Val_SELF_INTEN,

	RI_NUM_FEATURES

} ScoreModelFields_RI;



typedef enum ScoreModelFields_RNI {
	RNI_CONST,                          

	RNI_IND_DIS_FROM_MINMAX_LESS_50,	RNI_DIS_FROM_MINMAX0,
	RNI_IND_DIS_FROM_MINMAX_LESS_150,	RNI_DIS_FROM_MINMAX50,
	RNI_IND_DIS_FROM_MINMAX_LESS_250,	RNI_DIS_FROM_MINMAX150,
	RNI_IND_DIS_FROM_MINMAX_MORE,		RNI_DIS_FROM_MINMAX250,

	RNI_REL_POS0, RNI_REL_POS1, RNI_REL_POS2, RNI_REL_POS3, RNI_REL_POS4,
	RNI_REL_POS5, RNI_REL_POS6, RNI_REL_POS7, RNI_REL_POS8, RNI_REL_POS9,

	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_0,	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_1,
	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_2,	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_3,
	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_4,	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_5,
	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_6,	RNI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6,

	// exact combo for top 4 parents (if not viz, assume frag is present)
	RNI_IND_PARENT_COMBO_0,	RNI_IND_PARENT_COMBO_1, RNI_IND_PARENT_COMBO_2, RNI_IND_PARENT_COMBO_3,
	RNI_IND_PARENT_COMBO_4,	RNI_IND_PARENT_COMBO_5, RNI_IND_PARENT_COMBO_6, RNI_IND_PARENT_COMBO_7,
	RNI_IND_PARENT_COMBO_8,	RNI_IND_PARENT_COMBO_9, RNI_IND_PARENT_COMBO_10, RNI_IND_PARENT_COMBO_11,
	RNI_IND_PARENT_COMBO_12, RNI_IND_PARENT_COMBO_13, RNI_IND_PARENT_COMBO_14, RNI_IND_PARENT_COMBO_15,

	RNI_IND_GOT_BOTH_ORIS,	   RNI_IND_GOT_PREFIX,	   RNI_IND_GOT_SUFFIX,

	RNI_IND_PARENT1_NOT_VIZ, RNI_IND_PARENT1_NO_INTEN, RNI_PARENT1_LOG_INTEN, RNI_PARENT1_LOG_GLOBAL_RANK, RNI_PARENT1_ISO_LEVEL,	
	RNI_IND_PARENT2_NOT_VIZ, RNI_IND_PARENT2_NO_INTEN, RNI_PARENT2_LOG_INTEN, RNI_PARENT2_LOG_GLOBAL_RANK, RNI_PARENT2_ISO_LEVEL,
	RNI_IND_PARENT3_NOT_VIZ, RNI_IND_PARENT3_NO_INTEN, RNI_PARENT3_LOG_INTEN, RNI_PARENT3_LOG_GLOBAL_RANK, RNI_PARENT3_ISO_LEVEL,
	RNI_IND_PARENT4_NOT_VIZ, RNI_IND_PARENT4_NO_INTEN, RNI_PARENT4_LOG_INTEN, RNI_PARENT4_LOG_GLOBAL_RANK, RNI_PARENT4_ISO_LEVEL,
	RNI_IND_PARENT5_NOT_VIZ, RNI_IND_PARENT5_NO_INTEN, RNI_PARENT5_LOG_INTEN, RNI_PARENT5_LOG_GLOBAL_RANK, RNI_PARENT5_ISO_LEVEL,
	RNI_IND_PARENT6_NOT_VIZ, RNI_IND_PARENT6_NO_INTEN, RNI_PARENT6_LOG_INTEN, RNI_PARENT6_LOG_GLOBAL_RANK, RNI_PARENT6_ISO_LEVEL,
	RNI_IND_PARENT7_NOT_VIZ, RNI_IND_PARENT7_NO_INTEN, RNI_PARENT7_LOG_INTEN, RNI_PARENT7_LOG_GLOBAL_RANK, RNI_PARENT7_ISO_LEVEL,
	RNI_IND_PARENT8_NOT_VIZ, RNI_IND_PARENT8_NO_INTEN, RNI_PARENT8_LOG_INTEN, RNI_PARENT8_LOG_GLOBAL_RANK, RNI_PARENT8_ISO_LEVEL,

	// variable aa features
	RNI_IND_N_IS_GAP,	RNI_IND_C_IS_GAP,

	RNI_IND_N_HAS_N_TERM, RNI_IND_N_HAS_C_TERM, RNI_IND_N_HAS_Gap, RNI_IND_N_HAS_Xle, RNI_IND_N_HAS_Ala, RNI_IND_N_HAS_Arg, RNI_IND_N_HAS_Asn,
	RNI_IND_N_HAS_Asp,    RNI_IND_N_HAS_Cys,    RNI_IND_N_HAS_Gln, RNI_IND_N_HAS_Glu, RNI_IND_N_HAS_Gly, RNI_IND_N_HAS_His, RNI_IND_N_HAS_Ile,
	RNI_IND_N_HAS_Leu,    RNI_IND_N_HAS_Lys,    RNI_IND_N_HAS_Met, RNI_IND_N_HAS_Phe, RNI_IND_N_HAS_Pro, RNI_IND_N_HAS_Ser, RNI_IND_N_HAS_Thr,
	RNI_IND_N_HAS_Trp,    RNI_IND_N_HAS_Tyr,    RNI_IND_N_HAS_Val,

	RNI_IND_C_HAS_N_TERM, RNI_IND_C_HAS_C_TERM, RNI_IND_C_HAS_Gap, RNI_IND_C_HAS_Xle, RNI_IND_C_HAS_Ala, RNI_IND_C_HAS_Arg, RNI_IND_C_HAS_Asn,
	RNI_IND_C_HAS_Asp,    RNI_IND_C_HAS_Cys,    RNI_IND_C_HAS_Gln, RNI_IND_C_HAS_Glu, RNI_IND_C_HAS_Gly, RNI_IND_C_HAS_His, RNI_IND_C_HAS_Ile,
	RNI_IND_C_HAS_Leu,    RNI_IND_C_HAS_Lys,    RNI_IND_C_HAS_Met, RNI_IND_C_HAS_Phe, RNI_IND_C_HAS_Pro, RNI_IND_C_HAS_Ser, RNI_IND_C_HAS_Thr,
	RNI_IND_C_HAS_Trp,    RNI_IND_C_HAS_Tyr,    RNI_IND_C_HAS_Val,


	RNI_NUM_FEATURES

} ScoreModelFields_RNI;





/***************************************************************************
Virtual class.
Holds common elements to derived classes: strong fragment and regular fragment.
****************************************************************************/
class FragModel {
	friend class RegionalScoreModel;
public:
	FragModel() : config(NULL), tolerance(NEG_INF), frag_tolerance(NEG_INF), exact_peak_tolerance(NEG_INF), 
		one_over_tolerance(NEG_INF),  ind_has_models(0), model_frag_idx(NEG_INF), model_frag_charge(NEG_INF), 
		inten_log_scaling_factor(0), no_inten_log_scaling_factor(0) {}

	virtual void fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const =0;

	virtual void fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const =0;

	void set_config_and_tolerance(Config *_c)
	{
		config = _c;
		tolerance = config->get_tolerance();
		frag_tolerance = (tolerance<0.1 ? tolerance : tolerance * 0.6);
		exact_peak_tolerance = (tolerance<0.1 ? 0.5 * tolerance : 0.3 *tolerance);
		one_over_tolerance = 1.0 / tolerance;
	}


	virtual bool read_model(istream& is, bool silent_ind) =0;

	virtual bool write_model (ostream& os) const =0;

	int get_model_frag_idx() const { return model_frag_idx; }

protected:
	Config *config;
	mass_t tolerance;
	mass_t frag_tolerance;
	mass_t exact_peak_tolerance;
	mass_t one_over_tolerance;

	int ind_has_models;

	int model_frag_idx;
	int model_frag_charge;

	score_t inten_log_scaling_factor;
	score_t no_inten_log_scaling_factor;

	ME_Regression_Model inten_model;
	ME_Regression_Model no_inten_model;


};

/************************************************************************
Model for fragments designated "strong" in the regional fragments.
*************************************************************************/
class StrongFragModel : public FragModel {
	friend class RegionalScoreModel;
public:
	StrongFragModel() : mirror1_idx(NEG_INF), mirror2_idx(NEG_INF),
						parent1_idx(NEG_INF), parent2_idx(NEG_INF),
						mirror1_charge(0), mirror2_charge(0), parent1_charge(0), parent2_charge(0) {}



	bool read_model(istream& is, bool silent_ind);

	bool write_model (ostream& os) const;


	void fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const;

	void fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const;

private:
	int mirror1_idx,  mirror2_idx;
	int parent1_idx,  parent2_idx;

	int mirror1_charge, mirror2_charge, parent1_charge, parent2_charge;

	
	void fill_constant_vals(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage, 
							vector<fval>& f_vals) const;

	void fill_aa_variable_vals( Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const;
};


/************************************************************************
Model for fragments not designated "strong" in the regional fragments.
*************************************************************************/
class RegularFragModel : public FragModel {
	friend class RegionalScoreModel;
public:
	RegularFragModel() : num_parents(0), parent_idx_with_same_charge_ori(-1) {}
	
	bool read_model(istream& is, bool silent_ind);
	bool write_model (ostream& os) const;


	void fill_combo_vectors(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							const vector<BreakageInfo>& infos,
							vector< ME_Regression_Sample > & samples) const;

	void fill_single_frag_vector(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage,
							BreakageInfo& info,
							vector< fval >& f_vals) const;


private:
	int num_parents;
	int parent_idx_with_same_charge_ori;

	vector<int> parent_idxs;

	void fill_constant_vals(Spectrum *spec, 
							mass_t pm_with_19,  
							const Breakage *breakage, 
							vector<fval>& f_vals) const;

	void fill_aa_variable_vals( Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage,
							   const BreakageInfo* info,
							   vector<fval>& f_vals) const;
};


/*****************************************************************************
Scoring model for a specific region (charge / size_idx / region_idx).
Contains a FragModel for each fragment used for scoring in that region.
******************************************************************************/
class RegionalScoreModel
{
	friend class AdvancedScoreModel;
public:
	RegionalScoreModel() : charge(NEG_INF), size_idx(NEG_INF), region_idx(NEG_INF),
								   num_strong_frags(0), num_regular_frags(0), config(NULL),
								   was_initialized(false), has_all_breakage_models(false),
								   rand_prob(0), log_random(0), log_one_minus_random(0), 
								   missing_breakage_score(0) {};

	void init(Config* _c, int _charge, int _size_idx, int _region_idx);

	bool read_regional_score_model(const char *name, bool silent_ind);

	bool write_regional_score_model(const char *name) const;

	bool train_regional_score_model(Model *model, const char *name, const FileManager& fm);

	void calc_constant_element(Node& node,
							   Spectrum *spec, 
							   mass_t pm_with_19,  
							   const Breakage *breakage) const;

	score_t score_a_single_breakage_combo(
							   PrmGraph *prm,
							   Node& node, 
							   const Breakage *breakage,
							   BreakageInfo& info,
							   bool verbose=false) const;

	string make_model_file_name(const char *name) const;

	// Scores a breakage using a Dancik-like scoring model (independent frag probs vs. random)
	float basic_score_breakage(Breakage *breakage) const;

	const FragModel *get_frag_model(int frag_idx) const
	{
		int i;
		for (i=0; i<strong_models.size(); i++)
			if (strong_models[i].model_frag_idx == frag_idx)
				return &strong_models[i];

		for (i=0; i<regular_models.size(); i++)
			if (regular_models[i].model_frag_idx == frag_idx)
				return &regular_models[i];
		return NULL;
	}

	score_t get_frag_prob(int frag_type_idx) const
	{
		int i;
		for (i=0; i<frag_type_idxs.size(); i++)
			if (frag_type_idxs[i]==frag_type_idx)
				return frag_probs[i];
		return NEG_INF;
	}

private:
	int charge, size_idx, region_idx;
	int num_strong_frags;
	int num_regular_frags;
	Config *config;

	bool was_initialized;
	bool has_all_breakage_models;

	vector<StrongFragModel>  strong_models;
	vector<RegularFragModel> regular_models;

	// weights on how much to emphasize the FragModel vs. the Dancik prob 
	vector<float> strong_inten_weights,  strong_no_inten_weights;
	vector<float> regular_inten_weights, regular_no_inten_weights;

	// holds the (1-weight) * frag prob to be weight*model prob for the acutal probability
	vector<float> strong_inten_danc_part,  strong_no_inten_danc_part;
	vector<float> regular_inten_danc_part, regular_no_inten_danc_part;

	// basic score probs
	vector<int>     frag_type_idxs;
	vector<score_t> frag_probs;
	score_t		    rand_prob;
	score_t log_random; 
	score_t log_one_minus_random; 

	vector<score_t> frag_inten_scores;
	vector<score_t> frag_no_inten_scores;
	score_t			missing_breakage_score;

	void create_training_set(Model *model,
							 const FragModel& frag_model,
							 const FileManager& fm,
							 ME_Regression_DataSet& inten_ds,
							 ME_Regression_DataSet& no_inten_ds) const;

};




/*****************************************************************************
General scorin model.
Holds models of all regions.
This is the class used to interface with other classes (PrmGraph, etc.)
******************************************************************************/
class AdvancedScoreModel : public Model {
public:

    
	void train_score_model(const char *name, const FileManager& fm, 
						   int charge=-1, int size_idx=-1, int region_idx=-1);

	void score_all_node_combos(PrmGraph *prm) const;

	void initial_combos_score(PrmGraph *prm) const;

	void score_peptide_node_combos(PrmGraph *prm, const Peptide& peptide ) const;

	// performs scoring on demand (if the combo was not previously scored, calculates
	// values, otherwise returns hashed value
	score_t get_node_combo_score(PrmGraph *prm, int node_idx, 
								 int in_edge_idx, int in_var_idx, 
								 int out_edge_idx, int out_var_idx) const;
	// required Model functions
	void init_score_model() { return; }

	void init_model_for_scoring_spectrum(Spectrum *spec) { return; }

	//  simple score (Dancik style)
	void score_breakage(Spectrum *spec, Breakage *breakage, bool verbose=false) const;
	
	score_t get_missing_breakage_score(int charge, int size_idx, int region_idx) const
	{
		return this->regional_breakage_score_models[charge][size_idx][region_idx].missing_breakage_score;
	}


	void score_graph_edges(PrmGraph& prm) const;

	

	int get_max_score_model_charge() const;

	void *get_rank_model_ptr(int type) { return rank_models[type]; }

	void *get_rank_tag_model_ptr(int length) { return rank_tag_models[length]; }

	void read_rank_models(const char *name, bool silent_ind = false);

	void learn_prm_normalizer_values(const FileManager& fm);

	bool read_prm_normalizer_values();

	void write_prm_normalizer_values() const;

	void normalize_prm_scores(PrmGraph &prm) const;

protected:

	vector< vector< vector< RegionalScoreModel > > > regional_breakage_score_models;

	vector< vector< vector< score_t > > > regional_prm_normalizers;

	void *rank_models[3];

	void *rank_tag_models[10];

	void read_score_model(const char *name, bool silent_ind = false);
	
	void write_score_model(const char *name) const;


};



void predict_fragmentation(AdvancedScoreModel* model, const char* input_file, size_t num_peaks=20);


#endif


