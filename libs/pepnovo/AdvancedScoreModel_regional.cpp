#include "AdvancedScoreModel.h"
#include "auxfun.h"
#include "PrmGraph.h"

const char * ScoreModelFields_SI_names[]={
"SI_CONST",	"SI_IND_MIRROR1_VIZ",	"SI_IND_MIRROR1_NOT_VIZ",	"SI_IND_HAS_MIRROR1_INTEN",	"SI_MIRROR1_ISO_LEVEL",	"SI_IND_MIRROR1_HAS_MINUS_1",	"SI_MIRROR1_MINUS_1_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_MINUS_2",	"SI_MIRROR1_MINUS_2_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_PLUS_1",	"SI_MIRROR1_PLUS_1_INTEN_DIFF",	"SI_IND_MIRROR1_HAS_PLUS_2",	"SI_MIRROR1_PLUS_2_INTEN_DIFF",	"SI_MIRROR1_MASS_DIFF25",	
"SI_MIRROR1_MASS_DIFF75",	"SI_MIRROR1_MASS_DIFF_LARGE",	"SI_IND_MIRROR2_VIZ",	"SI_IND_MIRROR2_NOT_VIZ",	"SI_IND_HAS_MIRROR2_INTEN",	"SI_MIRROR2_ISO_LEVEL",	"SI_IND_MIRROR2_HAS_MINUS_1",	"SI_MIRROR2_MINUS_1_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_MINUS_2",	"SI_MIRROR2_MINUS_2_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_PLUS_1",	"SI_MIRROR2_PLUS_1_INTEN_DIFF",	"SI_IND_MIRROR2_HAS_PLUS_2",	"SI_MIRROR2_PLUS_2_INTEN_DIFF",	
"SI_MIRROR2_MASS_DIFF25",	"SI_MIRROR2_MASS_DIFF75",	"SI_MIRROR2_MASS_DIFF_LARGE",	"SI_IND_PARENT1_VIZ",	"SI_IND_PARENT1_NOT_VIZ",	"SI_IND_PARENT1_NO_INTEN",	"SI_PARENT1_ISO_LEVEL",	"SI_IND_PARENT1_LESS_THAN_100_MIN_MAX",	"SI_IND_PARENT1_LESS_THAN_200_MIN_MAX",	"SI_IND_PARENT1_INTEN_MORE",	"SI_PARENT1_INTEN_DIFF_MORE",	"SI_IND_PARENT1_INTEN_LESS",	"SI_PARENT1_INTEN_DIFF_LESS",	
"SI_IND_PARENT2_VIZ",	"SI_IND_PARENT2_NOT_VIZ",	"SI_IND_PARENT2_NO_INTEN",	"SI_PARENT2_ISO_LEVEL",	"SI_IND_PARENT2_LESS_THAN_100_MIN_MAX",	"SI_IND_PARENT2_LESS_THAN_200_MIN_MAX",	"SI_IND_PARENT2_INTEN_MORE",	"SI_PARENT2_INTEN_DIFF_MORE",	"SI_IND_PARENT2_INTEN_LESS",	"SI_PARENT2_INTEN_DIFF_LESS",	"SI_LOG_LOCAL_RANK",	"SI_LOG_GLOBAL_RANK",	"SI_ISO_LEVEL",	"SI_IND_HAS_MINUS_1",	
"SI_IND_HAS_MINUS_2",	"SI_IND_HAS_PLUS_1",	"SI_IND_HAS_PLUS_2",	"SI_IND_LOG_INTEN_LESS1",	"SI_LOG_INTEN_LESS1",	"SI_IND_LOG_INTEN_LESS2",	"SI_LOG_INTEN_LESS2",	"SI_IND_LOG_INTEN_LESS3",	"SI_LOG_INTEN_LESS3",	"SI_IND_LOG_INTEN_LESS4",	"SI_LOG_INTEN_LESS4",	"SI_IND_LOG_INTEN_MORE",	"SI_LOG_INTEN_MORE",	"SI_IND_DIS_FROM_MINMAX_LESS_50",	"SI_DIS_FROM_MINMAX0",	"SI_LOG_INTEN_DIS_50",	
"SI_IND_DIS_FROM_MINMAX_LESS_150",	"SI_DIS_FROM_MINMAX50",	"SI_LOG_INTEN_DIS_150",	"SI_IND_DIS_FROM_MINMAX_LESS_250",	"SI_DIS_FROM_MINMAX150",	"SI_LOG_INTEN_DIS_250",	"SI_IND_DIS_FROM_MINMAX_MORE",	"SI_DIS_FROM_MINMAX250",	"SI_LOG_INTEN_DIS_MORE",	"SI_REL_POS0",	"SI_REL_POS1",	"SI_REL_POS2",	"SI_REL_POS3",	"SI_REL_POS4",	"SI_REL_POS5",	"SI_REL_POS6",	"SI_REL_POS7",	"SI_REL_POS8",	
"SI_REL_POS9",	"SI_IND_HAS_PLUS_NH3",	"SI_PLUS_NH3_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2O",	"SI_PLUS_H2O_INTEN_DIFF",	"SI_IND_HAS_PLUS_CO",	"SI_PLUS_CO_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2ONH3",	"SI_PLUS_H2ONH3_INTEN_DIFF",	"SI_IND_HAS_PLUS_H2OH2O",	"SI_PLUS_H2OH2O_INTEN_DIFF",	"SI_IND_HAS_CHARGE_PLUS1",	"SI_CHARGE_PLUS1_INTEN_DIFF",	"SI_IND_HAS_CHARGE_MINUS1",	"SI_CHARGE_MINUS1_INTEN_DIFF",	
"SI_IND_N_IS_GAP",	"SI_IND_C_IS_GAP",	"SI_N_TERM_CAT20",	"SI_N_TERM_CAT18",	"SI_N_TERM_CAT12",	"SI_N_TERM_CAT8",	"SI_N_TERM_CAT4",	"SI_N_TERM_CAT2",	"SI_N_EDGE_CAT20",	"SI_N_EDGE_CAT18",	"SI_N_EDGE_CAT12",	"SI_N_EDGE_CAT8",	"SI_N_EDGE_CAT4",	"SI_N_EDGE_CAT2",	"SI_C_TERM_CAT20",	"SI_C_TERM_CAT18",	"SI_C_TERM_CAT12",	"SI_C_TERM_CAT8",	"SI_C_TERM_CAT4",	"SI_C_TERM_CAT2",	"SI_C_EDGE_CAT20",	
"SI_C_EDGE_CAT18",	"SI_C_EDGE_CAT12",	"SI_C_EDGE_CAT8",	"SI_C_EDGE_CAT4",	"SI_C_EDGE_CAT2",	"SI_SPAN_CAT20",	"SI_SPAN_CAT18",	"SI_SPAN_CAT12",	"SI_SPAN_CAT8",	"SI_SPAN_CAT4",	"SI_SPAN_CAT2",	"SI_ND_SPAN_CAT20",	"SI_ND_SPAN_CAT18",	"SI_ND_SPAN_CAT12",	"SI_ND_SPAN_CAT8",	"SI_ND_SPAN_CAT4",	"SI_ND_SPAN_CAT2",	"SI_CD_SPAN_CAT20",	"SI_CD_SPAN_CAT18",	"SI_CD_SPAN_CAT12",	"SI_CD_SPAN_CAT8",	
"SI_CD_SPAN_CAT4",	"SI_CD_SPAN_CAT2",	"SI_IND_CONNECTS_TO_N_TERM",	"SI_IND_PREFERRED_DIGEST_AA_N_TERM",	"SI_PREFERRED_DIGEST_AA_N_TERM_INTEN",	"SI_NON_PREFERRED_DIGEST_AA_N_TERM_INTEN",	"SI_IND_CONNECTS_TO_C_TERM",	"SI_IND_PREFERRED_DIGEST_AA_C_TERM",	"SI_PREFERRED_DIGEST_AA_C_TERM_INTEN",	"SI_NON_PREFERRED_DIGEST_AA_C_TERM_INTEN",	"SI_IND_NOT_CONNECTED_TO_TERMS",	"SI_IND_MISSED_CLEAVAGE",	
"SI_IND_N_FRAG_NOT_VIZ",	"SI_IND_C_FRAG_NOT_VIZ",	"SI_IND_N_INTEN",	"SI_IND_C_INTEN",	"SI_IND_N_NO_INTEN",	"SI_IND_C_NO_INTEN",	"SI_IND_N_N_TERM",	"SI_IND_N_C_TERM",	"SI_IND_N_Gap",	"SI_IND_N_Xle",	"SI_IND_N_Ala",	"SI_IND_N_Arg",	"SI_IND_N_Asn",	"SI_IND_N_Asp",	"SI_IND_N_Cys",	"SI_IND_N_Gln",	"SI_IND_N_Glu",	"SI_IND_N_Gly",	"SI_IND_N_His",	"SI_IND_N_Ile",	"SI_IND_N_Leu",	"SI_IND_N_Lys",	
"SI_IND_N_Met",	"SI_IND_N_Phe",	"SI_IND_N_Pro",	"SI_IND_N_Ser",	"SI_IND_N_Thr",	"SI_IND_N_Trp",	"SI_IND_N_Tyr",	"SI_IND_N_Val",	"SI_IND_N_N_TERM_INTEN",	"SI_IND_N_C_TERM_INTEN",	"SI_IND_N_Gap_INTEN",	"SI_IND_N_Xle_INTEN",	"SI_IND_N_Ala_INTEN",	"SI_IND_N_Arg_INTEN",	"SI_IND_N_Asn_INTEN",	"SI_IND_N_Asp_INTEN",	"SI_IND_N_Cys_INTEN",	"SI_IND_N_Gln_INTEN",	"SI_IND_N_Glu_INTEN",	"SI_IND_N_Gly_INTEN",	
"SI_IND_N_His_INTEN",	"SI_IND_N_Ile_INTEN",	"SI_IND_N_Leu_INTEN",	"SI_IND_N_Lys_INTEN",	"SI_IND_N_Met_INTEN",	"SI_IND_N_Phe_INTEN",	"SI_IND_N_Pro_INTEN",	"SI_IND_N_Ser_INTEN",	"SI_IND_N_Thr_INTEN",	"SI_IND_N_Trp_INTEN",	"SI_IND_N_Tyr_INTEN",	"SI_IND_N_Val_INTEN",	"SI_IND_N_SE",	"SI_SE_IND_HAS_N_FRAG_INTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	
"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_SE_IND_N_FRAG_DIFF_01",	"SI_SE_IND_N_FRAG_DIFF_05",	"SI_SE_IND_N_FRAG_DIFF_LARGE",	"SI_SE_IND_N_N_TERM_DIFF_INTEN",	"SI_SE_IND_N_C_TERM_DIFF_INTEN",	"SI_SE_IND_N_Gap_DIFF_INTEN",	"SI_SE_IND_N_Xle_DIFF_INTEN",	"SI_SE_IND_N_Ala_DIFF_INTEN",	"SI_SE_IND_N_Arg_DIFF_INTEN",	"SI_SE_IND_N_Asn_DIFF_INTEN",	"SI_SE_IND_N_Asp_DIFF_INTEN",	
"SI_SE_IND_N_Cys_DIFF_INTEN",	"SI_SE_IND_N_Gln_DIFF_INTEN",	"SI_SE_IND_N_Glu_DIFF_INTEN",	"SI_SE_IND_N_Gly_DIFF_INTEN",	"SI_SE_IND_N_His_DIFF_INTEN",	"SI_SE_IND_N_Ile_DIFF_INTEN",	"SI_SE_IND_N_Leu_DIFF_INTEN",	"SI_SE_IND_N_Lys_DIFF_INTEN",	"SI_SE_IND_N_Met_DIFF_INTEN",	"SI_SE_IND_N_Phe_DIFF_INTEN",	"SI_SE_IND_N_Pro_DIFF_INTEN",	"SI_SE_IND_N_Ser_DIFF_INTEN",	"SI_SE_IND_N_Thr_DIFF_INTEN",	
"SI_SE_IND_N_Trp_DIFF_INTEN",	"SI_SE_IND_N_Tyr_DIFF_INTEN",	"SI_SE_IND_N_Val_DIFF_INTEN",	"SI_SE_IND_HAS_NO_N_FRAG_INTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_N_ME",	"SI_ME_IND_HAS_N_FRAG_INTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	
"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_ME_IND_N_FRAG_DIFF_01",	"SI_ME_IND_N_FRAG_DIFF_05",	"SI_ME_IND_N_FRAG_DIFF_LARGE",	"SI_ME_IND_HAS_NO_N_FRAG_INTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_C_N_TERM",	"SI_IND_C_C_TERM",	"SI_IND_C_Gap",	"SI_IND_C_Xle",	
"SI_IND_C_Ala",	"SI_IND_C_Arg",	"SI_IND_C_Asn",	"SI_IND_C_Asp",	"SI_IND_C_Cys",	"SI_IND_C_Gln",	"SI_IND_C_Glu",	"SI_IND_C_Gly",	"SI_IND_C_His",	"SI_IND_C_Ile",	"SI_IND_C_Leu",	"SI_IND_C_Lys",	"SI_IND_C_Met",	"SI_IND_C_Phe",	"SI_IND_C_Pro",	"SI_IND_C_Ser",	"SI_IND_C_Thr",	"SI_IND_C_Trp",	"SI_IND_C_Tyr",	"SI_IND_C_Val",	"SI_IND_C_N_TERM_INTEN",	"SI_IND_C_C_TERM_INTEN",	"SI_IND_C_Gap_INTEN",	
"SI_IND_C_Xle_INTEN",	"SI_IND_C_Ala_INTEN",	"SI_IND_C_Arg_INTEN",	"SI_IND_C_Asn_INTEN",	"SI_IND_C_Asp_INTEN",	"SI_IND_C_Cys_INTEN",	"SI_IND_C_Gln_INTEN",	"SI_IND_C_Glu_INTEN",	"SI_IND_C_Gly_INTEN",	"SI_IND_C_His_INTEN",	"SI_IND_C_Ile_INTEN",	"SI_IND_C_Leu_INTEN",	"SI_IND_C_Lys_INTEN",	"SI_IND_C_Met_INTEN",	"SI_IND_C_Phe_INTEN",	"SI_IND_C_Pro_INTEN",	"SI_IND_C_Ser_INTEN",	"SI_IND_C_Thr_INTEN",	
"SI_IND_C_Trp_INTEN",	"SI_IND_C_Tyr_INTEN",	"SI_IND_C_Val_INTEN",	"SI_IND_C_SE",	"SI_SE_IND_HAS_C_FRAG_INTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_SE_IND_C_FRAG_DIFF_01",	"SI_SE_IND_C_FRAG_DIFF_05",	"SI_SE_IND_C_FRAG_DIFF_LARGE",	"SI_SE_IND_C_N_TERM_DIFF_INTEN",	"SI_SE_IND_C_C_TERM_DIFF_INTEN",	
"SI_SE_IND_C_Gap_DIFF_INTEN",	"SI_SE_IND_C_Xle_DIFF_INTEN",	"SI_SE_IND_C_Ala_DIFF_INTEN",	"SI_SE_IND_C_Arg_DIFF_INTEN",	"SI_SE_IND_C_Asn_DIFF_INTEN",	"SI_SE_IND_C_Asp_DIFF_INTEN",	"SI_SE_IND_C_Cys_DIFF_INTEN",	"SI_SE_IND_C_Gln_DIFF_INTEN",	"SI_SE_IND_C_Glu_DIFF_INTEN",	"SI_SE_IND_C_Gly_DIFF_INTEN",	"SI_SE_IND_C_His_DIFF_INTEN",	"SI_SE_IND_C_Ile_DIFF_INTEN",	"SI_SE_IND_C_Leu_DIFF_INTEN",	
"SI_SE_IND_C_Lys_DIFF_INTEN",	"SI_SE_IND_C_Met_DIFF_INTEN",	"SI_SE_IND_C_Phe_DIFF_INTEN",	"SI_SE_IND_C_Pro_DIFF_INTEN",	"SI_SE_IND_C_Ser_DIFF_INTEN",	"SI_SE_IND_C_Thr_DIFF_INTEN",	"SI_SE_IND_C_Trp_DIFF_INTEN",	"SI_SE_IND_C_Tyr_DIFF_INTEN",	"SI_SE_IND_C_Val_DIFF_INTEN",	"SI_SE_IND_HAS_NO_C_FRAG_INTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	
"SI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_IND_C_ME",	"SI_ME_IND_HAS_C_FRAG_INTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SI_ME_IND_C_FRAG_DIFF_01",	"SI_ME_IND_C_FRAG_DIFF_05",	"SI_ME_IND_C_FRAG_DIFF_LARGE",	"SI_ME_IND_HAS_NO_C_FRAG_INTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	
"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SI_NUM_FEATURES",	"ScoreModelFields_SI"};

const char * ScoreModelFields_SNI_names[]={
"SNI_CONST",	"SNI_IND_MIRROR1_VIZ",	"SNI_IND_MIRROR1_NOT_VIZ",	"SNI_IND_HAS_MIRROR1_INTEN",	"SNI_IND_MIRROR1_NO_INTEN",	"SNI_MIRROR1_ISO_LEVEL",	"SNI_IND_MIRROR1_HAS_MINUS_1",	"SNI_MIRROR1_MINUS_1_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_MINUS_2",	"SNI_MIRROR1_MINUS_2_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_PLUS_1",	"SNI_MIRROR1_PLUS_1_INTEN_DIFF",	"SNI_IND_MIRROR1_HAS_PLUS_2",	"SNI_MIRROR1_PLUS_2_INTEN_DIFF",	
"SNI_IND_MIRROR2_VIZ",	"SNI_IND_MIRROR2_NOT_VIZ",	"SNI_IND_HAS_MIRROR2_INTEN",	"SNI_IND_MIRROR2_NO_INTEN",	"SNI_MIRROR2_ISO_LEVEL",	"SNI_IND_MIRROR2_HAS_MINUS_1",	"SNI_MIRROR2_MINUS_1_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_MINUS_2",	"SNI_MIRROR2_MINUS_2_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_PLUS_1",	"SNI_MIRROR2_PLUS_1_INTEN_DIFF",	"SNI_IND_MIRROR2_HAS_PLUS_2",	"SNI_MIRROR2_PLUS_2_INTEN_DIFF",	
"SNI_IND_PARENT1_VIZ",	"SNI_IND_PARENT1_NOT_VIZ",	"SNI_IND_PARENT1_INTEN",	"SNI_IND_PARENT1_NO_INTEN",	"SNI_PARENT1_ISO_LEVEL",	"SNI_PARENT1_LOG_INTEN",	"SNI_PARENT1_LOG_GLOBAL_RANK",	"SNI_IND_PARENT2_VIZ",	"SNI_IND_PARENT2_NOT_VIZ",	"SNI_IND_PARENT2_INTEN",	"SNI_IND_PARENT2_NO_INTEN",	"SNI_PARENT2_ISO_LEVEL",	"SNI_PARENT2_LOG_INTEN",	"SNI_PARENT2_LOG_GLOBAL_RANK",	"SNI_IND_DIS_FROM_MINMAX_LESS_50",	
"SNI_DIS_FROM_MINMAX0",	"SNI_IND_DIS_FROM_MINMAX_LESS_150",	"SNI_DIS_FROM_MINMAX50",	"SNI_IND_DIS_FROM_MINMAX_LESS_250",	"SNI_DIS_FROM_MINMAX150",	"SNI_IND_DIS_FROM_MINMAX_MORE",	"SNI_DIS_FROM_MINMAX250",	"SNI_REL_POS0",	"SNI_REL_POS1",	"SNI_REL_POS2",	"SNI_REL_POS3",	"SNI_REL_POS4",	"SNI_REL_POS5",	"SNI_REL_POS6",	"SNI_REL_POS7",	"SNI_REL_POS8",	"SNI_REL_POS9",	"SNI_IND_N_IS_GAP",	
"SNI_IND_C_IS_GAP",	"SNI_N_TERM_CAT20",	"SNI_N_TERM_CAT18",	"SNI_N_TERM_CAT12",	"SNI_N_TERM_CAT8",	"SNI_N_TERM_CAT4",	"SNI_N_TERM_CAT2",	"SNI_N_EDGE_CAT20",	"SNI_N_EDGE_CAT18",	"SNI_N_EDGE_CAT12",	"SNI_N_EDGE_CAT8",	"SNI_N_EDGE_CAT4",	"SNI_N_EDGE_CAT2",	"SNI_C_TERM_CAT20",	"SNI_C_TERM_CAT18",	"SNI_C_TERM_CAT12",	"SNI_C_TERM_CAT8",	"SNI_C_TERM_CAT4",	"SNI_C_TERM_CAT2",	"SNI_C_EDGE_CAT20",	
"SNI_C_EDGE_CAT18",	"SNI_C_EDGE_CAT12",	"SNI_C_EDGE_CAT8",	"SNI_C_EDGE_CAT4",	"SNI_C_EDGE_CAT2",	"SNI_SPAN_CAT20",	"SNI_SPAN_CAT18",	"SNI_SPAN_CAT12",	"SNI_SPAN_CAT8",	"SNI_SPAN_CAT4",	"SNI_SPAN_CAT2",	"SNI_ND_SPAN_CAT20",	"SNI_ND_SPAN_CAT18",	"SNI_ND_SPAN_CAT12",	"SNI_ND_SPAN_CAT8",	"SNI_ND_SPAN_CAT4",	"SNI_ND_SPAN_CAT2",	"SNI_CD_SPAN_CAT20",	"SNI_CD_SPAN_CAT18",	"SNI_CD_SPAN_CAT12",	
"SNI_CD_SPAN_CAT8",	"SNI_CD_SPAN_CAT4",	"SNI_CD_SPAN_CAT2",	"SNI_IND_CONNECTS_TO_N_TERM",	"SNI_IND_CONNECTS_TO_C_TERM",	"SNI_IND_PREFERRED_DIGEST_AA_C_TERM",	"SNI_IND_PREFERRED_DIGEST_AA_N_TERM",	"SNI_IND_NOT_CONNECTED_TO_TERMS",	"SNI_IND_MISSED_CLEAVAGE",	"SNI_IND_N_NOT_VIZ",	"SNI_IND_C_NOT_VIZ",	"SNI_IND_N_INTEN",	"SNI_IND_N_NO_INTEN",	"SNI_IND_C_INTEN",	"SNI_IND_C_NO_INTEN",	
"SNI_IND_N_N_TERM",	"SNI_IND_N_C_TERM",	"SNI_IND_N_Gap",	"SNI_IND_N_Xle",	"SNI_IND_N_Ala",	"SNI_IND_N_Arg",	"SNI_IND_N_Asn",	"SNI_IND_N_Asp",	"SNI_IND_N_Cys",	"SNI_IND_N_Gln",	"SNI_IND_N_Glu",	"SNI_IND_N_Gly",	"SNI_IND_N_His",	"SNI_IND_N_Ile",	"SNI_IND_N_Leu",	"SNI_IND_N_Lys",	"SNI_IND_N_Met",	"SNI_IND_N_Phe",	"SNI_IND_N_Pro",	"SNI_IND_N_Ser",	"SNI_IND_N_Thr",	"SNI_IND_N_Trp",	
"SNI_IND_N_Tyr",	"SNI_IND_N_Val",	"SNI_IND_N_SE",	"SNI_SE_IND_HAS_N_FRAG_INTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_SE_IND_N_N_TERM_DIFF_INTEN",	"SNI_SE_IND_N_C_TERM_DIFF_INTEN",	"SNI_SE_IND_N_Gap_DIFF_INTEN",	"SNI_SE_IND_N_Xle_DIFF_INTEN",	"SNI_SE_IND_N_Ala_DIFF_INTEN",	
"SNI_SE_IND_N_Arg_DIFF_INTEN",	"SNI_SE_IND_N_Asn_DIFF_INTEN",	"SNI_SE_IND_N_Asp_DIFF_INTEN",	"SNI_SE_IND_N_Cys_DIFF_INTEN",	"SNI_SE_IND_N_Gln_DIFF_INTEN",	"SNI_SE_IND_N_Glu_DIFF_INTEN",	"SNI_SE_IND_N_Gly_DIFF_INTEN",	"SNI_SE_IND_N_His_DIFF_INTEN",	"SNI_SE_IND_N_Ile_DIFF_INTEN",	"SNI_SE_IND_N_Leu_DIFF_INTEN",	"SNI_SE_IND_N_Lys_DIFF_INTEN",	"SNI_SE_IND_N_Met_DIFF_INTEN",	"SNI_SE_IND_N_Phe_DIFF_INTEN",	
"SNI_SE_IND_N_Pro_DIFF_INTEN",	"SNI_SE_IND_N_Ser_DIFF_INTEN",	"SNI_SE_IND_N_Thr_DIFF_INTEN",	"SNI_SE_IND_N_Trp_DIFF_INTEN",	"SNI_SE_IND_N_Tyr_DIFF_INTEN",	"SNI_SE_IND_N_Val_DIFF_INTEN",	"SNI_SE_IND_HAS_NO_N_FRAG_INTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_SE_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_N_ME",	
"SNI_ME_IND_HAS_N_FRAG_INTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_ME_IND_HAS_NO_N_FRAG_INTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_ME_IND_N_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_C_N_TERM",	"SNI_IND_C_C_TERM",	
"SNI_IND_C_Gap",	"SNI_IND_C_Xle",	"SNI_IND_C_Ala",	"SNI_IND_C_Arg",	"SNI_IND_C_Asn",	"SNI_IND_C_Asp",	"SNI_IND_C_Cys",	"SNI_IND_C_Gln",	"SNI_IND_C_Glu",	"SNI_IND_C_Gly",	"SNI_IND_C_His",	"SNI_IND_C_Ile",	"SNI_IND_C_Leu",	"SNI_IND_C_Lys",	"SNI_IND_C_Met",	"SNI_IND_C_Phe",	"SNI_IND_C_Pro",	"SNI_IND_C_Ser",	"SNI_IND_C_Thr",	"SNI_IND_C_Trp",	"SNI_IND_C_Tyr",	"SNI_IND_C_Val",	"SNI_IND_C_SE",	
"SNI_SE_IND_HAS_C_FRAG_INTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_SE_IND_C_N_TERM_DIFF_INTEN",	"SNI_SE_IND_C_C_TERM_DIFF_INTEN",	"SNI_SE_IND_C_Gap_DIFF_INTEN",	"SNI_SE_IND_C_Xle_DIFF_INTEN",	"SNI_SE_IND_C_Ala_DIFF_INTEN",	"SNI_SE_IND_C_Arg_DIFF_INTEN",	"SNI_SE_IND_C_Asn_DIFF_INTEN",	
"SNI_SE_IND_C_Asp_DIFF_INTEN",	"SNI_SE_IND_C_Cys_DIFF_INTEN",	"SNI_SE_IND_C_Gln_DIFF_INTEN",	"SNI_SE_IND_C_Glu_DIFF_INTEN",	"SNI_SE_IND_C_Gly_DIFF_INTEN",	"SNI_SE_IND_C_His_DIFF_INTEN",	"SNI_SE_IND_C_Ile_DIFF_INTEN",	"SNI_SE_IND_C_Leu_DIFF_INTEN",	"SNI_SE_IND_C_Lys_DIFF_INTEN",	"SNI_SE_IND_C_Met_DIFF_INTEN",	"SNI_SE_IND_C_Phe_DIFF_INTEN",	"SNI_SE_IND_C_Pro_DIFF_INTEN",	"SNI_SE_IND_C_Ser_DIFF_INTEN",	
"SNI_SE_IND_C_Thr_DIFF_INTEN",	"SNI_SE_IND_C_Trp_DIFF_INTEN",	"SNI_SE_IND_C_Tyr_DIFF_INTEN",	"SNI_SE_IND_C_Val_DIFF_INTEN",	"SNI_SE_IND_HAS_NO_C_FRAG_INTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_SE_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_IND_C_ME",	"SNI_ME_IND_HAS_C_FRAG_INTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_WINTEN",	
"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_WINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_WINTEN",	"SNI_ME_IND_HAS_NO_C_FRAG_INTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_50_NOINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_150_NOINTEN",	"SNI_ME_IND_C_DIS_FROM_MINMAX_LESS_250_NOINTEN",	"SNI_NUM_FEATURES",	"ScoreModelFields_SNI"};

const char * ScoreModelFields_RI_names[]={
"RI_CONST",	"RI_LOG_LOCAL_RANK",	"RI_LOG_GLOBAL_RANK",	"RI_ISO_LEVEL",	"RI_IND_LOG_INTEN_LESS1",	"RI_LOG_INTEN_LESS1",	"RI_IND_LOG_INTEN_LESS2",	"RI_LOG_INTEN_LESS2",	"RI_IND_LOG_INTEN_LESS3",	"RI_LOG_INTEN_LESS3",	"RI_IND_LOG_INTEN_LESS4",	"RI_LOG_INTEN_LESS4",	"RI_IND_LOG_INTEN_MORE",	"RI_LOG_INTEN_MORE",	"RI_IND_DIS_FROM_MINMAX_LESS_50",	"RI_DIS_FROM_MINMAX0",	"RI_LOG_INTEN_DIS50",	
"RI_IND_DIS_FROM_MINMAX_LESS_150",	"RI_DIS_FROM_MINMAX50",	"RI_LOG_INTEN_DIS150",	"RI_IND_DIS_FROM_MINMAX_LESS_250",	"RI_DIS_FROM_MINMAX150",	"RI_LOG_INTEN_DIS250",	"RI_IND_DIS_FROM_MINMAX_MORE",	"RI_DIS_FROM_MINMAX250",	"RI_LOG_INTEN_DISMORE",	"RI_REL_POS0",	"RI_REL_POS1",	"RI_REL_POS2",	"RI_REL_POS3",	"RI_REL_POS4",	"RI_REL_POS5",	"RI_REL_POS6",	"RI_REL_POS7",	"RI_REL_POS8",	
"RI_REL_POS9",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_0",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_1",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_2",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_3",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_4",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_5",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_6",	"RI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6",	"RI_IND_PARENT_COMBO_0",	"RI_IND_PARENT_COMBO_1",	"RI_IND_PARENT_COMBO_2",	
"RI_IND_PARENT_COMBO_3",	"RI_IND_PARENT_COMBO_4",	"RI_IND_PARENT_COMBO_5",	"RI_IND_PARENT_COMBO_6",	"RI_IND_PARENT_COMBO_7",	"RI_IND_PARENT_COMBO_8",	"RI_IND_PARENT_COMBO_9",	"RI_IND_PARENT_COMBO_10",	"RI_IND_PARENT_COMBO_11",	"RI_IND_PARENT_COMBO_12",	"RI_IND_PARENT_COMBO_13",	"RI_IND_PARENT_COMBO_14",	"RI_IND_PARENT_COMBO_15",	"RI_IND_GOT_BOTH_ORIS",	"RI_IND_GOT_PREFIX",	"RI_IND_GOT_SUFFIX",	
"RI_IND_PARENT1_NOT_VIZ",	"RI_IND_PARENT1_NO_INTEN",	"RI_PARENT1_ISO_LEVEL",	"RI_IND_PARENT1_INTEN_MORE",	"RI_PARENT1_INTEN_DIFF_MORE",	"RI_IND_PARENT1_INTEN_LESS",	"RI_PARENT1_INTEN_DIFF_LESS",	"RI_IND_PARENT2_NOT_VIZ",	"RI_IND_PARENT2_NO_INTEN",	"RI_PARENT2_ISO_LEVEL",	"RI_IND_PARENT2_INTEN_MORE",	"RI_PARENT2_INTEN_DIFF_MORE",	"RI_IND_PARENT2_INTEN_LESS",	"RI_PARENT2_INTEN_DIFF_LESS",	
"RI_IND_PARENT3_NOT_VIZ",	"RI_IND_PARENT3_NO_INTEN",	"RI_PARENT3_ISO_LEVEL",	"RI_IND_PARENT3_INTEN_MORE",	"RI_PARENT3_INTEN_DIFF_MORE",	"RI_IND_PARENT3_INTEN_LESS",	"RI_PARENT3_INTEN_DIFF_LESS",	"RI_IND_PARENT4_NOT_VIZ",	"RI_IND_PARENT4_NO_INTEN",	"RI_PARENT4_ISO_LEVEL",	"RI_IND_PARENT4_INTEN_MORE",	"RI_PARENT4_INTEN_DIFF_MORE",	"RI_IND_PARENT4_INTEN_LESS",	"RI_PARENT4_INTEN_DIFF_LESS",	
"RI_IND_PARENT5_NOT_VIZ",	"RI_IND_PARENT5_NO_INTEN",	"RI_PARENT5_ISO_LEVEL",	"RI_IND_PARENT5_INTEN_MORE",	"RI_PARENT5_INTEN_DIFF_MORE",	"RI_IND_PARENT5_INTEN_LESS",	"RI_PARENT5_INTEN_DIFF_LESS",	"RI_IND_PARENT6_NOT_VIZ",	"RI_IND_PARENT6_NO_INTEN",	"RI_PARENT6_ISO_LEVEL",	"RI_IND_PARENT6_INTEN_MORE",	"RI_PARENT6_INTEN_DIFF_MORE",	"RI_IND_PARENT6_INTEN_LESS",	"RI_PARENT6_INTEN_DIFF_LESS",	
"RI_IND_PARENT7_NOT_VIZ",	"RI_IND_PARENT7_NO_INTEN",	"RI_PARENT7_ISO_LEVEL",	"RI_IND_PARENT7_INTEN_MORE",	"RI_PARENT7_INTEN_DIFF_MORE",	"RI_IND_PARENT7_INTEN_LESS",	"RI_PARENT7_INTEN_DIFF_LESS",	"RI_IND_PARENT8_NOT_VIZ",	"RI_IND_PARENT8_NO_INTEN",	"RI_PARENT8_ISO_LEVEL",	"RI_IND_PARENT8_INTEN_MORE",	"RI_PARENT8_INTEN_DIFF_MORE",	"RI_IND_PARENT8_INTEN_LESS",	"RI_PARENT8_INTEN_DIFF_LESS",	
"RI_IND_N_IS_GAP",	"RI_IND_C_IS_GAP",	"RI_IND_N_HAS_N_TERM",	"RI_IND_N_HAS_C_TERM",	"RI_IND_N_HAS_Gap",	"RI_IND_N_HAS_Xle",	"RI_IND_N_HAS_Ala",	"RI_IND_N_HAS_Arg",	"RI_IND_N_HAS_Asn",	"RI_IND_N_HAS_Asp",	"RI_IND_N_HAS_Cys",	"RI_IND_N_HAS_Gln",	"RI_IND_N_HAS_Glu",	"RI_IND_N_HAS_Gly",	"RI_IND_N_HAS_His",	"RI_IND_N_HAS_Ile",	"RI_IND_N_HAS_Leu",	"RI_IND_N_HAS_Lys",	"RI_IND_N_HAS_Met",	
"RI_IND_N_HAS_Phe",	"RI_IND_N_HAS_Pro",	"RI_IND_N_HAS_Ser",	"RI_IND_N_HAS_Thr",	"RI_IND_N_HAS_Trp",	"RI_IND_N_HAS_Tyr",	"RI_IND_N_HAS_Val",	"RI_N_N_TERM_SELF_INTEN",	"RI_N_C_TERM_SELF_INTEN",	"RI_N_Gap_SELF_INTEN",	"RI_N_Xle_SELF_INTEN",	"RI_N_Ala_SELF_INTEN",	"RI_N_Arg_SELF_INTEN",	"RI_N_Asn_SELF_INTEN",	"RI_N_Asp_SELF_INTEN",	"RI_N_Cys_SELF_INTEN",	"RI_N_Gln_SELF_INTEN",	"RI_N_Glu_SELF_INTEN",	
"RI_N_Gly_SELF_INTEN",	"RI_N_His_SELF_INTEN",	"RI_N_Ile_SELF_INTEN",	"RI_N_Leu_SELF_INTEN",	"RI_N_Lys_SELF_INTEN",	"RI_N_Met_SELF_INTEN",	"RI_N_Phe_SELF_INTEN",	"RI_N_Pro_SELF_INTEN",	"RI_N_Ser_SELF_INTEN",	"RI_N_Thr_SELF_INTEN",	"RI_N_Trp_SELF_INTEN",	"RI_N_Tyr_SELF_INTEN",	"RI_N_Val_SELF_INTEN",	"RI_IND_C_HAS_N_TERM",	"RI_IND_C_HAS_C_TERM",	"RI_IND_C_HAS_Gap",	"RI_IND_C_HAS_Xle",	
"RI_IND_C_HAS_Ala",	"RI_IND_C_HAS_Arg",	"RI_IND_C_HAS_Asn",	"RI_IND_C_HAS_Asp",	"RI_IND_C_HAS_Cys",	"RI_IND_C_HAS_Gln",	"RI_IND_C_HAS_Glu",	"RI_IND_C_HAS_Gly",	"RI_IND_C_HAS_His",	"RI_IND_C_HAS_Ile",	"RI_IND_C_HAS_Leu",	"RI_IND_C_HAS_Lys",	"RI_IND_C_HAS_Met",	"RI_IND_C_HAS_Phe",	"RI_IND_C_HAS_Pro",	"RI_IND_C_HAS_Ser",	"RI_IND_C_HAS_Thr",	"RI_IND_C_HAS_Trp",	"RI_IND_C_HAS_Tyr",	
"RI_IND_C_HAS_Val",	"RI_C_N_TERM_SELF_INTEN",	"RI_C_C_TERM_SELF_INTEN",	"RI_C_Gap_SELF_INTEN",	"RI_C_Xle_SELF_INTEN",	"RI_C_Ala_SELF_INTEN",	"RI_C_Arg_SELF_INTEN",	"RI_C_Asn_SELF_INTEN",	"RI_C_Asp_SELF_INTEN",	"RI_C_Cys_SELF_INTEN",	"RI_C_Gln_SELF_INTEN",	"RI_C_Glu_SELF_INTEN",	"RI_C_Gly_SELF_INTEN",	"RI_C_His_SELF_INTEN",	"RI_C_Ile_SELF_INTEN",	"RI_C_Leu_SELF_INTEN",	"RI_C_Lys_SELF_INTEN",	
"RI_C_Met_SELF_INTEN",	"RI_C_Phe_SELF_INTEN",	"RI_C_Pro_SELF_INTEN",	"RI_C_Ser_SELF_INTEN",	"RI_C_Thr_SELF_INTEN",	"RI_C_Trp_SELF_INTEN",	"RI_C_Tyr_SELF_INTEN",	"RI_C_Val_SELF_INTEN",	"RI_NUM_FEATURES",	"ScoreModelFields_RI"};

const char * ScoreModelFields_RNI_names[]={
"RNI_CONST",	"RNI_IND_DIS_FROM_MINMAX_LESS_50",	"RNI_DIS_FROM_MINMAX0",	"RNI_IND_DIS_FROM_MINMAX_LESS_150",	"RNI_DIS_FROM_MINMAX50",	"RNI_IND_DIS_FROM_MINMAX_LESS_250",	"RNI_DIS_FROM_MINMAX150",	"RNI_IND_DIS_FROM_MINMAX_MORE",	"RNI_DIS_FROM_MINMAX250",	"RNI_REL_POS0",	"RNI_REL_POS1",	"RNI_REL_POS2",	"RNI_REL_POS3",	"RNI_REL_POS4",	"RNI_REL_POS5",	"RNI_REL_POS6",	"RNI_REL_POS7",	
"RNI_REL_POS8",	"RNI_REL_POS9",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_0",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_1",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_2",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_3",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_4",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_5",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_6",	"RNI_IND_NUM_PARENTS_WITH_INTEN_IS_MORE6",	"RNI_IND_PARENT_COMBO_0",	"RNI_IND_PARENT_COMBO_1",	
"RNI_IND_PARENT_COMBO_2",	"RNI_IND_PARENT_COMBO_3",	"RNI_IND_PARENT_COMBO_4",	"RNI_IND_PARENT_COMBO_5",	"RNI_IND_PARENT_COMBO_6",	"RNI_IND_PARENT_COMBO_7",	"RNI_IND_PARENT_COMBO_8",	"RNI_IND_PARENT_COMBO_9",	"RNI_IND_PARENT_COMBO_10",	"RNI_IND_PARENT_COMBO_11",	"RNI_IND_PARENT_COMBO_12",	"RNI_IND_PARENT_COMBO_13",	"RNI_IND_PARENT_COMBO_14",	"RNI_IND_PARENT_COMBO_15",	"RNI_IND_GOT_BOTH_ORIS",	
"RNI_IND_GOT_PREFIX",	"RNI_IND_GOT_SUFFIX",	"RNI_IND_PARENT1_NOT_VIZ",	"RNI_IND_PARENT1_NO_INTEN",	"RNI_PARENT1_LOG_INTEN",	"RNI_PARENT1_LOG_GLOBAL_RANK",	"RNI_PARENT1_ISO_LEVEL",	"RNI_IND_PARENT2_NOT_VIZ",	"RNI_IND_PARENT2_NO_INTEN",	"RNI_PARENT2_LOG_INTEN",	"RNI_PARENT2_LOG_GLOBAL_RANK",	"RNI_PARENT2_ISO_LEVEL",	"RNI_IND_PARENT3_NOT_VIZ",	"RNI_IND_PARENT3_NO_INTEN",	"RNI_PARENT3_LOG_INTEN",	
"RNI_PARENT3_LOG_GLOBAL_RANK",	"RNI_PARENT3_ISO_LEVEL",	"RNI_IND_PARENT4_NOT_VIZ",	"RNI_IND_PARENT4_NO_INTEN",	"RNI_PARENT4_LOG_INTEN",	"RNI_PARENT4_LOG_GLOBAL_RANK",	"RNI_PARENT4_ISO_LEVEL",	"RNI_IND_PARENT5_NOT_VIZ",	"RNI_IND_PARENT5_NO_INTEN",	"RNI_PARENT5_LOG_INTEN",	"RNI_PARENT5_LOG_GLOBAL_RANK",	"RNI_PARENT5_ISO_LEVEL",	"RNI_IND_PARENT6_NOT_VIZ",	"RNI_IND_PARENT6_NO_INTEN",	
"RNI_PARENT6_LOG_INTEN",	"RNI_PARENT6_LOG_GLOBAL_RANK",	"RNI_PARENT6_ISO_LEVEL",	"RNI_IND_PARENT7_NOT_VIZ",	"RNI_IND_PARENT7_NO_INTEN",	"RNI_PARENT7_LOG_INTEN",	"RNI_PARENT7_LOG_GLOBAL_RANK",	"RNI_PARENT7_ISO_LEVEL",	"RNI_IND_PARENT8_NOT_VIZ",	"RNI_IND_PARENT8_NO_INTEN",	"RNI_PARENT8_LOG_INTEN",	"RNI_PARENT8_LOG_GLOBAL_RANK",	"RNI_PARENT8_ISO_LEVEL",	"RNI_IND_N_IS_GAP",	"RNI_IND_C_IS_GAP",	
"RNI_IND_N_HAS_N_TERM",	"RNI_IND_N_HAS_C_TERM",	"RNI_IND_N_HAS_Gap",	"RNI_IND_N_HAS_Xle",	"RNI_IND_N_HAS_Ala",	"RNI_IND_N_HAS_Arg",	"RNI_IND_N_HAS_Asn",	"RNI_IND_N_HAS_Asp",	"RNI_IND_N_HAS_Cys",	"RNI_IND_N_HAS_Gln",	"RNI_IND_N_HAS_Glu",	"RNI_IND_N_HAS_Gly",	"RNI_IND_N_HAS_His",	"RNI_IND_N_HAS_Ile",	"RNI_IND_N_HAS_Leu",	"RNI_IND_N_HAS_Lys",	"RNI_IND_N_HAS_Met",	"RNI_IND_N_HAS_Phe",	
"RNI_IND_N_HAS_Pro",	"RNI_IND_N_HAS_Ser",	"RNI_IND_N_HAS_Thr",	"RNI_IND_N_HAS_Trp",	"RNI_IND_N_HAS_Tyr",	"RNI_IND_N_HAS_Val",	"RNI_IND_C_HAS_N_TERM",	"RNI_IND_C_HAS_C_TERM",	"RNI_IND_C_HAS_Gap",	"RNI_IND_C_HAS_Xle",	"RNI_IND_C_HAS_Ala",	"RNI_IND_C_HAS_Arg",	"RNI_IND_C_HAS_Asn",	"RNI_IND_C_HAS_Asp",	"RNI_IND_C_HAS_Cys",	"RNI_IND_C_HAS_Gln",	"RNI_IND_C_HAS_Glu",	"RNI_IND_C_HAS_Gly",	
"RNI_IND_C_HAS_His",	"RNI_IND_C_HAS_Ile",	"RNI_IND_C_HAS_Leu",	"RNI_IND_C_HAS_Lys",	"RNI_IND_C_HAS_Met",	"RNI_IND_C_HAS_Phe",	"RNI_IND_C_HAS_Pro",	"RNI_IND_C_HAS_Ser",	"RNI_IND_C_HAS_Thr",	"RNI_IND_C_HAS_Trp",	"RNI_IND_C_HAS_Tyr",	"RNI_IND_C_HAS_Val",	"RNI_NUM_FEATURES",	"ScoreModelFields_RNI"};



void RegionalScoreModel::init(Config* _c, int _charge, int _size_idx, int _region_idx)
{
	config = _c;
	charge = _charge;
	size_idx = _size_idx;
	region_idx = _region_idx;

	const RegionalFragments& rf = config->get_regional_fragments(charge,size_idx,region_idx);

	frag_type_idxs = rf.get_frag_type_idxs();
	frag_probs     = rf.get_frag_probs();
	rand_prob	   = rf.get_rand_prob();
	log_random	   = log(rand_prob);
	log_one_minus_random = log(1.0-rand_prob);

	// init strong and regular models
	const vector<int>& strong_idxs = rf.get_strong_frag_type_idxs();
	const vector< vector<int> >& mirror_frag_idxs = rf.get_mirror_frag_idxs();
	
	strong_models.resize(strong_idxs.size());
	regular_models.resize(frag_type_idxs.size()-strong_idxs.size());

	num_strong_frags=0;
	num_regular_frags=0;
	int i;
	for (i=0; i<strong_idxs.size(); i++)
	{
		const int frag_idx = strong_idxs[i];

		StrongFragModel& strong_model=strong_models[num_strong_frags++];

		strong_model.model_frag_idx = frag_idx;
		strong_model.model_frag_charge = config->get_fragment(frag_idx).charge;
		strong_model.set_config_and_tolerance(config);

		strong_model.parent1_idx=-1;
		strong_model.parent2_idx=-1;
		strong_model.parent1_charge=-1;
		strong_model.parent2_charge=-1;

		// add parents
		int f;
		for (f=0; f<2 && f<frag_type_idxs.size(); f++)
		{
			if (frag_type_idxs[f] == frag_idx)
				break;

			const int parent_idx = frag_type_idxs[f];
			const int parent_charge = config->get_fragment(parent_idx).charge;
			if (f==0)
			{
				strong_model.parent1_idx = parent_idx;
				strong_model.parent1_charge = parent_charge;
			}
			else
			{
				strong_model.parent2_idx = parent_idx;
				strong_model.parent2_charge = parent_charge;
			}
		}

		// add mirrors
		for (f=0; f<2 && f<mirror_frag_idxs[frag_idx].size(); f++)
		{
			const int mirror_idx = mirror_frag_idxs[frag_idx][f];
			const int mirror_charge = config->get_fragment(mirror_idx).charge;
			if (f==0)
			{
				strong_model.mirror1_idx = mirror_idx;
				strong_model.mirror1_charge = mirror_charge;
			}
			else
			{
				strong_model.mirror2_idx = mirror_idx;
				strong_model.mirror2_charge = mirror_charge;
			}
		}
	}

	// add regular models
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		int j;
		for (j=0; j<strong_idxs.size(); j++)
			if (strong_idxs[j] == frag_type_idxs[i])
				break;

		if (j<strong_idxs.size())
			continue;

		const int reg_frag_idx = frag_type_idxs[i];

		RegularFragModel& regular_model=regular_models[num_regular_frags++];
		regular_model.model_frag_idx = reg_frag_idx;
		regular_model.model_frag_charge = config->get_fragment(reg_frag_idx).charge;
		regular_model.set_config_and_tolerance(config);
		regular_model.parent_idxs.clear();
		regular_model.parent_idx_with_same_charge_ori = -1;

		// add parents
		int f;
		for (f=0; f<frag_type_idxs.size() && f<8; f++)
		{
			if (frag_type_idxs[f] == reg_frag_idx)
				break;

			regular_model.parent_idxs.push_back(frag_type_idxs[f]);

			if (regular_model.parent_idx_with_same_charge_ori<0)
			{
				const FragmentType& parent_frag = config->get_fragment(frag_type_idxs[f]);
				if (parent_frag.charge == regular_model.model_frag_charge &&
					parent_frag.orientation == config->get_fragment(reg_frag_idx).orientation)
				{
					regular_model.parent_idx_with_same_charge_ori=frag_type_idxs[f];
				}
			}
		}
		regular_model.num_parents = regular_model.parent_idxs.size();
	}

	// set scores
	const int num_all_frags = config->get_all_fragments().size();
	frag_inten_scores.resize(num_all_frags,0);
	frag_no_inten_scores.resize(num_all_frags,0);
	missing_breakage_score=0;
	
	for (i=0; i<frag_type_idxs.size(); i++)
	{
		const int frag_idx=frag_type_idxs[i];
		const score_t missing_score = log(1.0 - frag_probs[i]) - log_one_minus_random;
		missing_breakage_score += missing_score;

		frag_no_inten_scores[frag_idx] = missing_score;
		frag_inten_scores[frag_idx] = log(frag_probs[i])- log_random;
	}

	//  init the weights 
//	vector<float> strong_inten_weights,  strong_no_inten_weights;
//	vector<float> regular_inten_weights, regular_no_inten_weights;

	const int num_strong = strong_models.size();
	strong_inten_weights.resize(num_strong,0.99);
	strong_no_inten_weights.resize(num_strong,0.98);

	strong_inten_danc_part.resize(num_strong);
	strong_no_inten_danc_part.resize(num_strong);

	for (i=0; i<num_strong; i++)
	{
		const float frag_prob = get_frag_prob(strong_models[i].model_frag_idx);
		strong_inten_danc_part[i]=(1.0-strong_inten_weights[i])*frag_prob;
		strong_no_inten_danc_part[i]=(1.0-strong_no_inten_weights[i])*(1.0-frag_prob);

//		if (size_idx == 1 && region_idx == 1)
//			cout << i << "\t" << frag_prob << "\t" << strong_inten_danc_part[i] << "\t" << strong_no_inten_danc_part[i] << endl;
	}


	// holds the (1-weight) * frag prob to be weight*model prob for the acutal probability
//	vector<float> strong_inten_danc_part,  strong_no_inten_danc_part;
//	vector<float> regular_inten_danc_part, regular_no_inten_danc_part;

	const int num_regular = regular_models.size();
	regular_inten_weights.resize(num_regular,0.97);
	regular_no_inten_weights.resize(num_regular,0.99);

	regular_inten_danc_part.resize(num_regular);
	regular_no_inten_danc_part.resize(num_regular);

	for (i=0; i<num_regular; i++)
	{
		regular_inten_weights[i] = 1.0 - 0.01*i;
		regular_no_inten_weights[i] = 0.5;
	}

	for (i=0; i<num_regular; i++)
	{
		const float frag_prob = get_frag_prob(regular_models[i].model_frag_idx);
		regular_inten_danc_part[i]=(1.0-regular_inten_weights[i])*frag_prob;
		regular_no_inten_danc_part[i]=(1.0-regular_no_inten_weights[i])*(1.0-frag_prob);

	//	if (size_idx == 1 && region_idx == 1)
	//		cout << i << "\t" << frag_prob << "\t" << regular_inten_danc_part[i] << "\t" << regular_no_inten_danc_part[i] << endl;
	}
}



/***************************************************************************
****************************************************************************/
void RegionalScoreModel::create_training_set(Model *model,
											 const FragModel& frag_model,
											 const FileManager& fm,
											 ME_Regression_DataSet& inten_ds,
											 ME_Regression_DataSet& no_inten_ds) const
{
	const mass_t min_mass = config->get_min_mass_for_size_idx(charge,size_idx);
	const mass_t max_mass = config->get_max_mass_for_size_idx(charge,size_idx);
	const mass_t tolerance_diff = config->get_tolerance()*0.5;
	const int frag_idx = frag_model.model_frag_idx;

	bool is_strong=false;
	int s;
	for (s=0; s<strong_models.size(); s++)
		if (frag_idx == strong_models[s].model_frag_idx)
			is_strong=true;

	FileSet fs;
	fs.select_files_in_mz_range(fm,min_mass/charge-2.0,max_mass/charge+2.0,charge);
	const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();

	cout << "Selected " << all_ssfs.size() << " headers..." << endl;
	cout << "Min m/z " << min_mass/charge << endl;
	cout << "Max m/z " << max_mass/charge << endl;

	PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model->get_pmcsqs_ptr();
	const bool use_pmcr = (pmcsqs && pmcsqs->get_ind_initialized_pmcr());

	if (use_pmcr)
		cout << "Using parent mass correction model to set PM.." << endl;

	BasicSpecReader bsr;
	vector<QCPeak> peaks;
	peaks.resize(10000);
	int i;
	for (i=0; i<all_ssfs.size(); i++)
	{
		BasicSpectrum bs;
		Spectrum s;
		PrmGraph prm;

		const mass_t true_mass_with_19 = all_ssfs[i]->peptide.get_mass_with_19();
		const int num_peaks = bsr.read_basic_spec(config,fm,all_ssfs[i],&peaks[0]);
		bs.num_peaks = num_peaks;
		bs.ssf = all_ssfs[i];
		bs.peaks = &peaks[0];

		// calc corrected pm_with_19, if it is good use it, otherwise, use a value with +- U[0,toleance/2]
		mass_t pm_with_19=NEG_INF;
	
		if (use_pmcr)
		{
			mass_t	mz1,mz2;
			int		charge1,charge2;
			float	prob1,prob2;
					

			// output m/z and prob values for the different charge states
			model->get_best_mz_charge(config,bs,&mz1,&charge1,&prob1,&mz2,&charge2,&prob2);

			const mass_t corr1_pm_with_19 = mz1*charge1 - MASS_PROTON*(charge1-1);
			const mass_t corr2_pm_with_19 = mz2*charge2 - MASS_PROTON*(charge2-1);
			if (fabs(corr2_pm_with_19-true_mass_with_19)<tolerance_diff)
				pm_with_19 = corr2_pm_with_19;
			if (fabs(corr1_pm_with_19-true_mass_with_19)<tolerance_diff)
				pm_with_19 = corr1_pm_with_19;
		}

		if (pm_with_19<0) // use a random value
		{
			double r=my_random();
			mass_t offset = r*r*tolerance_diff;
			if (my_random()<0.5)
				offset *= -1;
			pm_with_19 = true_mass_with_19 + offset;
	
		}

		int spec_size_idx = config->calc_size_idx(charge,pm_with_19);
		if (spec_size_idx != size_idx)
			continue;
	
		s.init_from_QCPeaks(config,&peaks[0],num_peaks,all_ssfs[i]);
		s.set_corrected_pm_with_19(pm_with_19);
		prm.create_graph_from_spectrum(model,&s,pm_with_19,all_ssfs[i]->charge);

	//	cout << i << endl;

	//	s.print_expected_by();
	//	s.print_spectrum();
	//	prm.print_with_combo_tables();
	//	prm.print();
	//	exit(0);

		vector<BreakageInfo> good_peak_examples, bad_peak_examples;
		prm.extract_breakage_infos_for_score_training(model,
													  frag_model.model_frag_idx, 
													  region_idx,
													  is_strong,
													  good_peak_examples, 
													  bad_peak_examples);

	//	cout << "GOOD:" << endl;
		int j;
		for (j=0; j<good_peak_examples.size(); j++)
		{
		//	good_peak_examples[j].print(config);
			ME_Regression_Sample sam;
			sam.label=0;

			frag_model.fill_single_frag_vector(&s, pm_with_19, good_peak_examples[j].breakage,
				good_peak_examples[j], sam.f_vals);

			if (good_peak_examples[j].breakage->get_position_of_frag_idx(frag_idx)>=0)
			{
				inten_ds.add_sample(sam);
			//	sam.print((is_strong ? ScoreModelFields_SI_names : ScoreModelFields_RI_names));
			}
			else
				no_inten_ds.add_sample(sam);
		}
		
	//	cout << "BAD:" << endl;
		for (j=0; j<bad_peak_examples.size(); j++)
		{
		//	bad_peak_examples[j].print(config);
			ME_Regression_Sample sam;
			sam.label=1;

			frag_model.fill_single_frag_vector(&s, pm_with_19, bad_peak_examples[j].breakage,
				bad_peak_examples[j], sam.f_vals);

			if (bad_peak_examples[j].breakage->get_position_of_frag_idx(frag_idx)>=0)
			{
				inten_ds.add_sample(sam);
			}
			else
				no_inten_ds.add_sample(sam);

	//		sam.print();
		}

		if (i>0 && i %1000 == 0)
		{
			cout << i << "/" << all_ssfs.size() << " ..." << endl;
		}
	}

	
	const score_t frag_prob  = this->get_frag_prob(frag_idx);
	cout << "Probability of observeing fragment: " << frag_prob << endl;

	inten_ds.num_classes=2;
	inten_ds.num_features = (is_strong ? (int)SI_NUM_FEATURES : (int)RI_NUM_FEATURES);
	inten_ds.calibrate_class_weights((is_strong ? 0.45 : 0.2));
	
	no_inten_ds.num_classes=2;
	no_inten_ds.num_features = (is_strong ? (int)SNI_NUM_FEATURES : (int)RNI_NUM_FEATURES);
	no_inten_ds.calibrate_class_weights((is_strong ? 0.2 : 0.5));

	if (is_strong)
	{
		vector<int> inten_features;
		inten_features.push_back(SI_IND_CONNECTS_TO_N_TERM);
		inten_features.push_back(SI_IND_N_IS_GAP);
		inten_features.push_back(SI_IND_C_IS_GAP);

		inten_ds.serial_scale(inten_features);

		vector<int> no_inten_features;
		no_inten_features.push_back(SNI_IND_CONNECTS_TO_N_TERM);
		no_inten_features.push_back(SNI_IND_N_IS_GAP);
		no_inten_features.push_back(SNI_IND_C_IS_GAP);

		no_inten_ds.serial_scale(no_inten_features);
	}
	else
	{
		vector<int> inten_features;
		inten_features.push_back(RI_IND_N_IS_GAP);
		inten_features.push_back(RI_IND_C_IS_GAP);

		inten_ds.serial_scale(inten_features);

		vector<int> no_inten_features;
		no_inten_features.push_back(RNI_IND_N_IS_GAP);
		no_inten_features.push_back(RNI_IND_C_IS_GAP);

		no_inten_ds.serial_scale(no_inten_features);
	}
}


bool RegionalScoreModel::train_regional_score_model(Model *model, const char *name, const FileManager& fm)
{
	const int min_num_of_samples_per_feature=6;
	const int num_me_rounds = 750;
	int i;
	for (i=0; i<strong_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = strong_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge << " SIZE " << size_idx << " REGION " << region_idx << " FRAGMENT " << i << " " <<
			config->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << get_random_seed() << endl;
			create_training_set(model, strong_models[i], fm, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_SI_names);
			if (! strong_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j << ")" << endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge << " size " << size_idx << " region " << region_idx << 
					" fragment " << frag_idx << " " << config->get_fragment(frag_idx).label <<  endl;
				strong_models[i].inten_model.print_ds_probs(inten_ds);
				strong_models[i].inten_log_scaling_factor = 
					log(strong_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.2*frag_prob));

				break;
			}
		}

		if (j==3)
		{
			cout << "Problem with models!!!" << endl;
			exit(1);
		}

		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
			charge << " SIZE " << size_idx << " REGION " << region_idx << " FRAGMENT " << i << " " <<
			config->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
		int num_bad=no_inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;

		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_SNI_names);
		if (! strong_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train no inten ME model, exiting!" << endl;
			exit(1);
		}
		else
		{
			cout << endl << "NO INTENSITY - Charge " << charge << " size " << size_idx << " region " << region_idx << 
				" fragment " << config->get_fragment(frag_idx).label << endl;
			strong_models[i].no_inten_model.print_ds_probs(no_inten_ds);
			strong_models[i].no_inten_log_scaling_factor = 
				log(strong_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,1.1*(1.0-frag_prob)));
		}

		strong_models[i].ind_has_models = true;
	}

	for (i=0; i<regular_models.size(); i++)
	{
		ME_Regression_DataSet inten_ds, no_inten_ds;
		const int frag_idx = regular_models[i].model_frag_idx;
		const float frag_prob = get_frag_prob(frag_idx);

		cout << endl << endl << "TRAINING INTENSITY MODEL FOR CHARGE " <<
			charge << " SIZE " << size_idx << " REGION " << region_idx << " FRAGMENT " << i << " " <<
			config->get_fragment(frag_idx).label << endl << endl;

		int j;
		for (j=0; j<3; j++)
		{
			cout << "SEED: " << get_random_seed() << endl;
			create_training_set(model, regular_models[i], fm, inten_ds, no_inten_ds);

			inten_ds.purge_low_count_features(min_num_of_samples_per_feature);
			int num_bad=inten_ds.check_samples(true);
			if (num_bad>0)
				cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


			inten_ds.print_summary();
			inten_ds.print_feature_summary(cout, ScoreModelFields_RI_names);
			if (! regular_models[i].inten_model.train_cg(inten_ds,num_me_rounds,2E-5))
			{
				cout << "Coudln't train ME model, setting all weights to 0! (" <<j<<")"<< endl;
			}
			else
			{
				cout << endl << "INTENSTY - Charge " << charge << " size " << size_idx << " region " << region_idx << 
					" fragment " << config->get_fragment(frag_idx).label << endl ;
				regular_models[i].inten_model.print_ds_probs(inten_ds);
				regular_models[i].inten_log_scaling_factor = 
					log(regular_models[i].inten_model.calc_log_scaling_constant(0,inten_ds,1.1*frag_prob));
				break;
			}
		}
		
		cout << endl << endl << "TRAINING NO INTENSITY MODEL FOR CHARGE " <<
				charge << " SIZE " << size_idx << " REGION " << region_idx << " FRAGMENT " << i << " " <<
				config->get_fragment(frag_idx).label << endl << endl;

		no_inten_ds.purge_low_count_features(min_num_of_samples_per_feature);

		int num_bad=no_inten_ds.check_samples(true);
		if (num_bad>0)
			cout << "Warning: had " << num_bad << " bad samples removed!" << endl;


		no_inten_ds.print_summary();
		no_inten_ds.print_feature_summary(cout, ScoreModelFields_RNI_names);
		if (! regular_models[i].no_inten_model.train_cg(no_inten_ds,num_me_rounds,2E-5))
		{
			cout << "Coudln't train ME model, setting all weights to 0!" << endl;
		}
		else
		{
			cout << endl << "NO INTENSTY - Charge " << charge << " size " << size_idx << " region " << region_idx << 
				" fragment " << config->get_fragment(frag_idx).label << endl;
				regular_models[i].no_inten_model.print_ds_probs(no_inten_ds);
				regular_models[i].no_inten_log_scaling_factor = 
					log(regular_models[i].no_inten_model.calc_log_scaling_constant(0,no_inten_ds,(1.0-frag_prob)));
		}

		regular_models[i].ind_has_models = true;
	}	


	was_initialized = true;

	write_regional_score_model(name);

	return true;
}




void BreakageInfo::print(Config *config) const
{
	const vector<string>& aa2label = config->get_aa2label();
	cout << setw(2) << this->type << "> ";
	cout << (this->connects_to_N_term ? '[' : ' ');
	cout << (this->n_edge_is_single ? '1' : '2');
	cout <<  " " << aa2label[this->n_aa] << " : " << aa2label[this->c_aa] << " ";
	cout << (this->c_edge_is_single ? '1' : '2');
	cout << (this->connects_to_C_term ? ']' : ' ');
	cout << " " << fixed << setprecision(3) << "  # " <<this->node_idx << " ";
	if (this->missed_cleavage)
		cout << "MC!";
	cout << endl;
}

void PrmGraph::fill_breakage_info(const Model *model, BreakageInfo *info, int node_idx,
								  int n_edge_idx, int n_variant_idx,
								  int c_edge_idx, int c_variant_idx, int type) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	const vector<int>& n_term_digest_aas = config->get_n_term_digest_aas();
	const vector<int>& c_term_digest_aas = config->get_c_term_digest_aas();
	
	info->node_idx = node_idx;
	info->breakage = &nodes[node_idx].breakage;
	info->type = type;

	info->n_edge_idx = n_edge_idx;
	info->n_var_idx  = n_variant_idx;
	info->c_edge_idx = c_edge_idx;
	info->c_var_idx  = c_variant_idx;

	int n_second_before_cut = NEG_INF;
	if (n_edge_idx>=0)
	{
		const MultiEdge& n_edge = multi_edges[n_edge_idx];
		const Node& n_node = nodes[n_edge.n_idx];
		int *v_ptr = n_edge.variant_ptrs[n_variant_idx];

		info->n_var_ptr = v_ptr;

		const int num_aa = *v_ptr++;
		int n_aa = v_ptr[num_aa-1];
		mass_t exp_edge_mass=0;
		int j;
		for (j=0; j<num_aa; j++)
			exp_edge_mass+=aa2mass[v_ptr[j]];

		if (n_aa == Ile || n_aa == Xle)
			n_aa = Leu;

		info->n_aa = n_aa;
		info->nn_aa = v_ptr[0];
		
		if (info->nn_aa == Ile || info->nn_aa == Xle)
			info->nn_aa = Leu;

		info->n_break = n_edge.n_break;
		info->exp_n_edge_mass = exp_edge_mass;
		info->n_edge_is_single = (num_aa ==1);

		if (! info->n_edge_is_single)
			n_second_before_cut = v_ptr[num_aa-2];
		if (n_second_before_cut == Ile || n_second_before_cut == Xle)
			n_second_before_cut = Leu;

		if (n_node.type == NODE_N_TERM)
		{
			info->connects_to_N_term=true;
			int j;
			for (j=0; j<n_term_digest_aas.size(); j++)
				if (n_term_digest_aas[j] == *v_ptr)
				{
					info->preferred_digest_aa_N_term=true;
					break;
				}
		}

		for (j=0; j<c_term_digest_aas.size(); j++)
			if (c_term_digest_aas[j] == v_ptr[num_aa-1])
			{
				info->missed_cleavage= true;
				break;
			}

		info->ind_n_edge_overlaps = (n_edge.ind_edge_overlaps);
		info->n_side_cat = model->get_aa_category(num_aa,v_ptr,info->connects_to_N_term,false);
	}
	else
	{
		info->n_aa=Gap;
		info->nn_aa=Gap;
	}

	if (c_edge_idx>=0)
	{
		const MultiEdge& c_edge = multi_edges[c_edge_idx];
		const Node& c_node = nodes[c_edge.c_idx];
		int *v_ptr = c_edge.variant_ptrs[c_variant_idx];

		info->c_var_ptr = v_ptr;

		const int num_aa = *v_ptr++;
		int c_aa = v_ptr[0];
		mass_t exp_edge_mass=0;
		int j;
		for (j=0; j<num_aa; j++)
			exp_edge_mass+=aa2mass[v_ptr[j]];

		if (c_aa == Ile || c_aa == Xle)
			c_aa = Leu;

		info->c_aa = c_aa;
		info->cc_aa = v_ptr[num_aa-1];

		if (info->cc_aa == Ile || info->cc_aa == Xle)
			info->cc_aa = Leu;

		info->c_break = c_edge.c_break;
		info->exp_c_edge_mass = exp_edge_mass;
		info->c_edge_is_single = (num_aa ==1);

		int c_second_after_cut=NEG_INF;
		if (! info->c_edge_is_single)
			c_second_after_cut = v_ptr[1];
		if (c_second_after_cut == Ile || c_second_after_cut == Xle)
			c_second_after_cut = Leu;

		if (c_node.type == NODE_C_TERM)
		{
			info->connects_to_C_term=true;
			int j;
			for (j=0; j<c_term_digest_aas.size(); j++)
				if (c_term_digest_aas[j] == v_ptr[num_aa-1])
				{
					info->preferred_digest_aa_C_term=true;
					break;
				}
		}

		for (j=0; j<n_term_digest_aas.size(); j++)
			if (n_term_digest_aas[j] == v_ptr[0])
			{
				info->missed_cleavage= true;
				break;
			}

		info->ind_c_edge_overlaps = (c_edge.ind_edge_overlaps);
		info->c_side_cat = model->get_aa_category(num_aa,v_ptr,false,info->connects_to_C_term);
		if (n_edge_idx>=0)
		{
			int aas[2]={info->n_aa,info->c_aa};
			info->span_cat = model->get_aa_category(2,aas,info->connects_to_N_term && info->n_edge_is_single,info->connects_to_C_term && info->c_edge_is_single);
			
			if (! info->n_edge_is_single)
			{
				if (n_second_before_cut<0)
				{
					cout << "Error: bad aa 2nd before n-side!" << endl;
					exit(1);
				}
				int aas[3]={n_second_before_cut,info->n_aa,info->c_aa};
				info->n_double_span_cat = model->get_aa_category(3,aas,info->connects_to_N_term ,info->connects_to_C_term && info->c_edge_is_single);
			}

			if (! info->c_edge_is_single)
			{
				if (c_second_after_cut<0)
				{
					cout << "Error: bad aa 2nd after c-side!" << endl;
					exit(1);
				}
				int aas[3]={info->n_aa,info->c_aa,c_second_after_cut};
				info->c_double_span_cat = model->get_aa_category(3,aas,info->connects_to_N_term && info->n_edge_is_single, info->connects_to_C_term);
			}
		}

	}
	else
	{
		info->c_aa=Gap;
		info->cc_aa=Gap;
	}
}

/************************************************************************************************
*************************************************************************************************/
void PrmGraph::extract_breakage_infos_for_score_training(Model *model,
														 int frag_idx,
														 int target_region_idx,
														 bool ind_strong_frag,
														 vector<BreakageInfo>& good_examples,
														 vector<BreakageInfo>& bad_examples) const
{

	const double Gap_ratio = 0.05;
	const mass_t tolerance = config->get_tolerance();

	const Peptide& true_pep = source_spectrum->get_peptide();
	const vector<int>& aas = true_pep.get_amino_acids();
	vector<mass_t> correct_break_masses;
	vector<int>    node_to_breakages, correct_node_idxs;
	vector<int>    correct_edge_variant_map;
	
	// find minimal and maximal nodes for which the frag is visible
	const mass_t min_mass = source_spectrum->get_min_peak_mass()-1.0;
	const mass_t max_mass = source_spectrum->get_max_peak_mass()+1.0;
	const FragmentType& frag = config->get_fragment(frag_idx);
	int min_viz_idx=nodes.size()+1;
	int max_viz_idx=0;
	int i;
	for (i=0; i<nodes.size(); i++)
	{
		mass_t exp_mass = frag.calc_expected_mass(nodes[i].mass,pm_with_19);
		if (exp_mass>min_mass && exp_mass<max_mass)
		{
			if (i<min_viz_idx)
				min_viz_idx=i;
			if (i>max_viz_idx)
				max_viz_idx=i;
		}
	}

	good_examples.clear();
	bad_examples.clear();

	true_pep.calc_expected_breakage_masses(config,correct_break_masses);
	node_to_breakages.resize(nodes.size(),NEG_INF);
	correct_node_idxs.clear();
	correct_edge_variant_map.resize(multi_edges.size(),NEG_INF);

	
	for (i=0; i<correct_break_masses.size(); i++)
	{
		const int max_node_idx = get_max_score_node(correct_break_masses[i],tolerance);
		if (max_node_idx>=0)
		{
			node_to_breakages[max_node_idx]=i;
			correct_node_idxs.push_back(max_node_idx);
		}
	}

//	cout << endl <<true_pep.as_string(config) << endl;
	for (i=0; i<multi_edges.size(); i++)
	{
		const MultiEdge& edge = multi_edges[i];
		if (node_to_breakages[edge.n_idx]>=0 && node_to_breakages[edge.c_idx]>=0)
		{
			const int brekage_idx = node_to_breakages[edge.n_idx];
			correct_edge_variant_map[i]=multi_edges[i].get_variant_idx(edge.num_aa,&aas[brekage_idx]);

			if (0 && correct_edge_variant_map[i]>=0)
			{
				cout << brekage_idx << "\t" << i << "\t" << correct_edge_variant_map[i] << "\t";
				int *var_ptr = edge.variant_ptrs[correct_edge_variant_map[i]];
				int num_aa = *var_ptr++;
				int j;
				for (j=0; j<num_aa; j++)
					cout << config->get_aa2label()[*var_ptr++];
				cout << endl;
			}
		}
	}

	vector<int> correct_in_edge_idxs;
	vector<int> correct_out_edge_idxs;
	correct_in_edge_idxs.resize(correct_node_idxs.size(),NEG_INF);
	correct_out_edge_idxs.resize(correct_node_idxs.size(),NEG_INF);

	bool had_good_connect_to_n_term = false;
	bool had_good_connect_to_c_term = false;

	// create infos for good peak samples
	for (i=0; i<correct_node_idxs.size(); i++)
	{
		const int node_idx = correct_node_idxs[i];
		const int node_region_idx = nodes[node_idx].breakage.region_idx;
		if (node_region_idx != target_region_idx || node_idx==0 || 
			node_idx==nodes.size()-1 || node_idx<min_viz_idx || node_idx>max_viz_idx)
			continue;

		const Node& node = nodes[node_idx];
		
		int n_edge_idx = NEG_INF;
		int c_edge_idx = NEG_INF;
		int n_varaint_idx = NEG_INF;
		int c_variant_idx = NEG_INF;

		int j;
		for (j=0; j<node.in_edge_idxs.size(); j++)
			if (correct_edge_variant_map[node.in_edge_idxs[j]]>=0)
			{
				n_edge_idx=node.in_edge_idxs[j];
				n_varaint_idx=correct_edge_variant_map[node.in_edge_idxs[j]];
				correct_in_edge_idxs[i]=n_edge_idx;
				break;
			}

		for (j=0; j<node.out_edge_idxs.size(); j++)
			if (correct_edge_variant_map[node.out_edge_idxs[j]]>=0)
			{
				c_edge_idx=node.out_edge_idxs[j];
				c_variant_idx=correct_edge_variant_map[node.out_edge_idxs[j]];
				correct_out_edge_idxs[i]=c_edge_idx;
				break;
			}

		BreakageInfo info;
		fill_breakage_info(model,&info,node_idx,n_edge_idx,n_varaint_idx,c_edge_idx,c_variant_idx,1);

		if (info.connects_to_N_term)
			had_good_connect_to_n_term=true;
		
		if (info.connects_to_C_term)
			had_good_connect_to_c_term=true;

		good_examples.push_back(info);
		if (n_edge_idx>=0 && c_edge_idx>=0 && my_random()<Gap_ratio*0.33)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,c_edge_idx,c_variant_idx,11);
			good_examples.push_back(gap_info);
		}

		if (n_edge_idx>=0 && c_edge_idx>=0 && my_random()<Gap_ratio*0.33)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,n_edge_idx,n_varaint_idx,NEG_INF,NEG_INF,11);
			good_examples.push_back(gap_info);
		}
	}

	// create for bad peak samples
	const int num_good = good_examples.size();

	// same as good nodes, but using bad edges
	// perform only for strong fragments!
	const int num_half_bad_examples = (ind_strong_frag ? int(num_good*0.333) : 0);
	for (i=0; i<num_half_bad_examples; i++)
	{
		const int node_idx = correct_node_idxs[i];
		const int node_region_idx = nodes[node_idx].breakage.region_idx;

		if (node_region_idx != target_region_idx || node_idx==0 || 
			node_idx==nodes.size()-1 || node_idx<min_viz_idx || node_idx>max_viz_idx)
			continue;

		const Node& node = nodes[node_idx];
		if (node.breakage.get_position_of_frag_idx(frag_idx)<0) // don't use such samples if there is no peak
			continue;

		if (node.in_edge_idxs.size()>1 && correct_in_edge_idxs[i]>=0 &&
			node.out_edge_idxs.size()>1 && correct_out_edge_idxs[i]>=0)
		{
			int j;
			vector<int> in_idxs,out_idxs;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;
				if (multi_edges[edge_idx].num_aa == 1)
				{
					in_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || my_random()<0.1))
					in_idxs.push_back(edge_idx);
			}
			
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;
				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<3 || my_random()<0.1))
					out_idxs.push_back(edge_idx);
			}

			if (in_idxs.size()==0 || out_idxs.size()==0)
				continue;

			const int bad_n_edge = in_idxs[int(my_random()*in_idxs.size())];
			const int bad_c_edge = out_idxs[int(my_random()*out_idxs.size())];
			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * my_random());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * my_random());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,2);
			bad_examples.push_back(info);
			if (bad_n_edge>=0 && bad_c_edge>=0 && my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,22);
				bad_examples.push_back(gap_info);
			}

			if (bad_n_edge>=0 && bad_c_edge>=0 && my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,22);
				bad_examples.push_back(gap_info);
			}
		}
	}


	// mostly single edges from the highest scoring incorrect nodes
	// gives advantage for aa combinations with high scoring nodes
	vector<score_pair> pairs;
	for (i=1; i<nodes.size()-1; i++)
		if (node_to_breakages[i]<0 && i>=min_viz_idx && i<=max_viz_idx &&
			nodes[i].breakage.region_idx==target_region_idx )
				pairs.push_back(score_pair(i,nodes[i].score));
	
	sort(pairs.begin(),pairs.end());

	vector<int> selected_idxs, random_idxs;
	int num_to_select = int(num_good * 3);
	for (i=0; i<num_to_select && i<pairs.size(); i++)
		selected_idxs.push_back(pairs[i].idx);

	if (pairs.size()< num_to_select)
		num_to_select = pairs.size();

	choose_k_from_n(num_to_select,pairs.size(),random_idxs);

	for (i=0; i<random_idxs.size(); i++)
		selected_idxs.push_back(pairs[random_idxs[i]].idx);

	for (i=0; i<selected_idxs.size(); i++)
	{
		const int node_idx = selected_idxs[i];
		const Node& node = nodes[node_idx];
		vector<int> in_idxs,out_idxs;
		score_t max_n_score=NEG_INF;
		score_t max_c_score=NEG_INF;
		int best_n_idx=-1;
		int best_c_idx=-1;
		if (node.in_edge_idxs.size()>0 && node.out_edge_idxs.size()>0)
		{
			int j;
			for (j=0; j<node.in_edge_idxs.size(); j++)
			{
				const int edge_idx = node.in_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					in_idxs.push_back(edge_idx);
					if (nodes[multi_edges[edge_idx].n_idx].score>max_n_score)
					{
						max_n_score=nodes[multi_edges[edge_idx].n_idx].score;
						best_n_idx=edge_idx;
					}
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<2 || my_random()<0.1))
					in_idxs.push_back(edge_idx);
			}

			if (in_idxs.size()==0)
				continue;
			
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
					if (nodes[multi_edges[edge_idx].c_idx].score>max_c_score)
					{
						max_c_score=nodes[multi_edges[edge_idx].c_idx].score;
						best_c_idx=edge_idx;
					}
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<2 || my_random()<0.1))
					out_idxs.push_back(edge_idx);
			}
		}
		if (in_idxs.size()==0 || out_idxs.size()==0)
			continue;

		// add entries with the highest scoring nodes to bias the results
		if (best_n_idx>=0)
		{
			int num_to_add = int(0.3 * in_idxs.size()+0.5);
			int k;
			for (k=0; k<num_to_add; k++)
			{
				in_idxs.push_back(best_n_idx);
				in_idxs.push_back(best_n_idx);
			}
		}

		if (best_c_idx>=0)
		{
			int num_to_add = int(0.3 * out_idxs.size()+0.5);
			int k;
			for (k=0; k<num_to_add; k++)
			{
				out_idxs.push_back(best_c_idx);
				out_idxs.push_back(best_c_idx);
			}
		}

		const int bad_n_edge = in_idxs[int(my_random()*in_idxs.size())];
		const int bad_c_edge = out_idxs[int(my_random()*out_idxs.size())];
		const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * my_random());
		const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * my_random());

		BreakageInfo info;
		fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,4);
		bad_examples.push_back(info);

		if (my_random()<Gap_ratio)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,44);
			bad_examples.push_back(gap_info);
		}

		if ( my_random()<Gap_ratio)
		{
			BreakageInfo gap_info;
			fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,44);
			bad_examples.push_back(gap_info);
		}
	}

	// add a few nodes connecting to the N-terminal
	if (had_good_connect_to_n_term && nodes[0].out_edge_idxs.size()>2)
	{
		vector<int> in_idxs;
		int j;
		for (j=0; j<nodes[0].out_edge_idxs.size(); j++)
		{
			const int edge_idx = nodes[0].out_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				in_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || my_random()<0.1))
				in_idxs.push_back(edge_idx);
		}

		vector<int> selected_in_idxs;
		if (in_idxs.size()>3)
		{
			vector<int> positions;
			choose_k_from_n(3,in_idxs.size(),positions);
			int j;
			for (j=0; j<3; j++)
				selected_in_idxs.push_back(in_idxs[positions[j]]);
		}
		else
			selected_in_idxs = in_idxs;
		
		int k;
		for (k=0; k<selected_in_idxs.size(); k++)
		{
			const int bad_n_edge = selected_in_idxs[k];
			const MultiEdge& in_edge = multi_edges[bad_n_edge];
			const int node_idx = in_edge.c_idx;
			const Node& node = nodes[node_idx];

			vector<int> out_idxs;
			int j;
			for (j=0; j<node.out_edge_idxs.size(); j++)
			{
				const int edge_idx = node.out_edge_idxs[j];
				if (correct_edge_variant_map[edge_idx]>=0)
					continue;

				if (multi_edges[edge_idx].num_aa == 1)
				{
					out_idxs.push_back(edge_idx);
				}
				else if (multi_edges[edge_idx].num_aa != 1 && (out_idxs.size()<3 || my_random()<0.1))
					out_idxs.push_back(edge_idx);
			}

			if (out_idxs.size() == 0)
				continue;

			const int bad_c_edge = out_idxs[int(my_random()*out_idxs.size())];

			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * my_random());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * my_random());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,5);
			bad_examples.push_back(info);

			if (my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,55);
				bad_examples.push_back(gap_info);
			}

			if ( my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,55);
				bad_examples.push_back(gap_info);
			}
		}
	}

	// add an example that connects to the C-terminal
	if (had_good_connect_to_c_term && nodes[nodes.size()-1].in_edge_idxs.size()>1)
	{
		vector<int> out_idxs;
		int j;
		for (j=0; j<nodes[nodes.size()-1].in_edge_idxs.size(); j++)
		{
			const int edge_idx = nodes[nodes.size()-1].in_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				out_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && my_random()<0.025)
				out_idxs.push_back(edge_idx);
		}

		if (out_idxs.size() == 0)
			return;

		const int bad_c_edge = out_idxs[int(my_random()*out_idxs.size())];
		
		const MultiEdge& out_edge = multi_edges[bad_c_edge];
		const int node_idx = out_edge.n_idx;
		const Node& node = nodes[node_idx];

		vector<int> in_idxs;
		for (j=0; j<node.in_edge_idxs.size(); j++)
		{
			const int edge_idx = node.in_edge_idxs[j];
			if (correct_edge_variant_map[edge_idx]>=0)
				continue;

			if (multi_edges[edge_idx].num_aa == 1)
			{
				in_idxs.push_back(edge_idx);
			}
			else if (multi_edges[edge_idx].num_aa != 1 && (in_idxs.size()<3 || my_random()<0.1))
				in_idxs.push_back(edge_idx);
		}

		if (in_idxs.size()>0)
		{
			const int bad_n_edge = in_idxs[int(my_random()*in_idxs.size())];

			const int n_var_idx = int(multi_edges[bad_n_edge].variant_ptrs.size() * my_random());
			const int c_var_idx = int(multi_edges[bad_c_edge].variant_ptrs.size() * my_random());

			BreakageInfo info;
			fill_breakage_info(model,&info,node_idx,bad_n_edge,n_var_idx,bad_c_edge,c_var_idx,6);
			bad_examples.push_back(info);

			if (my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,bad_n_edge,n_var_idx,NEG_INF,NEG_INF,66);
				bad_examples.push_back(gap_info);
			}

			if ( my_random()<Gap_ratio)
			{
				BreakageInfo gap_info;
				fill_breakage_info(model,&gap_info,node_idx,NEG_INF,NEG_INF,bad_c_edge,c_var_idx,66);
				bad_examples.push_back(gap_info);
			}
		}
		
	}
}

