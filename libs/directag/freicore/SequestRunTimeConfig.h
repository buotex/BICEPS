#ifndef _MYRIMATCHCONFIG_H
#define _MYRIMATCHCONFIG_H

#include "stdafx.h"
#include "freicore.h"
#include "BaseRunTimeConfig.h"

using namespace freicore;

#define SEQUEST_RUNTIME_CONFIG \
	RTCONFIG_VARIABLE( string,	database_name,						""											) \
	RTCONFIG_VARIABLE( float,	peptide_mass_tolerance,				2.5f										) \
	RTCONFIG_VARIABLE( int,		create_output_files,				1											) \
	RTCONFIG_VARIABLE( string,	ion_series,							"0 1 1 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0"	) \
	RTCONFIG_VARIABLE( float,	fragment_ion_tolerance,				0.0f										) \
	RTCONFIG_VARIABLE( int,		num_output_lines,					5											) \
	RTCONFIG_VARIABLE( int,		num_description_lines,				0											) \
	RTCONFIG_VARIABLE( int,		num_results,						500											) \
	RTCONFIG_VARIABLE( int,		show_fragment_ions,					0											) \
	RTCONFIG_VARIABLE( int,		print_duplicate_references,			1											) \
	RTCONFIG_VARIABLE( int,		enzyme_number,						1											) \
	RTCONFIG_VARIABLE( string,	diff_search_options,				""											) \
	RTCONFIG_VARIABLE( string,	term_diff_search_options,			"0.0 0.0"									) \
	RTCONFIG_VARIABLE( int,		max_num_differential_AA_per_mod,	2											) \
	RTCONFIG_VARIABLE( int,		nucleotide_reading_frame,			0											) \
	RTCONFIG_VARIABLE( int,		mass_type_parent,					0											) \
	RTCONFIG_VARIABLE( int,		mass_type_fragment,					1											) \
	RTCONFIG_VARIABLE( int,		remove_precursor_peak,				0											) \
	RTCONFIG_VARIABLE( float,	ion_cutoff_percentage,				0.0f										) \
	RTCONFIG_VARIABLE( string,	protein_mass_filter,				"0 0"										) \
	RTCONFIG_VARIABLE( int,		max_num_internal_cleavage_sites,	10											) \
	RTCONFIG_VARIABLE( int,		match_peak_count,					0											) \
	RTCONFIG_VARIABLE( int,		match_peak_allowed_error,			1											) \
	RTCONFIG_VARIABLE( float,	match_peak_tolerance,				1.0f										) \
	RTCONFIG_VARIABLE( string,	partial_sequence,					""											) \
	RTCONFIG_VARIABLE( string,	sequence_header_filter,				""											) \
	RTCONFIG_VARIABLE( float,	add_Cterm_protein,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Nterm_protein,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Cterm_peptide,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Nterm_peptide,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_G_Glycine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_A_Alanine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_S_Serine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_P_Proline,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_V_Valine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_T_Threonine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_C_Cysteine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_L_Leucine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_I_Isoleucine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_X_LorI,							0.0f										) \
	RTCONFIG_VARIABLE( float,	add_N_Asparagine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_O_Ornithine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_B_avg_NandD,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_D_Aspartic_Acid,				0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Q_Glutamine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_K_Lysine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Z_avg_QandE,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_E_Glutamic_Acid,				0.0f										) \
	RTCONFIG_VARIABLE( float,	add_M_Methionine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_H_Histidine,					0.0f										) \
	RTCONFIG_VARIABLE( float,	add_F_Phenylalanine,				0.0f										) \
	RTCONFIG_VARIABLE( float,	add_R_Arginine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_Y_Tyrosine,						0.0f										) \
	RTCONFIG_VARIABLE( float,	add_W_Tryptophan,					0.0f										)

namespace freicore
{
namespace sequest
{
	static const char userModChars[] = "*#@^~$";

	/*
	[SEQUEST_ENZYME_INFO]
	0.  No_Enzyme              0   -           -
	1.  Trypsin_Strict         1   KR   -
	2.  Trypsin                1   KRLNH   -
	3.  Chymotrypsin           1   FWYL   -
	4.  Chymotrypsin_WYF       1   FWY   -
	5.  Clostripain            1   R   -
	6.  Cyanogen_Bromide       1   M   -
	7.  IodosoBenzoate         1   W   -
	8.  Proline_Endopept       1   P   -
	9.  Staph_Protease         1   E   -
	10.  Pepsin                1   FLWYAEQ   -
	11.  Trypsin_R             1   R   P
	12.  GluC                  1   ED   -
	13.  LysC                  1   K   -
	14.  AspN                  0   D   -
	15.  Elastase              1   ALIV   P
	16.  Elastase/Tryp/Chymo   1   ALIVKRWFY   P
	17.  Trypsin/Chymo         1   KRLFWYNQCHSM   -
	*/
	static const char* enzymeNames[] =
	{
		"No_Enzyme",
		"Trypsin_Strict",
		"Trypsin",
		"Chymotrypsin",
		"Chymotrypsin_WYF",
		"Clostripain",
		"Cyanogen_Bromide",
		"IodosoBenzoate",
		"Proline_Endopept",
		"Staph_Protease",
		"Pepsin",
		"Trypsin_R",
		"GluC",
		"LysC",
		"AspN",
		"Elastase",
		"Elastase/Tryp/Chymo",
		"Trypsin/Chymo"
	};

	static const char* enzymeRules[] =
	{
		". .",
		"[|K|R . . ]",
		"[|K|R|L|N|H . . ]",
		"[|F|W|Y|L . . ]",
		"[|F|W|Y . . ]",
		"[|R . . ]",
		"[|M . . ]",
		"[|W . . ]",
		"[|P . . ]",
		"[|E . . ]",
		"[|F|L|W|Y|A|E|Q . . ]",
		"[|R A|C|D|E|F|G|H|I|K|L|M|N|T|V|W|Y . ]",
		"[|E|D . . ]",
		"[|K . . ]",
		"[ . . D|]",
		"[|A|L|I|V A|C|D|E|F|G|H|I|K|L|M|N|T|V|W|Y . ]",
		"[|A|L|I|V|K|R|W|F|Y A|C|D|E|F|G|H|I|K|L|M|N|T|V|W|Y . ]",
		"[|K|R|L|F|W|Y|N|Q|C|H|S|M . . ]"
	};

	struct RunTimeConfig : public BaseRunTimeConfig
	{
	public:
		RTCONFIG_DEFINE_MEMBERS( RunTimeConfig, SEQUEST_RUNTIME_CONFIG, "\r\n\t ", "sequest.params", "\r\n;" )

		DynamicModSet dynamicMods;
		StaticModSet staticMods;

	private:

		void finalize()
		{
			static const boost::char_separator<char> delim(" ");
			tokenizer parser( diff_search_options.begin(), diff_search_options.begin() + diff_search_options.length(), delim );
			tokenizer::iterator itr = parser.begin();

			size_t numMotifs = 0;
			while( itr != parser.end() )
			{
				if( numMotifs > strlen( userModChars ) )
					throw runtime_error( "too many mods specified in diff_search_options" );

				float modMass = lexical_cast<float>( *itr );
				string motif = '[' + *(++itr) + ']';
				char userModChar = userModChars[numMotifs];
				++itr;
				if( modMass == 0 && motif == "[X]" )
					continue;
				++numMotifs;
				//cout << motif << " " << userModChar << " " << modMass << endl;
				dynamicMods.parseMotif( motif, userModChar, modMass );
			}

			vector<string> terminalModStrings;
			split( terminalModStrings, term_diff_search_options, boost::is_space() );

			vector<float> terminalMods( terminalModStrings.size() );
			std::transform( terminalModStrings.begin(), terminalModStrings.end(), terminalMods.begin(), lexical_cast<float, string> );

			if( terminalMods[1] > 0.0f )
				dynamicMods.parseMotif( PEPTIDE_N_TERMINUS_STRING, ']', terminalMods[1] );

			if( terminalMods[0] > 0.0f )
				dynamicMods.parseMotif( PEPTIDE_C_TERMINUS_STRING, '[', terminalMods[0] );

			if( add_Nterm_peptide != 0 )	staticMods.insert( StaticMod( PEPTIDE_N_TERMINUS_SYMBOL, add_Nterm_peptide ) );
			if( add_Cterm_peptide != 0 )	staticMods.insert( StaticMod( PEPTIDE_C_TERMINUS_SYMBOL, add_Cterm_peptide ) );
			if( add_G_Glycine != 0 )		staticMods.insert( StaticMod( 'G', add_G_Glycine ) );
			if( add_A_Alanine != 0 )		staticMods.insert( StaticMod( 'A', add_A_Alanine ) );
			if( add_S_Serine != 0 )			staticMods.insert( StaticMod( 'S', add_S_Serine ) );
			if( add_P_Proline != 0 )		staticMods.insert( StaticMod( 'P', add_P_Proline ) );
			if( add_V_Valine != 0 )			staticMods.insert( StaticMod( 'V', add_V_Valine ) );
			if( add_T_Threonine != 0 )		staticMods.insert( StaticMod( 'T', add_T_Threonine ) );
			if( add_C_Cysteine != 0 )		staticMods.insert( StaticMod( 'C', add_C_Cysteine ) );
			if( add_L_Leucine != 0 )		staticMods.insert( StaticMod( 'L', add_L_Leucine ) );
			if( add_I_Isoleucine != 0 )		staticMods.insert( StaticMod( 'I', add_I_Isoleucine ) );
			if( add_X_LorI != 0 )			staticMods.insert( StaticMod( 'X', add_X_LorI ) );
			if( add_N_Asparagine != 0 )		staticMods.insert( StaticMod( 'N', add_N_Asparagine ) );
			if( add_O_Ornithine != 0 )		staticMods.insert( StaticMod( 'O', add_O_Ornithine ) );
			if( add_B_avg_NandD != 0 )		staticMods.insert( StaticMod( 'B', add_B_avg_NandD ) );
			if( add_D_Aspartic_Acid != 0 )	staticMods.insert( StaticMod( 'D', add_D_Aspartic_Acid ) );
			if( add_Q_Glutamine != 0 )		staticMods.insert( StaticMod( 'Q', add_Q_Glutamine ) );
			if( add_K_Lysine != 0 )			staticMods.insert( StaticMod( 'K', add_K_Lysine ) );
			if( add_Z_avg_QandE != 0 )		staticMods.insert( StaticMod( 'Z', add_Z_avg_QandE ) );
			if( add_E_Glutamic_Acid != 0 )	staticMods.insert( StaticMod( 'E', add_E_Glutamic_Acid ) );
			if( add_M_Methionine != 0 )		staticMods.insert( StaticMod( 'M', add_M_Methionine ) );
			if( add_H_Histidine != 0 )		staticMods.insert( StaticMod( 'H', add_H_Histidine ) );
			if( add_F_Phenylalanine != 0 )	staticMods.insert( StaticMod( 'F', add_F_Phenylalanine ) );
			if( add_R_Arginine != 0 )		staticMods.insert( StaticMod( 'R', add_R_Arginine ) );
			if( add_Y_Tyrosine != 0 )		staticMods.insert( StaticMod( 'Y', add_Y_Tyrosine ) );
			if( add_W_Tryptophan != 0 )		staticMods.insert( StaticMod( 'W', add_W_Tryptophan ) );
		}
	};
}
}

#endif
