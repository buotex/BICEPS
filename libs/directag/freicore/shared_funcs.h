#ifndef _SHARED_FUNCS_H
#define _SHARED_FUNCS_H

#include "stdafx.h"
#include "shared_types.h"
#include "shared_defs.h"
#include "lnFactorialTable.h"
#include "ResidueMap.h"
#include "BaseRunTimeConfig.h"
#include "BaseSpectrum.h"

using namespace freicore;

namespace std
{
	ostream&		operator<< ( ostream& o, const ProteinLocusByIndex& rhs );
	ostream&		operator<< ( ostream& o, const ProteinLocusByName& rhs );

	
	ostream&		operator<< ( ostream& o, const CleavageRule& rhs );
	istream&		operator>> ( istream& i, CleavageRule& rhs );
	istream&		operator>> ( istream& i, CleavageRuleSet& rhs );

	template< class T >
	ostream&		operator<< ( ostream& o, const topset<T>& rhs )
	{
		return o << reinterpret_cast< const set<T>& >( rhs );
	}

	ostream& operator<< ( ostream& o, MvIntKey& rhs );
	ostream& operator<< ( ostream& o, MvhTable& rhs );
}

namespace freicore
{
	double			lnCombin( int n, int k, lnFactorialTable& lnTable = g_lnFactorialTable );
	float			lnOdds( float p );
	//endianType_t	GetMachineEndianType();

	void			MakePtmVariants( const CandidateSequenceInfo& candidate, CandidateSequenceList& ptmSequences, int maxPtmCount, ResidueMap* residueMap = g_residueMap );
	string			GetUnmodifiedSequence( const string& ptmSeq, DynamicModSet& dynamicMods );
	string			GetUnmodifiedSequence( const string& ptmSeq, ResidueMap* residueMap = g_residueMap );
	string			GetRawSequence( const string& ptmSeq );
	string			ConvertFreiPtmToSqtPtm( const string& ptmSeq, ResidueMap* residueMap = g_residueMap );
	string			ConvertSqtPtmToFreiPtm( const string& ptmSeq, ResidueMap* residueMap = g_residueMap );
	string			ConvertSqtPtmToFreiPtm( const string& ptmSeq, DynamicModSet& dynamicMods );

	void			CommitCommonDatatypes();

	float			ChiSquaredToPValue( float x, int df );

	int				ClassifyError( float error, vector< float > mzFidelityThresholds );

	void			LoadTagInstancesFromIndexFile( const string& tagIndexFilename, const string& tag, tagMetaIndex_t& tagMetaIndex, tagIndex_t& tagIndex );
	void			LoadIndexFromFile( const string& tagIndexFilename, tagMetaIndex_t& tagMetaIndex );

	const static bool defaultIonTypes[4] = { true, true, false, false };
	void			CalculateSequenceIons(	const string& seq,
											int maxIonCharge,
											vector< float >* sequenceIonMasses,
											bool useSmartPlusThreeModel = true,
											vector< string >* sequenceIonLabels = NULL,
											float precursorMass = 0.0f,
											const bool ionTypes[4] = defaultIonTypes,
											ResidueMap* residueMap = g_residueMap );

	void			CalculateSequenceIons2(	const string& seq,
											int maxIonCharge,
											vector< float >* sequenceIonMasses,
											bool useSmartPlusThreeModel = true,
											vector< string >* sequenceIonLabels = NULL,
											float precursorMass = 0.0f,
											const bool ionTypes[4] = defaultIonTypes,
											ResidueMap* residueMap = g_residueMap );

	void			CreateScoringTableMVH(	const int minValue,
											const int totalValue,
											const int numClasses,
											vector< int > classCounts,
											MvhTable& mvTable,
											lnFactorialTable& lnFT,
											bool NormalizeOnMode = false,
											bool adjustRareOutcomes = true,
											bool convertToPValues = false,
											const MvIntKey* minKey = NULL );

	void			CreateScoringTableMVB(	const int minValue,
											const int totalValue,
											const int numClasses,
											const vector< float >& classProbabilities,
											MvhTable& mvTable,
											lnFactorialTable& lnFT,
											bool NormalizeOnMode = false,
											bool adjustRareOutcomes = true );

	int				paramIndex( const string& param, const char** atts, int attsCount );

	void			FindFilesByMask( const string& mask, fileList_t& filenames );
}

#endif
