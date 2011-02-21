#ifndef _DIRECTAGCONFIG_H
#define _DIRECTAGCONFIG_H

#include "stdafx.h"
#include "freicore.h"
#include "BaseRunTimeConfig.h"

using namespace freicore;

#define TAG_ONLY_HITS	0
#define TAG_ALWAYS		1
#define TAG_ONLY_MISSES	2

#define DIRECTAG_RUNTIME_CONFIG \
	COMMON_RTCONFIG SPECTRUM_RTCONFIG SEQUENCE_RTCONFIG MULTITHREAD_RTCONFIG VALIDATION_RTCONFIG \
	RTCONFIG_VARIABLE( string,			OutputSuffix,				"-tags"			) \
	RTCONFIG_VARIABLE( string,			InlineValidationFile,		""				) \
	RTCONFIG_VARIABLE( int,				InlineValidationMode,		TAG_ONLY_HITS	) \
	RTCONFIG_VARIABLE( bool,			AlwaysKeepValidTags,		true			) \
	RTCONFIG_VARIABLE( int,				MaxTagCount,				20				) \
	RTCONFIG_VARIABLE( float,			MaxTagScore,				20.0f			) \
	RTCONFIG_VARIABLE( int,				NumIntensityClasses,		3				) \
	RTCONFIG_VARIABLE( int,				NumMzFidelityClasses,		3				) \
	RTCONFIG_VARIABLE( int,				TagLength,					5				) \
	RTCONFIG_VARIABLE( int,				Tags,                       20              ) \
	RTCONFIG_VARIABLE( int,				StartSpectraScanNum,		0				) \
	RTCONFIG_VARIABLE( int,				EndSpectraScanNum,			-1				) \
	RTCONFIG_VARIABLE( int,				NumBatches,					50				) \
	RTCONFIG_VARIABLE( float,			TicCutoffPercentage,		1.0f			) \
	RTCONFIG_VARIABLE( size_t,			MaxPeakCount,				100				) \
	RTCONFIG_VARIABLE( float,			ClassSizeMultiplier,		2.0f			) \
	RTCONFIG_VARIABLE( float,			MinPrecursorAdjustment,		-2.5f			) \
	RTCONFIG_VARIABLE( float,			MaxPrecursorAdjustment,		2.5f			) \
	RTCONFIG_VARIABLE( float,			PrecursorAdjustmentStep,	0.1f			) \
	RTCONFIG_VARIABLE( bool,			NormalizeOnMode,			true			) \
	RTCONFIG_VARIABLE( bool,			AdjustPrecursorMass,		false			) \
	RTCONFIG_VARIABLE( int,				DeisotopingMode,			1				) \
	RTCONFIG_VARIABLE( bool,			MakeSpectrumGraphs,			false			) \
	RTCONFIG_VARIABLE( int,				MzFidelityErrorBinsSize,	20				) \
	RTCONFIG_VARIABLE( int,				MzFidelityErrorBinsSamples,	10000		) \
	RTCONFIG_VARIABLE( float,			MzFidelityErrorBinsLogMin,	-5.0f			) \
	RTCONFIG_VARIABLE( float,			MzFidelityErrorBinsLogMax,	1.0f			) \
	\
	RTCONFIG_VARIABLE( float,			IntensityScoreWeight,		1.0f			) \
	RTCONFIG_VARIABLE( float,			MzFidelityScoreWeight,		1.0f			) \
	RTCONFIG_VARIABLE( float,			ComplementScoreWeight,		1.0f			) \
	RTCONFIG_VARIABLE( float,			ContextScoreWeight,			1.0f			) \
	RTCONFIG_VARIABLE( float,			SatelliteScoreWeight,		0.0f			) \
	RTCONFIG_VARIABLE( float,			IsotopeScoreWeight,			0.0f			) \
	RTCONFIG_VARIABLE( float,			CompositionScoreWeight,		0.0f			) \
	RTCONFIG_VARIABLE( float,			OnLongestPathScoreWeight,	0.0f			) \
	RTCONFIG_VARIABLE( float,			RandomScoreWeight,			0.0f			)

namespace freicore
{
namespace directag
{

	struct RunTimeConfig : public BaseRunTimeConfig
	{
	public:
		BOOST_PP_SEQ_FOR_EACH( RTCONFIG_DECLARE_VAR, ~, DIRECTAG_RUNTIME_CONFIG )

		RunTimeConfig() : BaseRunTimeConfig(), rng(0), RandomScoreRange( 0.0, 0.99999999999 ), GetRandomScore( rng, RandomScoreRange )
		{
			BOOST_PP_SEQ_FOR_EACH( RTCONFIG_INIT_DEFAULT_VAR, ~, DIRECTAG_RUNTIME_CONFIG )
		}

		void initializeFromBuffer( string& cfgStr, const string& delim = "\r\n\t " )
		{
			BaseRunTimeConfig::initializeFromBuffer( cfgStr, delim );
			string strVal;
			BOOST_PP_SEQ_FOR_EACH( RTCONFIG_PARSE_BUFFER, ~, DIRECTAG_RUNTIME_CONFIG )
			finalize();
		}

		RunTimeVariableMap getVariables( bool hideDefaultValues = false )
		{
			BaseRunTimeConfig::getVariables();
			BOOST_PP_SEQ_FOR_EACH( RTCONFIG_FILL_MAP, m_variables, DIRECTAG_RUNTIME_CONFIG )
			return m_variables;
		}

		void setVariables( RunTimeVariableMap& vars )
		{
			BaseRunTimeConfig::setVariables( vars );
			BOOST_PP_SEQ_FOR_EACH( RTCONFIG_READ_MAP, vars, DIRECTAG_RUNTIME_CONFIG )
			finalize();
		}

		int initializeFromFile( const string& rtConfigFilename = "directag.cfg" )
		{
			return BaseRunTimeConfig::initializeFromFile( rtConfigFilename );
		}

		vector< float >			scoreThresholds;

		ResidueMap				inlineValidationResidues;
		tagMetaIndex_t			tagMetaIndex;

		map< string, float >	compositionScoreMap;

		boost::mt19937 rng;
		boost::uniform_real<> RandomScoreRange;
		boost::variate_generator< boost::mt19937&, boost::uniform_real<> > GetRandomScore;

		int		SpectraBatchSize;
		int		ProteinBatchSize;
		int		minIntensityClassCount;
		int		minMzFidelityClassCount;
		int		tagPeakCount;
		float	MzFidelityErrorBinsScaling;
		float	MzFidelityErrorBinsOffset;

	protected:
		void finalize()
		{
			BaseRunTimeConfig::finalize();

			string cwd;
			cwd.resize( MAX_PATH );
			getcwd( &cwd[0], MAX_PATH );
			WorkingDirectory = cwd.c_str();

			if( TicCutoffPercentage > 1.0f )
			{
				TicCutoffPercentage /= 100.0f;
				if( g_pid == 0 )
					cerr << g_hostString << ": TicCutoffPercentage > 1.0 (100%) corrected, now at: " << TicCutoffPercentage << endl;
			}

			if( !DynamicMods.empty() )
			{
				DynamicMods = TrimWhitespace( DynamicMods );
				g_residueMap->setDynamicMods( DynamicMods );
			}

			if( !StaticMods.empty() )
			{
				StaticMods = TrimWhitespace( StaticMods );
				g_residueMap->setStaticMods( StaticMods );
			}

			//if( g_residueMap )
			//	inlineValidationResidues = *g_residueMap;

			tagPeakCount = TagLength + 1;

			float m = ClassSizeMultiplier;
			if( m > 1 )
			{
				minIntensityClassCount = int( ( pow( m, NumIntensityClasses ) - 1 ) / ( m - 1 ) );
				minMzFidelityClassCount = int( ( pow( m, NumMzFidelityClasses ) - 1 ) / ( m - 1 ) );
			} else
			{
				minIntensityClassCount = NumIntensityClasses;
				minMzFidelityClassCount = NumMzFidelityClasses;
			}
		}
	};

	extern RunTimeConfig*					g_rtConfig;
}
}

#endif
