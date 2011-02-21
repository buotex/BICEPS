#ifndef _PEPXMLTYPES_H
#define _PEPXMLTYPES_H

#include "shared_defs.h"
#include "shared_funcs.h"
#include "searchResult.h"
#include <limits>
#include "expat_xml.h"

#define HAS_ATTR(name) (paramIndex(name,atts,attsCount) > -1)
#define GET_ATTR_AS(name, type) getAttributeAs<type>(name,atts,attsCount)
#define GET_ATTR(name) GET_ATTR_AS(name,std::string)

namespace freicore
{
	template< class T >
	T getAttributeAs( const string& name, const char** atts, int attsCount )
	{
		if( !HAS_ATTR(name) )
			throw out_of_range( "required attribute \"" + name + "\" not found" );
		return lexical_cast<T>( atts[paramIndex(name, atts, attsCount)+1] );
	}

	struct GenericSearchResult : public BaseSearchResult
	{
		GenericSearchResult( const float s = 0.0f, const MvIntKey& k = MvIntKey(), const string& seq = "" )
			:	BaseSearchResult( k, seq )
		{}

		SearchScoreList scoreList;
		float fdr;

		float getTotalScore() const
		{
			return scoreList[0].second;
		}

		SearchScoreList getScoreList() const
		{
			return scoreList;
		}

		bool operator< ( const GenericSearchResult& rhs ) const
		{
			if( getTotalScore() == rhs.getTotalScore() )
				if( mod == rhs.mod )
					return GetUnmodifiedSequence( sequence ) < GetUnmodifiedSequence( rhs.sequence );
				else
					return mod > rhs.mod;
			else
				return getTotalScore() < rhs.getTotalScore();
		}

		bool operator== ( const GenericSearchResult& rhs ) const
		{
			return ( getTotalScore() == rhs.getTotalScore() && sequence == rhs.sequence );
		}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & boost::serialization::base_object< BaseSearchResult >( *this );
			ar & scoreList & fdr;
		}
	};

	//typedef GenericSearchResult SearchResult;

	static const float MOD_MASS_EPSILON = 0.01f;

	/*struct EnzymeSpecificityInfo
	{
		char sense;
		int min_spacing;
		string cut;
		string no_cut;
	};

	struct EnzymeInfo
	{
		EnzymeInfo( const string& CleavageRules, int NumMinTerminiCleavages )
		{
			stringstream ruleStream(CleavageRules);
			CleavageRuleSet rules;
			ruleStream >> rules;

			name = "unknown";
			desc = "";
			switch( NumMinTerminiCleavages )
			{
				case 0:
					fidelity = "nonspecific";
					break;
				case 1:
					fidelity = "semispecific";
					break;
				default:
				case 2:
					fidelity = "specific";
					break;
			}
			independent = true;
		}

		string name;
		string desc;
		string fidelity;
		bool independent;
		vector< EnzymeSpecificityInfo > specifityList;
	};*/

	struct ModificationInstance
	{
		ModificationInstance( size_t position, float mass, float modMass )
			: position(position), mass(mass), modMass(modMass)
		{}

		size_t position;
		float mass;
		float modMass;
	};

	struct ModificationInfo
	{
		ModificationInfo( const string& modified_sequence, const ResidueMap& residueMap )
		{
			for( size_t i=0; i < modified_sequence.length(); ++i )
			{
				for( DynamicModSet::const_iterator itr = residueMap.dynamicMods.begin(); itr != residueMap.dynamicMods.end(); ++itr )
					if( itr->uniqueModChar == modified_sequence[i] )
						mods.push_back( ModificationInstance( i, residueMap.getMonoMassByName( itr->uniqueModChar ), itr->modMass ) );

				for( StaticModSet::const_iterator itr = residueMap.staticMods.begin(); itr != residueMap.staticMods.end(); ++itr )
					if( itr->name == modified_sequence[i] )
						mods.push_back( ModificationInstance( i, residueMap.getMonoMassByName( itr->name ), itr->mass ) );
			}
		}

		vector< ModificationInstance > mods;
	};

	struct ScoreInfo
	{
		ScoreInfo( bool higherIsBetter = true, const string& range = "[0,infinity)" ) : higherIsBetter( higherIsBetter )
		{
			string varRange = TrimWhitespace(range);
			minValueInclusive = ( *varRange.begin() == '[' );
			maxValueInclusive = ( *varRange.rbegin() == ']' );
			varRange = varRange.substr( 1, varRange.length()-2 );

			vector<string> values;
			split( values, varRange, boost::is_any_of(",") );

			if( TrimWhitespace( values[0] ) == "-infinity" )
				minValue = -std::numeric_limits<float>::infinity();
			else
				minValue = lexical_cast<float>( values[0] );

			if( TrimWhitespace( values[1] ) == "infinity" )
				maxValue = std::numeric_limits<float>::infinity();
			else
				maxValue = lexical_cast<float>( values[1] );
		}

		int testScore( float score )
		{
			if( ( minValueInclusive && score < minValue ) || ( !minValueInclusive && score <= minValue ) )
				return -1;
			if( ( maxValueInclusive && score > maxValue ) || ( !maxValueInclusive && score >= maxValue ) )
				return 1;
			return 0;
		}

		bool	higherIsBetter;
		float	minValue;
		bool	minValueInclusive;
		float	maxValue;
		bool	maxValueInclusive;
	};

	template< class SpectrumType, class SpectraListType >
	class PepXmlReader
	{
	public:

		typedef typename SpectraListType::iterator				SpectraListIterator;
		typedef PepXmlReader< SpectrumType, SpectraListType >	ReaderType;

		PepXmlReader( SpectraListType* pSpectra, const string& filename = "" )
			: m_pSpectra( pSpectra ), m_pParser( NULL ), m_maxResultRank(0)
		{
			/*m_absoluteScores["xcorr"]		= ScoreInfo( true, "[0,infinity)" );
			m_absoluteScores["mvh"]			= ScoreInfo( true, "[0,infinity)" );
			m_absoluteScores["hyperscore"]	= ScoreInfo( true, "[0,infinity)" );
			m_relativeScores["deltacn"]		= ScoreInfo( true, "[0,1]" );
			m_relativeScores["expect"]		= ScoreInfo( false, "(0,infinity)" );*/

			if( !filename.empty() )
				open( filename );
		}

		~PepXmlReader()
		{
			close();
		}

		void open( const string& filename )
		{
			close();

			m_vars["NumChargeStates"] = "3";
			m_vars["StaticMods"] = "";
			m_vars["DynamicMods"] = "";
			m_vars["ProteinDatabase"] = "";
			m_vars["CleavageRules"] = "K|R|[ . . ]";
			m_vars["UseAvgMassOfSequences"] = "1";
			m_vars["NumMaxMissedCleavages"] = "10";
			m_vars["NumMinTerminiCleavages"] = "2";

			m_curResultRank = 0;

			m_InputFileName = filename;
			m_InputFile.open( m_InputFileName.c_str(), std::ios::binary );
			if( !m_InputFile.is_open() )
				throw invalid_argument( string( "unable to open pepXML file \"" ) + filename + "\"" );

			m_InputFile.clear();
			m_pParser = XML_ParserCreate( NULL );
			XML_SetUserData( m_pParser, this );

			XML_SetElementHandler( m_pParser, StartElement, EndElement );
		}

		void close()
		{
			m_InputFile.close();

			if( m_pParser )
			{
				XML_ParserFree( m_pParser );
				m_pParser = NULL;
			}
		}

		int ReadSpectra( int nCount, int maxTotalPeakCount )
		{
				unsigned int done, bytesRead;
				char* buf = new char[READ_BUFFER_SIZE];

				do
				{
					m_InputFile.read( buf, READ_BUFFER_SIZE );
					bytesRead = m_InputFile.gcount();
					done = bytesRead < sizeof(buf);

					try
					{
						if( !XML_Parse( m_pParser, buf, bytesRead, done ) )
						{
							throw runtime_error( XML_ErrorString( XML_GetErrorCode( m_pParser ) ) );
						}
					} catch( exception& e )
					{
						throw runtime_error( string( e.what() ) + " at line " + lexical_cast<string>( XML_GetCurrentLineNumber( m_pParser ) ) );
					}

				} while( !done );

				delete buf;

				return m_nCount;
		}

		RunTimeVariableMap ReadSpectra( int maxResultRank = 1 )
		{
			m_maxResultRank = maxResultRank;
			ReadSpectra( 0, 0 );
			return m_vars;
		}

	private:
		ifstream						m_InputFile;
		SpectraListType					m_overflowList;

		SpectraListType*				m_pSpectra;
		SpectrumType*					m_pSpectrum;
		SearchResultSet<GenericSearchResult>	m_resultSet;
		GenericSearchResult				m_result;

		XML_Parser						m_pParser;

		int								m_nCount;
		string							m_InputFileName;
		string							m_scanName;
		int								m_maxResultRank;
		int								m_curResultRank;

		map<string,ScoreInfo>			m_absoluteScores;
		map<string,ScoreInfo>			m_relativeScores;

		RunTimeVariableMap				m_vars;
		ResidueMap						m_residueMap;

		static void StartElement( void *userData, const char *name, const char **atts )
		{
			ReaderType* pInstance = static_cast< ReaderType* >( userData );
			SpectrumType*& s = pInstance->m_pSpectrum;
			GenericSearchResult& result = pInstance->m_result;
			ResidueMap& residueMap = pInstance->m_residueMap;
			RunTimeVariableMap& vars = pInstance->m_vars;

			string tag(name);
			int attsCount = XML_GetSpecifiedAttributeCount( pInstance->m_pParser );

			try
			{
				if( tag == "search_hit" )
				{
					pInstance->m_curResultRank = GET_ATTR_AS("hit_rank", int);
					if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
						return;

					result.lociByName.clear();
					result.scoreList.clear();
					result.rank = pInstance->m_curResultRank;
					result.sequence = string( PEPTIDE_N_TERMINUS_STRING ) + GET_ATTR("peptide") + PEPTIDE_C_TERMINUS_STRING;
					result.lociByName.insert( ProteinLocusByName( GET_ATTR("protein") ) );
					result.key.resize(2);
					result.key[0] = HAS_ATTR("num_matched_ions") ? GET_ATTR_AS("num_matched_ions", int) : 0;
					result.key[1] = HAS_ATTR("tot_num_ions") ? GET_ATTR_AS("tot_num_ions", int) : result.key[0];
					result.mass = GET_ATTR_AS("calc_neutral_pep_mass", float);
					result.numTerminiCleavages = HAS_ATTR("num_tol_term") ? GET_ATTR_AS("num_tol_term", int) : 2;
					result.numMissedCleavages = HAS_ATTR("num_missed_cleavages") ? GET_ATTR_AS("num_missed_cleavages", int) : 0;

				} else if( tag == "search_score" )
				{
					if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
						return;

					string name = to_lower_copy( GET_ATTR("name") );
					try
					{
						result.scoreList.push_back( SearchScoreInfo( name, GET_ATTR_AS("value", float) ) );
					} catch( bad_lexical_cast& e )
					{
						e;
						// ignore non-numeric scores
					} catch( exception& e )
					{
						e;
						// ignore scores without values
					}

				} else if( tag == "alternative_protein" )
				{
					if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
						return;

					result.lociByName.insert( ProteinLocusByName( GET_ATTR("protein") ) );

				} else if( tag == "mod_aminoacid_mass" )
				{
					if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
						return;

					float mass = GET_ATTR_AS("mass", float);
					size_t position = GET_ATTR_AS("position", size_t);
					char unmodChar = result.sequence[position];
					float lowestMassError = residueMap.beginMonoMasses()->first;
					for( DynamicModSet::iterator itr = residueMap.dynamicMods.begin(); itr != residueMap.dynamicMods.end(); ++itr )
						if( itr->unmodChar == unmodChar )
						{
							float moddedMass = residueMap.getMonoMassByName(itr->uniqueModChar);
							float massError = fabs( moddedMass - mass );
							if( massError < lowestMassError )
							{
								result.sequence[position] = itr->uniqueModChar;
								lowestMassError = massError;
							}
						}

				} else if( tag == "modification_info" )
				{
					if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
						return;

					int attrPos;

					if( ( attrPos = paramIndex("mod_nterm_mass", atts, attsCount)+1 ) > 0 )
					{
						float nMass = lexical_cast<float>( atts[attrPos] );
						float lowestMassError = residueMap.beginMonoMasses()->first;
						for( DynamicModSet::iterator itr = residueMap.dynamicMods.begin(); itr != residueMap.dynamicMods.end(); ++itr )
							if( itr->unmodChar == PEPTIDE_N_TERMINUS_SYMBOL )
							{
								float moddedMass = residueMap.getMonoMassByName(itr->uniqueModChar);
								float massError = fabs( moddedMass - nMass );
								if( massError < lowestMassError )
								{
									result.sequence[0] = itr->uniqueModChar;
									lowestMassError = massError;
								}
							}
					}

					if( ( attrPos = paramIndex("mod_cterm_mass", atts, attsCount)+1 ) > 0 )
					{
						float cMass = lexical_cast<float>( atts[attrPos] );
						float lowestMassError = residueMap.beginMonoMasses()->first;
						for( DynamicModSet::iterator itr = residueMap.dynamicMods.begin(); itr != residueMap.dynamicMods.end(); ++itr )
							if( itr->unmodChar == PEPTIDE_C_TERMINUS_SYMBOL )
							{
								float moddedMass = residueMap.getMonoMassByName(itr->uniqueModChar);
								float massError = fabs( moddedMass - cMass );
								if( massError < lowestMassError )
								{
									*result.sequence.rbegin() = itr->uniqueModChar;
									lowestMassError = massError;
								}
							}
					}

				} else if( tag == "spectrum_query" )
				{
					// Read the spectrum scanNum
					int scan = GET_ATTR_AS("start_scan", int);
					int charge = GET_ATTR_AS("assumed_charge", int);

					// Create a new spectrum
					s = new SpectrumType;
					s->id.set( pInstance->m_scanName, scan, charge );

					s->mOfPrecursor = GET_ATTR_AS("precursor_neutral_mass", float);
					s->retentionTime = HAS_ATTR("retention_time_sec") ? ( GET_ATTR_AS("retention_time_sec", float) / 60.0f ) : 0;

					// Set the spectrum filename and scanName
					s->fileName = pInstance->m_InputFileName;

				} else if( tag == "search_result" )
				{
					s->numSequenceComparisons = HAS_ATTR("num_comparisons") ? GET_ATTR_AS("num_comparisons", int) : 0;

				} else if( tag == "msms_run_summary" )
				{
					pInstance->m_scanName = path( GET_ATTR("base_name") ).leaf();

				} else if( tag == "sample_enzyme" )
				{
					// not supported yet

				} else if( tag == "search_summary" )
				{
					bool useAvgMass = GET_ATTR("precursor_mass_type") == "average";
					vars["UseAvgMassOfSequences"] = useAvgMass ? "1" : "0";

					vars["SearchEngine: Name"] = GET_ATTR("search_engine");
					vars["SearchEngine: Version"] = "unknown";

				} else if( tag == "search_database" )
				{
					vars["ProteinDatabase"] = GET_ATTR("local_path");

				} else if( tag == "enzymatic_search_constraint" )
				{
					vars["NumMaxMissedCleavages"] = GET_ATTR("max_num_internal_cleavages");
					vars["NumMinTerminiCleavages"] = GET_ATTR("min_number_termini");

				} else if( tag == "aminoacid_modification" )
				{
					bool isDynamic = ( GET_ATTR("variable")[0] == 'Y' );
					if( isDynamic )
					{
						string unmodChar = GET_ATTR("aminoacid");
						string modMass = GET_ATTR("massdiff");
						string symbol = HAS_ATTR("symbol") ? GET_ATTR("symbol") : "*";
						vars["DynamicMods"] += unmodChar + " " + symbol + " " + modMass + " ";
					} else
					{
						string unmodChar = GET_ATTR("aminoacid");
						string modMass = GET_ATTR("massdiff");
						vars["StaticMods"] += unmodChar + " " + modMass + " ";
					}

				} else if( tag == "terminal_modification" )
				{
					bool isDynamic = ( GET_ATTR("variable")[0] == 'Y' );
					if( isDynamic )
					{
						string terminus = GET_ATTR("terminus");
						const char* unmodChar = STR_EQUAL( terminus, "n" ) ? PEPTIDE_N_TERMINUS_STRING : PEPTIDE_C_TERMINUS_STRING;
						string modMass = GET_ATTR("massdiff");
						string symbol = HAS_ATTR("symbol") ? GET_ATTR("symbol") : "*";
						vars["DynamicMods"] += string( unmodChar ) + " " + symbol + " " + modMass + " ";
					} else
					{
						string terminus = GET_ATTR("terminus");
						const char* unmodChar = STR_EQUAL( terminus, "n" ) ? PEPTIDE_N_TERMINUS_STRING : PEPTIDE_C_TERMINUS_STRING;
						string modMass = GET_ATTR("massdiff");
						vars["StaticMods"] += string( unmodChar ) + " " + modMass + " ";
					}

				} else if( tag == "parameter" )
				{
					string varName = GET_ATTR("name");
					if( varName.find("Config: ") == 0 )
						varName = varName.substr(8);
					vars[varName] = GET_ATTR("value");
				}
			} catch( exception& e )
			{
				throw runtime_error( string( "error parsing element \"" ) + tag + "\": " + e.what() );
			}
		}

		static void EndElement( void *userData, const char *name )
		{
			ReaderType* pInstance = static_cast< ReaderType* >( userData );
			SpectrumType*& s = pInstance->m_pSpectrum;
			GenericSearchResult& result = pInstance->m_result;
			ResidueMap& residueMap = pInstance->m_residueMap;

			string tag(name);

			if( tag == "search_hit" )
			{
				if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
					return;
				//cout << result.sequence << '\n';

				s->resultSet.insert( result );

			} else if( tag == "modification_info" )
			{
				if( pInstance->m_maxResultRank != -1 && pInstance->m_curResultRank > pInstance->m_maxResultRank )
					return;

			} else if( tag == "search_result" )
			{
				//s->resultSet.calculateRelativeScores();

			} else if( tag == "spectrum_query" )
			{
				if( s->resultSet.empty() )
					s->numTerminiCleavages = 2;
				else
					s->numTerminiCleavages = s->resultSet.rbegin()->numTerminiCleavages;

				// Push current spectrum onto the list
				pInstance->m_pSpectra->push_back( s );
				++pInstance->m_nCount;

			} else if( tag == "search_summary" )
			{
				residueMap.setDynamicMods( pInstance->m_vars["DynamicMods"] );
				residueMap.setStaticMods( pInstance->m_vars["StaticMods"] );
				//residueMap.dump();
			}
		}

	};
}

#endif
