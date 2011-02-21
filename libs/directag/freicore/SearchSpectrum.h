#ifndef _SQTFILE_H
#define _SQTFILE_H

#include "shared_defs.h"
#include "shared_funcs.h"
#include "BaseSpectrum.h"
#include "Profiler.h"
#include "searchResult.h"
#include "pepXmlTypes.h"
#include "SimpleXMLWriter.h"

using namespace freicore;

namespace freicore
{
	template< class SearchResultT >
	struct SearchSpectrum : public virtual BaseSpectrum
	{
		typedef SearchResultT						SearchResultType;
		typedef SearchResultSet< SearchResultT >	SearchResultSetType;

		SearchSpectrum() : BaseSpectrum() {}

		SearchSpectrum( const SearchSpectrum& old )
			:	BaseSpectrum( old ), resultSet( old.resultSet ), decoyState(0)
		{}

		virtual void	ScoreSequenceVsSpectrum(	SearchResultT& result,
													const string& seq,
													const vector<float>& seqIons,
													int NumIntensityClasses,
													float FragmentMzTolerance )
		{
		}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & resultSet;
		}

		SearchResultSetType	resultSet;
		char decoyState;
		char numTerminiCleavages; // 0, 1, or 2 termini
	};

	template< class SpectrumType = SearchSpectrum< GenericSearchResult > >
	struct SearchSpectraListSortByTotalScore
	{
		bool operator() ( const SpectrumType* lhs, const SpectrumType* rhs )
		{
			float lhsScore = lhs->resultSet.getBestTotalScore();
			float rhsScore = rhs->resultSet.getBestTotalScore();

			if( lhsScore == rhsScore )
			{
				return spectraSortByID()( lhs, rhs );
			}

			return lhsScore > rhsScore;
		}
	};

	template< class SpectrumType, class SpectraListType >
	struct SearchSpectraList : public virtual BaseSpectraList< SpectrumType, SpectraListType >
	{
		bool spectraDecoyStatesSet;

		SearchSpectraList()
			:	BaseSpectraList< SpectrumType, SpectraListType >(), spectraDecoyStatesSet(false)
		{}

//		typedef SearchSpectraList< SpectrumType, ListType >			ListType;
		typedef BaseSpectraList< SpectrumType, SpectraListType >	SearchBaseList;
		typedef typename SearchBaseList::ListConstIterator			ListConstIterator;
		typedef typename SearchBaseList::ListIterator				ListIterator;

		void filterByChargeStateAndTerminiCleavages(	int chargeState,
														int terminiCleavageCount,
														SpectraListType* passingSpectra = NULL,
														SpectraListType* failingSpectra = NULL )
		{
			vector< int > chargeStates( 1, chargeState );
			vector< int > terminiCleavageCounts( 1, terminiCleavageCount );
			filterByChargeStateAndTerminiCleavages( chargeStates, terminiCleavageCounts, passingSpectra, failingSpectra );
		}

		void filterByChargeStateAndTerminiCleavages(	const vector< int >& chargeStates,
														const vector< int >& terminiCleavageCounts,
														SpectraListType* passingSpectra = NULL,
														SpectraListType* failingSpectra = NULL )
		{
			set<int> passByChargeState;
			for( size_t i=0; i < chargeStates.size(); ++i )
				passByChargeState.insert( chargeStates[i] );

			set<int> passByTerminiCleavages;
			for( size_t i=0; i < terminiCleavageCounts.size(); ++i )
				passByTerminiCleavages.insert( terminiCleavageCounts[i] );

			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				if( passByChargeState.find( s->id.charge ) != passByChargeState.end() &&
					passByTerminiCleavages.find( s->numTerminiCleavages ) != passByTerminiCleavages.end() )
				{
					if( passingSpectra )
						passingSpectra->push_back( s );
				} else if( failingSpectra )
					failingSpectra->push_back( s );
			}
		}

		string readSQT(	const string& filepath,
						bool skipUnindexedSpectra = false,
						bool useNameInScanIds = false,
						size_t maxResultRank = 1,
						const string& proteinNameDelimiter = " " )
		{
			RunTimeVariableMap dummy;
			return readSQT( filepath, skipUnindexedSpectra, useNameInScanIds, maxResultRank, proteinNameDelimiter, dummy );
		}

		string readSQT(	const string& filepath,
						bool skipUnindexedSpectra,
						bool useNameInScanIds,
						size_t maxResultRank,
						const string& proteinNameDelimiter,
						RunTimeVariableMap& parameters )
		{
			ifstream fileStream( filepath.c_str() );
			if( !fileStream.is_open() )
				throw invalid_argument( string( "unable to open SQT file \"" ) + filepath + "\"" );

			Timer readTime(true);
			size_t fileSize = (size_t) GetFileSize( filepath );
			string fileStr;
			fileStr.resize( fileSize );
			fileStream.read( &fileStr[0], (streamsize) fileSize );
			fileStream.close();
			//cout << g_hostString << " finished reading " << fileSize << " bytes; " << readTime.End() << " seconds elapsed." << endl;

			return parseSQT( fileStr, filepath, skipUnindexedSpectra, useNameInScanIds, maxResultRank, proteinNameDelimiter, parameters );
		}

		string parseSQT(	const string& fileStr,
							const string& filepath,
							bool skipUnindexedSpectra,
							bool useNameInScanIds,
							size_t maxResultRank,
							const string& proteinNameDelimiter,
							RunTimeVariableMap& parameters )
		{
			size_t headerStart = 0;
			size_t headerEnd = fileStr.find( "S\t", headerStart );

			if( headerEnd == string::npos )
				throw runtime_error( string( "no S lines found in SQT file \"" ) + filepath + "\"" );

			string header = fileStr.substr( headerStart, headerEnd - headerStart );

			ResidueMap fileResidueMap( *g_residueMap );
			string realProteinNameDelimiter = proteinNameDelimiter;
			if( proteinNameDelimiter.find('\t') == string::npos )
				realProteinNameDelimiter += '\t';
			if( proteinNameDelimiter.find('\r') == string::npos )
				realProteinNameDelimiter += '\r';
			if( proteinNameDelimiter.find('\n') == string::npos )
				realProteinNameDelimiter += '\n';
			static const string SQTParametersToken = "SQTParameters";
			static const string inputFileToken = "InputFile";
			size_t startIdx, endIdx = 0;

			string sourceFilepath;
			if( ( startIdx = fileStr.find( inputFileToken, headerStart ) ) != string::npos )
			{
				startIdx += inputFileToken.length() + 1;
				endIdx = fileStr.find_first_of( "\r\n", startIdx );
				sourceFilepath = fileStr.substr( startIdx, endIdx - startIdx );
			}
			if( sourceFilepath.empty() )
				sourceFilepath = filepath.substr( 0, filepath.find_last_of('.') ) + ".xml";

			int numParams = 0;
			if( ( startIdx = fileStr.find( SQTParametersToken, headerStart ) ) != string::npos )
			{
				startIdx += SQTParametersToken.length();
				endIdx = fileStr.find_first_of( "\r\n", startIdx );
				numParams = atoi( fileStr.substr( startIdx, endIdx - startIdx ).c_str() );
			}

			for( int i=0; i < numParams; ++i )
			{
				endIdx = fileStr.find( ":", startIdx );
				startIdx = fileStr.find_last_of( ", \t", endIdx ) + 1;
				string paramName = fileStr.substr( startIdx, endIdx - startIdx );
				startIdx = fileStr.find_first_not_of( " \t", endIdx + 1 );
				endIdx = fileStr.find_first_of( ",\r\n", startIdx );
				map< string, string >::iterator param = parameters.find( paramName );
				if( param != parameters.end() )
					param->second = UnquoteString( fileStr.substr( startIdx, endIdx - startIdx ) );
			}

			if( parameters.count( "DynamicMods" ) )
				fileResidueMap.setDynamicMods( parameters["DynamicMods"] );
			if( parameters.count( "StaticMods" ) )
				fileResidueMap.setStaticMods( parameters["StaticMods"] );
			//cout << parameters << endl;

			string scanName;
			if( true )//useNameInScanIds )
			{
				scanName = GetFilenameFromFilepath( sourceFilepath );
				scanName = scanName.substr( 0, scanName.find_last_of('.') );
			}

			size_t tokenStart = headerEnd, tokenEnd;
			size_t fileSize = fileStr.length();
			SpectrumType* s;

			while( tokenStart < fileSize && fileStr[tokenStart] == 'S' )
			{
				//cout << "\"" << (*itr) << "\" ";

				tokenStart += 2; tokenEnd = fileStr.find( '\t', tokenStart+1 ); // skip S and \t
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
				int num = lexical_cast<int>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );

				tokenStart = fileStr.find( '\t', tokenEnd+1 )+1; tokenEnd = fileStr.find( '\t', tokenStart+1 ); // skip second scan number token and \t
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
				int charge = lexical_cast<int>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
				SpectrumId id( scanName, num, charge );
				if( !useNameInScanIds )
					id.setSource( "" );

				typename SpectraListType::SearchBaseList::ListIndexIterator indexItr = this->index.find( id );
				if( indexItr != this->index.end() )
				{
					s = *indexItr->second;
				} else if( skipUnindexedSpectra )
				{
					tokenStart = fileStr.find( "S\t", tokenEnd+13 ); // skip at least the next 13 characters
					continue;
				} else
				{
					s = new SpectrumType;
					s->id = id;
				}

				tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
				s->processingTime = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
				tokenStart = fileStr.find( '\t', tokenEnd+1 )+1; tokenEnd = fileStr.find( '\t', tokenStart+1 ); // skip hostname and \t
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
				s->mOfPrecursor = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) - PROTON;
				tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
				s->totalIonCurrent = exp( lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) );
				tokenStart = fileStr.find( '\t', tokenEnd+1 )+1; tokenEnd = fileStr.find_first_of( "\r\n", tokenStart+1 ); // skip score and \t
				//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << endl;
				s->numSequenceComparisons = lexical_cast<int>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );

				tokenStart = fileStr.find( '\n', tokenEnd ); // move to next line
				if( tokenStart >= fileSize )
					break;
				tokenStart += 1;

				while( fileStr[tokenStart] == 'M' )
				{
					tokenStart += 2; // skip M and \t

					typename SpectrumType::SearchResultType r;

					tokenEnd = fileStr.find( '\t', tokenStart+1 );
					size_t rank = lexical_cast<size_t>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
					if( maxResultRank > 0 && rank > maxResultRank )
					{
						tokenStart = fileStr.find( "S\t", tokenEnd+15 ); // skip at least the next 15 characters
						//cout << endl;
						break;
					}

					tokenStart = fileStr.find( '\t', tokenEnd+1 )+1; // skip rank2 and \t
					tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					r.mass = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) - PROTON;
					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					r.deltCn = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					r.score = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );

					if( r.score == 0 && r.deltCn == 0 )
					{
						tokenStart = fileStr.find( "S\t", tokenEnd+15 ); // skip at least the next 15 characters
						//cout << endl;
						break;
					}

					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					r.mod = lexical_cast<float>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ); // using sp to store modification mass
					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					int matched = lexical_cast<int>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << " ";
					int predicted = lexical_cast<int>( fileStr.substr( tokenStart, tokenEnd-tokenStart ) );
					tokenStart = tokenEnd+1; tokenEnd = fileStr.find( '\t', tokenStart+1 );
					//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << endl;
					r.sequence = fileStr.substr( tokenStart, tokenEnd-tokenStart );

					r.sequence = r.sequence.substr( 2, r.sequence.length() - 4 ); // trim flanking residue notation
					r.sequence = ConvertSqtPtmToFreiPtm( r.sequence, &fileResidueMap );
					//cout << r.sequence << endl;
					r.key = MvIntKey( matched, predicted );
					tokenStart += 2; tokenStart = fileStr.find( '\n', tokenStart )+1; // skip state and move to next line

					while( fileStr[tokenStart] == 'L' )
					{
						tokenStart += 2; // skip L and \t
						tokenEnd = fileStr.find_first_of( realProteinNameDelimiter, tokenStart+1 );
						//cout << QuoteString( fileStr.substr( tokenStart, tokenEnd-tokenStart ) ) << endl;
						r.lociByName.insert( ProteinLocusByName( fileStr.substr( tokenStart, tokenEnd-tokenStart ), 0 ) );

						tokenStart = fileStr.find( '\n', tokenEnd ); // move to next line
						if( tokenStart >= fileSize )
							break;
						tokenStart += 1;
					}
					s->resultSet.insert(r);
					if( tokenStart >= fileSize )
						break;
				}

				//s->resultSet.calculateRelativeScores();

				if( indexItr == this->index.end() )
				{
					this->push_back(s);
				}
			}

			return header;
		}

		void writeSQT(	const string& sourceFilepath,
						const string& filenameSuffix = "",
						const string& header = "" ) const
		{
			RunTimeVariableMap dummy;
			return writeSQT( sourceFilepath, filenameSuffix, header, dummy );
		}

		void writeSQT(	const string& sourceFilepath,
						const string& filenameSuffix,
						const string& header,
						const RunTimeVariableMap& vars ) const
		{
			string scanName = sourceFilepath.substr( sourceFilepath.find_last_of( SYS_PATH_SEPARATOR )+1,
													sourceFilepath.find_last_of( '.' ) - sourceFilepath.find_last_of( SYS_PATH_SEPARATOR )-1 );
			string filename = scanName + filenameSuffix + ".sqt";

			ofstream fileStream( filename.c_str() );
			if( !fileStream.is_open() )
				throw invalid_argument( string( "unable to write SQT file \"" ) + filename + "\"" );

			fileStream << header;

			fileStream << "H\tInputFile\t" << sourceFilepath << '\n';

			fileStream << showpoint << boolalpha;
			int n = 0;
			fileStream << "H\tSQTParameters\t" << vars.size();
			for( RunTimeVariableMap::const_iterator vItr = vars.begin(); vItr != vars.end(); ++vItr, ++n )
			{
				if( !(n % 4) )
					fileStream << "\nH\t";
				else
					fileStream << ", ";

				fileStream << vItr->first << ": " << vItr->second;
			}
			fileStream << "\n\n" << noshowpoint << noboolalpha;

			if( SearchBaseList::empty() )
				return;

			string hostname = GetHostname();
			hostname = hostname.substr( 0, hostname.find_first_of('.') );

			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				fileStream	<< "S\t"
							<< s->id.index << '\t'
							<< s->id.index << '\t'
							<< s->id.charge << '\t'
							<< s->processingTime << '\t'
							<< hostname << '\t'
							<< s->mOfPrecursor + PROTON << '\t'
							<< log( s->totalIonCurrent ) << '\t'
							<< 0 << '\t'
							<< s->numSequenceComparisons << '\n';

				s->resultSet.calculateRelativeScores();

				int n = 1;
				float lastScore = 0;
				for( typename SpectrumType::SearchResultSetType::const_reverse_iterator rItr = s->resultSet.rbegin(); rItr != s->resultSet.rend(); ++rItr )
				{
					if( rItr->score < lastScore )
						++n;
					lastScore = rItr->score;

					int fragmentsPredicted = accumulate( rItr->key.begin(), rItr->key.end(), 0 );
					int fragmentsFound = fragmentsPredicted - rItr->key.back();

					fileStream	<< "M\t"
								<< n << '\t'
								<< n << '\t'
								<< round( rItr->mass + PROTON, 4 )<< '\t'
								<< round( rItr->deltCn, 4 ) << '\t'
								<< round( rItr->score, 4 ) << '\t'
								<< round( rItr->mod, 2 ) << '\t'
								<< fragmentsFound << '\t'
								<< fragmentsPredicted << '\t'
								<< "-." << ConvertFreiPtmToSqtPtm( rItr->sequence ) << ".-" << "\tU\n";

					for( ProteinLociByName::const_iterator lItr = rItr->lociByName.begin(); lItr != rItr->lociByName.end(); ++lItr )
						fileStream << "L\t" << lItr->name << "\t" << lItr->offset << '\n';
				}
			}

			fileStream.close();
		}

		void readPepXml(	const string& filepath,
							size_t maxResultRank,
							RunTimeVariableMap& parameters )
		{
			PepXmlReader< SpectrumType, SpectraListType > pepXmlReader( reinterpret_cast< SpectraListType* >( this ) );
			pepXmlReader.open( filepath );
			parameters = pepXmlReader.ReadSpectra( maxResultRank );
			pepXmlReader.close();
		}

		void writePepXml(	const string& sourceFilepath,
							const string& filenameSuffix,
							const string& searchEngine,
							const string& searchDatabase,
							const RunTimeVariableMap& vars ) const
		{
			string scanName = basename( MAKE_PATH_FOR_BOOST(sourceFilepath) );
			string filename = scanName + filenameSuffix + ".pepXML";

			ResidueMap fileResidueMap;

			ofstream xmlStream( filename.c_str() );
			if( !xmlStream.is_open() )
				throw invalid_argument( string( "unable to write pepXML file \"" ) + filename + "\"" );

			SimpleXMLWriter xmlWriter;
			xmlWriter.condenseAttr_ = true;
			xmlWriter.setOutputStream( xmlStream );
			xmlWriter.startDocument();

			xmlWriter.open( "msms_pipeline_analysis" );
			xmlWriter.attr( "name", scanName );
			xmlWriter.attr( "date", GetDateTime() );
			xmlWriter.attr( "summary_xml", filename );
			xmlWriter.noattr();

			//xmlStream << xmlWriter.getIndentStr() << "<msms_run_index_offset file_offset=\""; // 2^64 is at most 20 decimal characters: 18446744073709551616

			//ofstream::pos_type msmsRunSummaryFileOffset = xmlWriter.getCurFileLength();

			xmlWriter.open( "msms_run_summary" );
			xmlWriter.attr( "base_name", scanName );
			xmlWriter.attr( "raw_data_type", "unknown" );
			xmlWriter.attr( "raw_data", "unknown" );

			vector<string> cleavageRuleNames;
			if( vars.count("Config: CleavageRuleNames") )
				split( cleavageRuleNames, vars.find("Config: CleavageRuleNames")->second, boost::is_space() );
			else
				cleavageRuleNames.resize( 1, "trypsin" );

			stringstream cleavageRulesStream;
			if( vars.count("Config: CleavageRules") )
				cleavageRulesStream.str( vars.find("Config: CleavageRules")->second );
			else
				cleavageRulesStream.str( "K|R A|C|D|E|F|G|H|I|K|L|M|N|T|V|W|Y" );

			CleavageRuleSet cleavageRules;
			cleavageRulesStream >> cleavageRules;
			//cout << cleavageRules << endl;
			size_t numCleavageRules = 0;

			for( size_t i=0; i < cleavageRules.size(); ++i )
			{
				CleavageHalfRule& n_terminal_residues = cleavageRules[i].first;
				CleavageHalfRule& c_terminal_residues = cleavageRules[i].second;
				if( n_terminal_residues.longestCleavageCandidate > 1 ||
					c_terminal_residues.longestCleavageCandidate > 1 )
				{
					cerr << "Warning: motif-style cleavage rule \"" << cleavageRules[i] << "\" not supported by pepXML!" << endl;
					//continue;
				}

				string cut_residues;
				for( CleavageHalfRule::iterator itr = n_terminal_residues.begin(); itr != n_terminal_residues.end(); ++itr )
					if( *itr != PROTEIN_N_TERMINUS_STRING && *itr != PROTEIN_C_TERMINUS_STRING )
						cut_residues += *itr->rbegin(); // copy only the most C terminal residue of each half rule

				string no_cut_residues;
				const set<char>& allResidues = fileResidueMap.getResidues();
				for( set<char>::const_iterator itr = allResidues.begin(); itr != allResidues.end(); ++itr )
					if( c_terminal_residues.find( lexical_cast<string>(*itr) ) != c_terminal_residues.end() )
						no_cut_residues += *itr;

				if( cut_residues.empty() && no_cut_residues.empty() )
					continue;

				++numCleavageRules;
				string cleavageRuleName = ( cleavageRuleNames.size() > i ? cleavageRuleNames[i] : "unknown" );

				xmlWriter.open( "sample_enzyme" );
				xmlWriter.attr( "name", cleavageRuleName );
				{
					xmlWriter.open( "specificity" );
					xmlWriter.attr( "cut", cut_residues );
					xmlWriter.attr( "no_cut", no_cut_residues );
					xmlWriter.attr( "sense", 'C' );
					xmlWriter.close(); // specificity
				}
				xmlWriter.close(); // sample_enzyme
			}

			xmlWriter.open( "search_summary" );
			xmlWriter.attr( "base_name", scanName );
			xmlWriter.attr( "search_engine", searchEngine );
			string precursor_mass_type = ( vars.count("Config: UseAvgMassOfSequences") && vars.find("Config: UseAvgMassOfSequences")->second == "1" ? "average" : "monoisotopic" );
			xmlWriter.attr( "precursor_mass_type", precursor_mass_type );
			xmlWriter.attr( "fragment_mass_type", "monoisotopic" );
			xmlWriter.attr( "out_data_type", "n/a" );
			xmlWriter.attr( "out_data", "n/a" );
			xmlWriter.attr( "search_id", 1 );
			{
				xmlWriter.open( "search_database" );
				xmlWriter.attr( "local_path", searchDatabase );
				xmlWriter.attr( "type", "AA" );
				xmlWriter.close(); // search_database

				for( size_t i=0; i < numCleavageRules; ++i )
				{
					string cleavageRuleName = ( cleavageRuleNames.size() > i ? cleavageRuleNames[i] : "unknown" );
					xmlWriter.open( "enzymatic_search_constraint" );
					xmlWriter.attr( "enzyme", cleavageRuleName );
					string max_num_internal_cleavages = ( vars.count("Config: NumMaxMissedCleavages") ? vars.find("Config: NumMaxMissedCleavages")->second : "10" );
					xmlWriter.attr( "max_num_internal_cleavages", max_num_internal_cleavages );
					string min_number_termini = ( vars.count("Config: NumMinTerminiCleavages") ? vars.find("Config: NumMinTerminiCleavages")->second : "2" );
					xmlWriter.attr( "min_number_termini", min_number_termini );
					xmlWriter.close(); // enzymatic_search_constraint
				}

				if( vars.count("Config: DynamicMods") )
				{
					DynamicModSet mods( vars.find("Config: DynamicMods")->second );
					fileResidueMap.setDynamicMods( vars.find("Config: DynamicMods")->second );
					//xmlStream << showpos;
					for( DynamicModSet::iterator itr = mods.begin(); itr != mods.end(); ++itr )
						switch( itr->unmodChar )
						{
							case PEPTIDE_N_TERMINUS_SYMBOL:
								xmlWriter.open( "terminal_modification" );
								xmlWriter.attr( "terminus", 'n' );
								xmlWriter.attr( "massdiff", itr->modMass );
								xmlWriter.attr( "mass", round( HYDROGEN_MONO + itr->modMass, 4 ) );
								xmlWriter.attr( "variable", 'Y' );
								xmlWriter.attr( "protein_terminus", "" );
								xmlWriter.attr( "symbol", itr->userModChar );
								xmlWriter.close(); // terminal_modification
								break;
							case PEPTIDE_C_TERMINUS_SYMBOL:
								xmlWriter.open( "terminal_modification" );
								xmlWriter.attr( "terminus", 'c' );
								xmlWriter.attr( "massdiff", itr->modMass );
								xmlWriter.attr( "mass", round( HYDROGEN_MONO + OXYGEN_MONO + itr->modMass, 4 ) );
								xmlWriter.attr( "variable", 'Y' );
								xmlWriter.attr( "protein_terminus", "" );
								xmlWriter.attr( "symbol", itr->userModChar );
								xmlWriter.close(); // terminal_modification
								break;
							default:
								xmlWriter.open( "aminoacid_modification" );
								xmlWriter.attr( "aminoacid", itr->unmodChar );
								xmlWriter.attr( "massdiff", itr->modMass );
								xmlWriter.attr( "mass", round( fileResidueMap.getMonoMassByName( itr->unmodChar ) + itr->modMass, 4 ) );
								xmlWriter.attr( "variable", 'Y' );
								xmlWriter.attr( "symbol", itr->userModChar );
								xmlWriter.close(); // aminoacid_modification
								break;
						}
					//xmlStream << noshowpos;
				}

				if( vars.count("Config: StaticMods") )
				{
					StaticModSet mods( vars.find("Config: StaticMods")->second );
					fileResidueMap.setStaticMods( vars.find("Config: StaticMods")->second );
					//xmlStream << showpos;
					for( StaticModSet::iterator itr = mods.begin(); itr != mods.end(); ++itr )
					{
						xmlWriter.open( "aminoacid_modification" );
						xmlWriter.attr( "aminoacid", itr->name );
						xmlWriter.attr( "massdiff", itr->mass );
						xmlWriter.attr( "mass", round( fileResidueMap.getMonoMassByName( itr->name ), 4 ) );
						xmlWriter.attr( "variable", 'N' );
						xmlWriter.close(); // aminoacid_modification
					}
					//xmlStream << noshowpos;
				}

				//xmlStream << showpoint << boolalpha;
				for( RunTimeVariableMap::const_iterator vItr = vars.begin(); vItr != vars.end(); ++vItr )
				{
					xmlWriter.open( "parameter" );
					xmlWriter.attr( "name", vItr->first );
					xmlWriter.attr( "value", vItr->second );
					xmlWriter.close(); // parameter
				}
				//xmlStream << noshowpoint << noboolalpha;
			}
			xmlWriter.close(); // search_summary

			map< SpectrumId, ofstream::pos_type > scanFileOffsets;

			size_t spectrumIndex = 0;
			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				//xmlStream << "\t\t";
				//xmlStream.flush();
				//scanFileOffsets[ s->id ] = xmlStream.tellp(); // store element offset

				xmlWriter.open( "spectrum_query" );
                stringstream spectrumQueryId;
                spectrumQueryId << "ID: \"" << s->stringID << "\"; Index: " << s->id.index;
                if( !s->nativeID.empty() )
                    spectrumQueryId << "; NativeID: \"" << s->nativeID << "\"";
                xmlWriter.attr( "spectrum", spectrumQueryId.str() );
				xmlWriter.attr( "start_scan", s->id.index );
				xmlWriter.attr( "end_scan", s->id.index );
				xmlWriter.attr( "precursor_neutral_mass", s->mOfPrecursor );
				xmlWriter.attr( "assumed_charge", s->id.charge );
				xmlWriter.attr( "index", ++spectrumIndex );
				xmlWriter.attr( "retention_time_sec", s->retentionTime * 60.0f );

				//s->resultSet.calculateRelativeScores();
				s->resultSet.calculateRanks();

				if( !s->resultSet.empty() )
				{
					xmlWriter.open( "search_result" );
					xmlWriter.attr( "num_comparisons", s->numSequenceComparisons );

					for( typename SpectrumType::SearchResultSetType::const_reverse_iterator rItr = s->resultSet.rbegin(); rItr != s->resultSet.rend(); ++rItr )
					{
						int fragmentsPredicted = accumulate( rItr->key.begin(), rItr->key.end(), 0 );
						int fragmentsFound = fragmentsPredicted - rItr->key.back();

						xmlWriter.open( "search_hit" );
						xmlWriter.attr( "hit_rank", rItr->rank );
						xmlWriter.attr( "peptide", GetUnmodifiedSequence( GetRawSequence( rItr->sequence ), &fileResidueMap ) );
						xmlWriter.attr( "peptide_prev_aa", '-' );
						xmlWriter.attr( "peptide_next_aa", '-' );
						xmlWriter.attr( "protein", rItr->lociByName.begin()->name );
						xmlWriter.attr( "peptide_offset", rItr->lociByName.begin()->offset );
						xmlWriter.attr( "num_tot_proteins", rItr->lociByName.size() );
						xmlWriter.attr( "num_matched_ions", fragmentsFound );
						xmlWriter.attr( "tot_num_ions", fragmentsPredicted );
						xmlWriter.attr( "calc_neutral_pep_mass", round( rItr->mass, 4 ) );
						xmlWriter.attr( "massdiff", round( s->mOfPrecursor - rItr->mass, 4 ) );
						xmlWriter.attr( "num_tol_term", rItr->numTerminiCleavages );
						xmlWriter.attr( "num_missed_cleavages", rItr->numMissedCleavages );
						//xmlWriter.attr( "is_rejected", 0 );
						{
							ModificationInfo modInfo( rItr->sequence, fileResidueMap );
							if( !modInfo.mods.empty() )
							{
								xmlWriter.open( "modification_info" );
								if( modInfo.mods.front().position == 0 )
								{
									xmlWriter.attr( "mod_nterm_mass", modInfo.mods.front().mass );
									modInfo.mods.erase( modInfo.mods.begin() );
								}
								if( !modInfo.mods.empty() && modInfo.mods.back().position == rItr->sequence.length()-1 )
								{
									xmlWriter.attr( "mod_cterm_mass", modInfo.mods.back().mass );
									modInfo.mods.pop_back();
								}
								for( vector< ModificationInstance >::iterator modItr = modInfo.mods.begin(); modItr != modInfo.mods.end(); ++modItr )
								{
									xmlWriter.open( "mod_aminoacid_mass" );
									xmlWriter.attr( "position", modItr->position );
									xmlWriter.attr( "mass", modItr->mass );
									xmlWriter.close(); // mod_aminoacid_mass
								}
								xmlWriter.close(); // modification_info
							}

							ProteinLociByName::const_iterator lItr = rItr->lociByName.begin();
							++lItr;
							for( ; lItr != rItr->lociByName.end(); ++lItr )
							{
								xmlWriter.open( "alternative_protein" );
								xmlWriter.attr( "protein", lItr->name );
								xmlWriter.attr( "peptide_offset", lItr->offset );
								xmlWriter.close(); // alternative_protein
							}

							SearchScoreList scoreList = rItr->getScoreList();
							for( SearchScoreList::const_iterator scoreItr = scoreList.begin(); scoreItr != scoreList.end(); ++scoreItr )
							{
								xmlWriter.open( "search_score" );
								xmlWriter.attr( "name", scoreItr->first );
								xmlWriter.attr( "value", round( scoreItr->second, 4 ) );
								xmlWriter.close(); // search_score
							}
						}
						xmlWriter.close(); // search_hit
					}
					xmlWriter.close(); // search_result
				}
				xmlWriter.close(); // spectrum_query
			}
			xmlWriter.close(); // msms_run_summary

			//xmlStream << "\t";
			/*xmlStream.flush();
			ofstream::pos_type indexFileOffset = xmlStream.tellp();
			xmlWriter.open( "msms_run_index" );
			{
				xmlWriter.open( "msms_run_summary_offset" );
				xmlWriter.attr( "base_name", scanName );
				xmlWriter.attr( "file_offset", msmsRunSummaryFileOffset );
				xmlWriter.close(); // msms_run_summary_offset

				for( map< SpectrumId, ofstream::pos_type >::iterator itr = scanFileOffsets.begin(); itr != scanFileOffsets.end(); ++itr )
				{
					xmlWriter.open( "spectrum_query_offset" );
					stringstream spectrumId;
					spectrumId << itr->first.source << '.' << itr->first.index << '.' << itr->first.index << '.' << itr->first.charge;
					xmlWriter.attr( "spectrum", spectrumId.str() );
					xmlWriter.attr( "file_offset", itr->second );
					xmlWriter.close(); // spectrum_query_offset
				}
			}
			xmlWriter.close();*/ // msms_run_index
			
			xmlWriter.close(); // msms_pipeline_analysis

			xmlStream.close();
		}

		vector< int >	calculateValidationThresholds(	vector< float >& scoreThresholds,
														int numChargeStates,
														float confidence,
														float realToDecoyRatio,
														const string& decoyPrefix,
														const int mode )
		{
			scoreThresholds.resize( numChargeStates, 0.0f );

			vector< int > avgComparisons;
			for( int z=0; z < numChargeStates; ++z )
			{
				float comparisonCount = 0;
				int chargeCount = 0;
				for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
				{
					if( (*sItr)->resultSet.empty() || (*sItr)->id.charge != z+1 )
						continue;

					comparisonCount += (*sItr)->numSequenceComparisons;
					++chargeCount;
				}

				avgComparisons.push_back( (int) round( comparisonCount / ( chargeCount > 0 ? chargeCount : 1 ), 0 ) );
			}

			// For every match with a deltCn of 0, consider its score to be a threshold
			// Let #R be the number of real loci above the threshold
			// Let #D be the number of decoy loci above the threshold
			// confidence is (#R - (#D * ratio of real to decoy entries in database)) / (#R + #D)

			// Sort the SQT entries in descending order by the chosen score
			SearchBaseList::sort( SearchSpectraListSortByTotalScore< SpectrumType >() );

			int numAmbiIds;

			map< SpectrumId, char > spectraDecoyStates;

			for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				if( (*sItr)->resultSet.empty() )
				{
					spectraDecoyStates[ (*sItr)->id ] = '-';
					continue;
				}

				char state;
				if( (*sItr)->resultSet.rbegin()->lociByName.begin()->name.find( decoyPrefix ) == 0 )
					state = 'D'; // initial locus is a decoy
				else
					state = 'R'; // initial locus is real

				for(	typename SpectrumType::SearchResultSetType::reverse_iterator resultItr = (*sItr)->resultSet.rbegin();
						state != 'B' && resultItr != (*sItr)->resultSet.rend();
						++resultItr )
				{
					if( resultItr->deltCn > 0.0f )
						break;

					for(	ProteinLociByName::iterator locusItr = resultItr->lociByName.begin();
							state != 'B' && locusItr != resultItr->lociByName.end();
							++locusItr )
					{
						if( locusItr->name.find( decoyPrefix ) == 0 ) // found a decoy locus
						{
							if( state == 'R' ) // previous loci were real and now one is a decoy, so state is mixed
								state = 'B';

						} else // found a real locus
						{
							if( state == 'D' ) // previous loci were decoys and now one is real, so state is mixed
								state = 'B';
						}
					}
				}

				spectraDecoyStates[ (*sItr)->id ] = state;
			}

			ofstream debugOut( "debug-out.txt" );
			for( int z=0; z < numChargeStates; ++z )
			{
				for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
				{
					if( (*sItr)->id.charge != z+1 )
						continue;

					debugOut << (*sItr)->id.id << "\t" <<
								(*sItr)->resultSet.getBestScore(mode) << "\t" <<
								spectraDecoyStates[ (*sItr)->id ] << "\n";
				}
			}
			debugOut.close();


			for( int z=0; z < numChargeStates; ++z )
			{
				float lastThreshold = 0;
				float curThreshold = 0;
				float curConfidence = 1.0f;

				int numRealIds = 0;
				int numDecoyIds = 0;
				numAmbiIds = 0;

				vector< pair< ListIterator, pair< float, float > > > allConfidences;

				for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
				{
					if( (*sItr)->resultSet.empty() || (*sItr)->id.charge != z+1 )
						continue;

					lastThreshold = curThreshold;
					curThreshold = (*sItr)->resultSet.getBestScore(mode);

					char state = spectraDecoyStates[ (*sItr)->id ];
					if( state == 'D' )
						++numDecoyIds;
					else if( state == 'R' )
						++numRealIds;
					else
						++numAmbiIds;

					if( numRealIds + numDecoyIds > 0 )
						curConfidence = max( float( numRealIds - ( numDecoyIds * realToDecoyRatio ) ) /
											 float( numRealIds + numDecoyIds ), 0.0f );
					else
						curConfidence = 0;

					//cout << curThreshold << " " << numRealIds << " " << numDecoyIds << " " << curConfidence << endl;
					scoreThresholds[z] = curThreshold - 0.00001f;//(curThreshold + lastThreshold) / 2.0f;
					allConfidences.push_back( pair< ListIterator, pair< float, float > >( sItr, pair< float, float >( curConfidence, scoreThresholds[z] ) ) );
				}

				size_t i;
				for( i = allConfidences.size(); i > 0; --i )
				{
					ListIterator sItr = allConfidences[i-1].first;
					char state = spectraDecoyStates[ (*sItr)->id ];
					if( allConfidences[i-1].second.first >= confidence && state == 'R' )
						break;
				}

				if( allConfidences.empty() || ( i == 0 && allConfidences[i].second.first < confidence ) )
				{
					scoreThresholds[z] = -1;
					cerr << "Warning: negative threshold for +" << z+1 << "s indicates not enough reverse matches for validation (check the decoy prefix)" << endl;
				} else
				{
					scoreThresholds[z] = allConfidences[i-1].second.second;
				}
			}

			return avgComparisons;
		}

		pair< SpectraListType, SpectraListType > filterByThresholds(	vector< float > scoreThresholds,
																		int numChargeStates, char mode,
																		vector< size_t >* potentialMatchCounts,
																		vector< size_t >* validMatchCounts )
		{
			SpectraListType passingSpectra;
			SpectraListType failingSpectra;

			potentialMatchCounts->resize( numChargeStates, 0 );
			validMatchCounts->resize( numChargeStates, 0 );

			for( int z=0; z < numChargeStates; ++z )
			{
				for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
				{
					if( (*sItr)->resultSet.empty() || (*sItr)->id.charge != z+1 )
						continue;

					++ potentialMatchCounts->at(z);
					if( scoreThresholds[z] < 0 )
						continue;

					float score = (*sItr)->resultSet.getBestScore( mode );
					if( score >= scoreThresholds[z] )
					{
						++ validMatchCounts->at(z);
						passingSpectra.push_back( (*sItr) );
					} else
						failingSpectra.push_back( (*sItr ) );
				}
			}

			return pair< SpectraListType, SpectraListType >( passingSpectra, failingSpectra );
		}


		void calculateFDRs(	int numChargeStates,
							float realToDecoyRatio,
							const string& decoyPrefix,
                            ostream* pQonversionDetailsStream = NULL )
		{
			if( SearchBaseList::empty() )
				return;

			if( !spectraDecoyStatesSet )
			{
				for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
				{
					SpectrumType* s = *sItr;

					if( s->resultSet.empty() )
					{
						s->decoyState = '-';
						continue;
					}

					s->resultSet.calculateRanks();

					char state;
					if( s->resultSet.rbegin()->lociByName.begin()->name.find( decoyPrefix ) == 0 )
						state = 'D'; // initial locus is a decoy
					else
						state = 'R'; // initial locus is real

					for(	typename SpectrumType::SearchResultSetType::reverse_iterator resultItr = s->resultSet.rbegin();
							state != 'B' && resultItr != s->resultSet.rend();
							++resultItr )
					{
						if( resultItr->rank > 1 )
							break;

						for(	ProteinLociByName::iterator locusItr = resultItr->lociByName.begin();
								state != 'B' && locusItr != resultItr->lociByName.end();
								++locusItr )
						{
							if( locusItr->name.find( decoyPrefix ) == 0 ) // found a decoy locus
							{
								if( state == 'R' ) // previous loci were real and now one is a decoy, so state is mixed
									state = 'B';

							} else // found a real locus
							{
								if( state == 'D' ) // previous loci were decoys and now one is real, so state is mixed
									state = 'B';
							}
						}
					}

					s->decoyState = state;
				}

				spectraDecoyStatesSet = true;
			}

			START_PROFILER(8);
			SearchBaseList::sort( SearchSpectraListSortByTotalScore< SpectrumType >() );
			STOP_PROFILER(8);
			for( int z=0; z < numChargeStates; ++z )
			{
				int numRealIds = 0;
				int numDecoyIds = 0;
				int numAmbiIds = 0;
				int rank = 1;

				for( ListIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr, ++rank )
				{
					SpectrumType* s = *sItr;

					if( s->id.charge != z+1 )
						continue;

					if( s->resultSet.empty() )
					{
						if( pQonversionDetailsStream )
							(*pQonversionDetailsStream) << rank << '\t' << s->id.index << '\t' << s->id.charge <<
										"\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\n";
						continue;
					}

					if( s->decoyState == 'D' )
						++numDecoyIds;
					else if( s->decoyState == 'R' )
						++numRealIds;
					else
						++numAmbiIds;

					if( pQonversionDetailsStream )
						(*pQonversionDetailsStream) <<
                            rank << '\t' << s->id.index << '\t' << s->id.charge << '\t' <<
	                        s->decoyState << '\t' << numRealIds << '\t' << numDecoyIds << '\t' <<
	                        numAmbiIds << '\t' << s->resultSet.rbegin()->getScoreList() << '\t' <<
	                        s->resultSet.rbegin()->getTotalScore();

					for(	typename SpectrumType::SearchResultSetType::reverse_iterator resultItr = s->resultSet.rbegin();
							resultItr != s->resultSet.rend() && resultItr->rank == 1;
							++resultItr )
					{
						if( numRealIds + numDecoyIds > 0 )
						{
							const_cast< typename SpectrumType::SearchResultType& >( *resultItr ).fdr =
									1- max(	float( numRealIds - ( numDecoyIds * realToDecoyRatio ) ) /
											float( numRealIds + numDecoyIds ), 0.0f );
							if( pQonversionDetailsStream && resultItr == s->resultSet.rbegin() )
								(*pQonversionDetailsStream) << '\t' << resultItr->fdr;

						} else
							const_cast< typename SpectrumType::SearchResultType& >( *resultItr ).fdr = 0;
					}

					if( pQonversionDetailsStream )
						(*pQonversionDetailsStream) << '\n';

				}
			}
		}

		void filterByFDR(	float maxFDR = 1,
							SpectraListType* passingSpectra = NULL,
							SpectraListType* failingSpectra = NULL )
		{
			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				if( !s->resultSet.empty() && s->resultSet.rbegin()->fdr <= maxFDR )
				{
					if( passingSpectra )
						passingSpectra->push_back( s );
				} else if( failingSpectra )
					failingSpectra->push_back( s );
			}
		}

		size_t getPassingCountByFDR( float maxFDR )
		{
			size_t passingSpectraCount = 0;
			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				if( !s->resultSet.empty() && s->resultSet.rbegin()->fdr <= maxFDR )
					++passingSpectraCount;
			}
			return passingSpectraCount;
		}

        float getScoreThresholdByFDR( float maxFDR )
        {
            float worstScore = std::numeric_limits<float>::max();
			for( ListConstIterator sItr = SearchBaseList::begin(); sItr != SearchBaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				if( !s->resultSet.empty() && s->resultSet.rbegin()->fdr <= maxFDR )
					if( s->resultSet.rbegin()->getTotalScore() < worstScore )
                        worstScore = s->resultSet.rbegin()->getTotalScore();
			}
			return worstScore;
        }
	};
}

//BOOST_CLASS_IMPLEMENTATION( freicore::SearchSpectrum, boost::serialization::object_serializable );
//BOOST_CLASS_TRACKING( freicore::SearchSpectrum, boost::serialization::track_never )

#endif
