#ifndef _TAGSFILE_H
#define _TAGSFILE_H

#include "stdafx.h"
#include "shared_defs.h"

namespace freicore
{
	struct TagInfo
	{
		TagInfo( bool sort = false ) : valid(false), lowPeakMz(0), worstPeakRank(0), ascendingScores( sort ) {}
		TagInfo( const map< string, float >& newScores, const string& tag = "", float mz = 0, float nT = 0, float cT = 0, bool sort = false )
			: tag(tag), valid(false), scores( newScores ), totalScore(0), lowPeakMz( mz ), worstPeakRank(0),
			nTerminusMass( nT ), cTerminusMass( cT ), ascendingScores( sort ) {}
		TagInfo( const string& tag, float nT, float cT, bool sort = false )
			: tag(tag), valid(false), totalScore(0), lowPeakMz(0), worstPeakRank(0),
			nTerminusMass( nT ), cTerminusMass( cT ), ascendingScores( sort ) {}

		//void CalculateTotal( const map< string, float >& scoreWeights )
		void CalculateTotal( float cWeight, float iWeight, float mWeight )
		{
			/*totalScore = 0;
			for( map< string, float >::const_iterator itr = scores.begin(); itr != scores.end(); ++itr )
				totalScore += itr->second * scoreWeights.find( itr->first )->second;
			scores[ "total" ] = totalScore;*/

			/*totalScore = 1;
			for( map< string, float >::const_iterator itr = scores.begin(); itr != scores.end(); ++itr )
				totalScore *= pow( itr->second, scoreWeights.find( itr->first )->second );
			float totalWeight = 0;
			for( map< string, float >::const_iterator itr = scoreWeights.begin(); itr != scoreWeights.end(); ++itr )
				totalWeight += itr->second;
			totalScore = pow( totalScore, 1.0f / totalWeight );
			scores[ "total" ] = totalScore;*/

			/*totalScore = 0;
			for( map< string, float >::const_iterator itr = scores.begin(); itr != scores.end(); ++itr )
				totalScore += log(itr->second) * scoreWeights.find( itr->first )->second;
			totalScore *= -2;
			int numScores = 0;
			for( map< string, float >::const_iterator itr = scoreWeights.begin(); itr != scoreWeights.end(); ++itr )
				if( itr->second > 0 )
					++numScores;
			totalScore = ChiSquaredToPValue( totalScore, numScores*2 );
			scores[ "total" ] = totalScore;*/

			totalScore = ( log(complementScore)*cWeight + log(intensityScore)*iWeight + log(mzFidelityScore)*mWeight ) * -2;
			totalScore = ChiSquaredToPValue( totalScore, (int) floor( cWeight + iWeight + mWeight )*2 );
			//scores[ "total" ] = totalScore;
		}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & tag & valid & scores & lowPeakMz & nTerminusMass & cTerminusMass & worstPeakRank;
			ar & complementScore & intensityScore & mzFidelityScore & totalScore;
			//totalScore = scores[ "total" ];
		}

		bool operator< ( const TagInfo& rhs ) const
		{
			if( totalScore == rhs.totalScore )
				if( tag == rhs.tag )
					if( lowPeakMz == rhs.lowPeakMz )
						if( nTerminusMass == rhs.nTerminusMass )
							return cTerminusMass < rhs.cTerminusMass;
						else
							return nTerminusMass < rhs.nTerminusMass;
					else
						return lowPeakMz < rhs.lowPeakMz;
				else
					return tag < rhs.tag;
			else if( ascendingScores )
				return totalScore < rhs.totalScore;
			else
				return totalScore > rhs.totalScore;
		}

		string					tag;
		bool					valid;
		map< string, float >	scores;
		float					complementScore;
		float					intensityScore;
		float					mzFidelityScore;
		float					totalScore;
		float					lowPeakMz;
		int						worstPeakRank;
		float					nTerminusMass;		// the mass of residues missing from the N terminus end of a tag
		float					cTerminusMass;		// the mass of residues missing from the C terminus end of a tag
		int						ranksum;

	private:
		bool					ascendingScores;
	};

	struct TagList : public topset< TagInfo >
	{
		bool ascendingScores;

		void tagExploder( const TagInfo& tag )
		{
			insert( tag, true );
			tagExploder_R( tag, 0 );
		}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & boost::serialization::base_object< topset< TagInfo > >( *this );
		}

	private:
		void tagExploder_R( const TagInfo& tag, size_t idx )
		{
			if( idx == tag.tag.length() )
			{
				insert( tag, true );
				return;
			}

			if( tag.tag[idx] == 'I' )
			{
				TagInfo newTag( tag );
				newTag.tag[idx] = 'L';
				tagExploder_R( newTag, idx+1 );
			}

			if( tag.tag[idx] == 'L' )
			{
				TagInfo newTag( tag );
				newTag.tag[idx] = 'I';
				tagExploder_R( newTag, idx+1 );
			}

			tagExploder_R( tag, idx+1 );
		}
	};

	struct TaggingSpectrum : public virtual BaseSpectrum
	{
		TaggingSpectrum()
			:	BaseSpectrum(), tagList(), tagCount(0)
		{}

		TaggingSpectrum( const TaggingSpectrum& old )
			:	BaseSpectrum( old ), tagList( old.tagList ), tagCount( old.tagCount )
		{}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & tagList & tagCount;
		}

		TagList			tagList;
		size_t			tagCount;
	};

	template< class SpectrumType, class SpectraListType >
	struct TaggingSpectraList : public virtual BaseSpectraList< SpectrumType, SpectraListType >
	{
		//typedef TaggingSpectraList< SpectrumType >				ListType;
		typedef BaseSpectraList< SpectrumType, SpectraListType >	TaggingBaseList;
		typedef typename TaggingBaseList::ListConstIterator			ListConstIterator;
		typedef typename TaggingBaseList::ListIterator				ListIterator;

		int trimByTagCount( size_t minTagCount = 1, bool deleteTrimmedSpectra = true )
		{
			size_t sizeBefore = TaggingBaseList::size();

			vector< ListIterator > trimmedItrs;
			for( ListIterator sItr = SpectraListType::begin(); sItr != SpectraListType::end(); ++sItr )
			{
				if( (*sItr)->tagList.size() < minTagCount )
					trimmedItrs.push_back( sItr );
			}

			for( size_t i=0; i < trimmedItrs.size(); ++i )
			{
				erase( trimmedItrs[i], deleteTrimmedSpectra );
			}

			return int( sizeBefore - TaggingBaseList::size() );
		}

		string readTags(	const string& filename,
							bool skipUnindexedSpectra = false )
		{
			RunTimeVariableMap dummy;
			return readTags( filename, skipUnindexedSpectra, dummy );
		}

		string readTags(	const string& filename,
							bool skipUnindexedSpectra,
							RunTimeVariableMap& parameters )
		{
			ifstream fileStream( filename.c_str() );
			if( !fileStream.is_open() )
				throw invalid_argument( string( "unable to open tags file \"" ) + filename + "\"" );

			size_t fileSize = (size_t) GetFileSize( filename );
			string fileStr;
			fileStr.resize( fileSize );
			fileStream.read( &fileStr[0], (streamsize) fileSize );
			fileStream.close();

			size_t headerStart = 0;
			size_t headerEnd = fileStr.find( "S\t", headerStart );

			string tagsGenerator;
			string sourceFilepath;
			bool highScoresAreBetter;
			string sourceName = GetFilenameWithoutExtension( GetFilenameFromFilepath( filename ) );

			static const string GutenTagToken = "GutenTag";
			if( fileStr.find( GutenTagToken ) == string::npos )
			{
				static const string generatorToken = "TagsGenerator";
				static const string parametersToken = "TagsParameters";
				static const string inputFileToken = "InputFile";
				size_t startIdx, endIdx = 0;

				if( ( startIdx = fileStr.find( generatorToken, headerStart ) ) != string::npos )
				{
					startIdx += generatorToken.length() + 1;
					endIdx = fileStr.find_first_of( "\r\n", startIdx );
					tagsGenerator = fileStr.substr( startIdx, endIdx - startIdx );
				}

				if( tagsGenerator == "DirecTag" || tagsGenerator == "FreiTag" )
					highScoresAreBetter = false;
				else if( tagsGenerator == "Inspect" || tagsGenerator == "PepNovo" )
					highScoresAreBetter = true;
				else
					highScoresAreBetter = true;

				if( ( startIdx = fileStr.find( inputFileToken, headerStart ) ) != string::npos )
				{
					startIdx += inputFileToken.length() + 1;
					endIdx = fileStr.find_first_of( "\r\n", startIdx );
					sourceFilepath = fileStr.substr( startIdx, endIdx - startIdx );
				}

				int numParams = 0;
				if( ( startIdx = fileStr.find( parametersToken, headerStart ) ) != string::npos )
				{
					startIdx += parametersToken.length();
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
						param->second = fileStr.substr( startIdx, endIdx - startIdx );
				}
				//cout << sqtInfo.sqtParameters << endl;

				static const string headerTagEntryDefinitionToken = "H\tTag\tlowPeakMz\tnTerminusMass\tcTerminusMass\tworstPeakRank\tValid\tTotal";
				size_t lineStart = fileStr.find( headerTagEntryDefinitionToken, headerStart );
				size_t lineEnd = 0;
				if( lineStart != string::npos )
				{
					lineStart += headerTagEntryDefinitionToken.length();
					lineEnd = min( fileSize, fileStr.find_first_of( "\r\n", lineStart ) );
				} else
					throw invalid_argument( string( "unable to find score names in tags file \"" ) + filename + "\"" );

				stringstream headerTagEntryDefinition( fileStr.substr( lineStart, lineEnd - lineStart ) );
				vector< string > tagEntryScoreNames;
				string scoreName;
				while( headerTagEntryDefinition >> scoreName )
					tagEntryScoreNames.push_back( scoreName );
				size_t numScores = tagEntryScoreNames.size();

				static const boost::char_separator<char> delim(" \t\r\n");
				tokenizer parser( fileStr.begin() + headerEnd, fileStr.begin() + fileStr.length(), delim );
				tokenizer::iterator itr = parser.begin();

				while( itr != parser.end() )
				{
					//cout << "\"" << (*itr) << "\" ";

					if( *itr == "S" )
					{
						sourceFilepath = filename;

						SpectrumType* s;

						int num = lexical_cast<int>( *(++itr) );
						int charge = lexical_cast<int>( *(++itr) );
						SpectrumId id( sourceName, num, charge );

						typename TaggingBaseList::ListIndexIterator indexItr = this->index.find( id );
						if( indexItr != this->index.end() )
						{
							s = *indexItr->second;
						} else if( skipUnindexedSpectra )
						{
							continue;
						} else
						{
							s = new SpectrumType;
							s->id = id;
						}

						s->mOfPrecursor = lexical_cast<float>( *(++itr) );
						s->mzLowerBound = lexical_cast<float>( *(++itr) );
						s->mzUpperBound = lexical_cast<float>( *(++itr) );
						s->peakPreCount = lexical_cast<int>( *(++itr) );
						s->peakCount = lexical_cast<int>( *(++itr) );
						s->tagCount = lexical_cast<int>( *(++itr) );
						++itr;

						while( itr != parser.end() && *(++itr) == "T" )
						{
							TagInfo t( highScoresAreBetter );
							t.tag = *(++itr);
							t.lowPeakMz = lexical_cast<float>( *(++itr) );
							t.nTerminusMass = lexical_cast<float>( *(++itr) );
							t.cTerminusMass = lexical_cast<float>( *(++itr) );
							t.worstPeakRank = lexical_cast<int>( *(++itr) );
							t.valid = lexical_cast<bool>( *(++itr) );
							t.totalScore = lexical_cast<float>( *(++itr) );

							for( size_t i=0; i < numScores; ++i )
								t.scores[ tagEntryScoreNames[i] ] = lexical_cast<float>( *(++itr) );

							s->tagList.insert(t);
						}

						if( indexItr == this->index.end() )
						{
							TaggingBaseList::push_back(s);
						}
					} else
						++itr;
				}
			} else
			{
				// Read GutenTag format
				tagsGenerator = "GutenTag";
				highScoresAreBetter = true;
				cout << "Detected GutenTag format, using alternate parser." << endl;
				static const boost::char_separator<char> delim(" \t\r\n");
				tokenizer parser( fileStr.begin() + headerEnd, fileStr.begin() + fileStr.length(), delim );
				tokenizer::iterator itr = parser.begin();

				while( itr != parser.end() )
				{
					//cout << "\"" << (*itr) << "\" ";

					if( *itr == "S" )
					{
						SpectrumType* s;

						int num = lexical_cast<int>( *(++itr) );
						++itr; // skip second scan number
						int charge = lexical_cast<int>( *(++itr) );
						SpectrumId id( sourceName, num, charge );

						typename TaggingBaseList::ListIndexIterator indexItr = this->index.find( id );
						if( indexItr != this->index.end() )
						{
							s = *indexItr->second;
						} else if( skipUnindexedSpectra )
						{
							continue;
						} else
						{
							s = new SpectrumType;
							s->id = id;
						}

						++itr; // skip time
						s->peakPreCount = lexical_cast<int>( *(++itr) );
						s->peakCount = lexical_cast<int>( *(++itr) );
						++itr; // skip original precursor mass
						s->mOfPrecursor = lexical_cast<float>( *(++itr) ) - PROTON;
						++itr; // skip TIC
						
						while( itr != parser.end() && *(++itr) == "T" )
						{
							TagInfo t( highScoresAreBetter );
							t.totalScore = lexical_cast<float>( *(++itr) );
							t.tag = *(++itr);
							t.nTerminusMass = lexical_cast<float>( *(++itr) );
							t.cTerminusMass = lexical_cast<float>( *(++itr) );

							s->tagList.insert(t);
						}

						if( indexItr == this->index.end() )
						{
							TaggingBaseList::push_back(s);
						}
					} else
						++itr;
				}
			}

			for( ListIterator itr = this->begin(); itr != this->end(); ++itr )
				(*itr)->tagList.ascendingScores = highScoresAreBetter;

			return sourceFilepath;
		}

		void writeTags(	const string& sourceFilepath,
						const string& filenameSuffix = "",
						const string& header = "" ) const
		{
            
			RunTimeVariableMap dummy;   
			writeTags( sourceFilepath, filenameSuffix, header, dummy );
        }

        
        void writeTagsMemory(string & tags, vector<float> & lowPeakMz)  //BX
        {
            ListConstIterator sItr = TaggingBaseList::begin();
            SpectrumType* s = *sItr;
            for( TagList::const_reverse_iterator tItr = s->tagList.rbegin(); tItr != s->tagList.rend(); ++tItr )
            {
                //small convert function
                string buffer = tItr->tag;
                tags.append(buffer);
                lowPeakMz.push_back(tItr->lowPeakMz);
            }
            
        }
        
        
		void writeTags(	const string& sourceFilepath,
                       const string& filenameSuffix,
                       const string& header,
                       const RunTimeVariableMap& vars ) const
		{
			string scanName = sourceFilepath.substr( sourceFilepath.find_last_of( SYS_PATH_SEPARATOR )+1,
													sourceFilepath.find_last_of( '.' ) - sourceFilepath.find_last_of( SYS_PATH_SEPARATOR )-1 );
			string filename = scanName + filenameSuffix + ".tags";
            
			ofstream fileStream( filename.c_str() );
			if( !fileStream.is_open() )
				throw invalid_argument( string( "unable to write tags file \"" ) + filename + "\"" );
            
			fileStream << header;
            
			fileStream << "H\tInputFile\t" << sourceFilepath << '\n';
            
			fileStream << showpoint << boolalpha;
			int n = 0;
			fileStream << "H\tTagsParameters\t" << vars.size();
			for( RunTimeVariableMap::const_iterator vItr = vars.begin(); vItr != vars.end(); ++vItr, ++n )
			{
				if( !(n % 4) )
					fileStream << "\nH\t";
				else
					fileStream << ", ";
                
				fileStream << vItr->first << ": " << vItr->second;
			}
			fileStream << "\n\n" << noshowpoint << noboolalpha;
            
			if( TaggingBaseList::empty() )
				return;
            
            // TODO: this will break if stringID or nativeID contain a tab
            fileStream << "H(S)\tID\tNativeID\tIndex\tCharge\tPrecursorNeutralMass"
                       << "\tTIC\tNormalizedTIC"
                       << "\tComplementaryPeakCount\tComplementaryTIC"
                       << "\tTagGraphPeakCount\tTagGraphTIC"
                       << "\tMzLowerBound\tMzUpperBound"
                       << "\tOriginalPeakCount\tFilteredPeakCount"
                       << "\tGeneratedTagCount\tFinalTagCount\n";
			fileStream << "H(T)\tTag\tLowPeakMz\tnTerminusMass\tcTerminusMass\tWorstPeakRank\tValid\tTotal";
			fileStream << "\tComplement\tIntensity\tMzFidelity";

			ListConstIterator firstTagPlaceItr;
			for( firstTagPlaceItr = TaggingBaseList::begin(); firstTagPlaceItr != TaggingBaseList::end(); ++firstTagPlaceItr )
				if( (*firstTagPlaceItr)->tagList.size() > 0 )
					break;

			if( firstTagPlaceItr != TaggingBaseList::end() )
			{
				/*for(	map< string, float >::const_iterator itr = (*firstTagPlaceItr)->tagList.begin()->scores.begin();
						itr != (*firstTagPlaceItr)->tagList.begin()->scores.end();
						++itr )
				{
					if( itr->first != "total" )
						fileStream << '\t' << itr->first;
				}*/
				fileStream << '\n';

				for( ListConstIterator sItr = TaggingBaseList::begin(); sItr != TaggingBaseList::end(); ++sItr )
				{
					SpectrumType* s = *sItr;

					fileStream	<< "S\t"
                                << s->stringID << '\t'
                                << s->nativeID << '\t'
								<< s->id.index << '\t'
								<< s->id.charge << '\t'
								<< s->mOfPrecursor << '\t'
                                << s->totalIonCurrent << '\t'
                                << s->totalIonCurrent / s->basePeakIntensity << '\t'
                                << s->complementClassCounts[0] << '\t'
                                << s->complementaryTIC << '\t'
                                << s->tagGraphPeakCount << '\t'
                                << s->tagGraphTIC << '\t'
								<< s->mzLowerBound << '\t'
								<< s->mzUpperBound << '\t'
								<< s->peakPreCount << '\t'
								<< s->peakCount << '\t'
								<< s->tagCount << '\t'
								<< s->tagList.size() << '\n';
                    string tagfilename = filename;
                    tagfilename.erase(tagfilename.size()-9);
                    tagfilename.append("short.tags");
                    ofstream tagstream(tagfilename.c_str());//TODO: test, this should output the tags needed
					for( TagList::const_reverse_iterator tItr = s->tagList.rbegin(); tItr != s->tagList.rend(); ++tItr )
					{
                        //TODO was muss ich denn hier noch machen? ...
                        if (tagstream.is_open())
                        {
                                //small convert function
                               string buffer = tItr->tag;
                               string tagik(buffer);
                               for (int i = 0; i < buffer.size(); i++)
                               {
                                   if (buffer[i] == 'Q') 
                                   {
                                       tagik.replace(i,1,1,'K');
                                       continue;
                                   }
                                   if (buffer[i] == 'L')
                                   {
                                       tagik.replace(i,1,1,'I');
                                       continue;
                                   }
                                   

                               }

                           tagstream << buffer << tagik << '\t'  << tItr->lowPeakMz <<'\n';
                        



                        }

                        fileStream	<< "T\t"
                            << tItr->tag << '\t'
                            << tItr->lowPeakMz << '\t'
                            << tItr->nTerminusMass << '\t'
                            << tItr->cTerminusMass << '\t'
                            << tItr->worstPeakRank << '\t'
                            << tItr->valid << '\t'
                            << tItr->totalScore;

                        //for( map< string, float >::const_iterator itr = tItr->scores.begin(); itr != tItr->scores.end(); ++itr )
                        //	if( itr->first != "total" )
                        //		fileStream << '\t' << itr->second;
                        //fileStream << '\n';
                        fileStream << '\t' << tItr->complementScore << '\t' << tItr->intensityScore << '\t' << tItr->mzFidelityScore << '\n';
                    }
                    tagstream.close();
                }
            }

            fileStream.close();
        }
    };
}

namespace std
{
    ostream& operator<< ( ostream& o, const TagInfo& rhs );
}

    BOOST_CLASS_IMPLEMENTATION( freicore::TagInfo, boost::serialization::object_serializable )
    BOOST_CLASS_IMPLEMENTATION( freicore::TagList, boost::serialization::object_serializable )
    BOOST_CLASS_IMPLEMENTATION( freicore::TaggingSpectrum, boost::serialization::object_serializable )
    BOOST_CLASS_TRACKING( freicore::TagInfo, boost::serialization::track_never )
    BOOST_CLASS_TRACKING( freicore::TagList, boost::serialization::track_never )
BOOST_CLASS_TRACKING( freicore::TaggingSpectrum, boost::serialization::track_never )

#endif
