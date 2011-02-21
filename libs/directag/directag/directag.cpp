#include "stdafx.h"
#include "freicore.h"
#include "base64.h"
#include "directagSpectrum.h"
#include "simplethreads.h"
#include "tagsFile.h"
#include "directag.h"
#include "Histogram.h"

//#ifdef USE_MPI
//#define NO_ZLIB 0
//#define BOOST_ZLIB_BINARY zlib.lib
//#define BOOST_IOSTREAMS_SOURCE 1
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/zlib.hpp>
//#endif

// Changed by BXu
// Moved from directag.h, messes up while linking otherwise
using namespace freicore;

namespace freicore
{
	namespace directag
	{
#ifdef USE_MPI
		void TransmitConfigsToChildProcesses();
		void ReceiveConfigsFromRootProcess();
#endif
		
		extern double lnCombin( int a, int b );
		extern float GetMassOfResidues( const string& a, bool b = false );
		
		struct tagFinder
		{
			tagFinder( const string& tagName = "" ) { name = tagName; }
			bool operator() ( TagInfo& test )
			{
				return name == test.tag;
			}
			string name;
		};
		
		struct taggingStats
		{
			taggingStats() :
			numSpectraTagged(0), numResidueMassGaps(0), 
			numTagsGenerated(0), numTagsRetained(0) {}
			size_t numSpectraTagged;
			size_t numResidueMassGaps;
			size_t numTagsGenerated;
			size_t numTagsRetained;
			
			taggingStats operator+ ( const taggingStats& rhs )
			{
				taggingStats tmp;
				tmp.numSpectraTagged = numSpectraTagged + rhs.numSpectraTagged;
				tmp.numResidueMassGaps = numResidueMassGaps + rhs.numResidueMassGaps;
				tmp.numTagsGenerated = numTagsGenerated + rhs.numTagsGenerated;
				tmp.numTagsRetained = numTagsRetained + rhs.numTagsRetained;
				return tmp;
			}
		};
		
		struct WorkerInfo : public BaseWorkerInfo
		{
			WorkerInfo( int num, int start, int end ) : BaseWorkerInfo( num, start, end ) {}
			taggingStats stats;
		};
		SpectraList					spectra;
		map< char, float >			compositionInfo;
		
		RunTimeConfig*				g_rtConfig;
		
		simplethread_mutex_t		resourceMutex;
		
		int						InitProcess( argList_t& args );
		int						ProcessHandler( int argc, vector<string> & argv, string &tags, vector<float> & lowPeakMzs, string & cachename);
		void					MakeResultFiles();
		void					GenerateForegroundTables();
		
		gapMap_t::iterator		FindPeakNear( gapMap_t&, float, float );
	}
}


//end changes BXu

namespace freicore
{
namespace directag
{
	double lnCombin( int a, int b ) { return lnCombin( a, b, g_lnFactorialTable ); }
	float GetMassOfResidues( const string& a, bool b ) { return g_residueMap->GetMassOfResidues( a, b ); }

	void WriteTagsToTagsFile(	const string& inputFilename,
								string startTime,
								string startDate,
								float totalTaggingTime,
                                string & tags,
                                vector<float> & lowPeakMzs)
	{
        
		//cout << g_hostString << " is generating output of tags." << endl;

		string filenameAsScanName;
		filenameAsScanName =	inputFilename.substr( inputFilename.find_last_of( SYS_PATH_SEPARATOR )+1,
								inputFilename.find_last_of( '.' ) - inputFilename.find_last_of( SYS_PATH_SEPARATOR )-1 );
		string outputFilename = filenameAsScanName + g_rtConfig->OutputSuffix + ".tags";

		stringstream header;
		header <<	"H\tTagsGenerator\tDirecTag\n" <<
					"H\tTagsGeneratorVersion\t" << DIRECTAG_VERSION_STRING << " (" << DIRECTAG_BUILD_DATE << ")\n";

		string license( DIRECTAG_LICENSE );
		boost::char_separator<char> delim("\n");
		tokenizer parser( license.begin(), license.begin() + license.length(), delim );
		for( tokenizer::iterator token = parser.begin(); token != parser.end(); ++token )
			header <<"H\t" << *token << "\n";

		header <<	"H\tTagging started at " << startTime << " on " << startDate << ".\n" <<
					"H\tTagging finished at " << GetTimeString() << " on " << GetDateString() << ".\n" <<
					"H\tTotal tagging time: " << totalTaggingTime << " seconds.\n" <<
					"H\tUsed " << g_numProcesses << " processing " << ( g_numProcesses > 1 ? "nodes" : "node" ) << ".\n";/* <<
					"H\tMean original (filtered) peak count: " << opcs[5] << " (" << fpcs[5] << ")\n" <<
					"H\tMin/max original (filtered) peak count: " << opcs[0] << " (" << fpcs[0] << ") / " << opcs[1] << " (" << fpcs[1] << ")\n" <<
					"H\tOriginal (filtered) peak count at 1st/2nd/3rd quartiles: " <<	opcs[2] << " (" << fpcs[2] << "), " <<
																						opcs[3] << " (" << fpcs[3] << "), " <<
																						opcs[4] << " (" << fpcs[4] << ")\n";*/
		if( !g_rtConfig->InlineValidationFile.empty() )
		{
			ofstream classCountsFile( string( filenameAsScanName + g_rtConfig->OutputSuffix + "-class-counts.txt" ).c_str() );
			classCountsFile << "<SpectrumId>\t<MatchedIonIntensityRankHistogram>\n";

			//ofstream scoreHistogramsFile( string( filenameAsScanName + g_rtConfig->OutputSuffix + "-histograms.txt" ).c_str() );
			Histogram<int> totalMatchedIonRanks;
			Histogram<int> totalMatchedBIonRanks;
			Histogram<int> totalMatchedYIonRanks;
			Histogram<int> totalMatchedYWaterLossRanks;
			Histogram<int> totalMatchedBWaterLossRanks;
			Histogram<int> totalMatchedYAmmoniaLossRanks;
			Histogram<int> totalMatchedBAmmoniaLossRanks;
			Histogram<int> totalLongestPathRanks;
			Histogram<int> totalValidLongestPathRanks;
			Histogram<int> totalIntensityRanksums;
			Histogram<int> totalValidIntensityRanksums;
			Spectrum* s;
			map< float, int > scoreToRanksumMap;

			/*for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
			{
				s = *sItr;

				for( TagList::iterator tItr = s->tagList.begin(); tItr != s->tagList.end(); ++tItr )
				{
					++ totalIntensityRanksums[tItr->ranksum];// += itr->second;
					if( tItr->valid )
						++ totalValidIntensityRanksums[tItr->ranksum];// += itr->second;
				}

				if( !s->resultSet.empty() )
				{
					Histogram<int> matchedIonRanks;
					SearchResultSet::reverse_iterator itr = s->resultSet.rbegin();
					vector< float > ionMasses;
					vector< string > ionLabels;
					bool allIonTypes[4] = { true, true, true, true };
					CalculateSequenceIons( itr->sequence, s->id.charge, &ionMasses, g_rtConfig->UseSmartPlusThreeModel, &ionLabels, 0, allIonTypes, &g_rtConfig->inlineValidationResidues );

					for( size_t i=0; i < ionMasses.size(); ++i )
					{
						PeakData::iterator pItr = s->peakData.findNear( ionMasses[i], g_rtConfig->FragmentMzTolerance );
						if( pItr != s->peakData.end() )
						{
							++ matchedIonRanks[ pItr->second.intensityRank ];
							++ totalMatchedIonRanks[ pItr->second.intensityRank ];
							if( ionLabels[i][0] == 'y' )
							{
								if( ionLabels[i].find( "H2O" ) != string::npos )
									++totalMatchedYWaterLossRanks[ pItr->second.intensityRank ];
								else if( ionLabels[i].find( "NH3" ) != string::npos )
									++totalMatchedYAmmoniaLossRanks[ pItr->second.intensityRank ];
								else
									++totalMatchedYIonRanks[ pItr->second.intensityRank ];
							} else if( ionLabels[i][0] == 'b' )
							{
								if( ionLabels[i].find( "H2O" ) != string::npos )
									++totalMatchedBWaterLossRanks[ pItr->second.intensityRank ];
								else if( ionLabels[i].find( "NH3" ) != string::npos )
									++totalMatchedBAmmoniaLossRanks[ pItr->second.intensityRank ];
								else
									++totalMatchedBIonRanks[ pItr->second.intensityRank ];
							}
						}
					}
				}
			}*/
			classCountsFile << "total\t" << totalMatchedIonRanks << endl;
			classCountsFile << "Ys\t" << totalMatchedYIonRanks << endl;
			classCountsFile << "Bs\t" << totalMatchedBIonRanks << endl;
			classCountsFile << "Y-H2Os\t" << totalMatchedYWaterLossRanks << endl;
			classCountsFile << "B-H2Os\t" << totalMatchedBWaterLossRanks << endl;
			classCountsFile << "Y-NH3s\t" << totalMatchedYAmmoniaLossRanks << endl;
			classCountsFile << "B-NH3s\t" << totalMatchedBAmmoniaLossRanks << endl;
			classCountsFile << "AllIntensityRanksums" << totalIntensityRanksums << endl;
			classCountsFile << "ValidIntensityRanksums" << totalValidIntensityRanksums << endl;
			classCountsFile << "AllLongestPathRanks" << totalLongestPathRanks << endl;
			classCountsFile << "ValidLongestPathRanks" << totalValidLongestPathRanks << endl;
		}

		//cout << g_hostString << " is writing tags to \"" << outputFilename << "\"." << endl;
		//spectra.writeTags( inputFilename, g_rtConfig->OutputSuffix, header.str(), g_rtConfig->getVariables() );
        
        spectra.writeTagsMemory(tags,lowPeakMzs); //BX
        
        spectra.clear();
	}

	void PrepareSpectra()
	{
		Timer timer;
/*PeakPreData foo;Sleep(3000);
Profiler foot(true);
for( size_t i=0; i < 6000000; ++i )
	foo[sqrt((float)i)] = pow((float)i,1.5f);
cout << endl << foot.End() << endl;
Sleep(3000);*/
		//if( g_numChildren == 0 )
			//cout << g_hostString << " is parsing " << spectra.size() << " spectra." << endl; //BX

		timer.Begin();
		for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
		{
			try
			{
				(*sItr)->Parse();
			} catch( exception& e )
			{
				throw runtime_error( string( "parsing scan \"" ) + string( (*sItr)->id ) + "\": " + e.what() );
			} catch( ... )
			{
				throw runtime_error( string( "parsing scan " ) + string( (*sItr)->id ) );
			}
		}

		if( g_numChildren == 0 )
		{
			//cout << g_hostString << " finished parsing its spectra; " << timer.End() << " seconds elapsed." << endl; //BX
			//cout << g_hostString << " is trimming spectra with less than " << 10 << " peaks." << endl; //BX
		}

		int preTrimCount = 0;
		try
		{
			preTrimCount = spectra.filterByPeakCount( 10 );
		} catch( exception& e )
		{
			throw runtime_error( string( "trimming spectra: " ) + e.what() );
		} catch( ... )
		{
			throw runtime_error( "trimming spectra" );
		}

		//if( g_numChildren == 0 ) //BX
//		{
//			cout << g_hostString << " trimmed " << preTrimCount << " spectra for being too small before peak filtering." << endl;
//			cout << g_hostString << " is determining spectrum charge states from " << spectra.size() << " spectra." << endl;
//		}

		timer.Begin();
		SpectraList duplicates;
		for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
		{
			try
			{
				if( !g_rtConfig->UseChargeStateFromMS )
						spectra.setId( (*sItr)->id, SpectrumId( (*sItr)->id.index, 0 ) );

				if( (*sItr)->id.charge == 0 )
				{
					SpectrumId preChargeId( (*sItr)->id );
					(*sItr)->DetermineSpectrumChargeState();
					SpectrumId postChargeId( (*sItr)->id );

					if( postChargeId.charge == 0 )
					{
						postChargeId.setCharge(2);

						if( g_rtConfig->DuplicateSpectra )
						{
							for( int z = 3; z <= g_rtConfig->NumChargeStates; ++z )
							{
								Spectrum* s = new Spectrum( *(*sItr) );
								s->id.setCharge(z);
								duplicates.push_back(s);
							}
						}
					}

					spectra.setId( preChargeId, postChargeId );
				}

			} catch( exception& e )
			{
				throw runtime_error( string( "duplicating scan " ) + string( (*sItr)->id ) + ": " + e.what() );
			} catch( ... )
			{
				throw runtime_error( string( "duplicating scan " ) + string( (*sItr)->id ) );
			}
		}

		try
		{
			spectra.insert( duplicates.begin(), duplicates.end(), spectra.end() );
			duplicates.clear(false);
		} catch( exception& e )
		{
			throw runtime_error( string( "adding duplicated spectra: " ) + e.what() );
		} catch( ... )
		{
			throw runtime_error( "adding duplicated spectra" );
		}

	//	if( g_numChildren == 0 ) //BX
//		{
//			cout << g_hostString << " finished determining spectrum charge states; " << timer.End() << " seconds elapsed." << endl;
//			cout << g_hostString << " is filtering peaks in " << spectra.size() << " spectra." << endl;
//		}

		timer.Begin();
		for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
		{
			try
			{
				(*sItr)->FilterPeaks();

				//if( !g_rtConfig->MakeSpectrumGraphs )
				//	(*sItr)->peakPreData.clear();
			} catch( exception& e )
			{
				throw runtime_error( string( "filtering peaks in scan: " ) + string( (*sItr)->id ) + e.what() );
			} catch( ... )
			{
				throw runtime_error( "filtering peaks in scan" );
			}
		}
//BX
//		if( g_numChildren == 0 )
//			cout << g_hostString << " finished filtering peaks; " << timer.End() << " seconds elapsed." << endl;

		int postTrimCount = 0;
		postTrimCount = spectra.filterByPeakCount( g_rtConfig->minIntensityClassCount );
//BX
//		if( g_numChildren == 0 )
//			cout << g_hostString << " trimmed " << postTrimCount << " spectra for being too small after peak filtering." << endl;
	}

	vector< int > workerNumbers;
	int numSearched;

	simplethread_return_t ExecutePipelineThread( simplethread_arg_t threadArg )
	{
		simplethread_lock_mutex( &resourceMutex );
		simplethread_id_t threadId = simplethread_get_id();
		WorkerThreadMap* threadMap = (WorkerThreadMap*) threadArg;
		WorkerInfo* threadInfo = reinterpret_cast< WorkerInfo* >( threadMap->find( threadId )->second );
		int numThreads = (int) threadMap->size();
		//if( g_numChildren == 0 ) //BX
		//	cout << threadInfo->workerHostString << " is initialized." << endl;
		simplethread_unlock_mutex( &resourceMutex );

		bool done;
		Timer executionTime(true);
		float totalExecutionTime = 0;
		float lastUpdate = 0;

		while( true )
		{
			simplethread_lock_mutex( &resourceMutex );
			done = workerNumbers.empty();
			if( !done )
			{
				threadInfo->workerNum = workerNumbers.back();
				workerNumbers.pop_back();
			}
			simplethread_unlock_mutex( &resourceMutex );

			if( done )
				break;

			threadInfo->endIndex = ( spectra.size() / g_numWorkers )-1;

			//cout << threadInfo->workerHostString << " " << numProteins << " " << g_numWorkers << endl;

			Spectrum* s;
			SpectraList::iterator sItr = spectra.begin();
			for(	advance_to_bound( sItr, spectra.end(), threadInfo->workerNum );
					sItr != spectra.end();
					advance_to_bound( sItr, spectra.end(), g_numWorkers ) )
			{
				s = (*sItr);
				++ threadInfo->stats.numSpectraTagged;

				//s->DetermineSpectrumChargeState();
				START_PROFILER(0)
				s->Preprocess();
				STOP_PROFILER(0)

				if( (int) (*sItr)->peakPreData.size() < g_rtConfig->minIntensityClassCount )
					continue;

				START_PROFILER(1)
				threadInfo->stats.numResidueMassGaps += s->MakeTagGraph();
				STOP_PROFILER(1)

				s->MakeProbabilityTables();

				//s->tagGraphs.clear();
				//s->nodeSet.clear();
				deallocate(s->nodeSet);

				START_PROFILER(2)
				threadInfo->stats.numTagsGenerated += s->Score();
				STOP_PROFILER(2)

				threadInfo->stats.numTagsRetained += s->tagList.size();

				//s->gapMaps.clear();
				//s->tagGraphs.clear();
				deallocate(s->gapMaps);
				deallocate(s->tagGraphs);
				if( ( !g_rtConfig->MakeSpectrumGraphs && g_rtConfig->InlineValidationFile.empty() ) )
				{
					//s->peakPreData.clear();
					//s->peakData.clear();
					deallocate(s->peakPreData);
					deallocate(s->peakData);
				}

				if( g_numChildren == 0 )
					totalExecutionTime = executionTime.TimeElapsed();

				if( g_numChildren == 0 && ( ( totalExecutionTime - lastUpdate > g_rtConfig->StatusUpdateFrequency ) || s->id == ((*spectra.rbegin())->id) ) )
				{
					//int curSpectrum = ( i + 1 ) / g_numWorkers;
					float spectraPerSec = float( threadInfo->stats.numSpectraTagged ) / totalExecutionTime;
					float estimatedTimeRemaining = float( spectra.size() - threadInfo->stats.numSpectraTagged ) / spectraPerSec / numThreads;

					simplethread_lock_mutex( &resourceMutex );
//					cout << threadInfo->workerHostString << " has sequence tagged " << threadInfo->stats.numSpectraTagged << " of " << spectra.size() <<
//							" spectra; " << spectraPerSec << " per second, " << totalExecutionTime << " elapsed, " << estimatedTimeRemaining << " remaining." << endl;
//					cout << threadInfo->workerHostString << " stats: " << 1 << " / " <<
//							threadInfo->stats.numSpectraTagged << " / " <<	
//							threadInfo->stats.numResidueMassGaps << " / " <<
//							threadInfo->stats.numTagsGenerated << " / " <<
//							threadInfo->stats.numTagsRetained << endl;

					PRINT_PROFILERS( cout, threadInfo->workerHostString + " profiling" )

					simplethread_unlock_mutex( &resourceMutex );

					lastUpdate = totalExecutionTime;
				}
			}

		}

		return 0;
	}
}
}

namespace std {
	ostream& operator<< ( ostream& o, const list< freicore::directag::Spectrum* >::iterator& itr )
	{
		return o << "itr";//*itr;
	}
}

namespace freicore {
namespace directag {
	taggingStats ExecutePipeline()
	{
		WorkerThreadMap workerThreads;
		int numProcessors = g_numWorkers;

		//if( !g_rtConfig->UseMultipleProcessors )
		//	numProcessors = 1;

		int numSpectra = (int) spectra.size();

		float maxPeakSpace = 0;
		for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
			if( (*sItr)->totalPeakSpace > maxPeakSpace )
				maxPeakSpace = (*sItr)->totalPeakSpace;

		//cout << "Resizing lnTable to " << maxPeakSpace << endl;
		g_lnFactorialTable.resize( (int) ceil( maxPeakSpace ) );

		if( g_rtConfig->UseMultipleProcessors && g_numWorkers > 1 )
		{
			g_numWorkers = min( numSpectra, g_numWorkers * g_rtConfig->ThreadCountMultiplier );

			for( int i=0; i < g_numWorkers; ++i )
				workerNumbers.push_back(i);

			simplethread_handle_array_t workerHandles;

			simplethread_lock_mutex( &resourceMutex );
			for( int t = 0; t < numProcessors; ++t )
			{
				simplethread_id_t threadId;
				simplethread_handle_t threadHandle = simplethread_create_thread( &threadId, &ExecutePipelineThread, &workerThreads );
				workerThreads[ threadId ] = new WorkerInfo( t, 0, 0 );
				workerHandles.array.push_back( threadHandle );
			}
			simplethread_unlock_mutex( &resourceMutex );

			simplethread_join_all( &workerHandles );
			//cout << g_hostString << " tagged " << numSearched << " of " << spectra.size() << " spectra." << endl;

		} else
		{
			//cout << g_hostString << " is preparing " << numSpectra << " unprepared spectra." << endl;
			g_numWorkers = 1;
			workerNumbers.push_back(0);
			simplethread_id_t threadId = simplethread_get_id();
			workerThreads[ threadId ] = new WorkerInfo( 0, 0, 0 );
			ExecutePipelineThread( &workerThreads );
		}

		taggingStats stats;

		for( WorkerThreadMap::iterator itr = workerThreads.begin(); itr != workerThreads.end(); ++itr )
		{	
			stats = stats + reinterpret_cast< WorkerInfo* >( itr->second )->stats;
			delete itr->second;
		}
		g_numWorkers = numProcessors;

		return stats;
	}

	#ifdef USE_MPI

	void TransmitConfigsToChildProcesses()
	{
		for( int p=0; p < g_numChildren; ++p )
		{
			int len;

			len = (int) g_rtConfig->cfgStr.length();
			MPI_Send( &len,									1,		MPI_INT,			p+1,	0x00, MPI_COMM_WORLD );
			MPI_Send( (void*) g_rtConfig->cfgStr.c_str(),	len,	MPI_CHAR,			p+1,	0x01, MPI_COMM_WORLD );

			len = (int) g_residueMap->cfgStr.length();
			MPI_Send( &len,									1,		MPI_INT,			p+1,	0x02, MPI_COMM_WORLD );
			MPI_Send( (void*) g_residueMap->cfgStr.c_str(),	len,	MPI_CHAR,			p+1,	0x03, MPI_COMM_WORLD );
		}
	}

	void ReceiveConfigsFromRootProcess()
	{
		int len;

		MPI_Recv( &len,								1,		MPI_INT,			0,		0x00, MPI_COMM_WORLD, &st );
		g_rtConfig->cfgStr.resize( len );
		MPI_Recv( &g_rtConfig->cfgStr[0],			len,	MPI_CHAR,			0,		0x01, MPI_COMM_WORLD, &st );
		g_rtConfig->initializeFromBuffer( g_rtConfig->cfgStr, "\r\n#" );

		g_residueMap = new ResidueMap();
		MPI_Recv( &len,								1,		MPI_INT,			0,		0x02, MPI_COMM_WORLD, &st );
		g_residueMap->cfgStr.resize( len );
		MPI_Recv( &g_residueMap->cfgStr[0],			len,	MPI_CHAR,			0,		0x03, MPI_COMM_WORLD, &st );
		g_residueMap->initializeFromBuffer( g_residueMap->cfgStr );
	}

	int TransmitUnpreparedSpectraToChildProcesses()
	{
		int numSpectra = (int) spectra.size();

		int sourceProcess, batchSize;
		bool IsFinished = false;

		Timer PrepareTime( true );
		float totalPrepareTime = 0.01f;
		float lastUpdate = 0.0f;

		int i = 0;
		int numChildrenFinished = 0;
		while( numChildrenFinished < g_numChildren )
		{
			stringstream packStream;
			binary_oarchive packArchive( packStream );
			// For every batch, listen for a worker process that is ready to receive it

			#ifdef MPI_DEBUG
				cout << g_hostString << " is listening for a child process to offer to prepare some spectra." << endl;
			#endif

			if( i < numSpectra )
			{
				batchSize = min( numSpectra-i, g_rtConfig->SpectraBatchSize );

				try
				{
					packArchive & batchSize;

					SpectraList::iterator sItr = spectra.begin();
					advance( sItr, i );
					for( int j = i; j < i + batchSize; ++j, ++sItr )
					{
						packArchive & **sItr;
						//packArchive & (*sItr)->resultSet;
					}
				} catch( exception& e )
				{
					cerr << g_hostString << " had an error: " << e.what() << endl;
					exit(1);
				}

				i += batchSize;
			} else
			{
				batchSize = 0;
				packArchive & batchSize;

				#ifdef MPI_DEBUG
					cout << "Process #" << sourceProcess << " has been informed that preparation is complete." << endl;
				#endif

				++numChildrenFinished;
			}

			MPI_Recv( &sourceProcess,		1,		MPI_INT,	MPI_ANY_SOURCE,	0xFF, MPI_COMM_WORLD, &st );

			string pack = packStream.str();
			int len = (int) pack.length();

			MPI_Send( &len,					1,		MPI_INT,	sourceProcess,	0x00, MPI_COMM_WORLD );
			MPI_Send( (void*) pack.c_str(),	len,	MPI_CHAR,	sourceProcess,	0x01, MPI_COMM_WORLD );

			totalPrepareTime = PrepareTime.TimeElapsed();
			if( !IsFinished && ( ( totalPrepareTime - lastUpdate > g_rtConfig->StatusUpdateFrequency ) || i == numSpectra ) )
			{
				if( i == numSpectra )
					IsFinished = true;

				float spectraPerSec = float(i) / totalPrepareTime;
				float estimatedTimeRemaining = float(numSpectra-i) / spectraPerSec;
				//cout << g_hostString << " has prepared " << i << " of " << numSpectra << " spectra; " << spectraPerSec <<
				//		" per second, " << estimatedTimeRemaining << " seconds remaining." << endl; //BX
				
                lastUpdate = totalPrepareTime;
			}
		}

		return 0;
	}

	int ReceiveUnpreparedSpectraBatchFromRootProcess()
	{
		int batchSize;

		MPI_Ssend( &g_pid,		1,	MPI_INT,		0,	0xFF, MPI_COMM_WORLD );
		string pack;
		int len;

		#ifdef MPI_DEBUG
			cout << g_hostString << " is receiving a batch of unprepared spectra." << endl;
			Timer receiveTime(true);
		#endif
		MPI_Recv( &len,					1,			MPI_INT,	0,	0x00, MPI_COMM_WORLD, &st );
		pack.resize( len );
		MPI_Recv( (void*) pack.data(),	len,		MPI_CHAR,	0,	0x01, MPI_COMM_WORLD, &st );

		stringstream packStream( pack );
		binary_iarchive packArchive( packStream );

		try
		{
			packArchive & batchSize;

			if( !batchSize )
			{
				#ifdef MPI_DEBUG
					cout << g_hostString << " is informed that all spectra have been prepared." << endl;
				#endif

				return 0; // do not expect another batch
			}

			#ifdef MPI_DEBUG
				cout << g_hostString << " finished receiving a batch of " << batchSize << " unprepared spectra; " <<
						receiveTime.End() << " seconds elapsed." << endl;
			#endif

			for( int j=0; j < batchSize; ++j )
			{
				Spectrum* s = new Spectrum;
				packArchive & *s;
				spectra.push_back( s );

			}
		} catch( exception& e )
		{
			cerr << g_hostString << " had an error: " << e.what() << endl;
			exit(1);
		}

		return 1; // expect another batch
	}

	int TransmitPreparedSpectraToRootProcess( SpectraList& preparedSpectra )
	{
		int numSpectra = (int) preparedSpectra.size();

		stringstream packStream;
		binary_oarchive packArchive( packStream );

		//Timer packTime(true);
		//cout << g_hostString << " is packing " << numSpectra << " results." << endl;
		packArchive & numSpectra;
		for( SpectraList::iterator sItr = preparedSpectra.begin(); sItr != preparedSpectra.end(); ++sItr )
		{
			Spectrum* s = *sItr;
			packArchive & *s;
		}
		//cout << g_hostString << " finished packing results; " << packTime.End() << " seconds elapsed." << endl;

		MPI_Ssend( &g_pid,		1,	MPI_INT,		0,	0xEE, MPI_COMM_WORLD );

		#ifdef MPI_DEBUG
			cout << g_hostString << " is sending " << numSpectra << " prepared spectra." << endl;
			Timer sendTime(true);
		#endif

		stringstream compressedStream;
		boost::iostreams::filtering_ostream compressorStream;
		compressorStream.push( boost::iostreams::zlib_compressor() );
		compressorStream.push( compressedStream );
		boost::iostreams::copy( packStream, compressorStream );
		compressorStream.reset();

		string pack = compressedStream.str();
		int len = (int) pack.length();
		MPI_Send( &len,					1,		MPI_INT,	0,	0x00, MPI_COMM_WORLD );
		MPI_Send( (void*) pack.c_str(),	len,	MPI_CHAR,	0,	0x01, MPI_COMM_WORLD );

		#ifdef MPI_DEBUG
			cout << g_hostString << " finished sending " << numSpectra << " prepared spectra; " <<
					sendTime.End() << " seconds elapsed." << endl;
		#endif

		return 0;
	}

	int ReceivePreparedSpectraFromChildProcesses()
	{
		Timer receiveTime( true );
		float totalReceiveTime = 0.01f;
		float lastUpdate = 0.0f;
		int sourceProcess, numSpectra;
		for( int p=0; p < g_numChildren; ++p )
		{
			MPI_Recv( &sourceProcess,		1,	MPI_INT,	MPI_ANY_SOURCE,	0xEE, MPI_COMM_WORLD, &st );

			#ifdef MPI_DEBUG
				cout << g_hostString << " is receiving " << numSpectra << " prepared spectra." << endl;
				Timer receiveTime(true);
			#endif

			string pack;
			int len;

			MPI_Recv( &len,					1,			MPI_INT,	sourceProcess,	0x00, MPI_COMM_WORLD, &st );
			pack.resize( len );
			MPI_Recv( (void*) pack.data(),	len,		MPI_CHAR,	sourceProcess,	0x01, MPI_COMM_WORLD, &st );

			stringstream compressedStream( pack );
			stringstream packStream;
			boost::iostreams::filtering_ostream decompressorStream;
			decompressorStream.push( boost::iostreams::zlib_decompressor() );
			decompressorStream.push( packStream );
			boost::iostreams::copy( compressedStream, decompressorStream );
			decompressorStream.reset();

			binary_iarchive packArchive( packStream );

			try
			{
				packArchive & numSpectra;

				//cout << g_hostString << " is unpacking results for " << numSpectra << " spectra." << endl;
				for( int j=0; j < numSpectra; ++j )
				{
					Spectrum* s = new Spectrum;
					packArchive & *s;
					spectra.push_back( s );

				}
				//cout << g_hostString << " is finished unpacking results." << endl;
			} catch( exception& e )
			{
				cerr << g_hostString << " had an error: " << e.what() << endl;
				exit(1);
			}

			#ifdef MPI_DEBUG
				cout << g_hostString << " finished receiving " << numSpectra << " prepared spectra; " <<
						receiveTime.End() << " seconds elapsed." << endl;
			#endif

			totalReceiveTime = receiveTime.TimeElapsed();
			if( ( totalReceiveTime - lastUpdate > g_rtConfig->StatusUpdateFrequency ) || p+1 == g_numChildren )
			{
				float nodesPerSec = float(p+1) / totalReceiveTime;
				float estimatedTimeRemaining = float(g_numChildren-p-1) / nodesPerSec;
				cout << g_hostString << " has received prepared spectra from " << p+1 << " of " << g_numChildren << " worker nodes; " << nodesPerSec <<
						" per second, " << estimatedTimeRemaining << " seconds remaining." << endl;
				lastUpdate = totalReceiveTime;
			}
		}
		return 0;
	}

	int TransmitUntaggedSpectraToChildProcesses()
	{
		int numSpectra = (int) spectra.size();

		int sourceProcess, batchSize;
		bool IsFinished = false;

		Timer PrepareTime( true );
		float totalPrepareTime = 0.01f;
		float lastUpdate = 0.0f;

		int i = 0;
		int numChildrenFinished = 0;

		while( numChildrenFinished < g_numChildren )
		{
			stringstream packStream;
			binary_oarchive packArchive( packStream );

			if( i < numSpectra )
			{
				batchSize = min( numSpectra-i, g_rtConfig->SpectraBatchSize );

				try
				{
					packArchive & batchSize;

					SpectraList::iterator sItr = spectra.begin();
					advance( sItr, i );
					for( int j = i; j < i + batchSize; ++j, ++sItr )
					{
						packArchive & **sItr;
						delete *sItr;
					}
					//spectra.erase( sItr );
				} catch( exception& e )
				{
					cerr << g_hostString << " had an error: " << e.what() << endl;
					exit(1);
				}

				i += batchSize;
			} else
			{
				batchSize = 0;
				packArchive & batchSize;

				#ifdef MPI_DEBUG
					cout << "Process #" << sourceProcess << " has been informed that tagging is complete." << endl;
				#endif

				++numChildrenFinished;
			}

			// For every batch, listen for a worker process that is ready to receive it

			#ifdef MPI_DEBUG
				cout << g_hostString << " is listening for a child process to offer to sequence tag some spectra." << endl;
			#endif

			MPI_Recv( &sourceProcess,			1,		MPI_INT,	MPI_ANY_SOURCE,	0xFF, MPI_COMM_WORLD, &st );

			#ifdef MPI_DEBUG
				cout << g_hostString << " is sending a batch of " << batchSize << " unsequenced spectra." << endl;
						Timer sendTime(true);
			#endif

			string pack = packStream.str();
			int len = (int) pack.length();
			MPI_Send( &len,						1,		MPI_INT,	sourceProcess,	0x00, MPI_COMM_WORLD );
			MPI_Send( (void*) pack.c_str(),		len,	MPI_CHAR,	sourceProcess,	0x01, MPI_COMM_WORLD );

			#ifdef MPI_DEBUG
				cout << g_hostString << " finished sending a batch of spectra; " <<
						sendTime.End() << " seconds elapsed." << endl;
			#endif

			totalPrepareTime = PrepareTime.TimeElapsed();
			if( ( totalPrepareTime - lastUpdate > g_rtConfig->StatusUpdateFrequency ) || ( !IsFinished && i == numSpectra ) )
			{
				if( i == numSpectra )
					IsFinished = true;

				float spectraPerSec = float(i) / totalPrepareTime;
				float estimatedTimeRemaining = float(numSpectra-i) / spectraPerSec;
//              cout << g_hostString << " has sequence tagged " << i << " of " << numSpectra << " spectra; " << spectraPerSec <<
//						" per second, " << estimatedTimeRemaining << " seconds remaining." << endl; //BX

				lastUpdate = totalPrepareTime;
			}
		}
		spectra.clear(false);
		return 0;
	}

	int ReceiveUntaggedSpectraBatchFromRootProcess()
	{
		int batchSize;

		MPI_Ssend( &g_pid,		1,	MPI_INT,		0,	0xFF, MPI_COMM_WORLD );

		string pack;
		int len;

		#ifdef MPI_DEBUG
			cout << g_hostString << " is receiving a batch of unsequenced spectra." << endl;
			Timer receiveTime(true);
		#endif
		MPI_Recv( &len,					1,			MPI_INT,	0,	0x00, MPI_COMM_WORLD, &st );
		pack.resize( len );
		MPI_Recv( (void*) pack.data(),	len,		MPI_CHAR,	0,	0x01, MPI_COMM_WORLD, &st );

		stringstream packStream( pack );
		binary_iarchive packArchive( packStream );

		try
		{
			packArchive & batchSize;

			if( !batchSize )
			{
				#ifdef MPI_DEBUG
					cout << g_hostString << " is informed that all spectra have been sequence tagged." << endl;
				#endif

				return 0; // do not expect another batch
			}

			#ifdef MPI_DEBUG
				cout << g_hostString << " finished receiving a batch of " << batchSize << " unsequenced spectra; " <<
						receiveTime.End() << " seconds elapsed." << endl;
			#endif

			for( int j=0; j < batchSize; ++j )
			{
				Spectrum* s = new Spectrum;
				packArchive & *s;
				spectra.push_back( s );
			}
		} catch( exception& e )
		{
			cerr << g_hostString << " had an error: " << e.what() << endl;
			exit(1);
		}

		return 1; // expect another batch
	}

	int ReceiveTaggedSpectraFromChildProcesses()
	{
		int sourceProcess, numSpectra;
		for( int p=0; p < g_numChildren; ++p )
		{
			MPI_Recv( &sourceProcess,					1,	MPI_INT,	MPI_ANY_SOURCE,	0xEE, MPI_COMM_WORLD, &st );

			#ifdef MPI_DEBUG
				cout << g_hostString << " is receiving " << numSpectra << " prepared spectra." << endl;
				Timer receiveTime(true);
			#endif

			string pack;
			int len;

			MPI_Recv( &len,					1,			MPI_INT,	sourceProcess,	0x00, MPI_COMM_WORLD, &st );
			pack.resize( len );
			MPI_Recv( (void*) pack.data(),	len,		MPI_CHAR,	sourceProcess,	0x01, MPI_COMM_WORLD, &st );

			stringstream compressedStream( pack );
			stringstream packStream;
			boost::iostreams::filtering_ostream decompressorStream;
			decompressorStream.push( boost::iostreams::zlib_decompressor() );
			decompressorStream.push( packStream );
			boost::iostreams::copy( compressedStream, decompressorStream );
			decompressorStream.reset();

			binary_iarchive packArchive( packStream );

			try
			{
				packArchive & numSpectra;

				//cout << g_hostString << " is unpacking results for " << numSpectra << " spectra." << endl;
				for( int j=0; j < numSpectra; ++j )
				{
					Spectrum* s = new Spectrum;

					packArchive & *s;

					spectra.push_back( s );

				}
				//cout << g_hostString << " is finished unpacking results." << endl;
			} catch( exception& e )
			{
				cerr << g_hostString << " had an error: " << e.what() << endl;
				exit(1);
			}

			#ifdef MPI_DEBUG
				cout << g_hostString << " finished receiving " << numSpectra << " prepared spectra; " <<
						receiveTime.End() << " seconds elapsed." << endl;
			#endif
		}

		return 0;
	}

	int TransmitTaggedSpectraToRootProcess( SpectraList& preparedSpectra )
	{
		int numSpectra = (int) preparedSpectra.size();

		stringstream packStream;
		binary_oarchive packArchive( packStream );

		//Timer packTime(true);
		//cout << g_hostString << " is packing " << numSpectra << " results." << endl;
		packArchive & numSpectra;
		for( SpectraList::iterator sItr = preparedSpectra.begin(); sItr != preparedSpectra.end(); ++sItr )
		{
			Spectrum* s = (*sItr);

			packArchive & *s;
		}
		//cout << g_hostString << " finished packing results; " << packTime.End() << " seconds elapsed." << endl;

		MPI_Ssend( &g_pid,			1,				MPI_INT,	0,	0xEE, MPI_COMM_WORLD );

	#ifdef MPI_DEBUG
			cout << g_hostString << " is sending " << numSpectra << " packed results." << endl;
			Timer sendTime(true);
	#endif

		stringstream compressedStream;
		boost::iostreams::filtering_ostream compressorStream;
		compressorStream.push( boost::iostreams::zlib_compressor() );
		compressorStream.push( compressedStream );
		boost::iostreams::copy( packStream, compressorStream );
		compressorStream.reset();

		string pack = compressedStream.str();
		int len = (int) pack.length();
		MPI_Send( &len,					1,		MPI_INT,	0,	0x00, MPI_COMM_WORLD );
		MPI_Send( (void*) pack.c_str(),	len,	MPI_CHAR,	0,	0x01, MPI_COMM_WORLD );

	#ifdef MPI_DEBUG
			cout << g_hostString << " finished sending " << numSpectra << " packed results; " <<
					sendTime.End() << " seconds elapsed." << endl;
	#endif

		return 0;
	}

	#endif

	void EncodeSpectraForOutput()
	{
		// Re-encode each spectrum's peak data in base64
		for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
		{
			//if( (*sItr)->id.index == 678 )
			//	cout << (*sItr)->peakData.begin()->first << " - " << (*sItr)->peakData.rbegin()->first << endl;

			mzData_t mzData;
			iData_t iData;
			int outSize = (int) (*sItr)->peakData.size() * sizeof(float) * 4;
			outSize = (outSize / 3) + (4 - ((outSize / 3) % 4));

			(*sItr)->mzData.resize( outSize );
			//(*sItr)->iData = new char[ outSize+1 ];
			(*sItr)->iData.resize( (*sItr)->peakData.size() );

			for( PeakData::iterator itr = (*sItr)->peakData.begin(); itr != (*sItr)->peakData.end(); ++itr )
			{
				mzData.push_back( itr->first );
				//iData.push_back( itr->second.inten );
				iData.push_back( (float) itr->second.intensityRank );
			}

			int len;
			float* dataArray = new float[ (*sItr)->peakData.size() ];

			for( size_t j=0; j < (*sItr)->peakData.size(); ++j )
				dataArray[j] = (float) mzData[j];

			unsigned char* mzDataBuffer = (unsigned char*) dataArray;

			len = b64_encode( &(*sItr)->mzData[0], mzDataBuffer, (int) (*sItr)->peakData.size() * sizeof(float) );
			(*sItr)->mzData = (*sItr)->mzData.c_str();

			//cout << (*sItr)->mzData << endl;

			char tmp[2];
			for( size_t j=0; j < (*sItr)->peakData.size(); ++j )
			{
				//cout << i << ":" << (int) iData[j] << "\t";
				sprintf( tmp, "%d", (int) iData[j] );
				(*sItr)->iData[j] = tmp[0];
			}

			//cout << (*sItr)->iData << endl;

			/*for( size_t j=0; j < (*sItr)->peakData.size(); ++j )
				dblArray[j] = (double) iData[j];
			unsigned char* iDataBuffer = (unsigned char*) dblArray;
			len = b64_encode( (*sItr)->iData, iDataBuffer, (*sItr)->peakData.size() * sizeof(double) );
			(*sItr)->iData[len] = 0;*/
			(*sItr)->iData[ (*sItr)->peakData.size() ] = 0;

			(*sItr)->mzDataIsDoublePrecision = false;
			(*sItr)->mzDataEndianType = g_endianType;
			(*sItr)->iDataIsDoublePrecision = false;
			(*sItr)->iDataEndianType = g_endianType;


			/*if( (*sItr)->id.index == 678 )
			{
				vector< float > data;

				//cout << g_elementData << ": ";
				string blah = (*sItr)->mzData;
				base64_decode_sequence( &blah, data, (*sItr)->mzDataIsDoublePrecision,
								g_endianType, (*sItr)->mzDataEndianType );
				cout << data[0] << " - " << data[ data.size()-1 ] << endl;
			}*/

			delete dataArray;
		}
	}

	int InitProcess( argList_t& args )
	{
		//cout << g_hostString << " is initializing." << endl;
//BX	if( g_pid == 0)
//		{
//			cout << "DirecTag " << DIRECTAG_VERSION_STRING << " (" << DIRECTAG_BUILD_DATE << ")\n" <<
//					DIRECTAG_LICENSE << endl;
//		}

		proteinStore proteins;

		g_residueMap = new ResidueMap;
		g_rtConfig = new RunTimeConfig;
		g_rtSharedConfig = (BaseRunTimeConfig*) g_rtConfig;
		g_endianType = GetHostEndianType();
		g_numWorkers = GetNumProcessors();

		// First set the working directory, if provided
		for( size_t i=1; i < args.size(); ++i )
		{
			if( args[i] == "-workdir" && i+1 <= args.size() )
			{
				chdir( args[i+1].c_str() );
				args.erase( args.begin()+i );
			} else if( args[i] == "-cpus" && i+1 <= args.size() )
			{
				g_numWorkers = atoi( args[i+1].c_str() );
				args.erase( args.begin()+i );
			} else
				continue;
			args.erase( args.begin()+i );
			--i;
		}

		if( g_pid == 0 )
		{
			for( size_t i=1; i < args.size(); ++i )
			{
				if( args[i] == "-cfg" && i+1 <= args.size() )
				{
					if( g_rtConfig->initializeFromFile( args[i+1] ) )
					{
						cerr << g_hostString << " could not find runtime configuration at \"" << args[i+1] << "\"." << endl;
						return 1;
					}
					args.erase( args.begin() + i );

				} else if( args[i] == "-rescfg" && i+1 <= args.size() )
				{
					if( g_residueMap->initializeFromFile( args[i+1] ) )
					{
						cerr << g_hostString << " could not find residue masses at \"" << args[i+1] << "\"." << endl;
						return 1;
					}
					args.erase( args.begin() + i );
				} else
					continue;

				args.erase( args.begin() + i );
				--i;
			}

			if( args.size() < 2 )
			{
				cout << "Not enough arguments.\nUsage: " << args[0] << " <input spectra filemask 1> [input spectra filemask 2] ..." << endl;
				return 1;
			}

			if( !g_rtConfig->initialized() )
			{
				if( g_rtConfig->initializeFromFile() )
				{
					//cerr << g_hostString << " could not find the default configuration file (hard-coded defaults in use)." << endl; //BX
				}
				//return 1;
			}

			if( !g_residueMap->initialized() )
			{
				if( g_residueMap->initializeFromFile() )
				{
					//cerr << g_hostString << " could not find the default residue masses file (hard-coded defaults in use)." << endl; //BX
				}
			}

			#ifdef USE_MPI
				if( g_numChildren > 0 )
					TransmitConfigsToChildProcesses();
			#endif

		} else // child process
		{
			#ifdef USE_MPI
				ReceiveConfigsFromRootProcess();
			#endif
		}

		// Command line overrides happen after config file has been distributed but before PTM parsing
		RunTimeVariableMap vars = g_rtConfig->getVariables();
		for( RunTimeVariableMap::iterator itr = vars.begin(); itr != vars.end(); ++itr )
		{
			string varName;
			varName += "-" + itr->first;

			for( size_t i=1; i < args.size(); ++i )
			{
				if( args[i] == varName && i+1 <= args.size() )
				{
					//cout << varName << " " << itr->second << " " << args[i+1] << endl;
					itr->second = args[i+1];
					args.erase( args.begin() + i );
					args.erase( args.begin() + i );
					--i;
				}
			}
		}
		g_rtConfig->setVariables( vars );

		if( g_pid == 0 )
		{
			for( size_t i=1; i < args.size(); ++i )
			{
				if( args[i] == "-dump" )
				{
					g_rtConfig->dump();
					g_residueMap->dump();
					args.erase( args.begin() + i );
					--i;
				}
			}

			for( size_t i=1; i < args.size(); ++i )
			{
				if( args[i][0] == '-' )
				{
					cerr << "Warning: ignoring unrecognized parameter \"" << args[i] << "\"" << endl;
					args.erase( args.begin() + i );
					--i;
				}
			}
		}

		return 0;
	}

	void ReadInlineValidationFile()
	{
		/*if( !g_rtConfig->InlineValidationFile.empty() )
		{
			g_rtConfig->inlineValidationResidues = *g_residueMap;
			fileList_t sqtFilenames;

			if( g_pid == 0 ) cout << "Finding SQT files matching mask \"" << g_rtConfig->InlineValidationFile << "\"" << endl;
			FindFilesByMask( g_rtConfig->InlineValidationFile, sqtFilenames );

			if( sqtFilenames.empty() )
			{
				if( g_pid == 0 ) cerr << "No files found matching given filemasks." << endl;
			} else
			{
				// Read SQT files for validation
				RunTimeVariableMap varsFromFile( "NumChargeStates DynamicMods StaticMods UseAvgMassOfSequences" );
				for( fileList_t::iterator fItr = sqtFilenames.begin(); fItr != sqtFilenames.end(); ++fItr )
				{
					if( g_pid == 0 ) cout << "Reading peptide identifications from \"" << *fItr << "\"" << endl;
					spectra.readSQT( *fItr, true, true, g_rtConfig->ValidationMode, " ", varsFromFile );
					cout << "Setting DynamicMods and StaticMods from SQT file: " << varsFromFile["DynamicMods"] << "; " << varsFromFile["StaticMods"] << endl;
					g_rtConfig->inlineValidationResidues.setDynamicMods( varsFromFile["DynamicMods"] );
					g_rtConfig->inlineValidationResidues.setStaticMods( varsFromFile["StaticMods"] );
					varsFromFile.erase( "DynamicMods" );
					varsFromFile.erase( "StaticMods" );
					g_rtConfig->setVariables( varsFromFile );
				}

				if( spectra.empty() )
				{
					if( g_pid == 0 ) cout << "No identifications found." << endl;

				} else
				{
					if( g_pid == 0 ) cout << "Finished reading " << spectra.size() << " identifications, now calculating validation thresholds." << endl;

					if( g_rtConfig->StartSpectraScanNum == 0 && g_rtConfig->EndSpectraScanNum == -1 )
					{
						spectra.calculateValidationThresholds(	g_rtConfig->scoreThresholds,
																g_rtConfig->NumChargeStates,
																g_rtConfig->Confidence,
																g_rtConfig->DecoyRatio,
																g_rtConfig->DecoyPrefix,														
																g_rtConfig->ValidationMode );

						//SpectraList originalSpectra = spectra;
						vector< size_t > potentialMatchCounts, validMatchCounts;
						pair< SpectraList, SpectraList > filteredSpectra = spectra.filterByThresholds(	g_rtConfig->scoreThresholds,
																										g_rtConfig->NumChargeStates,
																										g_rtConfig->ValidationMode,
																										&potentialMatchCounts,
																										&validMatchCounts );

						for(	SpectraList::ListIndexIterator itr = spectra.index.begin();
								itr != spectra.index.end();
								++itr )
						{
								if( g_rtConfig->InlineValidationMode == TAG_ONLY_HITS &&
									filteredSpectra.first.index.find( itr->first ) == filteredSpectra.first.index.end() )
								{
									deallocate( (*itr->second)->peakPreData );
									deallocate( (*itr->second)->peakData );
								} else if(	g_rtConfig->InlineValidationMode == TAG_ONLY_MISSES &&
											filteredSpectra.second.index.find( itr->first ) == filteredSpectra.first.index.end() )
								{
									deallocate( (*itr->second)->peakPreData );
									deallocate( (*itr->second)->peakData );
								}
						}

						spectra.filterByPeakCount();
						//spectra = originalSpectra;

						for( int z=0; z < g_rtConfig->NumChargeStates; ++z )
						{
							if( g_pid == 0 ) cout << "Threshold for " << g_rtConfig->Confidence * 100.0f << "% confidence in +" << z+1 << " IDs: " <<
									g_rtConfig->scoreThresholds[z] << "; " << validMatchCounts[z] << " of " << potentialMatchCounts[z] << " IDs pass." << endl;
						}

						if( g_pid == 0 ) cout << "All charge states: " << accumulate( validMatchCounts.begin(), validMatchCounts.end(), 0 ) <<
								" of " << accumulate( potentialMatchCounts.begin(), potentialMatchCounts.end(), 0 ) << " IDs pass." << endl;
						//cout << spectra.size() << " identifications pass confidence filter." << endl;
					}
				}
			}
		}*/
	}

	int ProcessHandler( int argc, vector<string> & argv, string &tags, vector<float> & lowPeakMzs, string & cachename)
	{
		simplethread_create_mutex( &resourceMutex );

		vector< string > args;
		for( int i=0; i < argc; ++i )
			args.push_back( argv[i] );

		if( InitProcess( args ) )
			return 1;

		SpectraList::InitMzFEBins();
		SpectraList::InitCEBins();

		INIT_PROFILERS(10)

		if( g_pid == 0 )
		{
			// The root process parses the XML data file and distributes spectra to child processes
			g_inputFilenames.clear();
			for( size_t i=1; i < args.size(); ++i )
			{
				//cout << g_hostString << " is reading spectra from files matching mask \"" << args[i] << "\"" << endl; //BX
				FindFilesByMask( args[i], g_inputFilenames );
			}

			if( g_inputFilenames.empty() )
			{
				cerr << g_hostString << " did not find any spectra matching given filemasks." << endl;
				return 1;
			}

			Timer overallTime(true);
			fileList_t finishedFiles;
			fileList_t::iterator fItr;
			for( fItr = g_inputFilenames.begin(); fItr != g_inputFilenames.end(); ++fItr )
			{
				spectra.clear();

				Timer fileTime(true);

				//cout << g_hostString << " is reading spectra from file \"" << *fItr << "\"" << endl; //BX
				finishedFiles.insert( *fItr );

				Timer readTime(true);
				spectra.readPeaks( *fItr, g_rtConfig->StartSpectraScanNum, g_rtConfig->EndSpectraScanNum );
				readTime.End();

				int totalPeakCount = 0;
				int numSpectra = (int) spectra.size();
				for( SpectraList::iterator sItr = spectra.begin(); sItr != spectra.end(); ++sItr )
					totalPeakCount += (*sItr)->peakPreCount;

				//cout << g_hostString << " read " << numSpectra << " spectra with " << totalPeakCount << " peaks; " << readTime.TimeElapsed() << " seconds elapsed." << endl; //BX

				int skip = 0;
				if( numSpectra == 0 )
				{
					cout << g_hostString << " is skipping a file with no spectra." << endl;
					skip = 1;
				}
				#ifdef USE_MPI
					for( int p=0; p < g_numChildren; ++p )
						MPI_Ssend( &skip,		1,		MPI_INT,	p+1, 0x00, MPI_COMM_WORLD );
				#endif

				Timer taggingTime;
				string startDate;
				string startTime;
				vector< size_t > opcs; // original peak count statistics
				vector< size_t > fpcs; // filtered peak count statistics

				if( !skip )
				{
					if( g_numProcesses > 1 )
					{
						#ifdef USE_MPI
							if( g_numChildren > 0 )
							{
								g_rtConfig->SpectraBatchSize = (int) ceil( (float) numSpectra / (float) g_numChildren / (float) g_rtConfig->NumBatches );
								cout << g_hostString << " calculates dynamic spectra batch size is " << g_rtConfig->SpectraBatchSize << endl;
							}

							//std::random_shuffle( spectra.begin(), spectra.end() );

							cout << g_hostString << " is sending spectra to worker nodes to prepare them for search." << endl;
							Timer prepareTime(true);
							TransmitUnpreparedSpectraToChildProcesses();

							deallocate( spectra );

							ReceivePreparedSpectraFromChildProcesses();

							SpectraList::PrecacheIRBins( spectra, cachename );

							numSpectra = (int) spectra.size();

							skip = 0;
							if( numSpectra == 0 )
							{
								cout << g_hostString << " is skipping a file with no suitable spectra." << endl;
								skip = 1;
							}

							for( int p=0; p < g_numChildren; ++p )
								MPI_Ssend( &skip,		1,		MPI_INT,	p+1, 0x00, MPI_COMM_WORLD );

							if( !skip )
							{
                                //BX: less output, too much clutter before
								opcs = spectra.getOriginalPeakCountStatistics();
								fpcs = spectra.getFilteredPeakCountStatistics();
//								cout << g_hostString << ": mean original (filtered) peak count: " <<
//										opcs[5] << " (" << fpcs[5] << ")" << endl;
//								cout << g_hostString << ": min/max original (filtered) peak count: " <<
//										opcs[0] << " (" << fpcs[0] << ") / " << opcs[1] << " (" << fpcs[1] << ")" << endl;
//								cout << g_hostString << ": original (filtered) peak count at 1st/2nd/3rd quartiles: " <<
//										opcs[2] << " (" << fpcs[2] << "), " <<
//										opcs[3] << " (" << fpcs[3] << "), " <<
//										opcs[4] << " (" << fpcs[4] << ")" << endl;

								float filter = 1.0f - ( (float) fpcs[5] / (float) opcs[5] );
//                              cout << g_hostString << " filtered out " << filter * 100.0f << "% of peaks." << endl;
//
//								cout << g_hostString << " has " << numSpectra << " spectra prepared now; " << prepareTime.End() << " seconds elapsed." << endl;

								ReadInlineValidationFile();

//								cout << g_hostString << " is sending " << spectra.size() << " spectra to worker nodes for sequence tagging." << endl;
								startTime = GetTimeString(); startDate = GetDateString(); taggingTime.Begin();
								TransmitUntaggedSpectraToChildProcesses();
//								cout << g_hostString << " has sequence tagged all its spectra; " << taggingTime.End() << " seconds elapsed." << endl;

								deallocate( spectra );

//								cout << g_hostString << " is receiving tag results from worker nodes." << endl;
								Timer resultsTime(true);
								ReceiveTaggedSpectraFromChildProcesses();
								//cout << g_hostString << " finished receiving tag results for " << spectra.size() << " spectra; " << resultsTime.End() << " seconds elapsed." << endl;
							
                            }

						#endif
					} else
					{
						spectra.random_shuffle();

						PrepareSpectra();

						skip = 0;
						if( spectra.size() == 0 )
						{
							cout << g_hostString << " is skipping a file with no suitable spectra." << endl;
							skip = 1;
						}

						if( !skip )
						{
							opcs = spectra.getOriginalPeakCountStatistics();
							fpcs = spectra.getFilteredPeakCountStatistics();
//							cout << g_hostString << ": mean original (filtered) peak count: " <<
//									opcs[5] << " (" << fpcs[5] << ")" << endl;
//							cout << g_hostString << ": min/max original (filtered) peak count: " <<
//									opcs[0] << " (" << fpcs[0] << ") / " << opcs[1] << " (" << fpcs[1] << ")" << endl;
//							cout << g_hostString << ": original (filtered) peak count at 1st/2nd/3rd quartiles: " <<
//									opcs[2] << " (" << fpcs[2] << "), " <<
//									opcs[3] << " (" << fpcs[3] << "), " <<
//									opcs[4] << " (" << fpcs[4] << ")" << endl;
//BX
							float filter = 1.0f - ( (float) fpcs[5] / (float) opcs[5] );
//BX							cout << g_hostString << " filtered out " << filter * 100.0f << "% of peaks." << endl;

							ReadInlineValidationFile();

							SpectraList::PrecacheIRBins( spectra, cachename );

//BX							cout << g_hostString << " is sequence tagging " << spectra.size() << " spectra." << endl;
							startTime = GetTimeString(); startDate = GetDateString(); taggingTime.Begin();
							taggingStats sumTaggingStats = ExecutePipeline();
							//EncodeSpectraForOutput();

//BX							cout << g_hostString << " has sequence tagged all its spectra; " << taggingTime.End() << " seconds elapsed." << endl;

//BX							cout << g_hostString << " stats: " << 1 << " / " <<
//									sumTaggingStats.numSpectraTagged << " / " <<
//									sumTaggingStats.numResidueMassGaps << " / " <<
//									sumTaggingStats.numTagsGenerated << " / " <<
//									sumTaggingStats.numTagsRetained << endl;
						}
					}

					if( !skip )
					{
						try
						{
							spectra.sort( spectraSortByID() );
							WriteTagsToTagsFile( *fItr, startTime, startDate, taggingTime.End(), tags,lowPeakMzs );
						//	cout << g_hostString << " finished file \"" << *fItr << "\"; " << fileTime.End() << " seconds elapsed." << endl;
						} catch( ... )
						{
							cerr << "Error while sorting and writing XML output." << endl;
							exit(1);
						}
					}
				}
			}

			#ifdef USE_MPI
					int done = ( ( g_inputFilenames.size() - finishedFiles.size() ) == 0 ? 1 : 0 );
					for( int p=0; p < g_numChildren; ++p )
						MPI_Ssend( &done,		1,		MPI_INT,	p+1, 0x00, MPI_COMM_WORLD );
			#endif

//BX			cout << g_hostString << " sequence tagged spectra from " << g_inputFilenames.size() << " files; " << overallTime.End() << " seconds elapsed." << endl;
		}
		#ifdef USE_MPI
			else
			{
				int allDone = 0;

				while( !allDone )
				{
					int skip;
					MPI_Recv( &skip,	1,		MPI_INT,	0,	0x00, MPI_COMM_WORLD, &st );

					if( !skip )
					{
						SpectraList preparedSpectra;

						while( ReceiveUnpreparedSpectraBatchFromRootProcess() )
						{
							PrepareSpectra();
							preparedSpectra.insert( spectra.begin(), spectra.end(), preparedSpectra.end() );
							spectra.clear( false );
						}

						TransmitPreparedSpectraToRootProcess( preparedSpectra );
						deallocate( preparedSpectra );

						MPI_Recv( &skip,	1,		MPI_INT,	0,	0x00, MPI_COMM_WORLD, &st );

						if( !skip )
						{
							SpectraList taggedSpectra;
							SpectraList::PrecacheIRBins( spectra, cachename );

							int numBatches = 0;
							taggingStats sumTaggingStats;
							taggingStats lastTaggingStats;
							while( ReceiveUntaggedSpectraBatchFromRootProcess() )
							{
								++ numBatches;

								lastTaggingStats = ExecutePipeline();
								sumTaggingStats = sumTaggingStats + lastTaggingStats;

								taggedSpectra.insert( spectra.begin(), spectra.end(), taggedSpectra.end() );
								spectra.clear( false );
							}

							cout << g_hostString << " stats: " << numBatches << " / " <<
									sumTaggingStats.numSpectraTagged << " / " <<
									sumTaggingStats.numResidueMassGaps << " / " <<
									sumTaggingStats.numTagsGenerated << " / " <<
									sumTaggingStats.numTagsRetained << endl;

							TransmitTaggedSpectraToRootProcess( taggedSpectra );
							taggedSpectra.clear();
						}
					}

					MPI_Recv( &allDone,	1,		MPI_INT,	0,	0x00, MPI_COMM_WORLD, &st );
				}
			}
		#endif

		return 0;
	}

	gapMap_t::iterator FindPeakNear( gapMap_t& peakData, float mz, float tolerance )
	{
		gapMap_t::iterator cur, min, max, best;

		min = peakData.lower_bound( mz - tolerance );
		max = peakData.lower_bound( mz + tolerance );

		if( min == max )
			return peakData.end(); // no peaks

		// find the peak closest to the desired mz
		best = min;
		float minDiff = (float) fabs( mz - best->first );
		for( cur = min; cur != max; ++cur )
		{
			float curDiff = (float) fabs( mz - cur->first );
			if( curDiff < minDiff )
			{
				minDiff = curDiff;
				best = cur;
			}
		}

		return best;
	}
}
}
 
//begin changes BXu
//try to convert directag into something call-able
int directagFunc(int argc, vector<string>& argv, string & tags, vector<float> & lowPeakMzs, string & cachename){
	try{
//		char buf[256];
//		GetHostname( buf, sizeof(buf) );
//		g_numProcesses = 1;
//		g_pid = 0;
//		g_numChildren = g_numProcesses - 1;
//		
//		ostringstream str;
//		str << "Process #" << g_pid << " (" << buf << ")";
//		g_hostString = str.str();
		int result = directag::ProcessHandler( argc, argv , tags, lowPeakMzs, cachename);
		return result;
		
		
	}
	catch( exception& e)
	{
        throw e;
    }
	return 1;
	
}



//end changes Bxu

//int main( int argc, const char* argv[] )
//{
//	try
//	{
//		char buf[256];
//		GetHostname( buf, sizeof(buf) );
//
//		// Initialize the message passing interface for the parallel processing system
//		#ifdef MPI_DEBUG
//			cout << buf << " is initializing MPI... " << endl;
//		#endif
//
//		#ifdef USE_MPI
//			int threadLevel;
//			MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &threadLevel );
//			if( threadLevel < MPI_THREAD_SINGLE )
//			{
//				cerr << "MPI library is not thread compliant: " << threadLevel << " should be " << MPI_THREAD_MULTIPLE << endl;
//				return 1;
//			}
//			MPI_Buffer_attach( malloc( MPI_BUFFER_SIZE ), MPI_BUFFER_SIZE );
//			//CommitCommonDatatypes();
//		#endif
//
//		#ifdef MPI_DEBUG
//			cout << buf << " has initialized MPI... " << endl;
//		#endif
//
//		// Get information on the MPI environment
//		#ifdef MPI_DEBUG
//			cout << buf << " is gathering MPI information... " << endl;
//		#endif
//
//		#ifdef USE_MPI
//			MPI_Comm_size( MPI_COMM_WORLD, &g_numProcesses );
//			MPI_Comm_rank( MPI_COMM_WORLD, &g_pid );
//		#else
//			g_numProcesses = 1;
//			g_pid = 0;
//		#endif
//
//		g_numChildren = g_numProcesses - 1;
//
//		ostringstream str;
//		str << "Process #" << g_pid << " (" << buf << ")";
//		g_hostString = str.str();
//
//		#ifdef MPI_DEBUG
//			cout << g_hostString << " has gathered its MPI information." << endl;
//		#endif
//
//
//		// Process the data
//		#ifdef MPI_DEBUG
//			cout << g_hostString << " is starting." << endl;
//		#endif
//
//		int result = directag::ProcessHandler( argc, argv );
//
//		#ifdef USE_MPI
//			if( g_pid == 0 && g_numChildren > 0 && result > 0 )
//				MPI_Abort( MPI_COMM_WORLD, 1 );
//		#endif
//
//		#ifdef MPI_DEBUG
//			cout << g_hostString << " has finished." << endl;	
//		#endif
//
//		// Destroy the message passing interface
//		#ifdef MPI_DEBUG
//			cout << g_hostString << " is finalizing MPI... " << endl;
//		#endif
//
//		#ifdef USE_MPI
//			int size;
//			MPI_Buffer_detach( &g_mpiBuffer, &size );
//			free( g_mpiBuffer );
//			MPI_Finalize();
//		#endif
//
//		#ifdef MPI_DEBUG
//			cout << g_hostString << " is terminating." << endl;
//		#endif
//
//		return result;
//
//	} catch( exception& e )
//	{
//		cerr << e.what() << endl;
//	}
//
//	return 1;
//}

