#include "stdafx.h"
#include "directagSpectrum.h"

using namespace freicore;

namespace std
{
	ostream& operator<< ( ostream& o, const freicore::directag::PeakInfo& rhs )
	{
		return o << "( " << rhs.intensityRank << " )";
	}

	ostream& operator<< ( ostream& o, const freicore::directag::GapInfo& rhs )
	{
		return o << "(GapInfo: " << rhs.fromPeakItr->first << " " << rhs.peakItr->first << " " << rhs.gapMass << " " << rhs.gapRes << " " << rhs.error << " )";
	}

	ostream& operator<< ( ostream& o, const freicore::directag::gapVector_t& rhs )
	{
		o << "(GapVector:";
		for( freicore::directag::gapVector_t::const_iterator itr = rhs.begin(); itr != rhs.end(); ++itr )
			 o << " " << *itr;
		o << " )";

		return o;
	}

	ostream& operator<< ( ostream& o, const freicore::directag::gapMap_t& rhs )
	{
		o << "(GapMap:";
		for( freicore::directag::gapMap_t::const_iterator itr = rhs.begin(); itr != rhs.end(); ++itr )
			o << " " << itr->first << "->" << itr->second << "\n";
		o << " )";

		return o;
	}
}

namespace freicore
{
namespace directag
{
	MzFEBins SpectraList::mzFidelityErrorBins;
	CEBinsList SpectraList::complementErrorBinsList;
	IRBinsTable SpectraList::intensityRanksumBinsTable;

	Spectrum::Spectrum()
		:	BaseSpectrum(), PeakSpectrum< PeakInfo >(), TaggingSpectrum(), SearchSpectrum< SearchResult >()
	{
		tagList.max_size( g_rtConfig->MaxTagCount );

		scoreWeights[ "intensity" ] = intensityScoreWeight = g_rtConfig->IntensityScoreWeight;
		scoreWeights[ "mzFidelity" ] = mzFidelityScoreWeight = g_rtConfig->MzFidelityScoreWeight;
		scoreWeights[ "complement" ] = complementScoreWeight = g_rtConfig->ComplementScoreWeight;
		scoreWeights[ "random" ] = g_rtConfig->RandomScoreWeight;
	}

	Spectrum::Spectrum( const Spectrum& old )
		:	BaseSpectrum( old ), PeakSpectrum< PeakInfo >( old ), TaggingSpectrum( old ), SearchSpectrum< SearchResult >( old )
	{
		scoreWeights[ "intensity" ] = intensityScoreWeight = g_rtConfig->IntensityScoreWeight;
		scoreWeights[ "mzFidelity" ] = mzFidelityScoreWeight = g_rtConfig->MzFidelityScoreWeight;
		scoreWeights[ "complement" ] = complementScoreWeight = g_rtConfig->ComplementScoreWeight;
		scoreWeights[ "random" ] = g_rtConfig->RandomScoreWeight;
	}

	void Spectrum::Parse( bool intenAsClasses )
	{
		//cout << "Precursor Ion MZ: " <<mzOfPrecursor <<
		//		"; time in minutes from start of run: " << retentionTime << endl;
		//cout << "mzData: " << mzData << endl << "iData: " << iData << endl << endl;

		vector< float > mzDataVec, iDataVec;

		if( mzDataIsDoublePrecision )
			base64_decode_sequence< float, double >( mzData, mzDataVec, g_endianType != mzDataEndianType );
		else
			base64_decode_sequence< float, float >( mzData, mzDataVec, g_endianType != mzDataEndianType );

		deallocate( mzData );

		if( mzDataVec.empty() )
			return;

		//for( int p=0; p < min(10, (int) mzData.size()); ++p )
		//	cout << mzData[p] << " ";
		//cout << endl;

		if( !intenAsClasses )
		{
			if( iDataIsDoublePrecision )
				base64_decode_sequence< float, double >( iData, iDataVec, g_endianType != iDataEndianType );
			else
				base64_decode_sequence< float, float >( iData, iDataVec, g_endianType != iDataEndianType );

			//mzData.clear();
			//iData.clear();
			deallocate( iData );

			for( size_t j=0; j < mzDataVec.size(); ++j )
				peakPreData.insert( peakPreData.end(), pair< float, float >( mzDataVec[j], iDataVec[j] ) );

			totalPeakSpace = peakPreData.rbegin()->first - peakPreData.begin()->first;
		} else
		{
			char tmp[2];
			tmp[1] = 0;
			for( size_t j=0; j < mzDataVec.size(); ++j )
			{
				tmp[0] = iData[j];
				peakData[ mzDataVec[j] ].intensityRank = atoi(tmp);
			}

			totalPeakSpace = peakData.rbegin()->first - peakData.begin()->first;
		}
		//cout << id.index << ": " << totalPeakSpace << endl;
		//cout << "Number of peaks: " << (int) spectra[i]->peakData.size() << endl;
	}

	void Spectrum::ClassifyPeakIntensities()
	{
		// Sort peaks by intensity.
		// Use multimap because multiple peaks can have the same intensity.
		typedef multimap< float, float > IntenSortedPeakPreData;
		IntenSortedPeakPreData intenSortedPeakPreData;
		for( PeakPreData::iterator itr = peakPreData.begin(); itr != peakPreData.end(); ++itr )
		{
			IntenSortedPeakPreData::iterator iItr = intenSortedPeakPreData.insert( pair< float, float >( itr->second, itr->second ) );
			iItr->second = itr->first;
		}

		// Restore the sorting order to be based on MZ
		IntenSortedPeakPreData::reverse_iterator r_iItr = intenSortedPeakPreData.rbegin();
		peakPreData.clear();
		peakData.clear();

		for( int i=0; i < g_rtConfig->NumIntensityClasses; ++i )
		{
			int numFragments = (int) round( (float) ( pow( (float) g_rtConfig->ClassSizeMultiplier, i ) * intenSortedPeakPreData.size() ) / (float) g_rtConfig->minIntensityClassCount, 0 );
			for( int j=0; r_iItr != intenSortedPeakPreData.rend() && j < numFragments; ++j, ++r_iItr )
			{
				float mz = r_iItr->second;
				float inten = r_iItr->first;
				peakPreData.insert( peakPreData.end(), pair<float, float>( mz, inten ) );
				peakData[ mz ].intensityRank = peakPreData.size();
			}
		}
		intenSortedPeakPreData.clear();
	}

	// Attempts to find a complement for each peak in the spectrum
	// Returns the sum of products of the found complements' intensities
	float Spectrum::FindComplements( float complementMzTolerance )
	{
		float sumOfProducts = 0;
		for( PeakPreData::iterator itr = peakPreData.begin(); itr != peakPreData.end(); ++itr )
		{
			PeakInfo& peak = peakData[ itr->first ];// = PeakInfo();

			for( int z=0; z < numFragmentChargeStates; ++z )
			{
				float complementMz = CalculateComplementMz( itr->first, z+1 );
				PeakPreData::iterator complementItr = peakPreData.findNear( complementMz, complementMzTolerance, true );
				if( complementItr != peakPreData.end() )
				{
					sumOfProducts += itr->second * complementItr->second;
					peak.hasComplementAsCharge[z] = true;
				} else
					peak.hasComplementAsCharge[z] = false;
			}
		}

		return sumOfProducts;
	}

	size_t Spectrum::MakeTagGraph()
	{
		PeakData::iterator left;	// main iterator pointing to the first peak in a comparison
		PeakData::iterator cur;		// main iterator pointing to the peak currently being looked at
		m2n_t::const_iterator resItr;
		size_t numResidueMassGaps = 0;

		gapMaps.clear();
		tagGraphs.clear();
		nodeSet.clear();

		gapMaps.resize( numFragmentChargeStates );
		tagGraphs.resize( numFragmentChargeStates );

		for( int z=0; z < numFragmentChargeStates; ++z )
		{
			gapMap_t& gapMap = gapMaps[z];
			spectrumGraph& tagGraph = tagGraphs[z];

			for( left = peakData.begin(); left != peakData.end(); ++left )
			{
				for( resItr = g_residueMap->beginMonoMasses(); resItr != g_residueMap->endMonoMasses(); ++resItr )
				{
					if( resItr->second == PEPTIDE_N_TERMINUS_SYMBOL || resItr->second == PEPTIDE_C_TERMINUS_SYMBOL )
						continue;

					float mzGap = resItr->first / (float) (z+1);
					float expectedMZ = left->first + mzGap;
					cur = peakData.findNear( expectedMZ, g_rtConfig->FragmentMzTolerance );

					if( cur != peakData.end() )
					{
						// Calculate the error between the m/z of the actual peak and the m/z that was expected for it
						float error = (cur->first - left->first) - mzGap;
						if( fabs( error ) > g_rtConfig->FragmentMzTolerance )
							continue;

						gapMap_t::iterator nextGapInfo = gapMap.insert( gapMap_t::value_type( cur->first, gapVector_t() ) ).first;
						gapMap[ left->first ].push_back( GapInfo( left, cur, nextGapInfo, mzGap, resItr->second, error ) );
						++ numResidueMassGaps;

						GapInfo newEdge( left, cur, nextGapInfo, mzGap, resItr->second, error, left->first, cur->first );
						tagGraph[ left->first ].cEdges.push_back( newEdge );
						tagGraph[ cur->first ].nEdges.push_back( newEdge );
						tagGraph[ cur->first ].nPathSize = max(	tagGraph[ cur->first ].nPathSize,
																tagGraph[ left->first ].nPathSize + 1 );
						nodeSet.insert( left->first );
						nodeSet.insert( cur->first );
					}
				}
			}

			for( spectrumGraph::reverse_iterator itr = tagGraph.rbegin(); itr != tagGraph.rend(); ++itr )
			{
				for( size_t i=0; i < itr->second.nEdges.size(); ++i )
				{
					tagGraph[ itr->second.nEdges[i].nTermMz ].cPathSize = max(	tagGraph[ itr->second.nEdges[i].nTermMz ].cPathSize,
																				itr->second.cPathSize + 1 );
				}

				itr->second.longestPath = itr->second.cPathSize + itr->second.nPathSize;
			}
		}

        tagGraphTIC = 0;
        for( nodeSet_t::iterator itr = nodeSet.begin(); itr != nodeSet.end(); ++itr )
            tagGraphTIC += peakPreData[*itr];

		return numResidueMassGaps;
	}

	void Spectrum::FilterPeaks()
	{
		// Secondly, determine the neutral mass of the precursor (m/z * z - z)
		mOfPrecursor = mzOfPrecursor * id.charge - ( id.charge * PROTON );

		numFragmentChargeStates = max( 1, id.charge - 1 );

		if( peakPreData.empty() )
			return;

		// Eliminate peaks above the precursor's mass with a given tolerance
		float maxPeakMass = mOfPrecursor + PROTON + g_rtConfig->PrecursorMzTolerance;
		PeakPreData::iterator itr = peakPreData.upper_bound( maxPeakMass );
		peakPreData.erase( itr, peakPreData.end() );

		if( peakPreData.empty() )
			return;

		// Thirdly, store the bounds of the spectrum before eliminating any peaks
		mzLowerBound = peakPreData.begin()->first;
		mzUpperBound = peakPreData.rbegin()->first;
		totalPeakSpace = mzUpperBound - mzLowerBound;

		if( g_rtConfig->MakeSpectrumGraphs )
			writeToSvgFile( "-unprocessed" + g_rtConfig->OutputSuffix );

		if( g_rtConfig->DeisotopingMode > 0 || g_rtConfig->AdjustPrecursorMass )
		{
			Deisotope( g_rtConfig->IsotopeMzTolerance );

			if( g_rtConfig->MakeSpectrumGraphs )
				writeToSvgFile( "-deisotoped" + g_rtConfig->OutputSuffix );
		}

		FilterByTIC( g_rtConfig->TicCutoffPercentage );
		FilterByPeakCount( g_rtConfig->MaxPeakCount );

		if( g_rtConfig->MakeSpectrumGraphs )
			writeToSvgFile( "-filtered" + g_rtConfig->OutputSuffix );

		if( peakPreData.size() < (size_t) g_rtConfig->minIntensityClassCount )
		{
			peakPreData.clear();
			return;
		}

		// Create peak data from pre peak data
		/*ClassifyPeakIntensities();

		MakeTagGraph();
		//map< int, int > pathSizeHistogram;
		int maxLongestPath = 0;
		map< float, int > longestPathMap;

		for( int z=0; z < numFragmentChargeStates; ++z )
		{
			spectrumGraph& tagGraph = tagGraphs[z];
			for( spectrumGraph::reverse_iterator itr = tagGraph.rbegin(); itr != tagGraph.rend(); ++itr )
			{
				if( longestPathMap[ itr->first ] < itr->second.longestPath )
					longestPathMap[ itr->first ] = itr->second.longestPath;

				if( itr->second.longestPath > maxLongestPath )
					maxLongestPath = itr->second.longestPath;

				//++ pathSizeHistogram[ itr->second.longestPath ];
			}
			deallocate(tagGraph);
		}
		//cout << id << " peak path histogram:\n" << pathSizeHistogram << endl;

		vector< PeakData::iterator > junkPeaks;
		for( PeakData::iterator cur = peakData.begin(); cur != peakData.end(); ++cur )
		{
			if( longestPathMap[ cur->first ] < g_rtConfig->tagPeakCount )
				junkPeaks.push_back( cur );
		}

		for( size_t i=0; i < junkPeaks.size(); ++i )
		{
			peakData.erase( junkPeaks[i] );
			peakPreData.erase( junkPeaks[i]->first );
		}

		if( peakData.size() < (size_t) g_rtConfig->minIntensityClassCount )
		{
			deallocate(peakPreData);
			deallocate(peakData);
			return;
		}*/

		// Create peak data from pre peak data
		ClassifyPeakIntensities();

		totalPeakSpace = peakData.rbegin()->first - peakData.begin()->first;
		peakCount = (int) peakData.size();
	}

	void Spectrum::Preprocess()
	{
		PeakPreData::iterator itr;
		PeakPreData::reverse_iterator r_itr;
		PeakPreData::iterator findItr;

		if( mzOfPrecursor < 1 )
		{
			peakPreData.clear();
			return;
		}

		if( g_rtConfig->AdjustPrecursorMass )
		{
			float originalPrecursorMass = mOfPrecursor;
			float originalPrecursorMz = mzOfPrecursor;
			float bestPrecursorAdjustment = 0.0f;
			float maxSumOfProducts = 0.0f;
			map<float, float> AdjustmentResults;

			for( mOfPrecursor += g_rtConfig->MinPrecursorAdjustment;
				 mOfPrecursor <= originalPrecursorMass + g_rtConfig->MaxPrecursorAdjustment;
				 mOfPrecursor += g_rtConfig->PrecursorAdjustmentStep )
			{
				mzOfPrecursor = ( mOfPrecursor + ( id.charge * PROTON ) ) / id.charge;

				float sumOfProducts = FindComplements( g_rtConfig->ComplementMzTolerance );

				if( sumOfProducts > maxSumOfProducts )
				{
					maxSumOfProducts = sumOfProducts;
					bestPrecursorAdjustment = mOfPrecursor - originalPrecursorMass;
				}

				AdjustmentResults[ mOfPrecursor ] = sumOfProducts;
			}

			if( maxSumOfProducts > 0.0f )
			{
				mOfPrecursor = originalPrecursorMass + bestPrecursorAdjustment;
				mzOfPrecursor = ( mOfPrecursor + ( id.charge * PROTON ) ) / id.charge;
			} else
			{
				mOfPrecursor = originalPrecursorMass;
				mzOfPrecursor = originalPrecursorMz;
			}

			if( g_rtConfig->MakeSpectrumGraphs )
			{
				writeToSvgFile( "-adjusted" + g_rtConfig->OutputSuffix );
				cout << "Original precursor m/z: " << originalPrecursorMz << endl;
				cout << "Corrected precursor m/z: " << mzOfPrecursor << endl;
				cout << "Sum of complement products: " << maxSumOfProducts << endl;

				/*cout << "Best complement total: " << BestComplementTotal << endl;
				cout << oldPrecursor << " (" << spectrum->mOfPrecursorFixed << ") corrected by " << spectrum->mzOfPrecursor - oldPrecursor <<
						" to " << spectrum->mzOfPrecursor << " (" << spectrum->mOfPrecursor << ") " << endl;*/

				cout << AdjustmentResults << endl;
			}
		}

		// Initialize the spectrum info tables
		initialize( g_rtConfig->NumIntensityClasses, g_rtConfig->NumMzFidelityClasses );

		if( peakData.size() < (size_t) g_rtConfig->minIntensityClassCount )
		{
			peakPreData.clear();
			peakData.clear();
			return;
		}

		// Reclassify intensities and find complements based on fully processed spectrum
		ClassifyPeakIntensities();
		FindComplements( g_rtConfig->ComplementMzTolerance );

		totalPeakSpace = peakData.rbegin()->first - peakData.begin()->first;
		peakCount = (int) peakData.size();
	}

	void Spectrum::MakeProbabilityTables()
	{
		for( PeakData::iterator itr = peakData.begin(); itr != peakData.end(); ++itr )
		{
			itr->second.hasSomeComplement = accumulate( itr->second.hasComplementAsCharge.begin(),
														itr->second.hasComplementAsCharge.begin() + numFragmentChargeStates,
														0 );
			++ complementClassCounts[ ( itr->second.hasSomeComplement == 0 ? 1 : 0 ) ];
		}

		if( complementClassCounts[0] == 0 )
		{
			//scoreWeights["complement"] = 0;
			complementScoreWeight = 0;
		} else
			CreateScoringTableMVH( 0, g_rtConfig->tagPeakCount, 2, complementClassCounts, bgComplements, g_lnFactorialTable, false, false, false );
	}

	size_t Spectrum::Score()
	{
		START_PROFILER(3)
		size_t numTagsGenerated = findTags();
		STOP_PROFILER(3)

		for( TagList::iterator itr = interimTagList.begin(); itr != interimTagList.end(); ++itr )
		{
			TagInfo& tag = const_cast< TagInfo& >( *itr );
			//tag.CalculateTotal( scoreWeights );
			START_PROFILER(4)
			tag.CalculateTotal( complementScoreWeight, intensityScoreWeight, mzFidelityScoreWeight );
			STOP_PROFILER(4)
			tag.totalScore *= numTagsGenerated;
			//for( map< string, float >::iterator scoreItr = tag.scores.begin(); scoreItr != tag.scores.end(); ++scoreItr )
			//	scoreItr->second *= numTagsGenerated;
			START_PROFILER(5)
			if( g_rtConfig->MaxTagScore == 0 || tag.totalScore <= g_rtConfig->MaxTagScore )
				tagList.insert( tag );
			STOP_PROFILER(5)
		}

		for( TagList::iterator itr = validTagList.begin(); itr != validTagList.end(); ++itr )
		{
			TagInfo& tag = const_cast< TagInfo& >( *itr );
			//tag.CalculateTotal( scoreWeights );
			tag.CalculateTotal( complementScoreWeight, intensityScoreWeight, mzFidelityScoreWeight );
			tag.totalScore *= numTagsGenerated;
			//for( map< string, float >::iterator scoreItr = tag.scores.begin(); scoreItr != tag.scores.end(); ++scoreItr )
			//	scoreItr->second *= numTagsGenerated;
			if( !g_rtConfig->AlwaysKeepValidTags )
			{
				if( g_rtConfig->MaxTagScore == 0 || tag.totalScore <= g_rtConfig->MaxTagScore )
					tagList.insert( tag );
			} else
				tagList.insert( tag, true );
		}

		START_PROFILER(6)
		deallocate( validTagList );
		deallocate( interimTagList );

		deallocate( bgComplements );
		STOP_PROFILER(6)

		return numTagsGenerated;
	}

	// Takes a tag and recursively fills a list of strings of variations of that tag based on I/L substitutions
	void TagExploder_R( const string& tag, int idx, vector< string >& tagList )
	{
		if( idx == (int) tag.length() )
		{
			tagList.push_back( tag );
			return;
		}

		if( tag[idx] == 'I' )
		{
			string newTag( tag );
			newTag[idx] = 'L';
			TagExploder_R( newTag, idx+1, tagList );
		}

		TagExploder_R( tag, idx+1, tagList );
	}

	void TagExploder( const string& tag, vector< string >& tagList )
	{
		tagList.push_back( tag );
		TagExploder_R( tag, 0, tagList );
	}

	void Spectrum::findTags_R(	gapMap_t::iterator gapInfoItr,
								int tagIndex,
								string& tag,
								vector< float >& peakErrors,
								vector< PeakData::iterator >& peakList,
								int peakChargeState,
								size_t& numTagsGenerated,
								IRBins& irBins )
	{
		if( tagIndex == 0 )
		{
			++ numTagsGenerated;

			TagInfo newTag;

			MvIntKey intensityClassKey, complementClassKey, mzFidelityKey;
			intensityClassKey.resize( g_rtConfig->NumIntensityClasses, 0 );
			complementClassKey.resize( 2, 0 );

			vector<float> modelPeaks;
			vector<float> modelErrors;
			vector<float> modelSquaredErrors;

			float gapMass = 0.0f;
			modelPeaks.push_back( peakList[0]->first );
			for( int i=0; i < g_rtConfig->TagLength; ++i )
			{
				gapMass += g_residueMap->getMonoMassByName( tag[i] ) / (float) peakChargeState;
				modelPeaks.push_back( peakList[i+1]->first - gapMass );
			}
			float sum = accumulate( modelPeaks.begin(), modelPeaks.end(), 0.0f );
			float avg = sum / g_rtConfig->tagPeakCount;

			sum = 0.0f;
			for( int i=0; i < g_rtConfig->tagPeakCount; ++i )
			{
				modelErrors.push_back( fabs( modelPeaks[i] - avg ) );
				sum += pow( modelErrors[i], 2 );
				//cout << e1 << " " << e2 << " " << e3 << ": " << errors << endl;
			}

			MzFEBins::iterator binItr = SpectraList::mzFidelityErrorBins.upper_bound( sum );
			-- binItr;
			//newTag.scores[ "mzFidelity" ] = binItr->second;
			newTag.mzFidelityScore = binItr->second;

			gapMass = 0.0f;
			modelPeaks[0] = avg;
			for( int i=0; i < g_rtConfig->TagLength; ++i )
			{
				gapMass += g_residueMap->getMonoMassByName( tag[i] ) / (float) peakChargeState;
				modelPeaks[i+1] = avg + gapMass;
			}
			//cout << peakList << " " << modelPeaks << " " << modelErrors << endl;

			//int totalPathLength = 0;
			int totalIntensityRanks = 1;
			//int totalContextRanks = 1;
			vector< float > complementPairMasses;
			//spectrumGraph& tagGraph = tagGraphs[peakChargeState-1];
			for( int i=0; i < g_rtConfig->tagPeakCount; ++i )
			{
				PeakInfo& peak = peakList[i]->second;

				newTag.worstPeakRank = max( peak.intensityRank, newTag.worstPeakRank );
				totalIntensityRanks += peak.intensityRank;

				bool hasComplement = peak.hasComplementAsCharge[ peakChargeState-1 ];
				++ complementClassKey[ hasComplement ? 0 : 1 ];
				if( hasComplement )
				{
					float complementMz = CalculateComplementMz( peakList[i]->first, peakChargeState );
					PeakData::iterator complementItr = peakData.findNear( complementMz, g_rtConfig->ComplementMzTolerance );
					complementPairMasses.push_back( peakList[i]->first + complementItr->first );
				}
			}

			//newTag.scores[ "intensity" ] = irBins[ totalIntensityRanks ];
			newTag.intensityScore = irBins[ totalIntensityRanks ];
			newTag.ranksum = totalIntensityRanks;

			float complementClassScore = 0;
			if( complementClassCounts[0] > 0 )
			{
				CEBins::iterator binItr = SpectraList::complementErrorBinsList[2].begin();
				if( complementClassKey[0] > 1 )
				{
					float complementPairMean = arithmetic_mean<float>( complementPairMasses );
					for( size_t i=0; i < complementPairMasses.size(); ++i )
						complementPairMasses[i] = pow( complementPairMasses[i] - complementPairMean, 2.0f );
					float sse = accumulate( complementPairMasses.begin(), complementPairMasses.end(), 0.0f );

					binItr = SpectraList::complementErrorBinsList[complementClassKey[0]].upper_bound( sse );
					-- binItr;
					while(binItr->second == 0)
						++ binItr;
				}
				MvhTable::reverse_iterator itr;
				int i = g_rtConfig->tagPeakCount;
				for( itr = bgComplements.rbegin(); itr != bgComplements.rend() && i >= (int) complementPairMasses.size(); ++itr, --i )
					complementClassScore += (float) exp(itr->second);
				--itr;
				if( i > 1 )
					complementClassScore -= (float) exp(itr->second) * ( 1.0f - binItr->second );
				else
					complementClassScore -= (float) exp(itr->second) / 2.0f;
				//newTag.scores[ "complement" ] = complementClassScore;
				newTag.complementScore = complementClassScore;
				//cout << id.index << ": " << complementClassKey << " " << complementClassCounts << " " << complementClassScore << " " << complementPairMasses << " " << binItr->second << " " << i << " " << itr->second << endl;
			} else
				//newTag.scores[ "complement" ] = 1.0f;
				newTag.complementScore = 1.0f;

			//newNode.peakList = peakList;

			if( g_rtConfig->RandomScoreWeight != 0 )
				newTag.scores[ "random" ] = (float) g_rtConfig->GetRandomScore();

			newTag.lowPeakMz = peakList.front()->first;

			//----------------------------------------- lower y - water+proton 
			//newNode->cTerminusMass = max( 0.0f, *peakList.begin() - WATER + PROTON );
			newTag.cTerminusMass = modelPeaks.front() * peakChargeState - ( PROTON * peakChargeState );
			newTag.cTerminusMass = max( 0.0f, newTag.cTerminusMass - WATER_MONO );

			//---------------------------- neutral precursor - proton ----- upper y
			//newNode->nTerminusMass = max( 0.0f, mOfPrecursor + 1 - *peakList.rbegin() );
			newTag.nTerminusMass = modelPeaks.back() * peakChargeState - ( PROTON * peakChargeState );
			newTag.nTerminusMass = max( 0.0f, mOfPrecursor - newTag.nTerminusMass );

			string properTag = tag;
			std::reverse( properTag.begin(), properTag.end() );

			newTag.tag = properTag;
			newTag.totalScore = (float) tagCount;

			if( resultSet.empty() || g_rtConfig->InlineValidationFile.empty() )
			{
				interimTagList.insert( newTag );
				
				++ tagCount;
			} else
			{
				++ tagCount;

				for(	Spectrum::SearchResultSetType::reverse_iterator rItr = resultSet.rbegin();
						rItr != resultSet.rend() && rItr->rank == 1;
						++rItr )
				{
					string seq = rItr->sequence;

					// Make a list of tag combinations expanding 'I' to 'I' and 'L'
					vector< string > tagVariants;
					TagExploder( newTag.tag, tagVariants );
					//cout << seq << endl;
					bool anyValidVariant = false;
					for( size_t i=0; i < tagVariants.size(); ++i )
					{
						newTag.tag = tagVariants[i];
						newTag.valid = false;

						//cout << properTag << endl;
						size_t offset = 0;
						size_t length = tagVariants[i].length();

						if( ( offset = seq.find( newTag.tag ) ) != string::npos )
						{
							float nTerminusMass = g_rtConfig->inlineValidationResidues.GetMassOfResidues( seq.substr( 0, offset ), g_rtConfig->UseAvgMassOfSequences );
							float cTerminusMass = g_rtConfig->inlineValidationResidues.GetMassOfResidues( seq.substr( offset+length ), g_rtConfig->UseAvgMassOfSequences );
							//cout << nTerminusMass << " " << nT << endl;
							//cout << cTerminusMass << " " << cT << endl;
							if( fabs( nTerminusMass - newTag.nTerminusMass ) < g_rtConfig->NTerminusMassTolerance &&
								fabs( cTerminusMass - newTag.cTerminusMass ) < g_rtConfig->CTerminusMassTolerance )
							{
									//spectrumHasValidTag = true;
									newTag.valid = true;
									anyValidVariant = true;
							}
						}

						if( newTag.valid )
						{
							validTagList.insert( newTag );
						}
					}

					if( !anyValidVariant )
					{
						newTag.tag = tagVariants[0];
						interimTagList.insert( newTag );
					}

				}
			}

		} else
		{
			if( gapInfoItr == gapMaps[peakChargeState-1].end() )
				return;

			gapVector_t& gapVector = gapInfoItr->second;

			if( gapVector.empty() )
				return;

			peakList.push_back( gapVector.front().fromPeakItr );

			size_t gapCount = gapVector.size();
			for( size_t i=0; i < gapCount; ++i )
			{
				if( tagIndex-1 == 0 )
					peakList.push_back( gapVector[i].peakItr );

				tag.push_back( gapVector[i].gapRes );
				peakErrors.push_back( gapVector[i].error );
				findTags_R( gapVector[i].nextGapInfo,
							tagIndex-1,
							tag,
							peakErrors,
							peakList,
							peakChargeState,
							numTagsGenerated,
							irBins );
				peakErrors.pop_back();
				tag.erase( tag.length()-1 );

				if( tagIndex-1 == 0 )
					peakList.pop_back();
			}

			peakList.pop_back();
		}
	}

	size_t Spectrum::findTags()
	{
		size_t numTagsGenerated = 0;
		gapMap_t::iterator gapInfoItr;
		string tag;
		vector< float > peakErrors;
		vector< PeakData::iterator > peakList;
		IRBins& irBins = SpectraList::intensityRanksumBinsTable[ g_rtConfig->TagLength ][ peakData.size() ];
//cout << peakData.size() << irBins << endl;
		peakErrors.reserve( g_rtConfig->tagPeakCount );
		peakList.reserve( g_rtConfig->tagPeakCount );

		tagCount = 0;

		for( int z=0; z < numFragmentChargeStates; ++z )
		{
			gapMap_t& gapMap = gapMaps[z];
			for( gapInfoItr = gapMap.begin(); gapInfoItr != gapMap.end(); ++gapInfoItr )
				findTags_R( gapInfoItr, g_rtConfig->TagLength, tag, peakErrors, peakList, z+1, numTagsGenerated, irBins );
		}

		return numTagsGenerated;
	}

	void SpectraList::CalculateIRBins_R( IRBins& theseRanksumBins, int tagLength, int numPeaks, int curRanksum, int curRank, int loopDepth )
	{
		if( loopDepth > tagLength )
			++ theseRanksumBins[ curRanksum ];
		else
			for( int rank = curRank + 1; rank <= numPeaks; ++rank )
				CalculateIRBins_R( theseRanksumBins, tagLength, numPeaks, curRanksum + rank, rank, loopDepth+1 );
	}

	void SpectraList::CalculateIRBins( int tagLength, int numPeaks )
	{
		if( intensityRanksumBinsTable.size() <= (size_t) tagLength )
			intensityRanksumBinsTable.resize( tagLength+1, vector< IRBins >() );
		if( intensityRanksumBinsTable[ tagLength ].size() <= (size_t) numPeaks )
			intensityRanksumBinsTable[ tagLength ].resize( numPeaks+1, IRBins() );
		IRBins& theseRanksumBins = intensityRanksumBinsTable[ tagLength ][ numPeaks ];
		theseRanksumBins.resize( (tagLength+1) * numPeaks, 0 );
		CalculateIRBins_R( theseRanksumBins, tagLength, numPeaks, 0, 0, 0 );

		float totalRanksum = 0;
		for( IRBins::iterator itr = theseRanksumBins.begin(); itr != theseRanksumBins.end(); ++itr )
			totalRanksum += *itr;

		float tmpRanksum = 0;
		for( IRBins::iterator itr = theseRanksumBins.begin(); itr != theseRanksumBins.end(); ++itr )
		{
			tmpRanksum += *itr;
			*itr = tmpRanksum / totalRanksum;
		}
	}

	void SpectraList::PrecacheIRBins( SpectraList& instance, string & cachename )
	{
		intensityRanksumBinsTable.clear();

		//cout << g_hostString << " is reading intensity ranksum bins cache file." << endl;
		ifstream cacheInputFile( cachename.c_str() );
	try{
		if( cacheInputFile.is_open() )
		{
			text_iarchive cacheInputArchive( cacheInputFile );
			cacheInputArchive & intensityRanksumBinsTable;
		}
	}
	catch(exception & e)
	{
		cout << e.what() << "\n";
		cout << "Error while reading "<< cachename << "please delete this file and start again";
		throw e;
	}
		cacheInputFile.close();

		//cout << g_hostString << " is calculating uncached ranksum bins (this could take a while)." << endl;
		for( iterator itr = instance.begin(); itr != instance.end(); ++itr )
		{
			if( intensityRanksumBinsTable.size() <= (size_t) g_rtConfig->TagLength ||
				intensityRanksumBinsTable[ g_rtConfig->TagLength ].size() <= (*itr)->peakData.size() ||
				intensityRanksumBinsTable[ g_rtConfig->TagLength ][ (*itr)->peakData.size() ].empty() )
			{
				//cout << g_rtConfig->TagLength << "," << (*itr)->peakData.size() << endl;
				CalculateIRBins( g_rtConfig->TagLength, (*itr)->peakData.size() );
			}
		}

		if( g_pid == 0 )
		{
		//	cout << g_hostString << " is writing intensity ranksum bins cache file." << endl;
			ofstream cacheOutputFile( cachename.c_str() );
			text_oarchive cacheOutputArchive( cacheOutputFile );
			cacheOutputArchive & intensityRanksumBinsTable;
			cacheOutputFile.close();
		}
	}

	void SpectraList::InitMzFEBins()
	{
		int numPeaks = g_rtConfig->tagPeakCount;
		vector< float > peakErrors( numPeaks );
		float peakErrorSum = 0.0f;
		for( int i=0; i < numPeaks; ++i )
		{
			peakErrors[i] = g_rtConfig->FragmentMzTolerance * i;
			peakErrorSum += peakErrors[i];
		}

		float peakErrorAvg = peakErrorSum / numPeaks;
		for( int i=0; i < numPeaks; ++i )
			peakErrors[i] -= peakErrorAvg;
		//cout << peakErrors << endl;

		float maxError = 0.0f;
		for( int i=0; i < numPeaks; ++i )
			maxError += pow( peakErrors[i], 2 );
		//cout << maxError << endl;

		mzFidelityErrorBins[ 0.0f ] = 0.0f;

		peakErrors.clear();
		peakErrors.resize( numPeaks, 0.0f );
		vector< float > sumErrors( numPeaks, 0.0f );
		vector< float > adjustedSumErrors( numPeaks, 0.0f );

		// Random sampling permits longer tag lengths
		boost::mt19937 rng(0);
		boost::uniform_real<float> MzErrorRange( -g_rtConfig->FragmentMzTolerance, g_rtConfig->FragmentMzTolerance );
		boost::variate_generator< boost::mt19937&, boost::uniform_real<float> > RandomMzError( rng, MzErrorRange );
		for( int i=0; i < g_rtConfig->MzFidelityErrorBinsSamples; ++i )
		{
			for( int p=1; p < numPeaks; ++p )
			{
				float e = RandomMzError();
				peakErrors[p] = e;
				sumErrors[p] = accumulate( peakErrors.begin(), peakErrors.begin()+p, e );
			}
			//cout << sumErrors << endl;
			//float sum = accumulate( peakErrors.begin(), peakErrors.end(), 0.0f );
			//float avg = sum / (int) peakErrors.size();
			float sum = accumulate( sumErrors.begin(), sumErrors.end(), 0.0f );
			float avg = sum / (int) sumErrors.size();

			sum = 0.0f;
			for( size_t i=0; i < sumErrors.size(); ++i )
			{
				adjustedSumErrors[i] = sumErrors[i] - avg;
				sum += pow( adjustedSumErrors[i], 2 );
			}

			mzFidelityErrorBins[ sum ] = 0;
		}

		float n = 0.0f;
		float totalSize = (float) mzFidelityErrorBins.size();
		for( MzFEBins::iterator itr = mzFidelityErrorBins.begin(); itr != mzFidelityErrorBins.end(); ++itr )
		{
			n += 1.0f;
			itr->second = n / totalSize;
		}
		//cout << mzFidelityErrorBins << endl;
	}

	void SpectraList::InitCEBins()
	{
		boost::mt19937 rng(0);
		boost::uniform_real<float> ComplementErrorRange( -g_rtConfig->ComplementMzTolerance, g_rtConfig->ComplementMzTolerance );
		boost::variate_generator< boost::mt19937&, boost::uniform_real<float> > RandomComplementError( rng, ComplementErrorRange );
		complementErrorBinsList.resize( g_rtConfig->tagPeakCount+1, CEBins() );
		for( int numComplements = 2; numComplements <= g_rtConfig->tagPeakCount; ++numComplements )
		{
			CEBins& errorBins = complementErrorBinsList[numComplements];
			errorBins[0.0f] = 0.0f;
			for( int i=0; i < g_rtConfig->MzFidelityErrorBinsSamples; ++i )
			{
				vector< float > errors;
				for( int j=0; j < numComplements; ++j )
					errors.push_back( RandomComplementError() );
				float mean = arithmetic_mean<float>(errors);
				for( int j=0; j < numComplements; ++j )
					errors[j] = pow( errors[j] - mean, 2.0f );
				float sse = accumulate( errors.begin(), errors.end(), 0.0f );
				errorBins[sse] = 0;
			}
			float count = 0;
			for( map< float, float >::iterator itr = errorBins.begin(); itr != errorBins.end(); ++itr, ++count )
				itr->second = count / (float) errorBins.size();
			//cout << errorBins << endl << endl;
		}
	}
}
}
