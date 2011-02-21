#include "stdafx.h"
#include "BaseSpectrum.h"

using namespace freicore;

namespace freicore
{
	BaseSpectrum::BaseSpectrum()
		:	peakPreCount(0), peakCount(0), numSequenceComparisons(0), numFragmentChargeStates(0),
			mzOfPrecursor(0), mOfPrecursor(0), mOfUnadjustedPrecursor(0), retentionTime(0),
			processingTime(0), mzUpperBound(0), mzLowerBound(0), totalIonCurrent(0), totalPeakSpace(0)
	{}

	BaseSpectrum::BaseSpectrum( const BaseSpectrum& old )
	{
		id					    = old.id;
        stringID                = old.stringID;
        nativeID                = old.nativeID;
        fileName				= old.fileName;

		peakPreCount			= old.peakPreCount;
		peakCount				= old.peakCount;

		mzOfPrecursor			= old.mzOfPrecursor;
		mOfPrecursor			= old.mOfPrecursor;
		mOfUnadjustedPrecursor	= old.mOfUnadjustedPrecursor;
		retentionTime			= old.retentionTime;
		processingTime			= old.processingTime;
		numSequenceComparisons	= old.numSequenceComparisons;
		mzUpperBound			= old.mzUpperBound;
		mzLowerBound			= old.mzLowerBound;
		totalIonCurrent			= old.totalIonCurrent;
		totalPeakSpace			= old.totalPeakSpace;
	}

	BaseSpectrum::~BaseSpectrum()
	{}

	float BaseSpectrum::CalculateComplementMz( float mz, int z )
	{
		float chargedPrecursorMass = mOfPrecursor + ( id.charge * PROTON );
		float chargedFragmentMass = mz * z;
		float chargedComplementMass = chargedPrecursorMass - chargedFragmentMass;
		int complementCharge = id.charge - z;
		return chargedComplementMass / complementCharge;
	}

	size_t BaseSpectrum::size()
	{
		size_t mySize = sizeof( *this );
		mySize += fileName.length();
		return mySize;
	}
}
