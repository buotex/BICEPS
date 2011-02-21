/*
    File: ScanFilter.h
    Description: parsing for Thermo/Xcalibur "filter line".
    Date: July 25, 2007

    Copyright (C) 2007 Joshua Tasman, ISB Seattle

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/


#ifndef _SCANFILTER_H_
#define _SCANFILTER_H_

#ifdef RAWFILE_DYN_LINK
#ifdef RAWFILE_SOURCE
#define RAWFILE_API __declspec(dllexport)
#else
#define RAWFILE_API __declspec(dllimport)
#endif  // RAWFILE_SOURCE
#endif  // RAWFILE_DYN_LINK

// if RAWFILE_API isn't defined yet define it now:
#ifndef RAWFILE_API
#define RAWFILE_API
#endif

// disable warning "class needs to have dll-interface to be used by clients of class"
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif

#include "RawFileTypes.h"
#include <string>
#include <vector>


namespace pwiz {
namespace raw {


class RAWFILE_API ScanFilter
{
    public:

	ScanFilterMassAnalyzerType parseMassAnalyzerType(const std::string& word);
	PolarityType parsePolarityType(const std::string& word);
	DataPointType parseDataPointType(const std::string& word);
	IonizationType parseIonizationType(const std::string & word);
	ScanType parseScanType(const std::string& word);
	ActivationType parseActivationType(const std::string& word);
	AccurateMassType parseAccurateMassType(const std::string& word);

	ScanFilterMassAnalyzerType massAnalyzerType_;
	PolarityType polarityType_;
	DataPointType dataPointType_;
	IonizationType ionizationType_;
	TriBool coronaOn_;
	TriBool photoIonizationOn_;
	TriBool sourceCIDOn_;
	TriBool detectorSet_;
	TriBool turboScanOn_;
	TriBool dependentActive_; // t: data-dependent active; f: non active
	TriBool widebandOn_; // wideband activation
	AccurateMassType accurateMassType_;
	ScanType scanType_;
	int msLevel_; // n, in MSn: >0
	ActivationType activationType_;

	std::vector<double> cidParentMass_; // one entry per ms level for level >= 2
	std::vector<double> cidEnergy_; // relative units; one entry per ms level for level >= 2

	std::vector<double> scanRangeMin_;
	std::vector<double> scanRangeMax_;


	ScanFilter();
	~ScanFilter();

	void print();

    void initialize();
	bool parse(std::string filterLine);

};

} // namespace raw
} // namespace pwiz

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // _SCANFILTER_H_
