/*
    File: ScanFilter.cpp
    Description: parsing for Thermo/Xcalibur "filter line".
    Date: April 16, 2008

    Originally copyright (C) 2007 Joshua Tasman, ISB Seattle


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


#define RAWFILE_SOURCE

#include "ScanFilter.h"

#include <sstream>
#include <stack>
#include <iostream>
#include <cctype> // for toupper
#include <algorithm>
#include <vector>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

using namespace pwiz::raw;
using namespace std;
using boost::lexical_cast;

/*

FilterLine dictionary
--From Thermo


Analyzer:

ITMS		Ion Trap
TQMS		Triple Quad
SQMS		Single Quad
TOFMS		TOF
FTMS		ICR
Sector		Sector

Segment Scan Event   (Sectors only)

Polarity
-		Negative
+		Positive


Scan Data
c		centroid
p		profile


Ionization Mode
EI		Electron Impact
CI		Chemical Ionization
FAB		Fast Atom Bombardment
ESI		Electrospray
APCI		Atmospheric Pressure Chemical Ionization
NSI		Nanospray
TSP		Thermospray
FD		Field Desorption
MALDI	Matrix Assisted Laser Desorption Ionization
GD		Glow Discharge

Corona
corona			corona on
!corona		corona off

PhotoIoniziation
pi			photo ionization on
!pi			photo ionization off

Source CID
sid			source cid on
!sid			source cid off
sid=<x>     source cid on at <x> energy

Detector set
det			detector set
!det			detector not set

TurboScan
t			turbo scan on
!t			turob scan off

Enhanced			(Sectors only)
E			enhanced on
!E			enhanced off

Dependent Type
d			data dependent active
!d			data dependent not-active

Wideband
w			wideband activation on
!w			wideband activation off

Accurate Mass
!AM			accurate mass not active
AM			accurate mass active 
AMI			accurate mass with internal calibration
AME			accurate mass with external calibration

Ultra
u			ultra on
!u			ultra off

Scan Type:
full			full scan
SIM			single ion monitor
SRM			single reaction monitor
CRM
z			zoom scan
Q1MS			q1 mass spec scan
Q3MS			q3 mass spec scan 

Sector Scan			(Sectors only)
BSCAN		b scan
ESCAN		e scan


MSorder
MS2			MSn order
MS3
�
MS15

Activation Type
cid			collision induced dissociation
mpd
ecd			electron capture dissociation
pqd			pulsed q dissociation
etd			electron transfer dissociation
hcd			high energy collision dissociation
sa			supplemental cid
ptr			proton transfer reaction

Free Region			(Sectors only)
ffr1			field free region 1
ffr2			field free region 2

Mass range
[low mass � high mass]

*/

ScanFilterMassAnalyzerType
ScanFilter::parseMassAnalyzerType(const string& word)
{
	if (word == "ITMS")
		return ScanFilterMassAnalyzerType_ITMS;
	else if (word == "TQMS")
		return ScanFilterMassAnalyzerType_TQMS;
	else if (word == "SQMS")
		return ScanFilterMassAnalyzerType_SQMS;
	else if (word == "TOFMS")
		return ScanFilterMassAnalyzerType_TOFMS;
	else if (word == "FTMS")
		return ScanFilterMassAnalyzerType_FTMS;
	else if (word == "SECTOR")
		return ScanFilterMassAnalyzerType_Sector;
	else
		return ScanFilterMassAnalyzerType_Unknown;
}

PolarityType 
ScanFilter::parsePolarityType(const string& word)
{
	if (word == "+")
		return PolarityType_Positive;
	else if (word == "-")
		return PolarityType_Negative;
	else
		return PolarityType_Unknown;
}

DataPointType 
ScanFilter::parseDataPointType(const string& word)
{
	if (word == "C")
		return DataPointType_Centroid;
	else if (word == "P")
		return DataPointType_Profile;
	else
		return DataPointType_Unknown;
}

IonizationType 
ScanFilter::parseIonizationType(const string & word)
{
	if (word == "EI" )
		return IonizationType_EI;
	else if (word == "CI")
		return IonizationType_CI;
	else if (word == "FAB")
		return IonizationType_FAB;
	else if (word == "ESI")
		return IonizationType_ESI;
	else if (word == "APCI")
		return IonizationType_APCI;
	else if (word == "NSI")
		return IonizationType_NSI;
	else if (word == "TSP")
		return IonizationType_TSP;
	else if (word == "FD")
		return IonizationType_FD;
	else if (word == "MALDI")
		return IonizationType_MALDI;
	else if (word == "GD")
		return IonizationType_GD;
	else
		return IonizationType_Unknown;
}

AccurateMassType 
ScanFilter::parseAccurateMassType(const string& word)
{
	if (word == "!AM") {
		return AccurateMass_NotActive;
	}
	else if (word == "AM") {
		return AccurateMass_Active;
	}
	else if (word == "AMI") {
		return AccurateMass_ActiveWithInternalCalibration;
	}
	else if (word == "AME") {
		return AccurateMass_ActiveWithExternalCalibration;
	}
	else {
		return AccurateMass_Unknown;
	}
}

ScanType 
ScanFilter::parseScanType(const string& word)
{
	if (word == "FULL")
		return ScanType_Full;
	else if (word == "SIM")
		return ScanType_SIM;
	else if (word == "SRM")
		return ScanType_SRM;
	else if (word == "CRM")
		return ScanType_CRM;
	else if (word == "Z")
		return ScanType_Zoom;
	else if (word == "Q1MS")
		return ScanType_Q1MS;
	else if (word == "Q3MS")
		return ScanType_Q3MS;
	else
		return ScanType_Unknown;
}

ActivationType 
ScanFilter::parseActivationType(const string& word)
{
	if (word == "CID")
		return ActivationType_CID;
	else if (word == "MPD")
		return ActivationType_MPD;
	else if (word == "ECD")
		return ActivationType_ECD;
	else if (word == "PQD")
		return ActivationType_PQD;
	else if (word == "ETD")
		return ActivationType_ETD;
	else if (word == "HCD")
		return ActivationType_HCD;
	else if (word == "SA")
		return ActivationType_SA;
	else if (word == "PTR")
		return ActivationType_PTR;
	else
		return ActivationType_Unknown;
}


ScanFilter::ScanFilter()
{
    initialize();
};


ScanFilter::~ScanFilter()
{
}


void 
ScanFilter::print()
{
	if (massAnalyzerType_ > ScanFilterMassAnalyzerType_Unknown) {
        cout << "mass analyzer: " << massAnalyzerType_ << endl;
	}

	if (polarityType_ > PolarityType_Unknown) {
        cout << "polarity: " << polarityType_ << endl;
	}

	if (dataPointType_ > DataPointType_Unknown) {
        cout << "data point type: " << dataPointType_ << endl;
	}

	if (ionizationType_ > IonizationType_Unknown) {
        cout << "ionization type: " << ionizationType_ << endl;
	}

	if (coronaOn_ != TriBool_Unknown) {
		cout << "corona: " << coronaOn_ << endl;
	}

	if (photoIonizationOn_ != TriBool_Unknown) {
		cout << "photoionization: " << photoIonizationOn_ << endl;
	}

	if (sourceCIDOn_ != TriBool_Unknown) {
		cout << "source CID: " << sourceCIDOn_ << endl;
	}

	if (detectorSet_ != TriBool_Unknown) {
		cout << "detector set: " << detectorSet_ << endl;
	}

	if (turboScanOn_ != TriBool_Unknown) {
		cout << "turboscan: " << turboScanOn_ << endl;
	}

	if (dependentActive_ != TriBool_Unknown) {
		cout << "data dependent: " << dependentActive_ << endl;
	}

	if (widebandOn_ != TriBool_Unknown) {
		cout << "wideband: " << widebandOn_ << endl;
	}

	if (accurateMassType_ > AccurateMass_Unknown) {
		cout << "accurate mass: " << accurateMassType_ << endl;
	}

	if (scanType_ > ScanType_Unknown) {
		cout << "scan type: " << scanType_ << endl;
	}

	if (msLevel_ > 0 ) {
		cout << "MS level: " << msLevel_ << endl;
	}

	if (activationType_ > ActivationType_Unknown) {
		cout << "activation type: " << activationType_ << endl;
	}

	cout << endl << endl << endl;
}


void
ScanFilter::initialize()
{
    massAnalyzerType_ = ScanFilterMassAnalyzerType_Unknown;
    polarityType_ = PolarityType_Unknown;
    dataPointType_ = DataPointType_Unknown;
    ionizationType_ = IonizationType_Unknown;
    scanType_ = ScanType_Unknown;
    accurateMassType_ = AccurateMass_Unknown;
    activationType_ = ActivationType_Unknown;

    coronaOn_ = TriBool_Unknown;
    photoIonizationOn_ = TriBool_Unknown;
    sourceCIDOn_ = TriBool_Unknown;
    detectorSet_ = TriBool_Unknown;
    turboScanOn_ = TriBool_Unknown;
    dependentActive_ = TriBool_Unknown;
    widebandOn_ = TriBool_Unknown;

    msLevel_ = 0;
    cidParentMass_.clear();
	cidEnergy_.clear();
	scanRangeMin_.clear();
	scanRangeMax_.clear();
}


bool 
ScanFilter::parse(string filterLine)
{
    initialize();

	/**
	almost all of the fields are optional
	*/
	boost::to_upper(filterLine);
	stringstream s(filterLine);
	string w;

	if (s.eof()) {
		return 1; // ok, empty line
	}
	s >> w;

	massAnalyzerType_ = parseMassAnalyzerType(w);
	if (massAnalyzerType_ > ScanFilterMassAnalyzerType_Unknown) {
		// "analyzer" field was present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	polarityType_ = parsePolarityType(w);
	if (polarityType_ > PolarityType_Unknown) {
		// "polarity" field was present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	dataPointType_ = parseDataPointType(w);
	if (dataPointType_ > DataPointType_Unknown) {
		// "scan data type" field present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	ionizationType_ = parseIonizationType(w);
	if (ionizationType_ > IonizationType_Unknown) {
		// "ionization mode" field present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	bool advance = false;

	// corona
	if (w == "!CORONA") {
		coronaOn_ = TriBool_False;
		advance = true;
	}
	else if (w == "CORONA") {
		coronaOn_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}

	// photoIonization
	if (w == "!PI") {
		photoIonizationOn_ = TriBool_False;
		advance = true;
	}
	else if (w == "PI") {
		photoIonizationOn_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}

	// source CID
	if (w == "!SID") {
		sourceCIDOn_ = TriBool_False;
		advance = true;
	}
	else if (w.find("SID") == 0) { // handle cases where SID energy is explicit
		sourceCIDOn_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}


	// detector
	if (w == "!DET") {
		detectorSet_ = TriBool_False;
		advance = true;
	}
	else if (w == "DET") {
		detectorSet_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}


	// turboscan
	if (w == "!T") {
		turboScanOn_ = TriBool_False;
		advance = true;
	}
	else if (w == "T") {
		turboScanOn_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}


	// dependent type
	if (w == "!D") {
		dependentActive_ = TriBool_False;
		advance = true;
	}
	else if (w == "D") {
		dependentActive_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}


	// wideband
	if (w == "!W") {
		widebandOn_ = TriBool_False;
		advance = true;
	}
	else if (w == "W") {
		widebandOn_ = TriBool_True;
		advance = true;
	}
	if (advance) {
		if (s.eof()) {
			return 1;
		}
		s >> w;
		advance = false;
	}


	accurateMassType_ = parseAccurateMassType(w);
	if (accurateMassType_ > AccurateMass_Unknown) {
		// "accurate mass" field present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	scanType_ = parseScanType(w);
	if (scanType_ > ScanType_Unknown) {
		if (scanType_ == ScanType_Q1MS || scanType_ == ScanType_Q3MS)
			msLevel_ = 1;

		// "scan type" field present
		if (s.eof()) {
			return 1;
		}
		s >> w;
	}

	// MS order
	if ( (w.substr(0,2) == "MS") && (w.length() >= 2) ) {
		if (w.length() == 2) {
			msLevel_ = 1; // just "MS"
		} else {
			// MSn: extract int n
			//cout << "len: " << w.length() << endl;
			msLevel_ = lexical_cast<int>(w.substr(2)); // take number after "ms"
		}
		if (s.eof()) {
			return 1;
		}
		s >> w;      
	}


	// CID info
	// TODO: MSn for n > 2: comma-separated
	// if msLevel >=2 there should be mass@energy pairs for each level >= 2
	if (msLevel_ > 1) {
		int expectedPairs = msLevel_ - 1;
		for (int i=0; i< expectedPairs; ++i) {
			char c=w[0];
			size_t markerPos = w.find('@',0);
			// make sure this word starts with a numeric char, and the word contains "@"
			if( markerPos == string::npos || ! ( (c >= '0') && (c <= '9') ) )
				return false;
			size_t energyPos = markerPos+1;
			c = w[energyPos];
			if ( ! ( (c >= '0') && (c <= '9') ) ) {
				energyPos = w.find_first_of("1234567890-+", energyPos); // find first numeric character after the "@"
				if (energyPos != string::npos) {
					activationType_ = parseActivationType(w.substr(markerPos+1, energyPos-markerPos-1));
					if (activationType_ == ActivationType_Unknown)
						return false;
				} else
					return false;
			}

			string mass = w.substr(0, markerPos);
			string energy = w.substr(energyPos);
			// cout << "got mass " << mass << " at " << energy << " energy using activation " << (int) activationMethod_ << " (from " << w << ")" << endl;
			cidParentMass_.push_back(lexical_cast<double>(mass));
			cidEnergy_.push_back(lexical_cast<double>(energy));			

			// prematurely done?
			if (s.eof()) {
				return false;
			}
			s >> w;

		}
	}

	// try to get activation type if not already set
	if (activationType_ == ActivationType_Unknown) {
		activationType_ = parseActivationType(w);
		if (activationType_ > ActivationType_Unknown) {
			// "activation type" field present
			if (s.eof()) {
				return 1;
			}
			s >> w;
		}
	}


	// product masses or mass ranges
	// TODO: parse single values, for SIM, SRM, CRM
	// some test based on ms level?

	string w2;
	std::getline(s, w2); // get all tokens until closing bracket
	w.append(w2);
	boost::trim_if(w, boost::is_any_of("[ ]")); // trim flanking brackets and whitespace
	vector<string> massRangeStrs;
	boost::split(massRangeStrs, w, boost::is_any_of(","));
	for(size_t i=0; i < massRangeStrs.size(); ++i)
	{
		string& massRangeStr = massRangeStrs[i]; // "<rangeMin>-<rangeMax>"
		boost::trim(massRangeStr); // trim flanking whitespace
		vector<string> rangeMinMaxStrs;
		boost::split(rangeMinMaxStrs, massRangeStr, boost::is_any_of("-"));
		scanRangeMin_.push_back(lexical_cast<double>(rangeMinMaxStrs[0]));
		scanRangeMax_.push_back(lexical_cast<double>(rangeMinMaxStrs[1]));
	}

	if (s.eof()) {
		// cout << "done parsing" << endl;
		return true;
	}
	else {
		do {
			cout << "unparsed scan filter element: " << w << endl;
		} while (s >> w);
		return false;
	}
	//     while (!s.eof()) {
	//       string w;
	//       s >> w;
	//       cout << "word: " << w << endl;
	//     }
}
