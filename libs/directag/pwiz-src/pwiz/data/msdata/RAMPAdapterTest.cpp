//
// RAMPAdapterTest.cpp
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#include "RAMPAdapter.hpp"
#include "MSDataFile.hpp"
#include "examples.hpp"
#include "utility/misc/unit.hpp"
#include <boost/filesystem/operations.hpp>
#include <iostream>
#include <iterator>


using namespace std;
using namespace pwiz::msdata;
using namespace pwiz::util;


ostream* os_ = 0;

/*
RAMP will need to instantiate an msdata::RAMPAdapter object, and call it to fill
in the appropriate RAMP structures: 

    using namespace pwiz::msdata;

    RAMPAdapter adapter("something.mzML");

    unsigned int scanIndex = rampAdapter.index(19); // get index for scan 19

    ScanHeaderStruct temp;
    adapter.getScanHeader(scanIndex, temp);

    vector<double> buffer; 
    adapter.getScanPeaks(scanIndex, buffer);

    RunHeaderStruct temp2;
    adapter.getRunHeader(temp2);
    
    InstrumentHeaderStruct temp3;
    adapter.getInstrumentHeader(temp3);

Note that the MSData library can throw exceptions, so the RAMP code will have to be able
to handle these gracefully.  Another option is to have these RAMPAdapter functions catch
any exceptions and return an error code if something goes wrong.
*/


string writeTempFile()
{
    const string& filename = "temp.RAMPAdapterTest.tiny.mzML";
    MSData tiny; 
    examples::initializeTiny(tiny);
    MSDataFile::write(tiny, filename);
    return filename;
}


ostream& operator<<(ostream& os, const ScanHeaderStruct& header)
{
   os << "seqNum: " << header.seqNum << endl;
   os << "acquisitionNum: " << header.acquisitionNum << endl;
   os << "msLevel: " << header.msLevel << endl;
   os << "peaksCount: " << header.peaksCount << endl;
   os << "totIonCurrent: " << header.totIonCurrent << endl;
   os << "retentionTime: " << header.retentionTime << endl;
   os << "basePeakMZ: " << header.basePeakMZ << endl;
   os << "basePeakIntensity: " << header.basePeakIntensity << endl;
   os << "collisionEnergy: " << header.collisionEnergy << endl;
   os << "ionisationEnergy: " << header.ionisationEnergy << endl;
   os << "lowMZ: " << header.lowMZ << endl;
   os << "highMZ: " << header.highMZ << endl;
   os << "precursorScanNum: " << header.precursorScanNum << endl;
   os << "precursorMZ: " << header.precursorMZ << endl;
   os << "precursorCharge: " << header.precursorCharge << endl;
   os << "precursorIntensity: " << header.precursorIntensity << endl;
   os << "scanType: " << header.scanType << endl;
   os << "mergedScan: " << header.mergedScan << endl;
   os << "mergedResultScanNum: " << header.mergedResultScanNum << endl;
   os << "mergedResultStartScanNum: " << header.mergedResultStartScanNum << endl;
   os << "mergedResultEndScanNum: " << header.mergedResultEndScanNum << endl;
   os << "filePosition: " << header.filePosition << endl;
   return os;
}


ostream& operator<<(ostream& os, const RunHeaderStruct& header)
{
    os << "scanCount: " << header.scanCount << endl;
    os << "lowMZ: " << header.lowMZ << endl;
    os << "highMZ: " << header.highMZ << endl;
    os << "startMZ: " << header.startMZ << endl;
    os << "endMZ: " << header.endMZ << endl;
    os << "dStartTime: " << header.dStartTime << endl;
    os << "dEndTime: " << header.dEndTime << endl;
    return os;
}


ostream& operator<<(ostream& os, const InstrumentStruct& instrument)
{
    os << "manufacturer: " << instrument.manufacturer << endl;
    os << "model: " << instrument.model << endl;
    os << "ionisation: " << instrument.ionisation << endl;
    os << "analyzer: " << instrument.analyzer << endl;
    os << "detector: " << instrument.detector << endl;
    return os;
}


void test(const string& filename)
{
    RAMPAdapter adapter(filename);

    size_t scanCount = adapter.scanCount();
    if (os_) *os_ << "scanCount: " << scanCount << "\n\n";
    unit_assert(scanCount == 4);

    unit_assert(adapter.index(19) == 0);
    unit_assert(adapter.index(20) == 1);
    unit_assert(adapter.index(21) == 2);
    unit_assert(adapter.index(22) == 3);

    // first scan (scan number == 19)

    ScanHeaderStruct header1;
    adapter.getScanHeader(0, header1);
    if (os_) *os_ << header1;
    unit_assert(header1.seqNum == 1);
    unit_assert(header1.acquisitionNum == 19);
    unit_assert(header1.msLevel == 1);
    unit_assert(header1.peaksCount == 15);
    const double epsilon = 1e-8;
    unit_assert_equal(header1.totIonCurrent, 1.66755e7, epsilon);
    unit_assert_equal(header1.retentionTime, 353.43, epsilon);
    unit_assert_equal(header1.basePeakMZ, 445.347, epsilon);
    unit_assert_equal(header1.basePeakIntensity, 120053, epsilon);
    unit_assert_equal(header1.collisionEnergy, 0., epsilon);
    unit_assert_equal(header1.lowMZ, 400.39, epsilon);
    unit_assert_equal(header1.highMZ, 1795.56, epsilon);
    unit_assert(header1.precursorScanNum == 0);
    unit_assert(header1.scanType == string("Full"));

    vector<double> peaks;
    adapter.getScanPeaks(0, peaks);
    unit_assert(peaks.size() == 30);
    if (os_)
    {
        const MZIntensityPair* begin = reinterpret_cast<const MZIntensityPair*>(&peaks[0]);
        copy(begin, begin+15, ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }

    // second scan (scan number == 20)

    ScanHeaderStruct header2;
    adapter.getScanHeader(1, header2);
    if (os_) *os_ << header2;
    unit_assert(header2.seqNum == 2);
    unit_assert(header2.acquisitionNum == 20);
    unit_assert(header2.msLevel == 2);
    unit_assert(header2.peaksCount == 10);
    unit_assert_equal(header2.totIonCurrent, 1.66755e7, epsilon);
    unit_assert_equal(header2.retentionTime, 359.43, epsilon);
    unit_assert_equal(header2.basePeakMZ, 456.347, epsilon);
    unit_assert_equal(header2.basePeakIntensity, 23433, epsilon);
    unit_assert_equal(header2.collisionEnergy, 35, epsilon);
    unit_assert_equal(header2.lowMZ, 320.39, epsilon);
    unit_assert_equal(header2.highMZ, 1003.56, epsilon);
    unit_assert(header2.precursorScanNum == 19);
    unit_assert_equal(header2.precursorMZ, 445.34, epsilon);
    unit_assert(header2.precursorCharge == 2);
    unit_assert_equal(header2.precursorIntensity, 120053, epsilon);
    unit_assert(header2.scanType == string("Full"));

    adapter.getScanPeaks(1, peaks);
    unit_assert(peaks.size() == 20);
    if (os_)
    {
        const MZIntensityPair* begin = reinterpret_cast<const MZIntensityPair*>(&peaks[0]);
        copy(begin, begin+10, ostream_iterator<MZIntensityPair>(*os_, "\n"));
        *os_ << endl;
    }

    // RunHeader

    RunHeaderStruct runHeader;
    adapter.getRunHeader(runHeader);
    unit_assert(runHeader.scanCount == 4);
    unit_assert(runHeader.lowMZ == 0);
    unit_assert(runHeader.highMZ == 0);
    unit_assert(runHeader.startMZ == 0);
    unit_assert(runHeader.endMZ == 0);
    unit_assert_equal(runHeader.dStartTime, header1.retentionTime, 1e-6);
    unit_assert_equal(runHeader.dEndTime, header2.retentionTime, 1e-6);

    if (os_)
        *os_ << "RunHeader:\n" << runHeader << endl;

    // Instrument
    InstrumentStruct instrument;
    adapter.getInstrument(instrument);
    if (os_)
        *os_ << "Instrument:\n" << instrument << endl;

    unit_assert(!strcmp(instrument.manufacturer, "Thermo Finnigan"));
    unit_assert(!strcmp(instrument.model, "LCQ Deca"));
    unit_assert(!strcmp(instrument.ionisation, "nanoelectrospray"));
    unit_assert(!strcmp(instrument.analyzer, "quadrupole ion trap"));
    unit_assert(!strcmp(instrument.detector, "electron multiplier"));
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        string filename = writeTempFile();
        test(filename);
       boost::filesystem::remove(filename);
        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
        return 1;
    }
}


