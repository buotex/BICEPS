//
// Reader_Thermo_Test.cpp
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


#include "Reader_Thermo.hpp"
#include "data/msdata/TextWriter.hpp"
#include "utility/misc/unit.hpp"
#include <iostream>
#include <fstream>


using namespace std;
using namespace pwiz::util;
using namespace pwiz::msdata;


ostream* os_ = 0;


void testAccept(const string& filename)
{
    if (os_) *os_ << "testAccept(): " << filename << endl;

    ifstream is(filename.c_str());
    if (!is)
        throw runtime_error(("Unable to open file: " + filename).c_str());

    string head(512, '\0');
    is.read(&head[0], (std::streamsize)head.size());
    is.close();

    Reader_Thermo reader;
    bool accepted = reader.accept(filename, head);
    if (os_) *os_ << "accepted: " << boolalpha << accepted << endl;

    #ifdef _MSC_VER
    unit_assert(accepted);
    #else // _MSC_VER
    unit_assert(!accepted);
    #endif // _MSC_VER
}


void testRead(const string& filename)
{
    if (os_) *os_ << "testRead(): " << filename << endl;

    // read RAW file into MSData object

    Reader_Thermo reader;
    MSData msd;
    reader.read(filename, "dummy", msd);

    // make assertions about msd

    //if (os_) TextWriter(*os_,0)(msd); 

    unit_assert(msd.run.spectrumListPtr.get());
    SpectrumList& sl = *msd.run.spectrumListPtr;
    if (os_) *os_ << "spectrum list size: " << sl.size() << endl;
    unit_assert(sl.size() == 48);

    SpectrumPtr spectrum = sl.spectrum(0, true);
    if (os_) TextWriter(*os_,0)(*spectrum); 
    unit_assert(spectrum->index == 0);
    unit_assert(spectrum->id == "S1"); // derived from scan number
    unit_assert(spectrum->nativeID == "1"); // scan number
    unit_assert(sl.spectrumIdentity(0).index == 0);
    unit_assert(sl.spectrumIdentity(0).id == "S1");
    unit_assert(sl.spectrumIdentity(0).nativeID == "1");
    int scanEvent = spectrum->spectrumDescription.scan.cvParamChild(MS_preset_scan_configuration).valueAs<int>();
    unit_assert(scanEvent == 1);
    unit_assert(spectrum->cvParam(MS_ms_level).valueAs<int>() == 1); 
    //CVParam massAnalyzer = spectrum->spectrumDescription.scan.cvParamChild(MS_mass_analyzer_type);
    //unit_assert(massAnalyzer.cvid == MS_FT_ICR);
    vector<MZIntensityPair> data;
    spectrum->getMZIntensityPairs(data);
    unit_assert(data.size() == 19914);

    spectrum = sl.spectrum(1, true);
    if (os_) TextWriter(*os_,0)(*spectrum); 
    unit_assert(spectrum->index == 1);
    unit_assert(spectrum->id == "S2"); // derived from scan number
    unit_assert(spectrum->nativeID == "2"); // scan number
    unit_assert(sl.spectrumIdentity(1).index == 1);
    unit_assert(sl.spectrumIdentity(1).id == "S2");
    unit_assert(sl.spectrumIdentity(1).nativeID == "2");
    scanEvent = spectrum->spectrumDescription.scan.cvParamChild(MS_preset_scan_configuration).valueAs<int>();
    unit_assert(scanEvent == 2);
    unit_assert(spectrum->cvParam(MS_ms_level).valueAs<int>() == 1); 
    //massAnalyzer = spectrum->spectrumDescription.scan.cvParamChild(MS_mass_analyzer_type);
    //unit_assert(massAnalyzer.cvid == MS_ion_trap);
    spectrum->getMZIntensityPairs(data);
    unit_assert(data.size() == 19800);

    spectrum = sl.spectrum(2, true);
    if (os_) TextWriter(*os_,0)(*spectrum); 
    unit_assert(spectrum->index == 2);
    unit_assert(spectrum->id == "S3"); // scan number
    unit_assert(sl.spectrumIdentity(2).index == 2);
    unit_assert(sl.spectrumIdentity(2).id == "S3");
    unit_assert(sl.spectrumIdentity(2).nativeID == "3");
    scanEvent = spectrum->spectrumDescription.scan.cvParamChild(MS_preset_scan_configuration).valueAs<int>();
    unit_assert(scanEvent == 3);
    unit_assert(spectrum->cvParam(MS_ms_level).valueAs<int>() == 2); 
    //massAnalyzer = spectrum->spectrumDescription.scan.cvParamChild(MS_mass_analyzer_type);
    //unit_assert(massAnalyzer.cvid == MS_ion_trap);
    spectrum->getMZIntensityPairs(data);
    unit_assert(data.size() == 485);
    unit_assert(spectrum->spectrumDescription.precursors.size() == 1);
    const Precursor& precursor = spectrum->spectrumDescription.precursors[0];
    const SelectedIon& selectedIon = precursor.selectedIons[0];
    unit_assert(precursor.spectrumID == "S2"); // previous ms1 scan
    unit_assert_equal(selectedIon.cvParam(MS_m_z).valueAs<double>(), 810.79, 1e-15);
    unit_assert_equal(selectedIon.cvParam(MS_intensity).valueAs<double>(), 0, 1e-15);

    spectrum = sl.spectrum(5, true);
    unit_assert(sl.spectrumIdentity(5).index == 5);
    unit_assert(sl.spectrumIdentity(5).id == "S6");
    unit_assert(sl.spectrumIdentity(5).nativeID == "6");
    unit_assert(spectrum->spectrumDescription.precursors.size() == 1);
    const Precursor& precursor2 = spectrum->spectrumDescription.precursors[0];
    unit_assert(precursor2.spectrumID == "S2"); // previous ms1 scan

    // test chromatogram list
    unit_assert(msd.run.chromatogramListPtr.get());
    ChromatogramList& cl = *msd.run.chromatogramListPtr;
    if (os_) *os_ << "chromatogram list size: " << cl.size() << endl;
    unit_assert(cl.size() == 1);
    ChromatogramPtr chromatogram = cl.chromatogram(0, true);
    unit_assert(chromatogram->id == "TIC");

    // test file-level metadata 
    unit_assert(msd.fileDescription.fileContent.hasCVParam(MS_MSn_spectrum));
}


void test(const string& filename)
{
    testAccept(filename);
    
    #ifdef _MSC_VER
    testRead(filename);
    #else
    if (os_) *os_ << "Not MSVC -- nothing to do.\n";
    #endif // _MSC_VER
}


int main(int argc, char* argv[])
{
    try
    {
        vector<string> filenames;

        for (int i=1; i<argc; i++)
        {
            if (!strcmp(argv[i],"-v")) os_ = &cout;
            else filenames.push_back(argv[i]);
        }

        if (filenames.size()!=1)
            throw runtime_error("Usage: Reader_Thermo_Test [-v] filename"); 
            
        test(filenames[0]);
        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
    }
    catch (...)
    {
        cerr << "Caught unknown exception.\n";
    }
    
    return 1;
}

