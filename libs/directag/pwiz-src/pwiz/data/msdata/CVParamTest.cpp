//
// CVParamTest.cpp
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2007 Spielberg Family Center for Applied Proteomics
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


#include "CVParam.hpp"
#include "utility/misc/unit.hpp"
#include <iostream>
#include <iterator>
#include <algorithm>


using namespace std;
using namespace pwiz::util;
using namespace pwiz::msdata;


ostream* os_ = 0;


class WriteCVParam
{
    public:

    WriteCVParam(ostream& os) : os_(os) {}

    void operator()(const CVParam& param)
    {
        os_ << "<cvParam " 
            << "cvLabel=\"" << cvinfo(param.cvid).id.substr(0,2) << "\" "
            << "accession=\"" << cvinfo(param.cvid).id << "\" "
            << "name=\"" << cvinfo(param.cvid).name << "\" "
            << "value=\"" << param.value << "\"";

        if (param.units != CVID_Unknown)
        {
            os_ << " unitAccession=\"" << cvinfo(param.units).id << "\" "
                << "unitName=\"" << cvinfo(param.units).name << "\""; 
        }

        os_ << "/>\n";
    }

    private:
    ostream& os_;
};


const char* mzmlScanTime = 
    "<cvParam cvLabel=\"MS\" accession=\"MS:1000016\" name=\"scan time\" value=\"5.890500\" "
    "unitAccession=\"MS:1000038\" unitName=\"minute\"/>\n";

const char* mzmlCollisionEnergy = 
    "<cvParam cvLabel=\"MS\" accession=\"MS:1000045\" name=\"collision energy\" value=\"35.00\" "
    "unitAccession=\"MS:1000137\" unitName=\"electron volt\"/>\n";


void test()
{
    vector<CVParam> params;

    params.push_back(CVParam(MS_lowest_m_z_value, 420));
    params.push_back(CVParam(MS_highest_m_z_value, 2000.012345));
    params.push_back(CVParam(MS_m_z, "goober"));
    params.push_back(CVParam(MS_scan_time, 5.890500, MS_minute)); 
    params.push_back(CVParam(MS_collision_energy, 35.00, MS_electron_volt)); 
    params.push_back(CVParam(MS_deisotoping, true)); 
    params.push_back(CVParam(MS_peak_picking, false)); 

    if (os_)
    {
        *os_ << "params:\n";
        copy(params.begin(), params.end(), ostream_iterator<CVParam>(*os_, "\n")); 
        *os_ << endl;
    
        *os_ << "as mzML <cvParam> elements:\n";
        for_each(params.begin(), params.end(), WriteCVParam(*os_));
        *os_ << endl;

        *os_ << "value casting:\n";
        int temp = params[0].valueAs<int>();
        *os_ << temp << endl;
        float temp2 = params[1].valueAs<float>();
        *os_ << temp2 << endl;
        string temp3 = params[2].valueAs<string>();
        *os_ << temp3 << "\n\n";
    }

    // verify simple things
    unit_assert(420 == params[0].valueAs<int>());
    unit_assert(2000.012345 == params[1].valueAs<double>());
    unit_assert("goober" == params[2].value);
    unit_assert(5.890500 == params[3].valueAs<double>());
    unit_assert(35.00 == params[4].valueAs<double>());
    unit_assert(params[0] == CVParam(MS_lowest_m_z_value, 420));
    unit_assert(params[1] != CVParam(MS_lowest_m_z_value, 420));
    unit_assert(CVParam(MS_m_z) == MS_m_z);
    unit_assert(params[5].valueAs<bool>() == true);
    unit_assert(params[6].valueAs<bool>() == false);

    // verify manual mzml writing -- this is to verify that we have enough
    // info to write <cvParam> elements as required by mzML

    ostringstream ossScanTime;
    CVParam scanTime(MS_scan_time, "5.890500", MS_minute); 
    (WriteCVParam(ossScanTime))(scanTime);
    unit_assert(ossScanTime.str() == mzmlScanTime);
    if (os_) *os_ << "scan time in seconds: " << scanTime.timeInSeconds() << endl;
    unit_assert_equal(scanTime.timeInSeconds(), 5.8905 * 60, 1e-10);

    ostringstream ossCollisionEnergy;
    (WriteCVParam(ossCollisionEnergy))(CVParam(MS_collision_energy, "35.00", MS_electron_volt));
    unit_assert(ossCollisionEnergy.str() == mzmlCollisionEnergy);
}


void testIs()
{
    vector<CVParam> params;
    params.push_back(CVParam(MS_plasma_desorption));
    params.push_back(CVParam(MS_lowest_m_z_value, 420));
    params.push_back(CVParam(MS_collision_induced_dissociation));

    vector<CVParam>::const_iterator it = 
        find_if(params.begin(), params.end(), CVParamIs(MS_lowest_m_z_value));

    unit_assert(it->value == "420");
}


void testIsChildOf()
{
    // example of how to search through a collection of CVParams
    // to find the first one whose cvid IsA specified CVID

    vector<CVParam> params;
    params.push_back(CVParam(MS_lowest_m_z_value, 420));
    params.push_back(CVParam(MS_plasma_desorption));
    params.push_back(CVParam(MS_collision_induced_dissociation));
    params.push_back(CVParam(MS_electron_volt));
    params.push_back(CVParam(MS_highest_m_z_value, 2400.0));

    vector<CVParam>::const_iterator itDiss = 
        find_if(params.begin(), params.end(), CVParamIsChildOf(MS_dissociation_method));

    vector<CVParam>::const_iterator itUnit = 
        find_if(params.begin(), params.end(), CVParamIsChildOf(MS_unit));

    if (os_)
    {
        *os_ << "find dissociation method: " 
             << (itDiss!=params.end() ? cvinfo(itDiss->cvid).name : "not found")
             << endl;

        *os_ << "find unit: " 
             << (itUnit!=params.end() ? cvinfo(itUnit->cvid).name : "not found")
             << endl;

    }

    unit_assert(itDiss!=params.end() && itDiss->cvid==MS_plasma_desorption);
    unit_assert(itUnit!=params.end() && itUnit->cvid==MS_electron_volt);
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        test();
        testIs();
        testIsChildOf();
        return 0;
    }
    catch (exception& e)
    {
        cerr << e.what() << endl;
    }
    
    return 1;
}

