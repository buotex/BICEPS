//
// Serializer_mzML_Test.cpp
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


#include "Serializer_mzML.hpp"
#include "Diff.hpp"
#include "examples.hpp"
#include "utility/misc/unit.hpp"
#include "boost/iostreams/positioning.hpp"
#include <iostream>
#include <fstream>


using namespace std;
using namespace pwiz::util;
using namespace pwiz::msdata;
using boost::shared_ptr;


ostream* os_ = 0;


void testWriteRead(const MSData& msd, const Serializer_mzML::Config& config)
{
    if (os_) *os_ << "testWriteRead() " << config << endl;

    Serializer_mzML mzmlSerializer(config);

    ostringstream oss;
    mzmlSerializer.write(oss, msd);

    if (os_) *os_ << "oss:\n" << oss.str() << endl; 

    boost::shared_ptr<istringstream> iss(new istringstream(oss.str()));
    MSData msd2;
    mzmlSerializer.read(iss, msd2);

    Diff<MSData> diff(msd, msd2);
    if (os_ && diff) *os_ << diff << endl; 
    unit_assert(!diff);
}


void testWriteRead()
{
    MSData msd;
    examples::initializeTiny(msd);

    Serializer_mzML::Config config;
    unit_assert(config.binaryDataEncoderConfig.precision == BinaryDataEncoder::Precision_64);
    testWriteRead(msd, config);

    config.binaryDataEncoderConfig.precision = BinaryDataEncoder::Precision_32;
    testWriteRead(msd, config);

    config.indexed = false;
    testWriteRead(msd, config);

    // TODO: test with compression 
}


int main(int argc, char* argv[])
{
    try
    {
        if (argc>1 && !strcmp(argv[1],"-v")) os_ = &cout;
        testWriteRead();
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

