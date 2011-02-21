//
// write_examples.cpp 
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


#include "MSDataFile.hpp"
#include "examples.hpp"
#include <iostream>
#include <fstream>


using namespace std;
using namespace pwiz;
using namespace pwiz::msdata;
using boost::shared_ptr;



void writeTiny()
{
    // create the MSData object in memory
    MSData msd;
    examples::initializeTiny(msd); 

    // write out mzML 
    string filename = "tiny.pwiz.mzML";
    cout << "Writing file " << filename << endl;
    MSDataFile::write(msd, filename);

    // write out mzXML
    filename = "tiny.pwiz.mzXML";
    cout << "Writing file " << filename << endl;
    MSDataFile::write(msd, filename, MSDataFile::Format_mzXML);
}


void writeMIAPE()
{
    const string& inputFile = "small.pwiz.mzML";
    const string& outputFile = "small_miape.pwiz.mzML";

    try
    {
        MSDataFile msd(inputFile);
        examples::addMIAPEExampleMetadata(msd);
        cout << "Writing file " << outputFile << endl;
        msd.write(outputFile);
    }
    catch (exception& e)
    {
        cerr << "Error opening file " << inputFile << endl;
    }
}


int main()
{
    try
    {
        writeTiny();
        writeMIAPE();

        cout << "\nhttp://proteowizard.sourceforge.net\n"
             << "support@proteowizard.org\n";

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

