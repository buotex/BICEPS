//
// ChemistryData.hpp 
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2006 Louis Warschaw Prostate Cancer Center
//   Cedars Sinai Medical Center, Los Angeles, California  90048
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


#ifndef _CHEMISTRYDATA_HPP_
#define _CHEMISTRYDATA_HPP_


#include "utility/misc/Export.hpp"
#include "Chemistry.hpp"


namespace pwiz {
namespace proteome {
namespace ChemistryData {


typedef pwiz::proteome::Chemistry::Element::Type Type;


struct PWIZ_API_DECL Isotope 
{
    double mass; 
    double abundance;
};


struct PWIZ_API_DECL Element 
{
    Type type;
    const char* symbol;
    int atomicNumber;
    double atomicWeight;
    Isotope* isotopes;
    int isotopesSize;
};


PWIZ_API_DECL Element* elements();
PWIZ_API_DECL int elementsSize();


} // namespace ChemistryData
} // namespace proteome
} // namespace pwiz


#endif // _CHEMISTRYDATA_HPP_
