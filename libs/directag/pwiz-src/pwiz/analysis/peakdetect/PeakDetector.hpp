//
// PeakDetector.hpp
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


#ifndef _PEAKDETECTOR_HPP_
#define _PEAKDETECTOR_HPP_


#include "utility/misc/Export.hpp"
#include "data/msdata/MSData.hpp"
#include "data/misc/PeakData.hpp"


namespace pwiz {
namespace analysis {


///
/// interface for peak detection
/// 
class PWIZ_API_DECL PeakDetector
{
    public:

    typedef pwiz::msdata::MZIntensityPair MZIntensityPair;
    typedef pwiz::data::peakdata::Peak Peak;
    
    /// find peaks in a specified array of MZIntensityPair 
    virtual void detect(const MZIntensityPair* begin,
                        const MZIntensityPair* end,
                        std::vector<Peak>& result) = 0;

    /// convenience function -- equivalent to:
    ///   detect(&data[0], &data[0]+data.size(), result) 
    virtual void detect(const std::vector<MZIntensityPair>& data,
                        std::vector<Peak>& result);

    virtual ~PeakDetector() {} 
};


} // namespace analysis 
} // namespace pwiz


#endif // _PEAKDETECTOR_HPP_

