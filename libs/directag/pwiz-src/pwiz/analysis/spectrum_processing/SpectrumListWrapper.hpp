//
// SpectrumListWrapper.hpp
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


#ifndef _SPECTRUMLISTWRAPPER_HPP_ 
#define _SPECTRUMLISTWRAPPER_HPP_ 


#include "data/msdata/MSData.hpp"
#include <stdexcept>


namespace pwiz {
namespace analysis {


/// Inheritable pass-through implementation for wrapping a SpectrumList 
class PWIZ_API_DECL SpectrumListWrapper : public msdata::SpectrumList
{
    public:

    SpectrumListWrapper(const msdata::SpectrumListPtr& inner)
    :   inner_(inner)
    {
        if (!inner.get()) throw std::runtime_error("[SpectrumListWrapper] Null SpectrumListPtr.");
    }

    static bool accept(const msdata::SpectrumListPtr& inner) {return true;}

    virtual size_t size() const {return inner_->size();}
    virtual bool empty() const {return inner_->empty();}
    virtual const msdata::SpectrumIdentity& spectrumIdentity(size_t index) const {return inner_->spectrumIdentity(index);} 
    virtual size_t find(const std::string& id) const {return inner_->find(id);}
    virtual size_t findNative(const std::string& nativeID) const {return inner_->findNative(nativeID);}
    virtual msdata::SpectrumPtr spectrum(size_t index, bool getBinaryData = false) const {return inner_->spectrum(index, getBinaryData);}

    protected:

    msdata::SpectrumListPtr inner_;
};


} // namespace analysis 
} // namespace pwiz


#endif // _SPECTRUMLISTWRAPPER_HPP_ 

