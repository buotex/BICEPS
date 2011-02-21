//
// MSDataCache.cpp
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


#define PWIZ_SOURCE

#include "MSDataCache.hpp"
#include "data/msdata/MSDataFile.hpp"
#include <stdexcept>
#include <list>
#include <iostream>


namespace pwiz {
namespace analysis {


using boost::shared_ptr;
using boost::lexical_cast;
using namespace std;


//
// MSDataCache::Impl
//


struct MSDataCache::Impl
{
    Impl(const Config& _config) : config(_config) {}

    MSDataCache::Config config;

    typedef list<SpectrumInfo*> MRU; // most recently used
    MRU mru;

    SpectrumListPtr spectrumListPtr;
};


//
// MSDataCache
//


PWIZ_API_DECL MSDataCache::MSDataCache(const MSDataCache::Config& config)
:   impl_(new Impl(config)) 
{}


PWIZ_API_DECL void MSDataCache::open(const DataInfo& dataInfo)
{
    clear();

    impl_->mru.clear();

    if (dataInfo.msd.run.spectrumListPtr.get())
    {
        resize(dataInfo.msd.run.spectrumListPtr->size());
        impl_->spectrumListPtr = dataInfo.msd.run.spectrumListPtr;
    }
}


PWIZ_API_DECL
void MSDataCache::update(const DataInfo& dataInfo,
                         const Spectrum& spectrum)
{
    if (!dataInfo.msd.run.spectrumListPtr.get() ||
        size()!=dataInfo.msd.run.spectrumListPtr->size())
        throw runtime_error("[MSDataCache::update()] Usage error."); 

    SpectrumInfo& info = at(spectrum.index);
    info.update(spectrum);

    // MRU binary data caching
    if (impl_->config.binaryDataCacheSize>0 && !info.data.empty())
    {
        // find and erase if we're already on the list
        Impl::MRU::iterator it = find(impl_->mru.begin(), impl_->mru.end(), &info);
        if (it!=impl_->mru.end()) 
            impl_->mru.erase(it);

        // put us at the front of the list
        impl_->mru.push_front(&info);

        // free binary data from the least recently used SpectrumInfo (back of list)
        if (impl_->mru.size() > impl_->config.binaryDataCacheSize)
        {
            SpectrumInfo* lru = impl_->mru.back();
            lru->clearBinaryData();
            impl_->mru.pop_back();
        }
    }
}


PWIZ_API_DECL const SpectrumInfo& MSDataCache::spectrumInfo(size_t index)
{
    if (!impl_->spectrumListPtr.get() ||
        size()!=impl_->spectrumListPtr->size())
        throw runtime_error("[MSDataCache::spectrumInfo()] Usage error."); 

    SpectrumInfo& info = at(index);

    // update cache if necessary 
    if (info.index == (size_t)-1)
    {
        SpectrumPtr spectrum = impl_->spectrumListPtr->spectrum(index); 
        info.update(*spectrum);
    }

    return info;
}


} // namespace analysis 
} // namespace pwiz

