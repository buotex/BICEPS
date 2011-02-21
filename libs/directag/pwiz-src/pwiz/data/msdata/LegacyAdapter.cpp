//
// LegacyAdapter.cpp
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
//


#define PWIZ_SOURCE

#include "LegacyAdapter.hpp"
#include "CVTranslator.hpp"
#include "boost/lambda/lambda.hpp"
#include "boost/lambda/bind.hpp"
#include <iostream>
#include <stdexcept>


namespace pwiz {
namespace msdata {


using namespace std;
using namespace boost::lambda;


//
// LegacyAdapter_Instrument::Impl
//


struct LegacyAdapter_Instrument::Impl
{
    Impl(InstrumentConfiguration& _instrumentConfiguration, const CVTranslator& _cvTranslator) 
    :   instrumentConfiguration(_instrumentConfiguration), cvTranslator(_cvTranslator)
    {}

    InstrumentConfiguration& instrumentConfiguration;
    const CVTranslator& cvTranslator;

    string get(const ParamContainer& paramContainer, 
               CVID cvid, 
               const string& userParamName);

    void set(ParamContainer& paramContainer, 
             CVID cvid, 
             const string& userParamName, 
             const string& value);
};


namespace {


void removeCVParams(vector<CVParam>& cvParams, CVID cvid)
{
    cvParams.erase(
    remove_if(cvParams.begin(), 
              cvParams.end(), 
              CVParamIsChildOf(cvid)),
              cvParams.end());
}


void removeUserParams(vector<UserParam>& userParams, const string& name)
{
    userParams.erase(
    remove_if(userParams.begin(), 
              userParams.end(), 
              boost::lambda::bind(std::equal_to<string>(), boost::lambda::bind(&UserParam::name,_1), name) ),
              userParams.end());
}


} // namespace


string LegacyAdapter_Instrument::Impl::get(const ParamContainer& paramContainer, CVID cvid, const string& userParamName)
{
    // cvParam
    CVParam param = paramContainer.cvParamChild(cvid);
    if (param.cvid != CVID_Unknown) 
        return param.name(); 

    // userParam
    string result = paramContainer.userParam(userParamName).value;
    if (result.empty()) result = "Unknown"; 
    return result;
}


void LegacyAdapter_Instrument::Impl::set(ParamContainer& paramContainer, CVID cvid, const string& userParamName, const string& value)
{
    // remove existing params
    removeCVParams(paramContainer.cvParams, cvid);
    removeUserParams(paramContainer.userParams, userParamName);

    // try to translate to cvParam
    CVID result = cvTranslator.translate(value);
    if (cvIsA(result, cvid))
    {
        paramContainer.cvParams.push_back(result);
        return;
    }

    // otherwise encode as userParam
    paramContainer.userParams.push_back(UserParam(userParamName, value));
}


//
// LegacyAdapter_Instrument
//


PWIZ_API_DECL
LegacyAdapter_Instrument::LegacyAdapter_Instrument(InstrumentConfiguration& instrumentConfiguration,
                                                   const CVTranslator& cvTranslator)
:   impl_(new Impl(instrumentConfiguration, cvTranslator)) 
{}


PWIZ_API_DECL string LegacyAdapter_Instrument::manufacturer() const
{
    // look first for cvParam

    CVParam model = impl_->instrumentConfiguration.cvParamChild(MS_instrument_model);
    if (model.cvid != CVID_Unknown)
    {
        // get the parent term
        const CVInfo& modelInfo = cvinfo(model.cvid);
        if (modelInfo.parentsIsA.empty())
            throw runtime_error("[LegacyAdapter_Instrument::manufacturer()] Model has no parents.");

        // s/ instrument model//
        string result = cvinfo(modelInfo.parentsIsA[0]).name;
        string::size_type index_suffix = result.find(" instrument model");
        if (index_suffix != string::npos) result.erase(index_suffix);
        return result;
    }

    // then try userParam

    string result = impl_->instrumentConfiguration.userParam("msManufacturer").value;
    if (result.empty()) result = "Unknown";
    return result;
}


PWIZ_API_DECL string LegacyAdapter_Instrument::model() const
{
    return impl_->get(impl_->instrumentConfiguration,
                      MS_instrument_model, 
                      "msModel");
}


PWIZ_API_DECL
void LegacyAdapter_Instrument::manufacturerAndModel(const string& valueManufacturer,
                                                    const string& valueModel)
{
    // remove existing params
    removeCVParams(impl_->instrumentConfiguration.cvParams, MS_instrument_model);
    removeUserParams(impl_->instrumentConfiguration.userParams, "msManufacturer");
    removeUserParams(impl_->instrumentConfiguration.userParams, "msModel");

    // try to translate to cvParam
    CVID cvid = impl_->cvTranslator.translate(valueModel);
    if (cvIsA(cvid, MS_instrument_model))
    {
        impl_->instrumentConfiguration.cvParams.push_back(cvid);
        return;
    }

    // otherwise encode as userParam
    impl_->instrumentConfiguration.userParams.push_back(UserParam("msManufacturer", valueManufacturer));
    impl_->instrumentConfiguration.userParams.push_back(UserParam("msModel", valueModel));
}


PWIZ_API_DECL string LegacyAdapter_Instrument::ionisation() const
{
    return impl_->get(impl_->instrumentConfiguration.componentList.source(0), 
                      MS_ionization_type, 
                      "msIonisation");
}


PWIZ_API_DECL void LegacyAdapter_Instrument::ionisation(const string& value)
{
    impl_->set(impl_->instrumentConfiguration.componentList.source(0), 
               MS_ionization_type, 
               "msIonisation", 
               value);
}


PWIZ_API_DECL string LegacyAdapter_Instrument::analyzer() const
{
    return impl_->get(impl_->instrumentConfiguration.componentList.analyzer(0), 
                      MS_mass_analyzer_type, 
                      "msMassAnalyzer");
}


PWIZ_API_DECL void LegacyAdapter_Instrument::analyzer(const string& value)
{
    impl_->set(impl_->instrumentConfiguration.componentList.analyzer(0), 
               MS_mass_analyzer_type, 
               "msMassAnalyzer", 
               value);
}


PWIZ_API_DECL string LegacyAdapter_Instrument::detector() const
{
    return impl_->get(impl_->instrumentConfiguration.componentList.detector(0), 
                      MS_detector_type, 
                      "msDetector");
}


PWIZ_API_DECL void LegacyAdapter_Instrument::detector(const string& value)
{
    impl_->set(impl_->instrumentConfiguration.componentList.detector(0), 
               MS_detector_type, 
               "msDetector", 
               value);
}


//
// LegacyAdapter_Software::Impl
//


struct LegacyAdapter_Software::Impl
{
    SoftwarePtr software;
    MSData& msd;
    const CVTranslator& cvTranslator;

    Impl(SoftwarePtr _software, MSData& _msd, const CVTranslator& _cvTranslator)
    :   software(_software), msd(_msd), cvTranslator(_cvTranslator)
    {
        if (!software.get())
            throw runtime_error("[LegacyAdapter_Software] Null SoftwarePtr.");
    }
};


//
// LegacyAdapter_Software
//


namespace {


ProcessingMethod& getProcessingMethod(SoftwarePtr software, MSData& msd)
{
    DataProcessingPtr dp;

    // find DataProcessing associated with Software
    for (vector<DataProcessingPtr>::const_iterator it=msd.dataProcessingPtrs.begin();
         it!=msd.dataProcessingPtrs.end(); ++it)
        if (it->get() && (*it)->softwarePtr.get()==software.get())
            dp = *it;

    // create a new DataProcessing if we didn't find one
    if (!dp.get())
    {
        const string& softwareID = software->id;
        if (softwareID.empty())
            throw runtime_error("[LegacyAdapter_Software::getProcessingMethod()] "
                                "Software::id not set.");

        dp = DataProcessingPtr(new DataProcessing(softwareID + " processing"));
        dp->softwarePtr = software; 
        msd.dataProcessingPtrs.push_back(dp);
    }

    // get ProcessingMethod 0
    if (dp->processingMethods.empty()) dp->processingMethods.push_back(ProcessingMethod());
    return dp->processingMethods[0];
}


string getProcessingMethodUserParamValue(const string& name, 
                                         const SoftwarePtr& software, 
                                         const MSData& msd)
{
    for (vector<DataProcessingPtr>::const_iterator it=msd.dataProcessingPtrs.begin();
         it!=msd.dataProcessingPtrs.end(); ++it)
    {
        if (!it->get() || (*it)->softwarePtr.get()!=software.get()) continue;

        for (vector<ProcessingMethod>::const_iterator jt=(*it)->processingMethods.begin();            
             jt!=(*it)->processingMethods.end(); ++jt)
        {
            UserParam result = jt->userParam(name);
            if (!result.empty())
                return result.value;
        }
    }

    return string();
}


} // namespace


LegacyAdapter_Software::LegacyAdapter_Software(SoftwarePtr software, 
                                               MSData& msd, 
                                               const CVTranslator& cvTranslator)
:   impl_(new Impl(software, msd, cvTranslator))
{}


PWIZ_API_DECL string LegacyAdapter_Software::name() const
{
    if (impl_->software->softwareParam.cvid != CVID_Unknown)
        return impl_->software->softwareParam.name();

    string result = getProcessingMethodUserParamValue("name", impl_->software, impl_->msd);
    return !result.empty() ? result : "unknown software name";
}


PWIZ_API_DECL void LegacyAdapter_Software::name(const string& value)
{
    impl_->software->softwareParam = impl_->cvTranslator.translate(value);

    if (impl_->software->softwareParam.cvid == CVID_Unknown)
    {
        ProcessingMethod& pm = getProcessingMethod(impl_->software, impl_->msd);
        removeUserParams(pm.userParams, "name");
        pm.userParams.push_back(UserParam("name", value));
    }
}


PWIZ_API_DECL string LegacyAdapter_Software::version() const
{
    return impl_->software->softwareParamVersion;
}


PWIZ_API_DECL void LegacyAdapter_Software::version(const string& value)
{
    impl_->software->softwareParamVersion = value;
}


PWIZ_API_DECL string LegacyAdapter_Software::type() const
{
    string result = getProcessingMethodUserParamValue("type", impl_->software, impl_->msd);
    return !result.empty() ? result : "unknown software type";
}


PWIZ_API_DECL void LegacyAdapter_Software::type(const string& value)
{
    ProcessingMethod& pm = getProcessingMethod(impl_->software, impl_->msd);
    removeUserParams(pm.userParams, "type");
    pm.userParams.push_back(UserParam("type", value));
}


} // namespace msdata
} // namespace pwiz

