//
// CVTranslator.cpp
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

#include "CVTranslator.hpp"
#include "boost/lexical_cast.hpp"
#include <map>
#include <iostream>
#include <sstream>
#include <iterator>
#include <stdexcept>


namespace pwiz {
namespace msdata {


using namespace std;
using boost::lexical_cast;


//
// default extra translations
//


namespace {

struct ExtraEntry
{
    const char* text;
    CVID cvid;
};

ExtraEntry defaultExtraEntries_[] =
{
    {"ITMS", MS_ion_trap},
    {"FTMS", MS_FT_ICR},
};

size_t defaultExtraEntriesSize_ = sizeof(defaultExtraEntries_)/sizeof(ExtraEntry);

} // namespace


//
// CVTranslator::Impl
//


class CVTranslator::Impl
{
    public:

    Impl();
    void insert(const string& text, CVID cvid);
    CVID translate(const string& text) const;

    private:

    typedef map<string,CVID> Map;
    Map map_;

    void insertCVTerms();
    void insertDefaultExtraEntries();
};


CVTranslator::Impl::Impl()
{
    insertCVTerms();
    insertDefaultExtraEntries();
}


namespace {


inline char alnum_lower(char c)
{
    // c -> lower-case or whitespace 
    return isalnum(c) ? static_cast<char>(tolower(c)) : ' ';
}


string preprocess(const string& s)
{
    string result = s;
    transform(result.begin(), result.end(), result.begin(), alnum_lower);
    return result;
}


string canonicalize(const string& s)
{
    // remove non-alnum characters
    istringstream iss(preprocess(s));

    // remove whitespace around tokens
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

    // concatenate with underscores
    ostringstream oss; 
    copy(tokens.begin(), tokens.end(), ostream_iterator<string>(oss, "_"));

    return oss.str();
}


} // namespace


bool shouldIgnore(const string& key, CVID value, CVID cvid)
{
    return (key=="unit_" && value==MS_unit && cvid==UO_unit ||
            key=="mass_unit_" && value==MS_mass_unit && cvid==UO_mass_unit ||
            key=="time_unit_" && value==MS_time_unit && cvid==UO_time_unit ||
            key=="energy_unit_" && value==MS_energy_unit && cvid==UO_energy_unit ||
            key=="pi_" && value==MS_PI && cvid==UO_pi); // MS_PI==photoionization, UO_pi==3.14

}


bool shouldReplace(const string& key, CVID value, CVID cvid)
{
    return (key=="second_" && value==MS_second && cvid==UO_second ||
            key=="minute_" && value==MS_minute && cvid==UO_minute ||
            key=="dalton_" && value==MS_Dalton && cvid==UO_dalton);
}


void CVTranslator::Impl::insert(const string& text, CVID cvid)
{
    string key = canonicalize(text);

    if (map_.count(key))
    {
        if (shouldIgnore(key, map_[key], cvid))
            return;

        if (!shouldReplace(key, map_[key], cvid))
        {
            throw runtime_error("[CVTranslator::insert()] Collision: " + 
                                lexical_cast<string>(map_[key]) + " " +
                                lexical_cast<string>(cvid));
        }
    }

    map_[key] = cvid;
}


CVID CVTranslator::Impl::translate(const string& text) const
{
    Map::const_iterator it = map_.find(canonicalize(text));
    if (it != map_.end())
        return it->second; 
    return CVID_Unknown;
}


void CVTranslator::Impl::insertCVTerms()
{
    for (vector<CVID>::const_iterator cvid=cvids().begin(); cvid!=cvids().end(); ++cvid)
    {
        if (cvIsA(*cvid, MS_purgatory)) continue;

        const CVInfo& info = cvinfo(*cvid);

        // insert name
        insert(info.name, *cvid);

        // insert synonyms
        if (*cvid < 100000000) // prefix == "MS"
        {
            for (vector<string>::const_iterator syn=info.exactSynonyms.begin(); 
                 syn!=info.exactSynonyms.end(); ++syn)
                insert(*syn, *cvid);
        }
    }
}


void CVTranslator::Impl::insertDefaultExtraEntries()
{
    for (const ExtraEntry* it=defaultExtraEntries_; 
         it!=defaultExtraEntries_+defaultExtraEntriesSize_; ++it)
        insert(it->text, it->cvid);
}


//
// CVTranslator
//


PWIZ_API_DECL CVTranslator::CVTranslator() : impl_(new Impl) {}
PWIZ_API_DECL void CVTranslator::insert(const string& text, CVID cvid) {impl_->insert(text, cvid);}
PWIZ_API_DECL CVID CVTranslator::translate(const string& text) const {return impl_->translate(text);}


} // namespace msdata
} // namespace pwiz

