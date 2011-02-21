//
// String.hpp
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars Sinai Medical Center, Los Angeles, California  90048
// Copyright 2008 Vanderbilt University - Nashville, TN 37232
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

#ifndef _STRING_HPP_
#define _STRING_HPP_

#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include "utility/misc/optimized_lexical_cast.hpp"

using std::string;
using std::getline;
using std::stringstream;
using std::istringstream;
using std::ostringstream;

namespace bal = boost::algorithm;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::format;

#endif // _STRING_HPP_
