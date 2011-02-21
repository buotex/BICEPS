//
// unit.hpp
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


#ifndef _UNIT_HPP_
#define _UNIT_HPP_


#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>


namespace pwiz {
namespace util {


//
// These are assertion macros for unit testing.  They throw a runtime_error 
// exception on failure, instead of calling abort(), allowing the application
// to recover and return an appropriate error value to the shell.
//
// unit_assert(x):                             asserts x is true
// unit_assert_equal(x, y, epsilon):           asserts x==y, within epsilon
// unit_assert_matrices_equal(A, B, epsilon):  asserts A==B, within epsilon
//


inline std::string unit_assert_message(const char* filename, int line, const char* expression)
{
    std::ostringstream oss;
    oss << "[" << filename << ":" << line << "] Assertion failed: " << expression; 
    return oss.str();
}


#define unit_assert(x) \
    (!(x) ? throw std::runtime_error(unit_assert_message(__FILE__, __LINE__, #x)) : 0) 


#define unit_assert_equal(x, y, epsilon) \
    unit_assert(fabs((x)-(y)) < (epsilon))


#define unit_assert_matrices_equal(A, B, epsilon) \
    unit_assert(boost::numeric::ublas::norm_frobenius((A)-(B)) < (epsilon))


#define unit_assert_vectors_equal(A, B, epsilon) \
    unit_assert(boost::numeric::ublas::norm_2((A)-(B)) < (epsilon))


} // namespace util
} // namespace pwiz


#endif // _UNIT_HPP_

