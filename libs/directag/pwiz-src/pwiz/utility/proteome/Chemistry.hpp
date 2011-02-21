//
// Chemistry.hpp 
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


#ifndef _CHEMISTRY_HPP_
#define _CHEMISTRY_HPP_


#include "utility/misc/Export.hpp"
#include <iosfwd>
#include <string>
#include <memory>
#include <vector>
#include <map>


namespace pwiz {
namespace proteome {
namespace Chemistry {


/// struct for holding isotope information
struct PWIZ_API_DECL MassAbundance
{
    double mass;
    double abundance;

    MassAbundance(double m = 0, double a = 0)
    :   mass(m), abundance(a)
    {}

    bool operator==(const MassAbundance& that) const;
    bool operator!=(const MassAbundance& that) const;
};


/// struct for holding isotope distribution
typedef std::vector<MassAbundance> MassDistribution;


PWIZ_API_DECL std::ostream& operator<<(std::ostream& os, const MassAbundance& ma);
PWIZ_API_DECL std::ostream& operator<<(std::ostream& os, const MassDistribution& md);

 
/// scope for declarations related to elements
namespace Element {


/// enumeration of the elements
enum PWIZ_API_DECL Type
{
    H, He, Li, Be, B, C, N, O, F, Ne, 
    Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, 
    Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, 
    Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, 
    Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, 
    Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, 
    Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, 
    Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, 
    Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, 
    Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, 
    Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Uun, 
    Uuu, Uub, Uuq, Uuh
};


PWIZ_API_DECL std::ostream& operator<<(std::ostream& os, Type type);


/// class for obtaining information about elements
class PWIZ_API_DECL Info
{
    public:

    Info();
    ~Info();

    /// structure for holding info about an element
    struct Record
    {
        Type type;
        std::string symbol;
        int atomicNumber;
        double atomicWeight;
        MassDistribution isotopes;
    };

    /// retrieve the record for an element
    const Record& operator[](Type type) const;

    private:

    /// hidden implementation
    class Impl;
    std::auto_ptr<Impl> impl_;

    // don't allow copying
    Info(const Info&);
    const Info& operator=(const Info&);
};


PWIZ_API_DECL std::ostream& operator<<(std::ostream& os, const Info::Record& record);


} // namespace Element


/// class to represent a chemical formula
class PWIZ_API_DECL Formula
{
    public:

    /// formula string given by symbol/count pairs, e.g. water: "H2 O1" (whitespace optional)  
    Formula(const std::string& formula = "");
    Formula(const Formula& formula);
    const Formula& operator=(const Formula& formula);
    ~Formula();

    double monoisotopicMass() const;
    double molecularWeight() const;
    std::string formula() const;

    /// access to the Element's count in the formula
    int operator[](Element::Type e) const;
    int& operator[](Element::Type e);

    // direct access to the map, for iteration
    typedef std::map<Element::Type, int> Map;
    const Map& data() const;

    // operations
    Formula& operator+=(const Formula& that);    
    Formula& operator-=(const Formula& that);    
    Formula& operator*=(int scalar);
    
    private:
    class Impl;
    std::auto_ptr<Impl> impl_;
};


PWIZ_API_DECL Formula operator+(const Formula& a, const Formula& b);
PWIZ_API_DECL Formula operator-(const Formula& a, const Formula& b);
PWIZ_API_DECL Formula operator*(const Formula& a, int scalar);
PWIZ_API_DECL Formula operator*(int scalar, const Formula& a);


/// output a Formula
PWIZ_API_DECL std::ostream& operator<<(std::ostream& os, const Formula& formula);


} // namespace Chemistry 
} // namespace proteome
} // namespace pwiz


#endif // _CHEMISTRY_HPP_

