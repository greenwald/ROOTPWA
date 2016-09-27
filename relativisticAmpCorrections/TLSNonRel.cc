#include "TLSNonRel.h"

#include "ClebschGordanBox.h"

#include <iostream>
#include <stdexcept>

//-------------------------
TLSNonRel::TLSNonRel(const TLSContrib& C)
    : JLS(C), RelLS_(1, C),
      GnrPrefac_(TFracNum(2 * L() + 1, 2 * J() + 1)
                 * ClebschGordanBox::instance()->GetCG(J(), L(), S())[ClebschGordanBox::CGIndex(L(), 0, S(), C.delta())]
                 * C.spinCG())
{}

//-------------------------
void TLSNonRel::Add(const TLSContrib& C)
{
    if (static_cast<const JLS&>(*this) != static_cast<const JLS&>(C))
        throw std::runtime_error("TLSNonRel::Add not appropriate.");

    RelLS_.push_back(C);
}

//-------------------------
void TLSNonRel::Print() const
{
    std::cout << " [ " << GnrPrefac_.FracStringSqrt() << "  G_" << L() << S() << " ] ";
    for (const auto c : RelLS_)
        c.PrintNR();
    std::cout << std::endl;
}

//-------------------------
void TLSNonRel::PrintG() const
{
    std::cout << " [ G_" << L() << S() << " ] ";
    for (const auto& c : RelLS_)
        c.Print(invert(GnrPrefac_));
    std::cout << std::endl;
}
