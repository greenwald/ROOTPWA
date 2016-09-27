#include "TFhh.h"

#include "ClebschGordanBox.h"
#include "MathUtils.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

unsigned TFhh::debugLevel_ = 1;

//-------------------------
TFhh::TFhh(unsigned J, const std::array<unsigned, 2>& j, unsigned lambda, unsigned nu,
           const std::vector<TLSAmpl>& LSampl, bool even_contraction)
    : J_(J),
      Lambda_(lambda),
      Nu_(nu),
      EvenContraction_(even_contraction),
      Name_("F_" + std::to_string(Lambda_) + "_" + std::to_string(Nu_))
{
    for (const auto& A : LSampl) {
        if (A.delta() != abs(Lambda_ - Nu_))
            continue;

        auto SpinCouplFac = ClebschGordanBox::instance()->GetCG(A.S(), j[0], j[1])[ClebschGordanBox::CGIndex(j[0], Lambda_, j[1], -Nu_)];

        if (SpinCouplFac != TFracNum::Zero)
            LSt_.emplace_back(A, Lambda_ - Nu_, SpinCouplFac);
        else if (debugLevel_ >= 2)
            std::cout << "Clebsch-Gordan is zero" << std::endl;
    }
        
    if (debugLevel_ > 0)
        Print();
}

//-------------------------
TFhh::TFhh(const TFhh& sFhh, symmetry s)
{
    *this = sFhh;

    if (debugLevel_)
        std::cout << "Initializing from single Amplitude" << std::endl;

    if ((s == symmetry::nu_nu and is_odd(sFhh.J())) or (s == symmetry::nu_minus_nu and is_even(sFhh.J()))) {
        std::cout << sFhh.name() << "[symm] = 0" << std::endl;
        LSt_.clear();        
    }
    else {
        Name_ += "[symm]";
        if (debugLevel_)
            Print();
    }
}

//-------------------------
void enforce_cancellations(std::vector<TLSContrib>& V)
{
    // search for cancelations
    for (auto it = V.begin() + 1; it != V.end(); ) {
        // look for duplicate parameters before current iterator
        auto it2 = std::find(V.begin(), it, *it);
        // if found
        if (it2 != it) {
            // add to find
            it2->Add(*it, false);
            // and remove
            it = V.erase(it);
        } else
            ++it;
    }
}

//-------------------------
TFhh::TFhh(const TFhh& sFhh, const TFhh& xFhh)
{
    if (sFhh.J() != xFhh.J() or sFhh.evenContraction() != xFhh.evenContraction() or
        sFhh.lambda() != xFhh.nu() or sFhh.nu() != xFhh.lambda())
        throw std::runtime_error("TFhh::TFhh(const TFhh&, const TFhh&): parameter mismatch");
    
    *this = sFhh;
    NRLSt_.clear();
    Name_ += "[symm]";

    enforce_cancellations(LSt_);

    auto temp = xFhh.LSt_;
    enforce_cancellations(temp);
    LSt_.reserve(LSt_.size() + temp.size());
    size_t s = LSt_.size();
    for (const auto& t : temp) {
        // look for duplicate parameters
        auto it = std::find(LSt_.begin(), LSt_.begin() + s, t);
        // if found
        if (it != LSt_.end())
            it->Add(t, true);
        else
            LSt_.emplace_back(t, true);
    }

    LSt_.erase(std::remove_if(LSt_.begin(), LSt_.end(), [](const TLSContrib& c){return c.polynomialTerms().empty();}),
               LSt_.end());

    Print();
}

//-------------------------
void TFhh::NonRelLimit()
{
    for (const auto& t : LSt_) {
        if (t.pureRelativistic())
            continue;
        auto it = std::find(NRLSt_.begin(), NRLSt_.end(), t);
        if (it != NRLSt_.end())
            it->Add(t);
        else
            NRLSt_.emplace_back(t);
    }
    std::cout << Name_ << " (NR) = " << std::endl;
    for (const auto& t : NRLSt_) {
        std::cout << "LS=" << t.L() << t.S() << " => ";
        t.Print();
    }
}

//-------------------------
void TFhh::PrintNRG() const
{
    std::cout << Name_ << " (NR) = " << std::endl;
    for (const auto& t : NRLSt_) {
        std::cout << "LS=" << t.L() << t.S() << " => ";
        t.PrintG();
    }
}

//-------------------------
void TFhh::Print() const
{
    std::cout << Name_ << " =" << (EvenContraction_ ? "" : " (iw)") << std::endl;
    for (const auto& t : LSt_)
        t.Print();
    if (LSt_.empty())
        std::cout << " 0" << std::endl;
}
