#ifndef TFhh_h
#define TFhh_h
/*!
  \class TFhh
  \brief Relativistic helicity-coupling structure

  The class describes relativistic helicity-coupling amplitudes, as
  described in the paper [PRD78(2008)], equation (6.6). It collects a field
  of terms TLSContrib with different LS. It also provides the relation to the
  non-relativistic Zeemach amplitudes.

  \author Jan.Friedrich@ph.tum.de
  */

#include "TLSAmpl.h"
#include "TLSContrib.h"
#include "TLSNonRel.h"

#include <array>
#include <string>
#include <vector>

class TFhh
{

public:

    // constructor
    TFhh(unsigned J, const std::array<unsigned, 2>& j, unsigned lambda, unsigned nu,
         const std::vector<TLSAmpl>& LSampl, bool even_contraction);

    enum class symmetry { nu_nu, nu_minus_nu };
    
    TFhh(const TFhh& sFhh, symmetry s);
    TFhh(const TFhh& sFhh, const TFhh& xFhh);

    const unsigned lambda() const
    { return Lambda_; }

    const unsigned nu() const
    { return Nu_; }

    const unsigned J() const
    { return J_; }
    
    const bool evenContraction() const
    { return EvenContraction_; }
    
    const std::vector<TLSContrib>& LSt() const
    { return LSt_; }

    const std::string& name() const
    { return Name_; }

    void NonRelLimit();
    void PrintNRG() const;
    void Print()    const;

private:

    unsigned J_;
    unsigned Lambda_;
    unsigned Nu_;
    bool EvenContraction_;
    std::string Name_;
    std::vector<TLSContrib> LSt_;
    std::vector<TLSNonRel> NRLSt_;

    static unsigned debugLevel_;

};

inline const bool is_nu_nu(const TFhh& f)
{ return f.lambda() == f.nu(); }

inline const bool is_nu_minus_nu(const TFhh& f)
{ return f.lambda() == -f.nu(); }

inline const bool nu_lambda_partners(const TFhh& f, const TFhh& g)
{ return f.lambda() == g.nu() and f.nu() == g.lambda(); }

#endif
