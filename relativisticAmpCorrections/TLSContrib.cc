#include "TLSContrib.h"

#include <reportingUtils.hpp>

#include <iostream>

using namespace std;

bool TLSContrib::Debug_ = false;

//-------------------------
TLSContrib::TLSContrib(const TLSAmpl& A, int delta, const TFracNum& scfac)
    : JLS(A),
      ContractionNumber_(A.contractionNumber()),
      Delta_(delta),
      SpinCG_(scfac)
{
    if (Debug_)
        cout << "(" << J() << ")" << L() << S() << "[" << Delta_ << "]" << endl;

    bool signReversal = Delta_ < 0 and is_odd(L() + S() + J());

    PolynomialTerms_.reserve(A.tensorSum().terms().size());
    for (const auto& t : A.tensorSum().terms()) {

        if (Debug_)
            cout << scfac.FracStringSqrt() << ","
                 << t.GetPreFac().FracStringSqrt()
                 << " L" << L() << "S" << S() << "J" << J()
                 << "Ampldelta" << A.delta()
                 << " delta" << Delta_ << ", sigrev" << signReversal;

        PolynomialTerms_.push_back(scfac * t.GetPreFac(), {t.GetGamS(), t.GetGamSig()});
        if (signReversal)
            PolynomialTerms_.back().PrefactorSquared.FlipSign();
        
        NormFactor_ = NormFactor_.SumSignedRoots(PolynomialTerms_.back().PrefactorSquared);
        if (Debug_)
            cout << " -> Normfactor: " << NormFactor_.FracStringSqrt() << endl;
    }

    PureRelativistic_ = (NormFactor_ == TFracNum::Zero);
    
    if (PureRelativistic_)
        NormFactor_ = TFracNum::One;
    else {
        auto  NormInv = invert(NormFactor_);
        for (auto& p : PolynomialTerms_)
            p.PrefactorSquared *= NormInv;
    }
}

//-------------------------
void TLSContrib::Add(const TLSContrib& rhs, bool particleExchange)
{
    if (_J != rhs._J or _L != rhs._L or _S != rhs._S) {
        printErr << "TLSContrib::Add : Something is wrong, trying to add different"
                 << " (J;L,S): (" << _J << ";" << _L << "," << _S << ") != ("
                 << rhs._J << ";" << rhs._L << "," << rhs._S << ")" << endl;
        throw;
    }

    //
    // Include normalisation factor and the factor (1/2)**2 in the squared
    // representation of prefactors
    //

    for (size_t i = 0; i < _polynomialTerms.size(); i++) {
        _polynomialTerms[i].squareOfPrefactor *= TFracNum::Quarter * _NormFactor;
    }

    for (size_t ib = 0; ib < rhs.GetNterms(); ib++) {
        bool termSummed = false;
        for (size_t i = 0; i < _polynomialTerms.size(); i++) {
            if (not termSummed and _cNum == rhs._cNum                                            and
                    (
                        (particleExchange                                                              and
                         _polynomialTerms[i].exponentOfGammaS     == rhs.getPolynomialTerms()[ib].exponentOfGammaSigma and
                         _polynomialTerms[i].exponentOfGammaSigma == rhs.getPolynomialTerms()[ib].exponentOfGammaS
                        )                                                                              or
                        (
                            not particleExchange                                                          and
                            _polynomialTerms[i].exponentOfGammaS     == rhs.getPolynomialTerms()[ib].exponentOfGammaS     and
                            _polynomialTerms[i].exponentOfGammaSigma == rhs.getPolynomialTerms()[ib].exponentOfGammaSigma
                        )
                    ) ) {
                termSummed = true;
                TFracNum bterm = TFracNum::Quarter * rhs._NormFactor * rhs.getPolynomialTerms()[ib].squareOfPrefactor;

                if (_J % 2) {
                    bterm.FlipSign();
                }

                _polynomialTerms[i].squareOfPrefactor = bterm.SumSignedRoots(_polynomialTerms[i].squareOfPrefactor);

            }
        }
        if (not termSummed) {

            polynomialTerms factor(rhs.getPolynomialTerms()[ib]);
            factor.squareOfPrefactor *= TFracNum::Quarter * rhs._NormFactor;
            if (_J % 2) {
                if (not factor.squareOfPrefactor.FlipSign()) {
                    printErr << "Flipping sign failed. Aborting..." << endl;
                    throw;
                }
            }
            if (particleExchange) {
                factor.swapExponents();
            }
            _polynomialTerms.push_back(factor);
        }
    }

    //
    // Eliminate zero entries
    //
    for (size_t i = _polynomialTerms.size(); i > 0; --i) {
        if (_polynomialTerms[i - 1].squareOfPrefactor == TFracNum::Zero) {
            _polynomialTerms.erase(_polynomialTerms.begin() + (i - 1));
        }
    }

    //
    // Recalculate Normalization Factor
    //
    _NormFactor = TFracNum::Zero;
    for (size_t i = 0; i < _polynomialTerms.size(); i++) {
        _NormFactor = _NormFactor.SumSignedRoots(_polynomialTerms[i].squareOfPrefactor);
    }

    //
    // Apply normalization
    //
    if (_NormFactor == TFracNum::Zero) {
        _pureRelativistic = true;
        _NormFactor = TFracNum::One;
    } else {
        _pureRelativistic = false;
        TFracNum NormInv = _NormFactor;
        NormInv.Invert();
        for (size_t i = 0; i < _polynomialTerms.size(); i++) {
            _polynomialTerms[i].squareOfPrefactor *= NormInv;
        }
    }
}


void TLSContrib::Print() const
{
    if (_cNum == 1) {
        cout << "g" << "[" << _cNum << "] (";
    }
    if (_cNum == 2) {
        cout << "f" << "[" << _cNum << "] (";
    }
    if (_cNum >= 3) {
        cout << "h" << "[" << _cNum << "] (";
    }
    cout << _J << ")" << _L << _S << "( ";
    if (not _pureRelativistic) {
        cout << _NormFactor.FracStringSqrt() << " ) ( ";
    }
    if (_polynomialTerms.empty()) {
        cout << "0";
    }
    for (size_t iT = 0; iT < _polynomialTerms.size(); iT++) {
        const polynomialTerms& factor = _polynomialTerms[iT];
        cout << factor.squareOfPrefactor.FracStringSqrt() << " ";
        if (factor.exponentOfGammaS) {
            if (factor.exponentOfGammaS == 1) {
                cout << " gs";
            } else {
                cout << " gs^" << factor.exponentOfGammaS;
            }
        }
        if (factor.exponentOfGammaSigma) {
            if (factor.exponentOfGammaSigma == 1) {
                cout << " gsig";
            } else {
                cout << " gsig^" << factor.exponentOfGammaSigma;
            }
        }
    }
    cout << " )" << endl;
}


void TLSContrib::PrintNR() const
{
    cout << _NormFactor.FracStringSqrt();
    if (_cNum == 1) {
        cout << "*g";
    }
    if (_cNum == 2) {
        cout << "*f";
    }
    if (_cNum == 3) {
        cout << "*h";
    }
}


void TLSContrib::PrintNRG(TFracNum m) const
{
    cout << (_NormFactor * m).FracStringSqrt();
    if (_cNum == 1) {
        cout << "*g";
    }
    if (_cNum == 2) {
        cout << "*f";
    }
    if (_cNum == 3) {
        cout << "*h";
    }
}
