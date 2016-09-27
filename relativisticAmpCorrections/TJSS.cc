#include "TJSS.h"

#include "MathUtils.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>

using namespace std;

unsigned int TJSS::_debugLevel = 9999;

//-------------------------
const std::string to_string(unsigned omega, unsigned epsilon, unsigned phi)
{
    return (phi == 0) ? std::string(omega, 'o') + std::string(epsilon, 'e')
        : std::string(phi, '#') + "(" + to_string(omega, epsilon, 0) + ")";
}


void TJSS::CalcAmpl()
{
    if (_debugLevel >= 2) {
        cout << "  Decay channel:   "
             << J_ << parity_to_string(P_)
             << " ->"
             << " " << j_[0] << parity_to_string(p_[0])
             << " " << j_[1] << parity_to_string(p_[1])
             << endl;
    }

    auto S = triangle(j_[0], j_[1]);

    if (_debugLevel >= 2) {
        cout << "possible S:";
        for (auto s : S)
            cout << " " << s;
        cout << endl;
    }

    auto intr_parity = P_ * p_[0] * p_[1];

    std::set<unsigned> L;
    for (auto s : S)
        for (auto l : triangle(J_, s))
            if (l % 2 == (intr_parity < 0))
                L.insert(l);
    
    if (_debugLevel >= 2) {
        cout << "possible L:";
        for (auto l : L)
            cout << " " << l;
        cout << endl;
    }

    const bool evenContraction = is_even(J_ + j_[0] + j_[1] + *L.begin());

    if (_debugLevel >= 2) {
        if (evenContraction)
            cout << "contraction only with g~" << endl;
        else
            cout << "contraction with eps*p and g~" << endl;
    }

    auto MaxPsiInternal = std::min(j_[0], j_[1]);
    auto MaxPsiPhi = std::min(j_[0] + j_[1], J_);

    for (auto l : L) {
         
        auto MaxPsiChi = std::min(j_[0] + j_[1], l);
        auto MaxChiPhi = std::min(J_, l);

        auto S_L = triangle(J_, l);
        S_L.erase(std::remove_if(S_L.begin(), S_L.end(), [&](unsigned s){return s < S.front() or s > S.back();}), S_L.end());

        for (auto s : S_L) {

            JLSIndices K(J_, l, s);

            if (_debugLevel >= 2)
                cout << "Amplitudes for L=" << K.L() << " S=" << K.S()
                     << "  Rank scheme [ " << j_[0] + j_[1] << " " << K.L() << " " << K.J() << "]" << endl;

            auto totalRank = K.J() + j_[0] + j_[1] + K.L();
            auto MaxDelta = std::min(s, K.J());

            auto IndexContractions = (totalRank - (evenContraction ? 0 : 3)) / 2;
            if (_debugLevel >= 2)
                cout << IndexContractions << " Lorentz contractions." << endl;
            
            unsigned MaxContractionNumber = 0;

            for (K.psiInternal() = 0; K.psiInternal() <= MaxPsiInternal; ++K.psiInternal()) {
                for (K.chiPhi() = 0; K.chiPhi() <= MaxChiPhi; ++K.chiPhi()) {
                    for (K.psiChi() = 0; K.psiChi() <= MaxPsiChi; ++K.psiChi()) {
                        for (K.psiPhi() = 0; K.psiPhi() <= MaxPsiPhi; ++K.psiPhi()) {

                            /////////////////////////
                            // check selection rules
                            if (_debugLevel >= 3)
                                cout << "Checking " << K.psiInternal() << " " << K.psiChi() << " " << K.chiPhi() << " " << K.psiPhi();
    
                            if ( K.psiInternal() + K.psiChi() + K.chiPhi() + K.psiPhi() != IndexContractions) {
                                if (_debugLevel >= 3)
                                    cout << " C-" << endl;
                                continue;
                            }
                            if (_debugLevel >= 3)
                                cout << " C+";
                            
                            if (evenContraction) {
                                if (2 * K.psiInternal() + K.psiChi() + K.psiPhi() != j_[0] + j_[1]) {
                                    if (_debugLevel >= 3)
                                        cout << "S-" << endl;
                                    continue;
                                }
                                if (_debugLevel >= 3)
                                    cout << "S+";
                                
                                if (K.psiChi() + K.chiPhi() != K.L()) {
                                    if (_debugLevel >= 3)
                                        cout << "L-" << endl;
                                    continue;
                                }
                                if (_debugLevel >= 3)
                                    cout << "L+";
                                
                                if (K.chiPhi() + K.psiPhi() != K.J()) {
                                    if (_debugLevel >= 3)
                                        cout << "J-" << endl;
                                    continue;
                                }
                                if (_debugLevel >= 3)
                                    cout << "J+";
                            } else {
                                if (K.L() - K.psiChi() - K.chiPhi() > 1) {
                                    if (_debugLevel >= 3)
                                        cout << "L-" << endl;
                                    continue;
                                }
                                if (_debugLevel >= 3)
                                    cout << "L+";
                                
                                if (K.J() - K.psiPhi() - K.chiPhi() > 1) {
                                    if (_debugLevel >= 3)
                                        cout << "J-" << endl;
                                    continue;
                                }
                                if (_debugLevel >= 3)
                                    cout << "J+";
                            }
                            /////////////////////////

                            /////////////////////////
                            // Get contraction ranks
                            std::array<unsigned, 2> ranks = {j_[0] - K.psiInternal(), j_[1] - K.psiInternal()};
                            if (!evenContraction) {
                                auto PsiRest = j_[0] + j_[1] - 2 * K.psiInternal() - K.psiChi() - K.psiPhi();
                                if (PsiRest < 0 or PsiRest > 2) {
                                    if (_debugLevel >= 3)
                                        cout << "R-" << endl;
                                    throw;
                                }
                                if (_debugLevel >= 3)
                                    cout << "R+";
                                
                                if (PsiRest == 2) {
                                    if (ranks[0] == 0) {
                                        if (_debugLevel >= 3)
                                            cout << "O-"; 
                                        continue;
                                    }
                                    if (ranks[1] == 0) {
                                        if (_debugLevel >= 3)
                                            cout << "E-"; 
                                        continue;
                                    }
                                    std::for_each(ranks.begin(), ranks.end(), [](unsigned& r){--r;});
                                }
                            }
                            /////////////////////////
                            
                            if (_debugLevel >= 3)
                                cout << "{" << ranks[0] << "}";

                            for (K.chiOmega() = 0; K.chiOmega() <= ranks[0] and K.chiOmega() <= K.psiChi(); ++K.chiOmega()) {
                                for (K.phiOmega() = 0; K.phiOmega() <= ranks[0] - K.chiOmega() and K.phiOmega() <= K.psiPhi(); ++K.phiOmega()) {

                                    K.chiEps() = K.psiChi() - K.chiOmega();
                                    K.phiEps() = K.psiPhi() - K.phiOmega();

                                    if (_debugLevel >= 3)
                                        cout << "[" << K.phiEps() << K.chiEps() << ranks[1] << "]";

                                    if (K.phiEps() + K.chiEps() > ranks[1])
                                        continue;

                                    if (_debugLevel >= 3)
                                        cout << "E+ OK" << endl;

                                    if (_debugLevel >= 2)
                                        cout << "Checking PsiInternal=" << K.psiInternal()
                                             << " PsiChi=" << K.psiChi()
                                             << " ChiPhi=" << K.chiPhi()
                                             << " PsiPhi=" << K.psiPhi()
                                             << endl
                                             << " PhiOmega=" << K.phiOmega()
                                             << " ChiOmega=" << K.chiOmega()
                                             << " PhiEps=" << K.phiEps()
                                             << " ChiEps=" << K.chiEps()
                                             << endl
                                             << "Contraction pattern "
                                             << to_string(K.phiOmega(), K.phiEps(), K.psiPhi())
                                             << " " << j_[0]  << j_[1] << std::string(K.psiInternal(), '\'') << " "
                                             << to_string(K.chiOmega(), K.chiEps(), K.psiChi())
                                             << " " << K.L() << " " << std::string(K.chiPhi(), '#') << " " << K.J() << " " 
                                             << to_string(K.phiOmega(), K.phiEps(), K.psiPhi())
                                             << endl;

                                    for (K.delta() = 0; K.delta() <= MaxDelta; ++K.delta()) {

                                        if (_debugLevel >= 2)
                                            cout << " Constructing LS-Amplitude " << LSAmplitudes_.size() << endl;

                                        auto it = std::find(LSAmplitudes_.begin(), LSAmplitudes_.end(), K);
                                        K.contractionNumber() = (it != LSAmplitudes_.end()) ? it->contractionNumber() : MaxContractionNumber + 1;

                                        LSAmplitudes_.emplace_back(K, j_);
                                        if (LSAmplitudes_.back().tensorSum().GetNterms() > 0)
                                            MaxContractionNumber = std::max(LSAmplitudes_.back().contractionNumber(), MaxContractionNumber);
                                        else
                                            LSAmplitudes_.pop_back();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (_debugLevel) {
        cout << LSAmplitudes_.size() << " LS-Amplitudes found to be non-zero." << endl;
        cout << "++++++++++++++++++++++++++++++++++++" << endl;
        cout << "+++ Helicity-coupling amplitudes +++" << endl;
        cout << "++++++++++++++++++++++++++++++++++++" << endl;
    }

    if (!FhhAmpl_.empty())
        throw std::runtime_error("FhhAmpl_ not empty.");

    for (unsigned lambda = 0; lambda <= j_[0]; ++lambda)
        for (unsigned nu = (lambda == 0) ? 0 : -j_[1]; nu <= j_[1]; ++nu)
            FhhAmpl_.emplace_back(J_, j_, lambda, nu, LSAmplitudes_, evenContraction);

    // remove amps with 0 terms
    FhhAmpl_.erase(std::remove_if(FhhAmpl_.begin(), FhhAmpl_.end(), [](const TFhh& f){return f.LSt().empty();}),
                   FhhAmpl_.end());
    
    // WHAT IS THE POINT OF THIS?
    // if identical particles:
    if (j_[0] == j_[1] and p_[0] == p_[1]) {

        if (_debugLevel)
            cout << endl << " for identical-particle decay:" << endl;

        FhhIdAmpl_.reserve(FhhAmpl_.size());
        for (const auto& f : FhhAmpl_)
            if (is_nu_nu(f))
                FhhIdAmpl_.emplace_back(f, TFhh::symmetry::nu_nu);
            else if (is_nu_minus_nu(f))
                FhhIdAmpl_.emplace_back(f, TFhh::symmetry::nu_minus_nu);
            else {
                auto it = std::find_if(FhhAmpl_.begin(), FhhAmpl_.end(), std::bind(nu_lambda_partners, f, std::placeholders::_1));
                if (it != FhhAmpl_.end())
                    FhhIdAmpl_.emplace_back(f, *it);
                else 
                    cerr << "No partner for amplitude " << f.name() << endl;
            }
    }
    
    cout << FhhAmpl_.size() << " amplitudes: non-relativistic limit" << endl;
    for (auto& f : FhhAmpl_)
        f.NonRelLimit();

    cout << "Check non-relativistic G's" << endl;
    for (const auto& f : FhhAmpl_)
        f.PrintNRG();
}

#if(0)
long TJSS::PrintHFILE()
{
    char DecayName[10];
    sprintf(DecayName, "%ld%ld%ld%c%c", J_, j_[0], j_[1],
            P_ * p_[0] * p_[1] == -1 ? 'n' : 'p', ' ');
    char ofname[20];
    sprintf(ofname, "CalcAmpl-%s.h", DecayName);
    ofstream ofs(ofname);
    ofs << "// CalcAmpl output for " << DecayName << endl;
    ofs << "const int FhhAmpl_" << DecayName << "[] = { " << endl;
    ofs << "  " << _NFhhAmpl << ",               // number of Fhh amplitudes"
        << endl;
    for (int i = 0; i < _NFhhAmpl; i++) {
        ofs << "  " << _FhhAmpl[i]->GetJ() << ", " << _FhhAmpl[i]->GetLambda()
            << ", " << _FhhAmpl[i]->GetNu() << ", "
            << _FhhAmpl[i]->GetEvenContr() << ",      // " << "F"
            << _FhhAmpl[i]->GetLambda() << _FhhAmpl[i]->GetNu()
            << ": J, lambda, nu, even_contr" << endl;
        ofs << "    " << _FhhAmpl[i]->GetNterms()
            << ",             // number of contributions" << endl;
        for (int j = 0; j < _FhhAmpl[i]->GetNterms(); j++) {
            TLSContrib* lsa = _FhhAmpl[i]->GetLStPtr()[j];
            ofs << "    " << lsa->GetJ() << ", " << lsa->GetL() << ", "
                << lsa->GetS() << ", " << lsa->GetDelta() << ", "
                << lsa->GetRunningNumber() << ",    // contr. F"
                << _FhhAmpl[i]->GetLambda() << _FhhAmpl[i]->GetNu() << "-"
                << j << ": J, L, S, delta, #[g,f,h,...]" << endl;
            ofs << "      " << lsa->GetNterms() << ", "
                << lsa->GetNormFactor()->GetSign() << ", "
                << lsa->GetNormFactor()->GetNumerator() << ", "
                << lsa->GetNormFactor()->GetDenominator()
                << ",     // number of terms, squared norm factor sign/N/D"
                << endl;
            for (int k = 0; k < lsa->GetNterms(); k++) {
                ofs << "      " << lsa->GetTermFracNum()[k].GetSign() << ", "
                    << lsa->GetTermFracNum()[k].GetNumerator() << ", "
                    << lsa->GetTermFracNum()[k].GetDenominator() << ", "
                    << lsa->GetTermg1pot()[k] << ", "
                    << lsa->GetTermg2pot()[k]
                    << ",  // squared sign/N/D, exponents of g_s and g_sigma"
                    << endl;
            }
        }
    }
    ofs << "};" << endl;

    return 0;
}
#endif
