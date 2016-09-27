#include "TLSAmpl.h"

#include "TSpinWaveFunction.h"

#include <reportingUtils.hpp>

#include <iostream>
#include <string>

using namespace std;

unsigned int TLSAmpl::_debugLevel = 1;

//-------------------------
TLSAmpl::TLSAmpl(const JLSIndices& jlsi,const std::array<unsigned, 2>& j)
    : JLSIndices(jlsi)
{
    const long totalR = j[0] + j[1] + L() + J();
    const long contractions = 2 * (psiInternal() + psiChi() + chiPhi() + psiPhi());
    const bool evenContraction = (not (totalR % 2));

    if ( (evenContraction and contractions != totalR) or
         (!evenContraction and contractions != totalR - 3)) {
        printErr << "Invalid contraction occurred. Aborting..." << endl;
        throw;
    }

    if (_debugLevel) {
        cout << "LSAmpl: " << j[0] << " " << j[1] << " L=" << L() << " "
             << J() << " d=" << delta() << " S=" << S() << " c: " << psiInternal()
             << " " << psiChi() << " " << chiPhi() << " " << psiPhi() << " s: "
             << phiOmega() << " " << chiOmega() << " " << phiEps() << " "
             << chiEps();
        if (not evenContraction) {
            cout << " (iw)";
        }
        cout << endl;
    }
    
    TSpinWaveFunction WFS1(j[0], 's');
    TSpinWaveFunction WFS2(j[1], 's');

    TTensorSum TSS = WFS1.GetSpinCoupledTensorSum(WFS2, delta(), S());
    if (_debugLevel >= 2)
        TSS.Print();

    if (psiInternal()) {
        long ires = TSS.SpinInnerContraction(psiInternal());
        if (ires == 0) {
            if (_debugLevel)
                cout << "Inner contraction is zero." << endl;
            return;
        }
        if (_debugLevel >= 2)
            TSS.Print();
    }

    TSpinWaveFunction WFL(L(), 'l');
    TTensorSum TSL = WFL.GetTensorSum('c', 0);
    if (_debugLevel >= 2)
        TSL.Print();

    TTensorSum TSLS = TSS.LSContraction(TSL, psiChi(), chiOmega(), chiEps(), 'c');
    if (TSLS.GetNterms() == 0) {
        if (_debugLevel)
            cout << "LS contraction is zero." << endl;
        return;
    }
    if (_debugLevel >= 2)
        TSLS.Print();

    TSpinWaveFunction WFJ(J(), 'c');
    TTensorSum TSJ = WFJ.GetTensorSum('p', delta());
    if (_debugLevel >= 2)
        TSJ.Print();

    TTensorSum TSLSJ = TSLS.LSContraction(TSJ, psiPhi(), phiOmega(), phiEps(), 'p');

    if (TSLSJ.GetNterms() == 0) {
        if (_debugLevel)
            cout << "JS contraction is zero." << endl;
        return;
    }

    if (_debugLevel >= 2)
        TSLSJ.Print();

    TSScalar_ = TSLSJ.LJContraction(chiPhi(), evenContraction);

    if (_debugLevel >= 2)
        TSLSJ.Print();

    if (TSScalar_.GetNterms() == 0) {
        if (_debugLevel)
            cout << "LJ contraction is zero." << endl;
        return;
    }

    if (_debugLevel)
        TSScalar_.Print('s');
}
