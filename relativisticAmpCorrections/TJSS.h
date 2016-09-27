#ifndef TJSS_HH
#define TJSS_HH

#include "TFhh.h"
#include "TLSAmpl.h"

#include <array>
#include <vector>

class TJSS
{

public:

    // contsructor
    TJSS(unsigned J, int eJ, unsigned S1, int e1, unsigned S2, int e2)
        : J_(J), P_(eJ), j_({S1, S2}), p_({e1, e2}) {}

    std::vector<TFhh>& fhh()
    { return FhhAmpl_; }

    void CalcAmpl();

private:

    unsigned J_; // mother spin
    int P_; // mother parity;
    std::array<unsigned, 2> j_; // daughter spins
    std::array<int, 2> p_; // daughter parities

    std::vector<TLSAmpl> LSAmplitudes_;
    std::vector<TFhh>    FhhAmpl_;
    std::vector<TFhh>    FhhIdAmpl_;

    static unsigned int _debugLevel;

};

#endif
