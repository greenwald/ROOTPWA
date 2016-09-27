#include "MathUtils.h"

//-------------------------
const std::vector<unsigned> triangle(unsigned s1, unsigned s2)
{
    unsigned smin = abs(s1 - s2);
    unsigned smax = s1 + s2;
    std::vector<unsigned> S;
    S.reserve(smax - smin + 1);
    for (auto s = smin; s <= smax; ++s)
        S.push_back(s);
    return S;
}
