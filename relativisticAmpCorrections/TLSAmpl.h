#ifndef TLSAmpl_h
#define TLSAmpl_h

#include "JLSIndices.h"
#include "TJwfTensor.h"

#include <array>

/// \class TLSAmpl
/// \brief Relativistic LS-coupling amplitudes
/// \author Jan.Friedrich@ph.tum.de, Daniel Greenwald
class TLSAmpl : public JLSIndices
{
public:

    TLSAmpl(const JLSIndices& jlsi, const std::array<unsigned,2>& j);
    
    const TTensorSum& tensorSum() const
    { return TSScalar_; }

    /* const TTensorTerm& GetTerm(const long& i) const */
    /* { return _TSScalar.GetTerm(i); } */

private:

    TTensorSum TSScalar_;

    static unsigned int _debugLevel;
};

#endif
