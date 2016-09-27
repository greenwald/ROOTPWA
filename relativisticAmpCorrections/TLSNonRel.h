#ifndef TLSNONREL_HH
#define TLSNONREL_HH

class TLSContrib;

#include "JLS.h"

#include <vector>

/// \class TLSNonRel
/// \brief Non-relativistic LS-coupling contributions
/// \author Jan.Friedrich@ph.tum.de, Daniel Greenwald
class TLSNonRel : public JLS
{
public:
    /// Constructor
    TLSNonRel(const TLSContrib& C);

    void Add(const TLSContrib& C);
    void Print()  const;
    void PrintG() const;

private:

    std::vector<TLSContrib> RelLS_;
    TFracNum GnrPrefac_;

};

#endif
