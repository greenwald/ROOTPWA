#ifndef __JLSIndices__h
#define __JLSIndices__h

#include "JLS.h"

/// \class JLSIndices
/// \brief extends JLS with additional indices
/// \author Daniel Greenwald
class JLSIndices : public JLS
{
public:
    /// constructor
    JLSIndices(unsigned j, unsigned l, unsigned s)
        : JLS(j, l, s) {} 

    const unsigned psiInternal() const { return PsiInternal_; }
    unsigned& psiInternal() { return PsiInternal_; }

    const unsigned chiPhi() const { return ChiPhi_; }
    unsigned& chiPhi() { return ChiPhi_; }
    
    const unsigned psiChi() const { return PsiChi_; }
    unsigned& psiChi() { return PsiChi_; }
    
    const unsigned psiPhi() const { return PsiPhi_; }
    unsigned& psiPhi() { return PsiPhi_; }

    const unsigned chiOmega() const { return ChiOmega_; }
    unsigned& chiOmega() { return ChiOmega_; }
    
    const unsigned phiOmega() const { return PhiOmega_; }
    unsigned& phiOmega() { return PhiOmega_; }

    const unsigned phiEps() const { return PhiEps_; }
    unsigned& phiEps() { return PhiEps_; }
    
    const unsigned chiEps() const { return ChiEps_; }
    unsigned& chiEps() { return ChiEps_; }

    const unsigned contractionNumber() const { return ContractionNumber_; }
    unsigned& contractionNumber() {return ContractionNumber_; }

    const unsigned delta() const { return Delta_; }
    unsigned& delta() { return Delta_; }
    
private:

    /// \name Compared
    /// @{
    unsigned PsiInternal_{0};
    unsigned ChiPhi_{0};
    unsigned PsiChi_{0};
    unsigned PsiPhi_{0};
    unsigned ChiOmega_{0};
    unsigned PhiOmega_{0};
    unsigned PhiEps_{0};
    unsigned ChiEps_{0};
    /// @}
    
    /// \name Not compared
    /// @{
    unsigned ContractionNumber_{0};
    unsigned Delta_{0};
    /// @}

};

/// equality operator
inline const bool operator==(const JLSIndices& lhs, const JLSIndices& rhs)
{
    return static_cast<const JLS&>(lhs) == static_cast<const JLS&>(rhs)
    and lhs.psiInternal() == rhs.psiInternal()
    and lhs.chiPhi() == rhs.chiPhi()
    and lhs.psiChi() == rhs.psiChi()
    and lhs.psiPhi() == rhs.psiPhi()
    and lhs.chiOmega() == rhs.chiOmega()
    and lhs.phiOmega() == rhs.phiOmega()
    and lhs.phiEps() == rhs.phiEps()
    and lhs.chiEps() == rhs.chiEps();
}

#endif
