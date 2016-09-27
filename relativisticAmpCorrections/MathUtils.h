#ifndef __MathUtils__h
#define __MathUtils__h

#include <string>
#include <vector>

const bool is_odd(int v)
{ return v & 0x1; }

const bool is_even(int v)
{ return !is_odd(v); }

inline std::string parity_to_string(int p)
{ return (p > 0) ? "+" : "-"; }

// \return vector of spins from |s1-s2| to (s1 + s2)
const std::vector<unsigned> triangle(unsigned s1, unsigned s2);


#endif
