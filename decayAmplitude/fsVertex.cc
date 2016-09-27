///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      class that describes final state vertex decay topology
//      class is just used for internal book keeping
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "conversionUtils.hpp"
#include "fsVertex.h"


using namespace std;
using namespace rpwa;


bool fsVertex::_debug = false;


fsVertex::fsVertex(const particlePtr& fsParticle)
    : interactionVertex()
{
    if (not fsParticle) {
        printErr << "null pointer to final state particle. Aborting..." << endl;
        throw;
    }
    interactionVertex::addInParticle(fsParticle);
    if (_debug)
        printDebug << "constructed " << *this << endl;
}


fsVertex::fsVertex(const fsVertex& vert)
    : interactionVertex()
{
    *this = vert;
}


fsVertex::~fsVertex()
{ }


fsVertex*
fsVertex::doClone(const bool cloneInParticles,
                  const bool) const
{
    fsVertex* vertexClone = new fsVertex(*this);
    if (cloneInParticles)
        vertexClone->cloneInParticles();
    if (_debug)
        printDebug << "cloned " << *this << "; " << this << " -> " << vertexClone << " "
                   << ((cloneInParticles ) ? "in" : "ex") << "cluding incoming particles" << std::endl;
    return vertexClone;
}


bool
fsVertex::addInParticle(const particlePtr&)
{
    if (_debug)
        printWarn << "cannot add incoming particle to " << *this << endl;
    return false;
}


bool
fsVertex::addOutParticle(const particlePtr&)
{
    if (_debug)
        printWarn << "cannot add outgoing particle to " << *this << endl;
    return false;
}


ostream&
fsVertex::print(ostream& out) const
{
    out << name() << ": " << fsParticle()->qnSummary();
    return out;
}


ostream&
fsVertex::dump(ostream& out) const
{
    out << name() << ":" << endl
        << "    final state particle: " << *fsParticle() << endl;
    return out;
}


ostream&
fsVertex::printPointers(ostream& out) const
{
    out << name() << " " << this << ": final state particle: " << fsParticle() << endl;
    return out;
}
