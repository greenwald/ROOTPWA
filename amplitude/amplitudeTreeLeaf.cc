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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      TTree leaf persistency storage class for amplitude information
//      needed by fit program
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "amplitudeTreeLeaf.h"

	
#if AMPLITUDETREELEAF_ENABLED


using namespace std;
using namespace rpwa;


ClassImp(amplitudeTreeLeaf);


amplitudeTreeLeaf::amplitudeTreeLeaf()
	: TObject()
{ }


amplitudeTreeLeaf::~amplitudeTreeLeaf()
{ }


// void amplitudeTreeLeaf::clear()
// {
// 	_subAmps.clear();
// }


#endif  // AMPLITUDETREELEAF_ENABLED
