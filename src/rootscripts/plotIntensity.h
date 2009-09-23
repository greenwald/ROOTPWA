///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Draws intensity graph for single wave from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef PLOTINTENSITY_HH
#define PLOTINTENSITY_HH


#include <string>

#include "TGraphErrors.h"
#include "TTree.h"


// signature with wave name
TGraphErrors*
plotIntensity(TTree*             tree,                    // TFitBin tree
	      const std::string& waveName,                // wave name
	      const std::string& selectExpr    = "",      // TTree::Draw() selection expression
	      const std::string& graphTitle    = "",      // name and title of graph (default is waveId)
	      const char*        drawOption    = "APZ",   // draw option for graph
	      const double       normalization = 1,       // scale factor for intensities
	      const int          graphColor    = kBlack,  // color of line and marker
	      const bool         saveEps       = false);  // if set, EPS file with name waveId is created

// signature with wave index
TGraphErrors*
plotIntensity(TTree*             tree,                    // TFitBin tree
	      const int          waveIndex,               // wave index
	      const std::string& selectExpr    = "",      // TTree::Draw() selection expression
	      const std::string& graphTitle    = "",      // name and title of graph (default is waveId)
	      const char*        drawOption    = "APZ",   // draw option for graph
	      const double       normalization = 1,       // scale factor for intensities
	      const int          graphColor    = kBlack,  // color of line and marker
	      const bool         saveEps       = false);  // if set, EPS file with name waveId is created


#endif  //PLOTINTENSITY_HH
