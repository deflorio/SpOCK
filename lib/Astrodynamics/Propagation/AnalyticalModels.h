//==========================================================================
/*
*    Copyright 2020 Sergio De Florio
*    All rigths reserved
*
*    This file is part of SpOCK
* 
*    SpOCK is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation version 3
* 
*    SpOCK is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
* 
*    You should have received a copy of the GNU General Public License
*    along with SpOCK. If not, see <https://www.gnu.org/licenses/>.
*
*/
//==========================================================================

#ifndef ANALYTICALMOD_H_
#define ANALYTICALMOD_H_

#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>
      
using namespace std;
using namespace math;
 
//------------------------------------------------------------------------------
// Group of functions for propagation with analytical models
//------------------------------------------------------------------------------

// Non-spherical Earth gravitational model up to J4 spherical zonal coefficient
Vector6d J4model(double a, double e, double inc);

#endif // ANALYTICALMOD_H_