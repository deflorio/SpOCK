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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <Transformations.h>
#include <Constants.h>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>

extern "C"
      {
       #include <SpiceUsr.h>
      }

using namespace std;
using namespace math;
using namespace Eigen;
using namespace constants;
using namespace mathconst;
using namespace astro;

//------------------------------------------------------------------------------
// Vector6d J4model(double h, double inc, double e)
//------------------------------------------------------------------------------
/**
 * Non-spherical Earth gravitational model up to J4 zonal coefficient
 *
 * @param a     Mean semi-major axis
 * @param e     Mean eccentricity
 * @param inc   Mean inclination
 *
 * @return   6-D vector containing n, dn/dt, dOmega/dt, domega/dt, Ta, Td with
 *           n       Mean motion n = sqrt(mu/a^3) [rad/s]
 *           dn/dt   First derivative of n [rad/s^2]
 *           dOmega/dt   First derivative of right ascension of ascending node (RAAN) [rad/s]
 *           domega/dt   First derivative of argument of latitude [rad/s]
 *           Ta          Anomalistic period [rad/s]
 *           Td          Draconic period [rad/s]
 */
//------------------------------------------------------------------------------                      
Vector6d J4model(double a,
                 double e,
                 double inc)
                {
                Vector6d output;
                
                double n0, T0, p, e1, e2, e4, si, ci, Rp2, Rp4, J2_2, si2, si4;
                double n, dOMdt, domdt, delta_n, Ta, Td;
                
                n0 = sqrt( GM_EARTH/(a*a*a) );
                T0 = 2.0*PI/n0;

                // First derivatives of orbital elements

                p = a*(1.0 - e*e);
                e1 = sqrt(1.0 - e*e);
                
                si = sin(inc);
                ci = cos(inc);
                
                Rp2 = (R_EARTH/p)*(R_EARTH/p); Rp4 = (R_EARTH/p)*(R_EARTH/p)*(R_EARTH/p)*(R_EARTH/p);
                J2_2 = J2*J2;
                si2 = si*si; si4 = si*si*si*si;
                e2 = e*e; e4 = e*e*e*e;
                
                dOMdt = n0*(   J2*Rp2*ci*(-3.0/2.0)
                             + J2_2*Rp4*ci*( ( -45.0/8.0 + (3.0/4.0)*e2 + (9.0/32.0)*e4 ) + ( 57.0/8.0 - (69.0/32.0)*e2 - (27.0/64.0)*e4 )*si2 )
                             + J4*Rp4*ci*( (15.0/4.0) - (105.0/16.0)*si2 )*( 1.0 + (3.0/2.0)*e2 )
                           );
                        
                        
                        
                domdt = n0*(   J2*Rp2*( 3.0 - (15.0/4.0)*si2 )
                             + J2_2*Rp4*( ( 27.0/2.0 - (15.0/16.0)*e2 - (9.0/16.0)*e4 ) + ( -507.0/16.0 + (171.0/31.0)*e2 + (99.0/64.0)*e4 )*si2 + ( 1185.0/64.0 - (675.0/128.0)*e2 - (135.0/128.0)*e4 )*si4 )
                             + J4*Rp4*( ( -(3.0/8.0) + (15.0/8.0)*si2 - (105.0/64.0)*si4 )*( 10.0 + (15.0/2.0)*e2 ) + ( -(15.0/4.0) + (165.0/16.0)*si2 - (105.0/16.0)*si4 )*( 1.0 + (3.0/2.0)*e2 ) )
                           );        
                
                
                delta_n = n0*(   J2*Rp2*e1*(3.0/4.0)*( 2.0 - 3.0*si2)*( 1.0 + J2*Rp2*(1.0/8.0)*( ( 10.0 + 5.0*e2 + 8.0*e1 ) - ( 65.0/6.0 - (25.0/12.0)*e2 + 12.0*e1 )*si2 ) )
                               - J2_2*Rp4*e1*(5.0/64.0)*( 2.0 - e2 )*si2
                               - J4*Rp4*e1*(45.0/128.0)*e2*( 8.0 - 40.0*si2 + 35.0*si4 )
                             );
                
                n = n0 + delta_n;

                // Anomalistic period
                Ta = T0/( 1 + delta_n/n0 );
                // Draconic period
                Td = Ta/( 1 + domdt/n );
                
                // Put in output vector
                
                output << n, delta_n, dOMdt, domdt, Ta, Td;
            
                return(output);
                };
