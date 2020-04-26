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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <map>
#include <string>

#include <VarTypes.h>


namespace constants
           {
            namespace astro
                    {
                    // Astronomical Unit [m] (Standish, 1998)
                    const static double AU = 1.49597870691E11;
                    // Speed of light in meters per second (Standish, 1995)
                    const static double C_LIGHT = 299792458.0;
                    // Solar flux at a distance of 1 AU from the Sun [W/m^2]
                    const static double C_SUN = 1367.0;
                    // Solar radiation pressure constant [N/m^2] IERS 1996 (McCarthy 1996)
                    const static double P = 4.560E-6;
                    // Flattening WGS-84
                    const static double F_EARTH = 1.0/298.257223563;
                    // Gravitational coefficient Earth [m^3/s^2]; JGM3
                    const static double GM_EARTH = 398600.4415E+9;
                    // Gravitational coefficient Moon [m^3/s^2]; DE200
                    const static double GM_MOON = GM_EARTH/81.300587;
                    // Gravitational coefficient Sun [m^3/s^2]; IAU 1976
                    const static double GM_SUN = 1.32712438E+20;
                    // Standard gravitational acceleration at sea-level
                    const static double GRAV_ACC = 9.80665;
                    // Gravitational constant [m^3*kg*s^2] (Standish, 1995)
                    const static double GRAV_CONST = 6.67259E-11;
                    // Precomputed inverse-square of speed of light
                    const static double INV_C_LIGHT2 = 1.0/(C_LIGHT*C_LIGHT);
                    // Precomputed inverse 3rd power of speed of light
                    const static double INV_C_LIGHT3 = 1.0/(C_LIGHT*C_LIGHT*C_LIGHT);
                    // Second-order zonal coefficient
                    const static double J2 = 0.00108263;
                    // Third-order zonal coefficient
                    const static double J3 = -2.53615069E-6;
                    // Fourth-order zonal coefficient
                    const static double J4 = -1.61936355E-6;
                    // Earth rotation (derivative of GMST at J2000; differs from inertial period by precession) [rad/s] Aoki 1982, NIMA 1997
                    const static double OMEGA_EARTH = 7.2921158553E-5;
                    // Radius of Earth [m] WGS-84
                    const static double R_EARTH = 6378.137E3;
                    // Radius of Moon [m]
                    const static double R_MOON = 1738.0E3;          
                    // Radius of Sun [m]; Seidelmann 1992
                    const static double R_SUN = 696000.0E3;
                    }
            
            namespace timescales
                    {
                    // Julian year in Julian days (NASA, 2012)
                    const static double JULIAN_YEAR_DAYS = 365.25;
                    // Julian day [s] (NASA, 2012)
                    const static double JULIAN_DAY = 86400.0;
                    // Julian year [s]. Result of JULIAN_YEAR_IN_DAYS*JULIAN_DAY
                    const static double JULIAN_YEAR = 3.15576E7;
                    // J2000 epoch in GPS time (1 Jan 2000 11:59:08.816)
                    //const static double J2000_GPSSECS = 630763161.816;
                    // GPS seconds at 1 Jan 2000 11:58:55.816 UTC
                    const static double J2000_GPSSECS = 630763148.816;//630763213 + 19.0 + 32.184;
                    
                    // GPS seconds at 1 Jan 2000 12:00:00 UTC
                    //const static double J2000_UTC_GPSSECS = 630763213;
                    // GPS seconds at 1 Jan 2000 12:00:00 UTC - (19.0 + 32.184)
                    //const static double J2000_GPSSECS = 630763161.816;
                    
                    // Relative time rate difference between TCB and TDB time scales
                    const static double LB_TIME_RATE_TERM = 1.550519768E-8;
                    // Relative time rate difference between TCG and TT time scales
                    const static double LG_TIME_RATE_TERM = 6.969290134E-10;
                    // Modified Julian Date of J2000.0
                    const static double MJD_J2000 = 51544.5;
                    // Sidereal day [s] (NASA, 2012)
                    const static double SIDEREAL_DAY = 86164.09054;
                    // Sidereal year in quasar reference frame. Result of SIDEREAL_YEAR_IN_DAYS*JULIAN_DAY
                    const static double SIDEREAL_YEAR = 3.1558149504E7;
                    // Sidereal year in Julian days in quasar reference frame (NASA, 2012)
                    const static double SIDEREAL_YEAR_IN_DAYS = 365.25636;
                    // Conversion factor from seconds to minutes
                    const static double SEC2MIN = 60.0;
                    // Conversion factor from seconds to hours
                    const static double SEC2H = 3600.0;
                    // Conversion factor from seconds to day
                    const static double SEC2DAY = 86400.0;
                    }
                    
            namespace physics
                    {
                    // The Boltzmann constant (gas constant per particle) [m^2*kg/( s^2 * K)] (NIST, 2013)
                    const static double BOLTZMANN_CONST = 1.3806488E-23;
                    // Molar gas constant. The molar gas constant [J/{mol*K)] (NIST: http://physics.nist.gov/cgi-bin/cuu/Value?r, 2016)
                    const static double MOLAR_GAS_CONSTANT = 8.3144598;
                    // Planck constant [m^2*kg/s] (NIST, 2013)
                    const static double PLANCK_CONST = 6.62606957E-34;
                    // The specific gas constant of air [J/(kg*K)] (Anderson, 2006)
                    const static double SPEC_GAS_CONST_AIR = 2.87E2;
                    }
                    
            namespace mathconst
                    {
                    const static double PI = 3.141592653589793238;
                    const static long double LONG_PI = 3.14159265358979323846264338328L;
                    const static double PI2 = 2.0*PI;
                    const static double DEG2RAD = PI/180.0;
                    const static double RAD2DEG = 180.0/PI;
                    const static double RAD2ARCS = 3600.0*180.0/PI;
                    }
                    
            namespace materials
                      {
                      // The three coefficients are in order specular reflectivity, diffuse reflectivity and transmitted portions of incoming photons coefficients
                      const Vec3d SolArrays_optics(0.1, 0.01, 0.0);
                      const Vec3d MLI_optics(0.3, 0.3, 0.0);
                      const Vec3d Alluminium_optics(0.8, 1.1, 0.0);
                      // Map containing optical features of different materials
                      const std::map<const std::string, const Vec3d> optical_properties =
                          {
                          
                          {"Solar Array", SolArrays_optics},
                          {"MLI", MLI_optics},
                          {"Alluminium", Alluminium_optics}
                          
                          };
                        
                      const Vec3d SC_aerodynamics(0.8, 0.8, 0.05); 
                      }
                    
            namespace spacecraft
                    {
                    const Vec3d n_Xplus(1.0, 0.0, 0.0);    
                    const Vec3d n_Xminus = -n_Xplus;
                    const Vec3d n_Yplus(0.0, 1.0, 0.0);    
                    const Vec3d n_Yminus = -n_Yplus;
                    const Vec3d n_Zplus(0.0, 0.0, 1.0);    
                    const Vec3d n_Zminus = -n_Zplus;
                    // Map containing normal unit vectors to the SC faces (SC body-fixed frame)
                    const std::map<const std::string, const Vec3d> Normals =
                                {
                                
                                {"+X", n_Xplus},
                                {"-X", n_Xminus},
                                {"+Y", n_Yplus},
                                {"-Y", n_Yminus},
                                {"+Z", n_Zplus},
                                {"-Z", n_Zminus}
                                
                                };
                    }
                    
            } // namespace constants

#endif
