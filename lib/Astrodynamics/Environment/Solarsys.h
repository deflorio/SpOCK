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

#ifndef SOLARSYS_H
#define SOLARSYS_H

#include <SpaceEnvironment.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>
extern "C"
      {
       #include "extlib/cspice/include/SpiceUsr.h"
      }

using namespace spaceenvironment;
using namespace math;
using namespace std;
    
namespace solarsystem
   {
    //------------------------------------------------------------------------------
    //! Class SOLSYS
    //------------------------------------------------------------------------------
    /*!
       Class which implements functions related to the relative positioning of spacecraft
       and solar system objects
     */
    //------------------------------------------------------------------------------ 
    class SOLSYS
        {
        public:
        //! Constructor.
        SOLSYS();
        //! Destructor.
        ~SOLSYS();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------  
        // Load SPICE planet ephemeris file. This method is to be used in case the ephemeris are not loaded in the main file.*/
        void load_ephemeris_file(const string ephe_filepath);
        // Return the position of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.*/
        Vec3d OBJ_pos(const char* OBJtarget_name, const char* OBJobserver_name, const char* refframe, const char* aberration_corr, const double time);
        // Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.*/
        Vector6d OBJ_posvel(const char* OBJtarget_name, const char* OBJobserver_name, const char* refframe, const char* aberration_corr, const double time);
        // Return the position of the Sun in rectangular coordinates (ECI) from JPL ephemerides.*/
        Vec3d sunposREC(double GPStime);
        // Return the position of the Sun in rectangular coordinates (ECI) from low-precision (LP) Solar coordinates.*/
        Vec3d LP_sunposREC(double GPStime);
        // Return the position of the Moon in rectangular coordinates (ECI) from JPL ephemerides.*/
        Vec3d moonposREC(double GPStime);
        // Return the position of the Moon in rectangular coordinates (ECI) from low-precision (LP) Lunar coordinates.*/
        Vec3d LP_moonposREC(double GPStime);
        // Return right ascension, declination and range of the sun (ECI).*/
        Vec3d sunposRAD(double GPStime);
        // Return the angle between the spacecraft position vector(ECI) and the Sun position vector (ECI).*/
        double xi_angle(double GPStime, const Vec3d& SC_pos);
        // Return the eclipse condition of the spacecraft at a certain time.*/
        void eclipse(double GPStime, const Vec3d& SC_pos, bool& umbra, bool& penumbra);
        // Return the unit vector (ECI) of sun direction from the spacecraft's center of mass.*/
        Vec3d sundir_u(double GPStime, const Vec3d& SC_pos);
        // Return the dot product angle between the sun-spacecraft direction unit vector and a generic unit vector.*/
        double sun_angle(double GPStime, Vec3d& SC_pos, Vec3d& v3D_u);
        //double xi_angle;  // Angle between the spacecraft position vector (Earth centered) and the position vector of the Sun (Earth centered)
        };

   }; // End of namespace solarsystem

#endif
