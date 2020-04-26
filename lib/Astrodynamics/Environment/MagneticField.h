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

#ifndef MAGFIELD_H
#define MAGFIELD_H

#include <SpaceEnvironment.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>
//#include <GeographicLib/include/MagneticModel.hpp>

//#include <emm_sph_point.h>


extern "C"
      {
      #include "extlib/cspice/include/SpiceUsr.h"
      //#include "extlib/MagneticField/IGRF/C/geomag70.h"
      }

using namespace spaceenvironment;
using namespace math;
using namespace std;
//using namespace GeographicLib;

namespace magnetic
   {
    //------------------------------------------------------------------------------
    //! Class MAGFIELD
    //------------------------------------------------------------------------------
    /*!
       Class derived from SPACEENV which implements the Earth's magnetic field model
     */
    //------------------------------------------------------------------------------ 
    class MAGFIELD : public SPACEENV
          {
      
          public:
          //! Constructor.
          /*!
              Using class SPACEENV constructor
            */
          using SPACEENV::SPACEENV;
          //! Destructor.
          ~MAGFIELD();
          //------------------------------------------------------------------------------
          //
          // Class methods specification
          //
          //------------------------------------------------------------------------------
          // Earth's magnetic field vector.*/
          Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate);
          
          void getmodel_coeff();
          
          public:
          
          /** Magnetic field vector.*/  
          static Vec3d magfield;
          /** Initial epoch to read magnetic field model file. */
          double init_epoch;
          
          private:
          /** WMM magnet model indices file.*/
          static char c_wmm_file[200];
          /** Initial year to read magnetic field model file. */
          static float init_year;
          };

   }; // End of namespace magnetic

#endif
