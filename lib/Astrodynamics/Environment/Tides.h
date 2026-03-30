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

#ifndef TIDES_H
#define TIDES_H

#include <SpaceEnvironment.h>
#include <Solarsys.h>
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
using namespace solarsystem;

namespace tides
   {
    //------------------------------------------------------------------------------
    //! Class TIDES
    //------------------------------------------------------------------------------
    /*!
       Class derived from SPACEENV which implements the solar radiation pressure model
     */
    //------------------------------------------------------------------------------ 
    class TIDES : public SPACEENV
        {
        public:
        //! Constructor.
        /*!
            Using class SPACEENV constructor
          */
        using SPACEENV::SPACEENV;
        //! Destructor.
        ~TIDES();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Compute acceleration from simplified model of solid Earth tides
        Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate);
        // Contribution of tides to gravity field model coefficients
        void getmodel_coeff();
        
        private:
        /** Object of class SOLSYS type.*/
        SOLSYS Solar;
        };

   }; // End of namespace tides

#endif
