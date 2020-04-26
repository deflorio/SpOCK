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

#ifndef SOLARPAN_H
#define SOLARPAN_H

#include <random>

#include <Subsystem.h>
#include <Solarsys.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>

using namespace subsystem;
using namespace math;
using namespace std;
using namespace solarsystem;
//using boost::math::normal;

namespace solarpan
   {
    //------------------------------------------------------------------------------
    //! Class SOLRAD
    //------------------------------------------------------------------------------
    /*!
       Class derived from SUBSYS which implements coarse sun sensor, solar camera and
       nadir camera models
     */
    //------------------------------------------------------------------------------ 
    class SOLARPAN : public SUBSYS
        {
        public:
        //! Constructor.
        /*!
            Using class SUBSYS constructor
          */
        using SUBSYS::SUBSYS;
        //! Destructor.
        ~SOLARPAN();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Generated power and current
        VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate);
        // Store in clas members constant parameters contained in base class' member vector ConstPrm
        void Init();
        
        public:
        /** Unit vector from the spacecraft to the Sun (ECI).*/
        Vec3d s_uvec;
        /** Unit vector from the spacecraft to the Sun (SC body-fixed frame).*/
        Vec3d s_uvec_SC;
        
        private:
        // Evaluate the operability status (subsystem on or off)
        void status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) ;
        
        private:
        /** Object of class SOLSYS type.*/
        SOLSYS Solar;
        /** True when the spacecraft is in eclipse.*/
        static bool eclipse;
        /** Solar panelsurface [m^2].*/
        static double Surface;
        /** Voltage of a string of solar cells.*/
        static double V_string;
        /** Solar cells efficiency.*/
        static double epsilon; // Signal noise
        /** Scalar product of sun unit vector and solar panel unit normal vector.*/
        static double costheta;
        };

   }; // End of namespace magnetic

#endif
