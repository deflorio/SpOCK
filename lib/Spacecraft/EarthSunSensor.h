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

#ifndef EARTHSUNSENSOR_H
#define EARTHSUNSENSOR_H

#include <random>

#include <Subsystem.h>
#include <Solarsys.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>
//#include <boost/math/distributions/normal.hpp>

using namespace subsystem;
using namespace math;
using namespace std;
using namespace solarsystem;
//using boost::math::normal;

namespace earthsun
   {
    //------------------------------------------------------------------------------
    //! Class SOLRAD
    //------------------------------------------------------------------------------
    /*!
       Class derived from SUBSYS which implements coarse sun sensor, solar camera and
       nadir camera models
     */
    //------------------------------------------------------------------------------ 
    class EARTHSUNSENS : public SUBSYS
        {
        public:
        //! Constructor.
        /*!
            Using class SUBSYS constructor
          */
        using SUBSYS::SUBSYS;
        //! Destructor.
        ~EARTHSUNSENS();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Evaluate the operability status (subsystem on or off)
        VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate);
        // Store in clas members constant parameters contained in base class' member vector ConstPrm
        void Init();
        
        public:
        /** Unit vector from the spacecraft to the Sun (ECI).*/
        Vec3d s_uvec;
        /** Unit vector from the spacecraft to the Sun (SC body-fixed frame).*/
        Vec3d s_uvec_SC;
        /** Unit vector from the spacecraft to the Earth center (ECI).*/
        Vec3d e_uvec;
        /** Unit vector from the spacecraft to the Sun (SC body-fixed frame).*/
        Vec3d e_uvec_SC;
        /** Sensor capture status.*/
        double capture;
        /** Sensor detection result.*/
        double detection;
        /** Angle between Sun direction and main sensor axis.*/
        double theta;
        /** cos(theta)*/
        double costheta;
        /** sin(theta)*/
        double sintheta;
        /** Sensor field of view semi-angle.*/
        double fov;
        
        
        private:
        // Evaluate the operability status (subsystem on or off)
        void status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) ;
        
        private:
        /** Object of class SOLSYS type.*/
        SOLSYS Solar;
        //SOLSYS* Solar_ptr = &Solar;
        /** Coarse Sun sensor nominal current at zero Sun incidence angle.*/
        double I0;
        /** Coarse Sun sensor signal noise.*/
        double nu; // Signal noise
        /** Normal distribution error of measurements.*/
        //normal_distribution<double> MeasError;//(0.0,1.0); //boost::math::normal::normal SunMeasError;
        ///** Normal distribution error of Earth camera measurements.*/
        //normal_distribution<double> EarthMeasError(0.0,1.0); //boost::math::normal::normal EarthMeasError;
        ///** Normal distribution error of coarse Sun sensor measurements.*/
        //normal_distribution<double> CSS(0.0,1.0); //boost::math::normal::normal CSS;
        
        //random_device rd;
        //mt19937 generator(random_device());
        };

   }; // End of namespace magnetic

#endif
