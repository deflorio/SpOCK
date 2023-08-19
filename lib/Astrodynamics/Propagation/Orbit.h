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

#ifndef ORBIT_H
#define ORBIT_H

#include <Propagator.h>
#include <VarTypes.h>
#include <Gravity.h>
#include <ThirdBody.h>
#include <SolarRadiation.h>
#include <Atmosphere.h>
// External libraries: Eigen
#include <Eigen/Core>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

//#include <boost/numeric/odeint/stepper/XYZ.hpp>   // Include path for all steppers, XYZ is a placeholder for a stepper
//#include <boost/numeric/odeint/algebra/XYZ.hpp>   // All algebras
//#include <boost/numeric/odeint/util/XYZ.hpp>      // Utility functions like is_resizeable , same_size , or resize
//#include <boost/numeric/odeint/integrate/XYZ.hpp> // Integrate routines.
//#include <boost/numeric/odeint/iterator/XYZ.hpp>  // Range and iterator functions.
//#include <boost/numeric/odeint/external/XYZ.hpp>  // Any binders to external libraries

using namespace boost::numeric::odeint;
using namespace SC;
using namespace math;
using namespace propagator;
using namespace Eigen;
using namespace gravity;
using namespace thirdbody;
using namespace solradiation;
using namespace atmosphere;

typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

namespace orbit
   {
    //------------------------------------------------------------------------------
    //! Class ORB
    //------------------------------------------------------------------------------
    /*!
       Class derived from PROP which implements the orbit dynamics of a spacecraft
     */
    //------------------------------------------------------------------------------ 
    class ORB : public PROP
        {
        public:
        //! Constructor.
        /*!
            Using class PROP constructor
          */
        using PROP::PROP;
        //! Destructor.
        ~ORB();
        //------------------------------------------------------------------------------
        // Class methods specification
        //------------------------------------------------------------------------------
        // Set up file paths and parameters for environmental forces models
        void ForceModelsSetup();
        // Insert total actuators' torque
        void Maneuver(const Ref<const VectorXd>& maneuver);
        // Compute perturbation torques
        void ComputeAction(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& orb_state);
        // Integrate equations of motion
        void Integrate(double t, double step);
        // Setup numerical integrator parameters
        //void StepperSetup(double eps_abs, double eps_rel, double factor_x, double factor_dxdt);
        
        public:
    
        /** Environmental accelerations.*/
        static VectorNd<15> Acceleration_env;
        /** Vecor containing ECI orbital state and ECEF position.*/
        static VectorNd<9> orb_state_ECI_ECEF;
        /** Spacecraft orbit state in Euler angles form (phi,theta,psi,om_x,om_y,om_z) [rad,rad,rad,rad/s,rad/s,rad/s].*/
        static state_type x;
        
        protected:
        
        GRAV GravityField;
        THIRDB SunMoonPerturbation;
        SOLRAD SolarRadiation;
        ATMO Atmosphere;
    
        /** Acceleration vector.*/
        static Vec3d Acceleration;
        /** Thrusters acceleration vector.*/
        static Vec3d dv_CMD;
        
        private:
        // Implementation of dynamics model
        void DynModel(const state_type &x , state_type &dxdt , const double t);
    
        private:
            
        /** Rotation matrix from ECI to spacecraft body-fixed frame.*/
        Mat3x3d ECItoBody; // Transformation matrix ECI to body-fixed
        
        static Vector6d rv_vec;
        static Vec4d q_currentstate;
        static Vec3d pos_ECI;
        static Vec3d pos_ECEF;
        /** boost odeint stepper for integrator. @see Method DynModel.*/
        runge_kutta4< state_type > stepper;
        //runge_kutta_dopri5< state_type > stepper;
        //adams_bashforth_moulton< 4 , state_type > stepper;
        //controlled_stepper_type controlled_stepper;
        //bulirsch_stoer<state_type> bulirsch_stoer_stepper;
        };

   }; // End of namespace pod_trj

#endif
