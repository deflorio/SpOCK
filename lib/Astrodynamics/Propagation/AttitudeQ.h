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

#ifndef ATTITUDEQ_H
#define ATTITUDEQ_H

#include <Propagator.h>
#include <Attitude.h>
#include <VarTypes.h>
// External libraries: Eigen
#include <Eigen/Core>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

//#include <boost/numeric/odeint/stepper/XYZ.hpp>   // Include path for all steppers, XYZ is a placeholder for a stepper
//#include <boost/numeric/odeint/algebra/XYZ.hpp>   // All algebras
//#include <boost/numeric/odeint/util/XYZ.hpp>      // Utility functions like is_resizeable , same_size , or resize
//#include <boost/numeric/odeint/integrate/XYZ.hpp> // Integrate routines.
//#include <boost/numeric/odeint/iterator/XYZ.hpp>  // Range and iterator functions.
//#include <boost/numeric/odeint/external/XYZ.hpp>  // Any binders to external libraries

using namespace boost::numeric::odeint;
//typedef boost::array< double , 7 > state_typeQ;

typedef runge_kutta_cash_karp54< state_typeQ > error_stepper_typeQ;
typedef controlled_runge_kutta< error_stepper_typeQ > controlled_stepper_typeQ;

using namespace SC;
using namespace math;
using namespace propagator;
using namespace attitude;
using namespace Eigen;

namespace attitudeq
   {
    //------------------------------------------------------------------------------
    //! Class ATTQ
    //------------------------------------------------------------------------------
    /*!
       Class derived from ATT which implements the attitude dynamics of a spacecraft
       using quaternions
     */
    //------------------------------------------------------------------------------ 
    class ATTQ : public ATT
        {
        public:
        //! Constructor.
        /*!
            Using class PROP constructor
         */
        //using PROP::PROP;
        using ATT::ATT;
        /** Using variable Torque of class ATT.*/
        using ATT::Torque;
        //! Destructor.
        ~ATTQ();
        //------------------------------------------------------------------------------
        // Class methods specification
        //------------------------------------------------------------------------------
        // Translate in quaternion form the initial state given as Euler angles
        //void Eul2Quat_inistate(Vector6d& init_att_state);
        // Integrate equations of motion
        void Integrate(double t, double step);  // Execute fixed step numerical integration
        // Setup numerical integrator parameters
        void StepperSetup(double eps_abs, double eps_rel, double factor_x, double factor_dxdt);
        
        public:
    
        /** Spacecraft attitude state (q1,q2,q3,q4,om_x,om_y,om_z).*/
        //VectorNd<7> stateQ;
        /** Spacecraft attitude state (q1,q2,q3,q4,om_x,om_y,om_z) in format for odeint integrator.*/
        static state_typeQ x;
    
        private:
        // Implementation od dynamics model with quaternions
        void DynModel(const state_typeQ &x , state_typeQ &dxdt , const double t); // Definition of dynamic model
    
        private:
            
        /** Rotation matrix from ECI to spacecraft body-fixed frame.*/
        //Mat3x3d ECItoBody;
        /** Spacecraft attitude initial state (q1,q2,q3,q4,om_x,om_y,om_z).*/
        //VectorNd<7> initstateQ;
        /** boost odeint stepper for integrator. @see Method DynModel.*/
        //runge_kutta4< state_typeQ > stepper;
        //controlled_stepper_typeQ controlled_stepper;
        bulirsch_stoer<state_typeQ> bulirsch_stoer_stepperQ;
        
        //double abs_err = 1.0E-10;
        //double rel_err = 1.0E-6;
        //double a_x = 1.0;
        //double a_dxdt = 1.0;
        //controlled_stepper_typeQ controlled_stepper(default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
        //controlled_stepper_typeQ controlled_stepper;
        };

   }; // End of namespace attitudeq

#endif
