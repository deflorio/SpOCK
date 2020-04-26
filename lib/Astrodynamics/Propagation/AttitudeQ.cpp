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
#include <string>
#include <math.h>

#include <AttitudeQ.h>
#include <Constants.h>
#include <Transformations.h>

using namespace std;
using namespace math;
using namespace boost::numeric::odeint;
using namespace constants;
using namespace propagator;

namespace pl = std::placeholders;

namespace attitudeq
    {
    //------------------------------------------------------------------------------
    // Class ATTQ implementation
    //------------------------------------------------------------------------------
    state_typeQ ATTQ::x = {};
    // Destructor
    //------------------------------------------------------------------------------
    ATTQ::~ATTQ() {};
    //------------------------------------------------------------------------------
    // Method void Eul2Quat_inistate(Vector6d& init_att_state)
    //------------------------------------------------------------------------------
    /**
     * Translate in quaternion form the initial state given as Euler angles
     *
     * @param init_att_state      Input initial attitude state in Euler angles form (phi,theta,psi,om_x,om_y,om_z) [rad,rad,rad,rad/s,rad/s,rad/s]
     */
    //------------------------------------------------------------------------------   
    //void ATTQ::Eul2Quat_inistate(Vector6d& init_att_state)
    //              {
    //              initstateQ = VectorNd<7>::Zero();
    //              double phi, theta, psi;
    //              phi = init_att_state(0);
    //              theta = init_att_state(1);
    //              psi = init_att_state(2);
    //              
    //              ECItoBody = RotationMatrix321(phi, theta, psi);
    //              
    //              initstateQ.segment(0,4) = RotationMatrix2Quaternion(ECItoBody);
    //              
    //              initstateQ.segment(4,3) = init_att_state.segment(3,3);
    //              
    //              //cout << "InitstateQ: " << initstateQ(0) << "   " << initstateQ(1) << "   " << initstateQ(2) << "   " << initstateQ(3) << "   " << initstateQ(4)*180/3.14 << "   " << initstateQ(5)*180/3.14 << "   " << initstateQ(6)*180/3.14 << endl;
    //              };
    //------------------------------------------------------------------------------
    // Method static void DynModel(const state_type &x , state_type &dxdt , const double t)
    //------------------------------------------------------------------------------
    /**
     * Implementation of dynamics model with quaternions.
     * @note The dynamics of reaction/momentum wheels is also included.
     *
     * @param x      State vector
     * @param dxdt   First derivative of state vector
     * @param t      Step [s]
     */
    //------------------------------------------------------------------------------   
    void ATTQ::DynModel(
            const state_typeQ &x,
            state_typeQ &dxdt,
            const double t)
            {
            double abs_t = inittime + t;
            Vec4d state_quat;
            state_quat << x[0], x[1], x[2], x[3];
            ComputeAction(abs_t, state_quat, orbstate);
                
            dxdt[0] = 0.5*( x[6]*x[1] - x[5]*x[2] + x[4]*x[3] );
            
            dxdt[1] = 0.5*( -x[6]*x[0] + x[4]*x[2] + x[5]*x[3] );
            
            dxdt[2] = 0.5*( x[5]*x[0] - x[4]*x[1] + x[6]*x[3] );
            
            dxdt[3] = 0.5*( -x[4]*x[0] - x[5]*x[1] - x[6]*x[2] );
                    
            dxdt[4] = -( invMoI(0,0)*( ( MoI(2,0)*x[5] - MoI(1,0)*x[6] )*x[4] + ( MoI(2,1)*x[5] - MoI(1,1)*x[6] )*x[5] + ( MoI(2,2)*x[5] - MoI(1,2)*x[6] )*x[6] ) +
                         invMoI(0,1)*( ( MoI(0,0)*x[6] - MoI(2,0)*x[4] )*x[4] + ( MoI(0,1)*x[6] - MoI(2,1)*x[4] )*x[5] + ( MoI(0,2)*x[6] - MoI(2,2)*x[4] )*x[6] ) +
                         invMoI(0,2)*( ( MoI(1,0)*x[4] - MoI(0,0)*x[5] )*x[4] + ( MoI(1,1)*x[4] - MoI(0,1)*x[5] )*x[5] + ( MoI(1,2)*x[4] - MoI(0,2)*x[5] )*x[6] ) ) +
            
                      -( invMoI(0,0)*( -hw(1)*x[6] + hw(2)*x[5] ) + invMoI(0,1)*( hw(0)*x[6] - hw(2)*x[4] ) + invMoI(0,2)*( -hw(0)*x[5] + hw(1)*x[4] ) ) +
            
                         invMoI(0,0)*Torque(0) + invMoI(0,1)*Torque(1) + invMoI(0,2)*Torque(2);
                        
            dxdt[5] = -( invMoI(1,0)*( ( MoI(2,0)*x[5] - MoI(1,0)*x[6] )*x[4] + ( MoI(2,1)*x[5] - MoI(1,1)*x[6] )*x[5] + ( MoI(2,2)*x[5] - MoI(1,2)*x[6] )*x[6] ) +
                         invMoI(1,1)*( ( MoI(0,0)*x[6] - MoI(2,0)*x[4] )*x[4] + ( MoI(0,1)*x[6] - MoI(2,1)*x[4] )*x[5] + ( MoI(0,2)*x[6] - MoI(2,2)*x[4] )*x[6] ) +
                         invMoI(1,2)*( ( MoI(1,0)*x[4] - MoI(0,0)*x[5] )*x[4] + ( MoI(1,1)*x[4] - MoI(0,1)*x[5] )*x[5] + ( MoI(1,2)*x[4] - MoI(0,2)*x[5] )*x[6] ) ) +
            
                      -( invMoI(1,0)*( -hw(1)*x[6] + hw(2)*x[5] ) + invMoI(1,1)*( hw(0)*x[6] - hw(2)*x[4] ) + invMoI(1,2)*( -hw(0)*x[5] + hw(1)*x[4] ) ) +    
            
                         invMoI(1,0)*Torque(0) + invMoI(1,1)*Torque(1) + invMoI(1,2)*Torque(2);
                        
            dxdt[6] = -( invMoI(2,0)*( ( MoI(2,0)*x[5] - MoI(1,0)*x[6] )*x[4] + ( MoI(2,1)*x[5] - MoI(1,1)*x[6] )*x[5] + ( MoI(2,2)*x[5] - MoI(1,2)*x[6] )*x[6] ) +
                         invMoI(2,1)*( ( MoI(0,0)*x[6] - MoI(2,0)*x[4] )*x[4] + ( MoI(0,1)*x[6] - MoI(2,1)*x[4] )*x[5] + ( MoI(0,2)*x[6] - MoI(2,2)*x[4] )*x[6] ) +
                         invMoI(2,2)*( ( MoI(1,0)*x[4] - MoI(0,0)*x[5] )*x[4] + ( MoI(1,1)*x[4] - MoI(0,1)*x[5] )*x[5] + ( MoI(1,2)*x[4] - MoI(0,2)*x[5] )*x[6] ) ) +
            
                      -( invMoI(2,0)*( -hw(1)*x[6] + hw(2)*x[5] ) + invMoI(2,1)*( hw(0)*x[6] - hw(2)*x[4] ) + invMoI(2,2)*( -hw(0)*x[5] + hw(1)*x[4] ) ) +
            
                         invMoI(2,0)*Torque(0) + invMoI(2,1)*Torque(1) + invMoI(2,2)*Torque(2);
            };
    //------------------------------------------------------------------------------
    // Method void StepperSetup(double eps_abs, double eps_rel, double factor_x, double factor_dxdt)
    //------------------------------------------------------------------------------
    /**
     * Setup numerical integrator parameters
     *
     * @param eps_abs       Absolute tolerance level
     * @param eps_rel       Relative tolerance level
     * @param factor_x      Factor for the weight of the derivative
     * @param factor_dxdt   Factor for the weight of the state
     */
    //------------------------------------------------------------------------------   
    void ATTQ::StepperSetup(double eps_abs,
                            double eps_rel,
                            double factor_x,
                            double factor_dxdt)
                            {
                            bulirsch_stoer<state_typeQ> setup_stepper(eps_abs, eps_rel, factor_x, factor_dxdt);
                            bulirsch_stoer_stepperQ = setup_stepper;
                            };
    //------------------------------------------------------------------------------
    // Method void Integrate(double t, double step)
    //------------------------------------------------------------------------------
    /**
     * Integrate equations of motion in quaternion form
     *
     * @param t     Step [s]
     * @param step  Step length [s]
     */
    //------------------------------------------------------------------------------
    void ATTQ::Integrate(
            double t,
            double step)
            {
            if(integ_first_step)
               {
               for(int i = 0; i < 7; i++) x[i] = initstate(i);
               //cout << "Initial state: " << initstate(0) << "   " << initstate(1) << "   " << initstate(2) << "   " << initstate(3) << "   " << initstate(4)*180/3.14 << "   " << initstate(5)*180/3.14 << "   " << initstate(6)*180/3.14 << endl;
               state = initstate;
               integ_first_step = false;
               cout << "First integration step" << endl;
               }
            //q_attstate = Eigen::Quaterniond(state(0),state(1),state(2),state(3));
            Vec4d state_quat = state.segment(0,4);
            ECItoBody = Quaternion2RotationMatrix(state_quat);
            //stepper.do_step(DynModel, x, t, step);
            //stepper.do_step(std::bind(&ATTQ::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            
            bulirsch_stoer_stepperQ.try_step(std::bind(&ATTQ::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            //controlled_stepper.do_step(DynModel, x, t, step);
            //integrate_adaptive( make_controlled< error_stepper_typeQ >( 1.0e-10 , 1.0e-6 ) , DynModel, x, t, 2.0*86400 , 0.01 );
            
            //cout << t << "   " << x[0]*180/3.14 << "   " << x[1]*180/3.14 << "   " << x[2]*180/3.14 << "   " << x[3]*180/3.14 << "   " << x[4]*180/3.14 << "   " << x[5]*180/3.14 << endl;
            //cout << "Torque: " << t << "   " << Torque(0) << "   " << Torque(1) << "   " << Torque(2) << "\n" << endl;
            
            for(int i = 0; i < 7; i++) state(i) = x[i];
            
            //state_quat = state.segment(0,4);
            //state_quat = state_quat.normalized();
            //for(int i = 0; i < 4; i++) x[i] = state_quat(i);
            
            //cout << "state_quat.norm(): " << state_quat.norm() << endl;
            //cout << t << "   " << state(0)*180/3.14 << "   " << state(1)*180/3.14 << "   " << state(2)*180/3.14 << "   " << state(3)*180/3.14 << "   " << state(4)*180/3.14 << "   " << state(5)*180/3.14 << endl;
            
            //cout << "Torque in Integrate: " << Torque(0) << "   " << Torque(1) << "   " << Torque(2) << "\n" << endl;
            };
            
}; // End of namespace attitudeq



