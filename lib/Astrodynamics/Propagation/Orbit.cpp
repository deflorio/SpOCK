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

#include <Orbit.h>
#include <Constants.h>
#include <Transformations.h>

using namespace std;
using namespace math;
using namespace boost::numeric::odeint;
using namespace constants;
using namespace propagator;

namespace pl = std::placeholders;

namespace orbit
    {
    //------------------------------------------------------------------------------
    // Class ORB implementation
    //------------------------------------------------------------------------------
    VectorNd<15> ORB::Acceleration_env = VectorNd<15>::Zero();
    VectorNd<9> ORB::orb_state_ECI_ECEF = VectorNd<9>::Zero();
    Vector6d ORB::rv_vec = Vector6d::Zero();
    Vec4d ORB::q_currentstate = Vec4d::Zero();
    Vec3d ORB::pos_ECI = Vec3d::Zero();
    Vec3d ORB::pos_ECEF = Vec3d::Zero();
    Vec3d ORB::Acceleration = Vec3d::Zero();
    Vec3d ORB::dv_CMD = Vec3d::Zero();
    state_type ORB::x = {};
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    ORB::~ORB() {};
    //------------------------------------------------------------------------------
    // Method void ForceModelsSetup()
    //------------------------------------------------------------------------------
    /**
     * Set up file paths and parameters for environmental forces models
     */
    //------------------------------------------------------------------------------ 
    void ORB::ForceModelsSetup()
            {
            GravityField.setfilespath(gravityfield,datapath);
            GravityField.grav_epoch = inittime;
			//(&GravityField)->n_max = nMAX;
            //SPACEENV::n_max = nMAX;
            GravityField.getmodel_coeff();
            
            if(drag_on)
                {
                Atmosphere.init_epoch = inittime;
                Atmosphere.simduration = simdur;
                Atmosphere.Drag_Model = Drag_Model;
                Atmosphere.CF = CD;
                Atmosphere.Area = Area_D;
                Atmosphere.setfilespath(atmosphere,datapath);
                Atmosphere.getmodel_coeff();
                }
                
            if(srp_on)
                {
                SolarRadiation.SRP_Model = SRP_Model;
                SolarRadiation.CF = C_SRP;
                SolarRadiation.Area = Area_R;
                }
            };
    //------------------------------------------------------------------------------
    // Method static void DynModel(const state_type &x , state_type &dxdt , const double t)
    //------------------------------------------------------------------------------
    /**
     * Implementation of dynamics model with Euler angles.
     * @note The dynamics of reaction/momentum wheels is also included.
     *
     * @param x      State vector
     * @param dxdt   First derivative of state vector
     * @param t      Step [s]
     */
    //------------------------------------------------------------------------------         
    void ORB::DynModel(
                    const state_type &x,
                    state_type &dxdt,
                    const double t)
                    {
                    double abs_t = inittime + t;
                    
                    rv_vec << x[0], x[1], x[2], x[3], x[4], x[5];
                    ComputeAction(abs_t, state, rv_vec);
                    
                    dxdt[0] = x[3];
                    
                    dxdt[1] = x[4];
                    
                    dxdt[2] = x[5];
                            
                    dxdt[3] = Acceleration(0);
                                
                    dxdt[4] = Acceleration(1);
                                
                    dxdt[5] = Acceleration(2);
                    
                    //cout << "Acceleration(0): " << Acceleration(0) << "  Acceleration(1): " << Acceleration(1) << "  Acceleration(2): " << Acceleration(2) << endl;
                    //cout << "acc_ECEF: " << Acceleration.norm() << endl;
                    //Vec3d r_vec;
                    //r_vec << x[0], x[1], x[2];
                    //cout << fixed << "phi_acc: " << Acceleration.dot(r_vec)/(Acceleration.norm()*r_vec.norm()) << endl;
                    };
    //------------------------------------------------------------------------------
    // Abstract method void Maneuver(const Ref<const VectorXd>& maneuver)
    //------------------------------------------------------------------------------
    /**
      * Insert the total torque generated by the actuators
      *
      * @param maneuver   Total actuators' torque
      */
    //------------------------------------------------------------------------------  
    void ORB::Maneuver(const Ref<const VectorXd>& maneuver)
                {
                dv_CMD = maneuver;
                };
    //------------------------------------------------------------------------------
    // Method void ComputeAction(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& orb_state)
    //------------------------------------------------------------------------------
    /**
     * Compute perturbation torques
     *
     * @param epoch         GPS epoch (seconds) of the input states
     * @param currentstate  Orbit state at epoch GPStime
     * @param orb_state     12-D vector containing ECI orbital state (0-5) and ECEF orbital state (6-11) at epoch GPStime
     */
    //------------------------------------------------------------------------------       
    void ORB::ComputeAction(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& orb_state)
            {
            Acceleration = Vec3d::Zero();
            Acceleration_env = VectorNd<15>::Zero();
            orb_state_ECI_ECEF = VectorNd<9>::Zero();
            pos_ECI = Vec3d::Zero();
            pos_ECEF = Vec3d::Zero();
            //vel_ECI = Vec3d::Zero();
            //Vec3d pos_ECI, vel_ECI, pos_ECEF;
            //VectorNd<9> orb_state_ECI_ECEF;
            
            pos_ECI = orb_state.segment(0,3);
            pos_ECEF = v3D_transform(epoch, pos_ECI, "J2000", "ITRF93");
            //vel_ECI = orb_state.segment(3,3);
            
            orb_state_ECI_ECEF.segment(0,6) = orb_state;
            orb_state_ECI_ECEF.segment(6,3) = pos_ECEF;
            
            q_currentstate = currentstate;
            
            // Computation of gravity field acceleration (default acceleration)
            Acceleration = GravityField.field_vec(epoch, pos_ECEF);
            Acceleration_env.segment(0,3) = Acceleration;
            
            // Computation of third body acceleration (Sun and Moon)
            if(sunmoon_on)
                {
                Vec3d acc_sunmoon;
                acc_sunmoon = SunMoonPerturbation.field_vec(epoch, pos_ECI);
                
                Acceleration += acc_sunmoon;
                Acceleration_env.segment(3,3) = acc_sunmoon;
                }
            
            // Computation of solar radiation pressure acceleration   
            if(srp_on)
                {
                Vec3d F_SRP = Vec3d::Zero();
                
                if( SRP_Model.compare("Panels") == 0 )
                    {
                    Vec3d F_SRP_ECI = Vec3d::Zero();
                    SC::Face F;
                    
                    for(auto element : SC_Faces) // Ranged based for loop on map SC_Faces (member of Class PROP)
                       {
                       F = SC_Faces[element.first];
                       
                       SolarRadiation.SetSurfaceParameters(F, q_currentstate); // C_SRP is a member of Class PROP
                       F_SRP_ECI = SolarRadiation.field_vec(epoch, pos_ECI); // Solar radiation pressure force (ECI) on surface F
                       
                       F_SRP += F_SRP_ECI;
                       }
                   }
                else if( SRP_Model.compare("RefArea") == 0 )
                    {
                    F_SRP = SolarRadiation.field_vec(epoch, pos_ECI); // Solar radiation pressure force (ECI)
                    }
                
                Acceleration_env.segment(6,3) = F_SRP/SC_mass;;
                   
                Acceleration += F_SRP/SC_mass;
                }
                
            // Computation of atmospheric drag acceleration   
            if(drag_on)
                {
                Vec3d F_ATM = Vec3d::Zero();
                
                Atmosphere.AtmosphericDensity(epoch, orb_state_ECI_ECEF);
                
                //cout << setprecision(20) <<Atmosphere.rho_atm << endl;
                
                if( Drag_Model.compare("Panels") == 0 )
                    {
                    Vec3d F_ATM_ECI = Vec3d::Zero();
                    SC::Face F;
                    
                    for(auto element : SC_Faces) // Ranged based for loop on map SC_Faces (member of Class PROP)
                       {
                       F = SC_Faces[element.first];
                       
                       Atmosphere.SetSurfaceParameters(F, q_currentstate); // C_ATM is a member of Class PROP
                       F_ATM_ECI = Atmosphere.field_vec(epoch, orb_state_ECI_ECEF); // Atmospheric drag force (ECI) on surface F
                       
                       F_ATM += F_ATM_ECI;
                       }
                    }
                else if( Drag_Model.compare("RefArea") == 0 )
                    {
                    F_ATM = Atmosphere.field_vec(epoch, orb_state_ECI_ECEF); // Atmospheric drag force (ECI)
                    }
                
                Acceleration_env.segment(9,3) = F_ATM/SC_mass;
                
                Acceleration += F_ATM/SC_mass;
                }    
            
            //cout << "Acceleration(0): " << Acceleration(0) << "  Acceleration(1): " << Acceleration(1) << "  Acceleration(2): " << Acceleration(2) << endl;
            //cout << "Acceleration: " << sqrt(Acceleration(0)*Acceleration(0) +  Acceleration(1)*Acceleration(1) + Acceleration(2)*Acceleration(2)) << endl;
            
            // Put environmental torque values in Acceleration_env
            Acceleration_env.segment(12,3) = Acceleration;
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
    void ORB::StepperSetup(double eps_abs,
                            double eps_rel,
                            double factor_x,
                            double factor_dxdt)
                            {
                            bulirsch_stoer<state_type> setup_stepper(eps_abs, eps_rel, factor_x, factor_dxdt);
                            bulirsch_stoer_stepper = setup_stepper;
                            };
    //------------------------------------------------------------------------------
    // Method void Integrate(double t, double step)
    //------------------------------------------------------------------------------
    /**
     * Integrate equations of motion in Euler form
     *
     * @param t     Step [s]
     * @param step  Step length [s]
     */
    //------------------------------------------------------------------------------
    void ORB::Integrate(
            double t,
            double step)
            {
            if(integ_first_step)
               {
               for(int i = 0; i < 6; i++) x[i] = orbstate(i);
               
               integ_first_step = false;
               cout << "First integration step\n" << endl;
               }
            
            x[3] = x[3] + dv_CMD(0);
                    
            x[4] = x[4] + dv_CMD(1);
                    
            x[5] = x[5] + dv_CMD(2);
            
            stepper.do_step(std::bind(&ORB::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            //bulirsch_stoer_stepper.try_step(std::bind(&ORB::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            
            for(int i = 0; i < 6; i++) orbstate(i) = x[i];
            
            //cout << t << "   " << state(0)*180/3.14 << "   " << state(1)*180/3.14 << "   " << state(2)*180/3.14 << "   " << state(3)*180/3.14 << "   " << state(4)*180/3.14 << "   " << state(5)*180/3.14 << endl;
            //cout << "Acceleration: " << t << "   " << Acceleration(0) << "   " << Acceleration(1) << "   " << Acceleration(2) << "\n" << endl;
            };

}; // End of namespace orbit



