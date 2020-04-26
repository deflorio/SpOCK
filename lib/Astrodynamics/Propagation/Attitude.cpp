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

#include <Attitude.h>
#include <Constants.h>
#include <Transformations.h>

using namespace std;
using namespace math;
using namespace boost::numeric::odeint;
using namespace constants;
using namespace propagator;

namespace pl = std::placeholders;

namespace attitude
    {
    //------------------------------------------------------------------------------
    // Class ATT implementation
    //------------------------------------------------------------------------------
    Vec3d ATT::Torque = Vec3d::Zero();
    VectorNd<15> ATT::Torque_env = VectorNd<15>::Zero();
    Vec3d ATT::TorqueACT = Vec3d::Zero();
    Vec3d ATT::hw = Vec3d::Zero();
    state_type ATT::x = {};
    bool ATT::magnetometer_on = false;
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    ATT::~ATT() {};
    //------------------------------------------------------------------------------
    // Method void ForceModelsSetup()
    //------------------------------------------------------------------------------
    /**
     * Set up file paths and parameters for environmental forces models
     */
    //------------------------------------------------------------------------------ 
    void ATT::ForceModelsSetup()
            {
            if(drag_on)
                {
                Atmosphere.init_epoch = inittime;
                Atmosphere.simduration = simdur;
                Atmosphere.Drag_Model = "Panels";
                Atmosphere.CF = CD;
                //Atmosphere.Area = Area_D;
                Atmosphere.setfilespath(atmosphere,datapath);
                Atmosphere.getmodel_coeff();
                }
                
            if(srp_on)
                {
                SolarRadiation.SRP_Model = "Panels";
                SolarRadiation.CF = C_SRP;
                }
                
            if(mag_on || magnetometer_on)
                {
                MagneticField.init_epoch = inittime;
                MagneticField.setfilespath(magneticfield,datapath);
                MagneticField.getmodel_coeff();
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
    void ATT::DynModel(
            const state_type &x,
            state_type &dxdt,
            const double t)
            {
            double abs_t = inittime + t;
            Vec3d euler_ang;
            Vec4d state_quat;
            
            euler_ang << x[0], x[1], x[2];
            ECItoBody = RotationMatrix321(euler_ang(0),euler_ang(1),euler_ang(2));  
            // Compute torques
            state_quat = RotationMatrix2Quaternion(ECItoBody);
            
            ComputeAction(abs_t, state_quat, orbstate);
            
            dxdt[0] = x[4]*sin(x[2])/cos(x[1]) + x[5]*cos(x[2])/cos(x[1]);
            
            dxdt[1] = x[4]*cos(x[2]) - x[5]*sin(x[2]);
            
            dxdt[2] = x[3] + x[4]*sin(x[2])*sin(x[1])/cos(x[1]) + x[5]*cos(x[2])*sin(x[1])/cos(x[1]);
                    
            dxdt[3] = -( invMoI(0,0)*( ( MoI(2,0)*x[4] - MoI(1,0)*x[5] )*x[3] + ( MoI(2,1)*x[4] - MoI(1,1)*x[5] )*x[4] + ( MoI(2,2)*x[4] - MoI(1,2)*x[5] )*x[5] ) +
                         invMoI(0,1)*( ( MoI(0,0)*x[5] - MoI(2,0)*x[3] )*x[3] + ( MoI(0,1)*x[5] - MoI(2,1)*x[3] )*x[4] + ( MoI(0,2)*x[5] - MoI(2,2)*x[3] )*x[5] ) +
                         invMoI(0,2)*( ( MoI(1,0)*x[3] - MoI(0,0)*x[4] )*x[3] + ( MoI(1,1)*x[3] - MoI(0,1)*x[4] )*x[4] + ( MoI(1,2)*x[3] - MoI(0,2)*x[4] )*x[5] ) ) +
            
                      -( invMoI(0,0)*( -hw(1)*x[5] + hw(2)*x[4] ) + invMoI(0,1)*( hw(0)*x[5] - hw(2)*x[3] ) + invMoI(0,2)*( -hw(0)*x[4] + hw(1)*x[3] ) ) +
                           
                         invMoI(0,0)*Torque(0) + invMoI(0,1)*Torque(1) + invMoI(0,2)*Torque(2);
                        
            dxdt[4] = -( invMoI(1,0)*( ( MoI(2,0)*x[4] - MoI(1,0)*x[5] )*x[3] + ( MoI(2,1)*x[4] - MoI(1,1)*x[5] )*x[4] + ( MoI(2,2)*x[4] - MoI(1,2)*x[5] )*x[5] ) +
                         invMoI(1,1)*( ( MoI(0,0)*x[5] - MoI(2,0)*x[3] )*x[3] + ( MoI(0,1)*x[5] - MoI(2,1)*x[3] )*x[4] + ( MoI(0,2)*x[5] - MoI(2,2)*x[3] )*x[5] ) +
                         invMoI(1,2)*( ( MoI(1,0)*x[3] - MoI(0,0)*x[4] )*x[3] + ( MoI(1,1)*x[3] - MoI(0,1)*x[4] )*x[4] + ( MoI(1,2)*x[3] - MoI(0,2)*x[4] )*x[5] ) ) +
            
                      -( invMoI(1,0)*( -hw(1)*x[5] + hw(2)*x[4] ) + invMoI(1,1)*( hw(0)*x[5] - hw(2)*x[3] ) + invMoI(1,2)*( -hw(0)*x[4] + hw(1)*x[3] ) ) +
            
                         invMoI(1,0)*Torque(0) + invMoI(1,1)*Torque(1) + invMoI(1,2)*Torque(2);
                        
            dxdt[5] = -( invMoI(2,0)*( ( MoI(2,0)*x[4] - MoI(1,0)*x[5] )*x[3] + ( MoI(2,1)*x[4] - MoI(1,1)*x[5] )*x[4] + ( MoI(2,2)*x[4] - MoI(1,2)*x[5] )*x[5] ) +
                         invMoI(2,1)*( ( MoI(0,0)*x[5] - MoI(2,0)*x[3] )*x[3] + ( MoI(0,1)*x[5] - MoI(2,1)*x[3] )*x[4] + ( MoI(0,2)*x[5] - MoI(2,2)*x[3] )*x[5] ) +
                         invMoI(2,2)*( ( MoI(1,0)*x[3] - MoI(0,0)*x[4] )*x[3] + ( MoI(1,1)*x[3] - MoI(0,1)*x[4] )*x[4] + ( MoI(1,2)*x[3] - MoI(0,2)*x[4] )*x[5] ) ) +
            
                      -( invMoI(2,0)*( -hw(1)*x[5] + hw(2)*x[4] ) + invMoI(2,1)*( hw(0)*x[5] - hw(2)*x[3] ) + invMoI(2,2)*( -hw(0)*x[4] + hw(1)*x[3] ) ) +
            
                         invMoI(2,0)*Torque(0) + invMoI(2,1)*Torque(1) + invMoI(2,2)*Torque(2);
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
    void ATT::Maneuver(const Ref<const VectorXd>& maneuver)
                {
                TorqueACT = maneuver;
                };
    //------------------------------------------------------------------------------
    // Method void ComputeAction(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& orb_state)
    //------------------------------------------------------------------------------
    /**
     * Compute perturbation torques
     *
     * @param epoch         GPS epoch (seconds) of the input states
     * @param currentstate  Attitude state at epoch GPStime
     * @param orb_state     12-D vector containing ECI orbital state (0-5) and ECEF orbital state (6-11) at epoch GPStime
     */
    //------------------------------------------------------------------------------       
    void ATT::ComputeAction(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& orb_state)
            {
            Torque = Vec3d::Zero();
            Torque_env = VectorNd<15>::Zero();
            Vec3d pos_ECI, pos_ECEF, pos_Body, Tg_body;
            
            pos_ECI = orb_state.segment(0,3);
            pos_ECEF = orb_state.segment(6,3);
            
            Vec4d q_currentstate = currentstate;
            // Transformation ECI -> Body by means of quaternion
            pos_Body = TransbyQ(pos_ECI,q_currentstate);
            
            if(ggrad_on)
                {
                Vec3d Tg;
                double r = pos_ECI.norm();
                double r5 = r*r*r*r*r;
                
                //Tg(0) = (3.0*astro::GM_EARTH/r5)*( ( MoI(2,0)*pos_Body(1) - MoI(1,0)*pos_Body(2) )*pos_Body(0) +
                //                                   ( MoI(2,1)*pos_Body(1) - MoI(1,1)*pos_Body(2) )*pos_Body(1) +
                //                                   ( MoI(2,2)*pos_Body(1) - MoI(1,2)*pos_Body(2) )*pos_Body(2) );
                //
                //Tg(1) = (3.0*astro::GM_EARTH/r5)*( ( MoI(0,0)*pos_Body(2) - MoI(2,0)*pos_Body(0) )*pos_Body(0) +
                //                                   ( MoI(0,1)*pos_Body(2) - MoI(2,1)*pos_Body(0) )*pos_Body(1) +
                //                                   ( MoI(0,2)*pos_Body(2) - MoI(2,2)*pos_Body(0) )*pos_Body(2) );
                //
                //Tg(2) = (3.0*astro::GM_EARTH/r5)*( ( MoI(1,0)*pos_Body(0) - MoI(0,0)*pos_Body(1) )*pos_Body(0) +
                //                                   ( MoI(1,1)*pos_Body(0) - MoI(0,1)*pos_Body(1) )*pos_Body(1) +
                //                                   ( MoI(1,2)*pos_Body(0) - MoI(0,2)*pos_Body(1) )*pos_Body(2) );
                
                Vec3d Irpos = Vec3d::Zero();
                
                Irpos = MoI*pos_Body;
                
                Tg = (3.0*astro::GM_EARTH/r5)*pos_Body.cross(Irpos);
                
                Torque += Tg;
                
                Torque_env.segment(0,3) = Tg;
                }
            
            // Computation of magnetic torque   
            if(mag_on)
                {
                Vec3d Tm = Vec3d::Zero();
                Vec3d B_ECI = Vec3d::Zero();
                Vec3d B_Body = Vec3d::Zero();
                
                MagneticField.SetReferenceFrame("ECI");
                B_ECI = MagneticField.field_vec(epoch, pos_ECEF);
                
                B_Body = TransbyQ(B_ECI,q_currentstate);
                
                B_Body = B_Body*1e-9; // Conversion from nanotesla to tesla
                //cout << "Magnetic field vector (Body): " << B_Body(0) << "   " << B_Body(1) << "   " << B_Body(2) << "\n" << endl;
                Tm = Mdip.cross(B_Body);
                
                Torque += Tm;
                
                Torque_env.segment(3,3) = Tm;
               }
            // Computation of solar radiation pressure torque   
            if(srp_on)
                {
                Vec3d Ts = Vec3d::Zero();
                Vec3d F_SRP_ECI = Vec3d::Zero();
                Vec3d F_SRP_Body = Vec3d::Zero();
                Vec3d cP = Vec3d::Zero();
                
                SC::Face F;
                
                for(auto element : SC_Faces) // Ranged based for loop on map SC_Faces (member of Class PROP)
                   {
                   F = SC_Faces[element.first];
                   
                   cP = F.cP;
                   
                   SolarRadiation.SetSurfaceParameters(F, q_currentstate); // C_SRP is a member of Class PROP
                   F_SRP_ECI = SolarRadiation.field_vec(epoch, pos_ECI); // Solar radiation pressure force (ECI) on surface F
                   
                   F_SRP_Body = TransbyQ(F_SRP_ECI,q_currentstate);
                   
                   Ts += cP.cross(F_SRP_Body);
                   }
                
                //F_Xminus = SC_Faces["-X"];
                //F_Yplus = SC_Faces["+Y"];
                //F_Yminus = SC_Faces["-Y"];
                //F_Zplus = SC_Faces["+Z"];
                //F_Zminus = SC_Faces["-Z"];
                
                Torque += Ts;
                
                Torque_env.segment(6,3) = Ts;
                }
                
            // Computation of atmospheric drag torque   
            if(drag_on)
                {
                Vec3d Ts = Vec3d::Zero();
                Vec3d F_ATM_ECI = Vec3d::Zero();
                Vec3d F_ATM_Body = Vec3d::Zero();
                Vec3d cA = Vec3d::Zero();
                
                SC::Face F;
                
                for(auto element : SC_Faces) // Ranged based for loop on map SC_Faces (member of Class PROP)
                   {
                   F = SC_Faces[element.first];
                   
                   cA = F.cA;
                   
                   Atmosphere.SetSurfaceParameters(F, q_currentstate); // C_ATM is a member of Class PROP
                   F_ATM_ECI = Atmosphere.field_vec(epoch, orb_state); // Solar radiation pressure force (ECI) on surface F
                   
                   F_ATM_Body = TransbyQ(F_ATM_ECI,q_currentstate);
                   
                   Ts += cA.cross(F_ATM_Body);
                   }
                
                Torque += Ts;
                
                Torque_env.segment(9,3) = Ts;
                }    
            
            // Put environmental torque values in Torque_env
            Torque_env.segment(12,3) = Torque;
            // Add actuators torque to environmental torque   
            Torque += TorqueACT;    
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
    void ATT::Integrate(
            double t,
            double step)
            {
            if(integ_first_step)
               {
               for(int i = 0; i < 6; i++) x[i] = initstate(i);
               
               state = initstate;
               integ_first_step = false;
               cout << "First integration step" << endl;
               }
            // Execute one integration step
            //stepper.do_step(DynModel, x, t, step);
            //stepper.do_step(std::bind(&ATT::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            bulirsch_stoer_stepper.try_step(std::bind(&ATT::DynModel, *this , pl::_1 , pl::_2 , pl::_3), x, t, step);
            
            //cout << t << "   " << x[0]*180/3.14 << "   " << x[1]*180/3.14 << "   " << x[2]*180/3.14 << "   " << x[3]*180/3.14 << "   " << x[4]*180/3.14 << "   " << x[5]*180/3.14 << endl;
            //cout << "Torque: " << t << "   " << Torque(0) << "   " << Torque(1) << "   " << Torque(2) << "\n" << endl;
            
            for(int i = 0; i < 6; i++) state(i) = x[i];
            
            //cout << t << "   " << state(0)*180/3.14 << "   " << state(1)*180/3.14 << "   " << state(2)*180/3.14 << "   " << state(3)*180/3.14 << "   " << state(4)*180/3.14 << "   " << state(5)*180/3.14 << endl;
            //cout << "Torque: " << t << "   " << Torque(0) << "   " << Torque(1) << "   " << Torque(2) << "\n" << endl;
            };

}; // End of namespace attitude



