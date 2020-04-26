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

#include <SolarPanel.h>
#include <Constants.h>

//#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;

namespace solarpan
    {
    //------------------------------------------------------------------------------
    // SOLRAD implementation
    //------------------------------------------------------------------------------
    bool SOLARPAN::eclipse = false;
    double SOLARPAN::Surface = 0.0;
    double SOLARPAN::V_string = 0.0;
    double SOLARPAN::epsilon = 0.0;
    double SOLARPAN::costheta = 0.0;
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    SOLARPAN::~SOLARPAN() {};
    //------------------------------------------------------------------------------
    // Abstract method void Init(const Ref<const VectorXd>& ConstPrm)
    //------------------------------------------------------------------------------
    /**
     * Store in class members constant parameters contained in base class' members of
     * struct SYS_Parameters
     *
     * @param ConstPrm    Subsystem constant parameters
     * @see SYS_Parameters.ConstPrm in base class method Setup
     */
    //------------------------------------------------------------------------------      
    void SOLARPAN::Init()
                        {
                        Surface = ConstPrm(0);
                        V_string = ConstPrm(1);
                        epsilon = AuxPrm(0);
                        };
    //------------------------------------------------------------------------------
    // Method bool status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
    //------------------------------------------------------------------------------
    /**
     * Evaluate the operability status (subsystem on or off) and put its value in the class member subsystem_on
     *
     * @param epoch          GPS epoch (seconds) of the input state
     * @param currentstate   Spacecraft attitude state (quaternion)
     * @param auxstate       Spacecraft position vector (ECI)
     */
    //------------------------------------------------------------------------------  
    void SOLARPAN::status(double epoch,
                            const Ref<const VectorXd>& currentstate,
                            const Ref<const VectorXd>& auxstate)
                            {
                            subsystem_on = false;
                            s_uvec = Vec3d::Zero();
                            //e_uvec = Vec3d::Zero();
                            s_uvec = Solar.sundir_u(epoch, auxstate);          
                            //e_uvec = -auxstate.normalized();  // e_uvec is the unit vector from the spacecraft to the Earth center, hence the sign '-'
                            
                            // OPS_limits(0) = semi-FOV of Sun camera; OPS_limits(1) = semi-FOV of Earth camera;
                            Vec3d Z_ECI = Vec3d::Zero();
                            Vec4d q_currentstate_inv = Vec4d::Zero();
                            
                            Vec4d q_currentstate = currentstate;
                            q_currentstate_inv = q_inv(q_currentstate);
                            Z_ECI = TransbyQ(Z,q_currentstate_inv);
                            
                            bool penumbra;
                            Solar.eclipse(epoch, auxstate, eclipse, penumbra);
                            double sintheta, theta, fov;
                            
                            bool sun_on = false;
                            
                            costheta = s_uvec.dot(Z_ECI);
                            sintheta = sqrt(1.0 - costheta*costheta);
                            theta = atan2(sintheta, costheta);
                            theta = mod(theta, PI2);
                            fov = fabs(2.0*theta); // atan2 coputes theta in the interval [-pi,+pi] radians
                            //s_uvec = Solar.sundir_u(epoch, auxstate);
                            if( fov < OPS_limits(0) ) sun_on = true; // The Sun is inside the solar panel FOV
                                      
                            subsystem_on = !eclipse && sun_on;
                            };
    //------------------------------------------------------------------------------
    // Method VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
    //------------------------------------------------------------------------------
    /**
      * Compute the subsystem output (Output data) given the input data
      *
      * @param epoch          GPS epoch (seconds) of the input state
      * @param currentstate   Spacecraft attitude state (quaternion)
      * @param auxstate       Spacecraft position vector (ECI)
      *
      * @return 2D vector containing power [W] and current [A] generated by the solar panel
      */
    //------------------------------------------------------------------------------  
    VectorXd SOLARPAN::Output(double epoch,
                                const Ref<const VectorXd>& currentstate,
                                const Ref<const VectorXd>& auxstate)
                                {
                                VectorXd PI_out;
                                
                                Vec3d PI = Vec3d::Zero();
                                PI_out = PI;
                              
                                Vec4d q_currentstate = currentstate;
                                s_uvec_SC = TransbyQ(s_uvec,q_currentstate);
                                
                                status(epoch, q_currentstate, auxstate); // Give a value 'true' or 'false' to class member 'subsystem_on'
                                
                                if(On && subsystem_on)
                                    {
                                    PI_out(0) = Surface*astro::C_SUN*epsilon*costheta; // Generated power (costheta is computed in function status))
                                    PI_out(1) = PI_out(0)/V_string; // Generated current
                                    }
                                    
                                PI_out(2) = (double)eclipse;
                                
                                return(PI_out);     
                                  };
    
    }; // End of namespace solarpan



