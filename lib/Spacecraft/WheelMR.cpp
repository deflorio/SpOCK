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

#include <WheelMR.h>
#include <Constants.h>

//#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;
using namespace timescales;
//using boost::math::normal;
//using namespace spaceenvironment;

namespace mrwheel
    {
    //------------------------------------------------------------------------------
    // SOLRAD implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    MRWHEEL::~MRWHEEL() {};
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
    void MRWHEEL::Init()
                        {
                        Iw = ConstPrm(0);
                        wheelspeed_SYS(2) = AuxPrm(0);
                        wheelspeed = (SC2SYS.transpose())*wheelspeed_SYS;
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
    void MRWHEEL::status(double epoch,
                         const Ref<const VectorXd>& currentstate,
                         const Ref<const VectorXd>& auxstate)
                        {
                        subsystem_on = true;
                        };
    //------------------------------------------------------------------------------
    // Method VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
    //------------------------------------------------------------------------------
    /**
      * Compute the subsystem output (Output data) given the input data
      *
      * @param epoch          GPS epoch (seconds) of the input state
      * @param currentstate   Spacecraft attitude state (quaternion and angular rotation vector)
      * @param auxstate       Spacecraft orbital state vector (ECI)
      *
      * @return Readings of momentum wheels speed measurements or rate sensors or momentum wheels' torque vector depending on the value of class member 'Name'.
      *         Output momentum wheels speed measurements [rpm]
      *         Output of the rate sensor is the spacecraft estimated angular rotation vector [centdeg/s] (sensor frame)
      *         Output of the reaction wheels is the torque vector [Nm] produced by the actuator (SC body-fixed frame).
      */
    //------------------------------------------------------------------------------  
    VectorXd MRWHEEL::Output(double epoch,
                             const Ref<const VectorXd>& currentstate,
                             const Ref<const VectorXd>& auxstate)
                            {
                            VectorXd output;
                            
                            Vec3d omegaSC = currentstate.segment(4,3);
                            
                            if(On)
                              {
                              double err;
                              
                              if(Name.find("Wheel") != string::npos)
                                {
                                Tw_SC = Vec3d::Zero();
                                
                                Vec3d d_wheelspeed_SYS = Vec3d::Zero();
                                Vec3d d_wheelspeed = Vec3d::Zero();
                                // Wheel speed
                                if(wheelspeedCMD.norm() >= 1.0) d_wheelspeed_SYS = wheelspeedCMD - wheelspeed_SYS;
                                wheelspeed_SYS = wheelspeed_SYS + d_wheelspeed_SYS;
                                wheelspeed = (SC2SYS.transpose())*wheelspeed_SYS;
                                d_wheelspeed = (SC2SYS.transpose())*d_wheelspeed_SYS;
                                // Angular momentum and torque in spacecraft frame
                                hw_SC = Iw*wheelspeed/SEC2MIN; // hw is the wheel's angular momentum in SC frame
                                Tw_SC = -Iw*d_wheelspeed/SEC2MIN;
                                
                                for(int i = 0; i < 3; i++)
                                    {
                                    err = Error2(generator2); // Measurement error
                                    if(Tw_SC(i) != 0.0) Tw_SC(i) = Tw_SC(i) + err;
                                    }
                                //output = 0.0;
                                
                                output = wheelspeed_SYS;
                                
                                for(int i = 0; i < 3; i++)
                                    {
                                    err = Error1(generator1); // Measurement error
                                    output(i) = output(i) + err;
                                    }
                                
                                //output = wheelspeed;
                                //
                                //err = Error1(generator1); // Measurement error
                                //output(0) = output(0) + err;
                                
                                //cout << "wheelspeed_SYS: " << wheelspeed_SYS(0) << "   " << wheelspeed_SYS(1) << "   " << wheelspeed_SYS(2) << endl;
                                //cout << "wheelspeed: " << wheelspeed(0) << "   " << wheelspeed(1) << "   " << wheelspeed(2) << endl;
                                //cout << "wheelspeedCMD: " << wheelspeedCMD(0) << "   " << wheelspeedCMD(1) << "   " << wheelspeedCMD(2) << endl;
                                //cout << "d_wheelspeed_SYS: " << d_wheelspeed_SYS(0) << "   " << d_wheelspeed_SYS(1) << "   " << d_wheelspeed_SYS(2) << endl;
                                //cout << "d_wheelspeed: " << d_wheelspeed(0) << "   " << d_wheelspeed(1) << "   " << d_wheelspeed(2) << "\n" << endl;
                                }
                                
                              if(Name.compare("Rate Sensor") == 0)
                                {
                                output = SC2SYS*omegaSC;
                                
                                for(int i = 0; i < 3; i++)
                                    {
                                    err = Error1(generator1); // Measurement error
                                    output(i) = output(i) + err;
                                    }
                                
                                output = output*RAD2DEG;
                                }
                              }
                            else
							  {
                              if(Name.compare("Rate Sensor") == 0) output = Vec3d::Zero();
                                
							  if(Name.find("Wheel") != string::npos)
								{
								wheelspeed_SYS = Vec3d::Zero();
								output = wheelspeed_SYS;
								hw_SC = Vec3d::Zero();
								Tw_SC = Vec3d::Zero();
								}
							  }
                                
                            return(output);     
                            };
    
}; // End of namespace spaceenvironment



