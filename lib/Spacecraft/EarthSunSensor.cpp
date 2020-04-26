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

#include <EarthSunSensor.h>
#include <Constants.h>

//#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;
//using boost::math::normal;

namespace earthsun
    {
    //------------------------------------------------------------------------------
    // SOLRAD implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    EARTHSUNSENS::~EARTHSUNSENS() {};
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
    void EARTHSUNSENS::Init()
                            {
                            if(Name.find("CSS") != string::npos)
                                {
                                I0 = ConstPrm(0);
                                nu = ConstPrm(1);
                                }
                            s_uvec = Vec3d::Zero();
                            s_uvec_SC = Vec3d::Zero();
                            e_uvec = Vec3d::Zero();
                            e_uvec_SC = Vec3d::Zero();
                            capture = 0.0;
                            detection = 0.0;
                            theta = 0.0;
                            costheta = 0.0;
                            sintheta = 0.0;
                            fov = 0.0;
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
    void EARTHSUNSENS::status(double epoch,
                              const Ref<const VectorXd>& currentstate,
                              const Ref<const VectorXd>& auxstate)
                              {
                              subsystem_on = false;
                              bool eclipse, penumbra;
                              
                              s_uvec = Solar.sundir_u(epoch, auxstate);          
                              e_uvec = -auxstate.normalized();  // e_uvec is the unit vector from the spacecraft to the Earth center, hence the sign '-'
                              
                              // OPS_limits(0) = semi-FOV of Sun camera; OPS_limits(1) = semi-FOV of Earth camera;
                              Vec3d Z_ECI = Vec3d::Zero();
                              Vec4d q_currentstate_inv = Vec4d::Zero();
                              
                              Vec4d q_currentstate = currentstate;
                              q_currentstate_inv = q_inv(q_currentstate);
                              Z_ECI = TransbyQ(Z,q_currentstate_inv);
                            
                              Solar.eclipse(epoch, auxstate, eclipse, penumbra);
                              
                              if( (Name.compare("Sun Camera") == 0) || (Name.find("CSS") != string::npos) )
                                {
                                bool sun_on = false;
                                
                                costheta = s_uvec.dot(Z_ECI);
                                sintheta = sqrt(1.0 - costheta*costheta);
                                theta = atan2(sintheta, costheta);
                                theta = mod(theta, PI2);
                                fov = fabs(2.0*theta); // atan2 coputes theta in the interval [-pi,+pi] radians
                                
                                if( fov <= OPS_limits(0) ) sun_on = true; // The Sun is inside the camera FOV
                                          
                                subsystem_on = !eclipse && sun_on; 
                                }
                                  
                              if(Name.compare("Earth Camera") == 0)
                                {
                                bool earth_on = false;
                                
                                costheta = e_uvec.dot(Z_ECI);
                                sintheta = sqrt(1.0 - costheta*costheta);
                                theta = atan2(sintheta, costheta);
                                theta = mod(theta, PI2);
                                fov = fabs(2.0*theta); // atan2 coputes theta in the interval [-pi,+pi] radians
                                if( fov <= OPS_limits(0) ) earth_on = true; // The Earth is inside the camera FOV
                                        
                                subsystem_on = !eclipse && earth_on;
                                }
                                
                              //if(Name.compare("Coarse Sun Sensor") == 0) subsystem_on = !eclipse; // Since the coarse sun sensors are on all the six faces of the spacecraft, the only off-condition is when the spacecraft is in the eclipse part of the orbit
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
      * @return Readings of Sun camera or Earth camera or coarse Sun sensor depending on the value of class member 'Name'.
      *         Output of Sun and Earth cameras are Azimuth and Elevation angles in the camera frame [centideg].
      *         Output of the coarse Sun sensors are the currents produced by the sun arrays composing the sensor [A].
      */
    //------------------------------------------------------------------------------  
    VectorXd EARTHSUNSENS::Output(double epoch,
                                  const Ref<const VectorXd>& currentstate,
                                  const Ref<const VectorXd>& auxstate)
                                  {
                                  VectorXd SensorReadings;
                                  //normal_distribution<double> Error1 = Error[0];
                                  //mt19937 generator1 = generator[0];
                                
                                  Vec4d q_currentstate = currentstate;
                                  s_uvec_SC = TransbyQ(s_uvec,q_currentstate);
                                  e_uvec_SC = TransbyQ(e_uvec,q_currentstate);
                                
                                  status(epoch, q_currentstate, auxstate); // Give a value 'true' or 'false' to class member 'subsystem_on'
                                
                                  if(On && subsystem_on)
                                    {
                                    double alpha = 0.0, betha = 0.0;
                                    Vec4d alphabetha;
                                  
                                    if(Name.compare("Sun Camera") == 0)
                                      {
                                      double s_uvec_SC_X, s_uvec_SC_Y, s_uvec_SC_Z;
                                      
                                      double err_alpha = Error1(generator1); // Measurement error
                                      double err_betha = Error1(generator1); // Measurement error
                                      
                                      s_uvec_SC_X = s_uvec_SC.dot(X);
                                      s_uvec_SC_Y = s_uvec_SC.dot(Y);
                                      s_uvec_SC_Z = s_uvec_SC.dot(Z);
                                  
                                      alpha = atan2(s_uvec_SC_X, s_uvec_SC_Z) + err_alpha; // Azimuth
                                      betha = atan2(s_uvec_SC_Y, s_uvec_SC_Z) + err_betha; // Elevation
                                      capture = 2.0;
                                      detection = 7.0;
                                   
                                      alphabetha(0) = (alpha*RAD2DEG)*100.0;
                                      alphabetha(1) = (betha*RAD2DEG)*100.0;
                                      alphabetha(2) = capture;
                                      alphabetha(3) = detection;
                                     
                                      SensorReadings = alphabetha;
                                      }
                                  
                                    if(Name.compare("Earth Camera") == 0)
                                      {
                                      double e_uvec_SC_X, e_uvec_SC_Y, e_uvec_SC_Z;
                                      
                                      double err_alpha = Error1(generator1); // Measurement error
                                      double err_betha = Error1(generator1); // Measurement error
                                 
                                      e_uvec_SC_X = e_uvec_SC.dot(X);
                                      e_uvec_SC_Y = e_uvec_SC.dot(Y);
                                      e_uvec_SC_Z = e_uvec_SC.dot(Z);
                                 
                                      alpha = atan2(e_uvec_SC_X, e_uvec_SC_Z) + err_alpha; // Azimuth
                                      betha = atan2(e_uvec_SC_Y, e_uvec_SC_Z) + err_betha; // Elevation
                                      capture = 2.0;
                                      detection = 7.0;
                                    
                                      alphabetha(0) = (alpha*RAD2DEG)*100.0;
                                      alphabetha(1) = (betha*RAD2DEG)*100.0;
                                      alphabetha(2) = capture;
                                      alphabetha(3) = detection;
                                    
                                      SensorReadings = alphabetha;   
                                      }
                                    
                                    if(Name.find("CSS") != string::npos)
                                        {
                                        VectorNd<1> Current = VectorNd<1>::Zero();
                                        
                                        double err = Error1(generator1); // Measurement error 
                                        double meas_costheta = costheta*cos(err) - sintheta*sin(err);
                                        
                                        Current(0) = fabs(I0*meas_costheta) + nu;
                                         
                                        SensorReadings = Current;
                                        }
                                    }
                                    else
                                        {
                                        if(Name.compare("Sun Camera") == 0)
                                          {
                                          Vec4d alphabetha;
                                       
                                          capture = 4.0;
                                          detection = 6.0;
                                       
                                          alphabetha(0) = 0.0;
                                          alphabetha(1) = 0.0;
                                          alphabetha(2) = capture;
                                          alphabetha(3) = detection;
                                       
                                          SensorReadings = alphabetha;
                                          }
                                         
                                        if(Name.compare("Earth Camera") == 0 )
                                          {
                                          Vec4d alphabetha;
                                       
                                          capture = 4.0;
                                          detection = 4.0;
                                       
                                          alphabetha(0) = 0.0;
                                          alphabetha(1) = 0.0;
                                          alphabetha(2) = capture;
                                          alphabetha(3) = detection;
                                       
                                          SensorReadings = alphabetha;
                                          }
                                         
                                        if(Name.find("CSS") != string::npos)
                                          {
                                          //VectorNd<10> Currents = VectorNd<10>::Zero();
                                          VectorNd<1> Current = VectorNd<1>::Zero();
                                          SensorReadings = Current;
                                          }
                                        }
                                  
                                  return(SensorReadings);     
                                  };
    
    }; // End of namespace spaceenvironment



