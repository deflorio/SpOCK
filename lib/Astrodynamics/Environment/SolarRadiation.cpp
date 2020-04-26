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
#include <iomanip>

#include <SolarRadiation.h>
#include <Constants.h>

using namespace std;
using namespace math;
using namespace constants;

namespace solradiation
    {
    //------------------------------------------------------------------------------
    // SOLRAD implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    SOLRAD::~SOLRAD() {};
    //------------------------------------------------------------------------------
    // Method void SetSurfaceParameters(double A, double C_SRP, Vec3d& n, const string mat)
    //------------------------------------------------------------------------------
    /**
     * Set up surface parameters and orientation required to compute the solar radiation pressure
     *
     * @param A       Surface area [m^2]
     * @param C_SRP   Solar radiation pressure coefficient
     * @param n       3-dimensional unit vector normal to the surface (ECI)
     * @param mat     Material composing the surface (@see file Constants.h for otpions)
     *
     */
    //------------------------------------------------------------------------------
    void SOLRAD::SetSurfaceParameters(SC::Face face_prms,
                                      Vec4d& q_attitude)
                                     {
                                     Area = face_prms.Area;
                                     Vec3d n_body;
                                     n_body = face_prms.n;
                                     
                                     Vec4d q_attitude_inv;
                                     q_attitude_inv = q_inv(q_attitude); // Attitude state inverted quaternion
                                     n = TransbyQ(n_body, q_attitude_inv); // Surface normal in ECI frame
                                     
                                     string mat = face_prms.Material;
                                     //CF = C_SRP;
                                     
                                     Vec3d opt_coefficients;
                                     try{ opt_coefficients = materials::optical_properties.at(mat); } // See Constants.h
                                     catch(const std::out_of_range& oor)
                                          {
                                           cerr << "Exception of type " << oor.what() << " in function void SOLRAD::SetSurfaceParameters." << endl;
                                           cerr << "One of the surface materials selected is a NOT valid input (see Constants.h) or a space has introducted in the name." << endl;
                                           exit(EXIT_FAILURE);
                                          }
                                     
                                     rho1 = opt_coefficients(0); // Specular reflectivity coefficient
                                     rho2 = opt_coefficients(1); // Diffuse reflectivity coefficient
                                     rho3 = opt_coefficients(2); // Transmitted portions of incoming photons 
                                     };
    //------------------------------------------------------------------------------
    // Method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate)
    //------------------------------------------------------------------------------
    /**
     * Implementation of the abstract method of class SPACEENV
     * Compute the force applied on a surface by the solar radiation pressure using a multiple-surfaces model
     *
     * @param time        GPS epoch (seconds) of the input state
     * @param orbstate    Spacecraft position vector (ECI)
     *
     * @return 3-dimensional solar radiation pressure force on the surface
     * considered
     *
     * @note In this function coefficient CF (@see method SetSurfaceParameters)
     * represents a scale factor to be estimated and should have a nominal value of 1
     */
    //------------------------------------------------------------------------------
    Vec3d SOLRAD::field_vec(double time,
                            const Ref<const VectorXd>& orbstate)
                            {
                            Vec3d F_SRP = Vec3d::Zero();
                            Vec3d s = Vec3d::Zero();
                            //bool sun_on = false;
                            bool eclipse, penumbra;
                            
                            Solar.eclipse(time, orbstate, eclipse, penumbra);
                            s = Solar.sundir_u(time, orbstate);// Versor from the spacecraft to the Sun
                            
                            //Vec3d sunpos = Vec3d::Zero();
                            //Vec3d sunsc = Vec3d::Zero();
                            //
                            //sunpos = sunposREC(time);
                            //sunsc = (sunpos - orbstate);
                            //
                            //double d0divd = astro::AU/sunsc.norm();
                            
                            //if( s.dot(n) > 0.0 ) sun_on = true; // The surface is facing the sun
                            //bool force_on = !eclipse && sun_on;
                            //cout << eclipse << " " << time << endl;
                            //if(force_on)
                            if(!eclipse)
                                {
                                if( SRP_Model.compare("Panels") == 0 )
                                    {
                                    if( s.dot(n) > 0.0 ) // The surface is facing the sun
                                        {
                                        F_SRP = -CF*astro::P*Area*( s.dot(n) )*( ( 1.0 - rho1 - rho3 )*s + 2.0*( rho1*( s.dot(n) ) + (1.0/3.0)*rho2 )*n );
                                        //if(force_on) F_SRP = -CF*astro::P*Area*( s.dot(n) )*( ( 1.0 - rho1)*s + 2.0*rho1*( s.dot(n) )*n );
                                        
                                        //cout.precision(20);
                                        //if( fabs(F_SRP(0)) != 0.0 || fabs(F_SRP(1)) != 0.0 || fabs(F_SRP(2)) != 0.0) cout << "F_SRP.norm() = " << F_SRP.norm() << endl;
                                        //cout << "rho1 = " << rho1 << " rho2 = " << rho2 << " rho3 = " << rho3 << endl;
                                        }
                                    }
                                else if( SRP_Model.compare("RefArea") == 0 )
                                    {
                                    F_SRP = -CF*astro::P*Area*s;
                                    //cout << "s(0) = " << s(0) << " s(1) = " << s(1) << " s(2) = " << s(2) << endl;
                                    //cout.precision(20);
                                    //if( fabs(F_SRP(0)) != 0.0 || fabs(F_SRP(1)) != 0.0 || fabs(F_SRP(2)) != 0.0) cout << "F_SRP.norm() = " << F_SRP.norm() << endl;
                                    }
                            
                                }
                            
                            return(F_SRP);
                            };
    //------------------------------------------------------------------------------
    // Method getmodel_coeff()
    //------------------------------------------------------------------------------
    /**
     * Get surface optical coefficients
     *
     * @return 3-dimensional vector containing in order the specular reflectivity,
     * diffuse reflectivity and transmitted portions of incoming photons coefficients
     * currently in use
     */
    //------------------------------------------------------------------------------
    void SOLRAD::getmodel_coeff()
                                {
                                //Vec3d coefficients(rho1, rho2, rho3);
                                //    
                                //return(coefficients);
                                };

}; // End of namespace spaceenvironment



