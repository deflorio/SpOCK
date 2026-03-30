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
#include <omp.h>

#include <Tides.h>
#include <Constants.h>
#include <IO_utils.h>

using namespace std;
using namespace math;
using namespace constants;

namespace tides
    {
    //------------------------------------------------------------------------------
    // TIDES implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    TIDES::~TIDES() {};
    //------------------------------------------------------------------------------
    // Method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate)
    //------------------------------------------------------------------------------
    /**
     * Implementation of the abstract method of class SPACEENV
     * Compute the third-body perturbing acceleration vector due to the Sun and Moon gravity fields
     *
     * @param time        GPS epoch (seconds) of the input state
     * @param orbstate    Spacecraft position vector (ECI)
     *
     * @return 3-dimensional third-body perturbing acceleration vector (ECI)
     *
     */
    //------------------------------------------------------------------------------
    Vec3d TIDES::field_vec(double time,
                            const Ref<const VectorXd>& orbstate)
                            {
                            Vec3d acc = Vec3d::Zero();
                            Vec3d sun_acc = Vec3d::Zero();
                            Vec3d moon_acc = Vec3d::Zero();
                            
                            Vec3d SC_posECI = Vec3d::Zero();
                            Vec3d u_SC = Vec3d::Zero();
                            Vec3d sunpos = Vec3d::Zero();
                            Vec3d u_Sun = Vec3d::Zero();
                            Vec3d moonpos = Vec3d::Zero();
                            Vec3d u_Moon = Vec3d::Zero();
                            
                            //Vec3d SCtoSunpos = Vec3d::Zero();
                            //Vec3d SCtoMoonpos = Vec3d::Zero();
                            //double norm_SCtoSunpos, norm_SCtoMoonpos;
                            
                            double norm_sunpos, norm_sunpos_3, norm_moonpos, norm_moonpos_3, R_EARTH_5;
                            double SC_posECI_norm, SC_posECI_norm_4;
                            double costheta;
                            
                            R_EARTH_5 = astro::R_EARTH*astro::R_EARTH*astro::R_EARTH*astro::R_EARTH*astro::R_EARTH;
                            
                            // Position of spacecraft in ECI
                            SC_posECI = orbstate;
                            SC_posECI_norm = SC_posECI.norm();
                            SC_posECI_norm_4 = SC_posECI_norm*SC_posECI_norm*SC_posECI_norm*SC_posECI_norm;
                            u_SC = SC_posECI.normalized();
                            
                            // Position at epoch of Sun and Moon in ECI
                            if( modelname.compare("LP") == 0 )
                                {
                                sunpos = Solar.LP_sunposREC(time);
                                moonpos = Solar.LP_moonposREC(time);    
                                }
                            #ifdef USE_SPICE
                            if( modelname.compare("DE") == 0 )
                                {
                                sunpos = Solar.sunposREC(time);
                                moonpos = Solar.moonposREC(time);
                                }
                            #endif
                            
                            // Position of Sun and Moon relative to the spacecraft
                            //SCtoSunpos = sunpos - SC_posECI;
                            //SCtoMoonpos = moonpos - SC_posECI;
                            
                            // Position vectors norms
                            norm_sunpos = sunpos.norm();
                            norm_sunpos_3 = norm_sunpos*norm_sunpos*norm_sunpos;
                            
                            norm_moonpos = moonpos.norm();
                            norm_moonpos_3 = norm_moonpos*norm_moonpos*norm_moonpos;
                            
                            u_Sun = sunpos.normalized();
                            u_Moon = moonpos.normalized();
                            
                            //norm_SCtoSunpos = SCtoSunpos.norm();
                            //norm_SCtoMoonpos = SCtoMoonpos.norm();
                            // Accelerations
                            costheta = u_SC.dot(u_Sun);
                            sun_acc = (astro::Love_k2/2.0)*(astro::GM_SUN/norm_sunpos_3)*(R_EARTH_5/SC_posECI_norm_4)*( (3.0 - 15.0*costheta*costheta)*u_SC + 6.0*costheta*u_Sun );
                            
                            costheta = u_SC.dot(u_Moon);
                            moon_acc = (astro::Love_k2/2.0)*(astro::GM_MOON/norm_moonpos_3)*(R_EARTH_5/SC_posECI_norm_4)*( (3.0 - 15.0*costheta*costheta)*u_SC + 6.0*costheta*u_Moon );
                            
                            acc = sun_acc + moon_acc;
                            
                            return(acc);
                            };					
    //------------------------------------------------------------------------------
    // Method getmodel_coeff()
    //------------------------------------------------------------------------------
    /**
     * Get coefficients
     *
     */
    //------------------------------------------------------------------------------
    void TIDES::getmodel_coeff() {};                        
    
}; // End of namespace tides



