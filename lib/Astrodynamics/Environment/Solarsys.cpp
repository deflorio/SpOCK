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

#include <Solarsys.h>
#include <Constants.h>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;

namespace solarsystem
    {
    //------------------------------------------------------------------------------
    // SOLSYS implementation
    //------------------------------------------------------------------------------
    // Constructor
    //------------------------------------------------------------------------------
    SOLSYS::SOLSYS() {};
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    SOLSYS::~SOLSYS() {};
    //------------------------------------------------------------------------------
    // Method load_ephemeris_file(const string ephe_filepath)
    //------------------------------------------------------------------------------
    /**
     * Load SPICE planet ephemeris file. This method is to be used in case the ephemeris are not loaded in the main file
     *
     * @param ephe_filepath   Path of SPICE planet ephemeris file
     *
     */
    //------------------------------------------------------------------------------  
    void SOLSYS::load_ephemeris_file(const string ephe_filepath)
                                {
                                furnsh_c(ephe_filepath.c_str( ));
                                };
    //------------------------------------------------------------------------------
    // Method Vec3d OBJ_pos(const char* OBJtarget_name, const char* OBJobserver_name, const char* refframe, const char* aberration_corr, const double time)
    //------------------------------------------------------------------------------
    /**
     * Return the position of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration
     *
     * @note This function is based on the NASA SPICE library's function spkpos_c
     *
     * @param OBJtarget_name     Target body name
     * @param OBJobserver_name   Observing body name
     * @param refframe           Reference frame of output position vector
     * @param aberration_corr    Aberration correction flag
     * @param time               Observer epoch (seconds past J2000 TDB)
     *
     * @return Position vector at epoch of desired object in desired frame [m]
     */
    //------------------------------------------------------------------------------  
    Vec3d SOLSYS::OBJ_pos(const char* OBJtarget_name,
                            const char* OBJobserver_name,
                            const char* refframe,
                            const char* aberration_corr,
                            const double time)
                            {
                            Vec3d pos = Vec3d::Zero();
                            
                            double positionAtEpoch[3], lightTime; // Declare variables for cartesian position and light-time to be determined by SPICE
                            
                            spkpos_c(OBJtarget_name, time, refframe, aberration_corr, OBJobserver_name, positionAtEpoch, &lightTime );

                            for( int i = 0; i < 3 ; i++ ) pos(i) = positionAtEpoch[i]*1E3; // Put result in Eigen Vector

                            return(pos);    
                            };
    //------------------------------------------------------------------------------
    // Method Vector6d OBJ_posvel(const char* OBJtarget_name, const char* OBJobserver_name, const char* refframe, const char* aberration_corr, const double time)
    //------------------------------------------------------------------------------
    /**
     * Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration
     *
     * @note This function is based on the NASA SPICE library's function spkezr_c
     *
     * @param OBJtarget_name     Target body name
     * @param OBJobserver_name   Observing body name
     * @param refframe           Reference frame of output position vector
     * @param aberration_corr    Aberration correction flag
     * @param time               Observer epoch (seconds past J2000 TDB)
     *
     * @return State vector at epoch of desired object in desired frame [m], [m/s]
     */
    //------------------------------------------------------------------------------  
    Vector6d SOLSYS::OBJ_posvel(const char* OBJtarget_name,
                                  const char* OBJobserver_name,
                                  const char* refframe,
                                  const char* aberration_corr,
                                  const double time)
                                  {
                                   Vector6d posvel = Vector6d::Zero();
                                   // Declare variables for cartesian position and light-time to be determined by Spice.
                                   double stateAtEpoch[6];
                                   double lightTime;
    
                                   // Call Spice function to calculate position and light-time.
                                   spkezr_c(OBJtarget_name, time, refframe, aberration_corr, OBJobserver_name, stateAtEpoch, &lightTime );

                                   // Put result in Eigen Vector.
                                   for( int i = 0; i < 6 ; i++ ) posvel(i) = stateAtEpoch[i]*1E3;
    
                                   return(posvel);
                                   };
    //------------------------------------------------------------------------------
    // Method Vec3d sunposREC(double GPStime)
    //------------------------------------------------------------------------------
    /**
     * Return the position of the sun in rectangular coordinates (ECI) from JPL ephemerides
     *
     * @note This function is based on the NASA SPICE library's function spkpos_c
     *       and uses the JPL ephemerides
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     *
     * @return Position vector of the Sun in ECI Cartesian coordinates [m]
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::sunposREC(double GPStime)
                            {
                            Vec3d sunpos = Vec3d::Zero();
                            
                            const string OBJtarget = "SUN";
                            //const string OBJobserver = "EARTH BARYCENTER";
                            const string OBJobserver = "EARTH BARYCENTER";
                            const string reference_frame = "J2000";
                            const string abcorr = "NONE";
                            double time = GPS2ET(GPStime);
                            
                            //sunpos = OBJ_pos(OBJtarget.c_str( ), OBJobserver.c_str( ), reference_frame.c_str( ), abcorr.c_str( ), time);
                            
                            double positionAtEpoch[3];
                            double lightTime;
                            spkpos_c(OBJtarget.c_str( ), time, reference_frame.c_str(), abcorr.c_str(), OBJobserver.c_str(), positionAtEpoch, &lightTime );

                            // Put result in Eigen Vector.
                            for( int i = 0; i < 3 ; i++ ) sunpos(i) = positionAtEpoch[i]*1E3;
                            
                            return(sunpos);
                            };
    //------------------------------------------------------------------------------
    // Method Vec3d LP_sunposREC(double GPStime)
    //------------------------------------------------------------------------------
    /**
     * Return the position of the Sun in rectangular coordinates (ECI) from low-precision (LP) Solar coordinates
     *
     * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
     *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X.
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     *
     * @return Position vector of the Sun in ECI Cartesian coordinates [m]
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::LP_sunposREC(double GPStime)
                            {
                            Vec3d sunpos = Vec3d::Zero();
                            Mat3x3d Tx_eps = Mat3x3d::Zero();
                            Vec3d Vec_rL = Vec3d::Zero();
                            
                            double eps, JD, T, L_, L, intpart, M_, M, r;
                            // Obliquity of J2000 ecliptic
                            eps = astro::EPS_J2000*DEG2RAD;
                            // Modified Julian date
                            JD = GPS2JD(GPStime);
                            // Julian centuries since J2000
                            T = (JD - timescales::JD_J2000)/36525.0;
                            
                            // Mean anomaly, ecliptic longitude and radius
                            M_ = 0.9931267 + 99.9973583*T;
                            M = PI2*modf(M_, &intpart); // [rad]
                            
                            L_ = 0.7859444 + M/PI2 + ( 6892.0*sin(M) + 72.0*sin(2.0*M) )/1296.0e3;
                            L = PI2*modf(L_, &intpart); // [rad]
                            
                            L = L + 1.3972*DEG2RAD*T;
                            r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2*M); // [m]
                            
                            Tx_eps = Rot_x(-eps);
                            Vec_rL << r*cos(L), r*sin(L), 0.0;
                        
                            // Equatorial position vector
                            sunpos = Tx_eps*Vec_rL;
                            
                            return(sunpos);
                            };
    //------------------------------------------------------------------------------
    // Method Vec3d moonposREC(double GPStime)
    //------------------------------------------------------------------------------
    /**
     * Return the position of the Moon in rectangular coordinates (ECI) from JPL ephemerides
     *
     * @note This function is based on the NASA SPICE library's function spkpos_c
     *       and uses the JPL ephemerides
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     *
     * @return Position vector of the Moon in ECI Cartesian coordinates [m]
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::moonposREC(double GPStime)
                            {
                            Vec3d moopos = Vec3d::Zero();
                            
                            const string OBJtarget = "MOON";
                            const string OBJobserver = "EARTH BARYCENTER";
                            const string reference_frame = "J2000";
                            const string abcorr = "NONE";
                            double time = GPS2ET(GPStime);
                            
                            //moopos = OBJ_pos(OBJtarget.c_str( ), OBJobserver.c_str( ), reference_frame.c_str( ), abcorr.c_str( ), time);
                            
                            double positionAtEpoch[3];
                            double lightTime;
                            spkpos_c(OBJtarget.c_str( ), time, reference_frame.c_str( ), abcorr.c_str( ), OBJobserver.c_str( ), positionAtEpoch, &lightTime );

                            // Put result in Eigen Vector.
                            for( int i = 0; i < 3 ; i++ ) moopos(i) = positionAtEpoch[i]*1E3;
                            
                            return(moopos);
                            };
    //------------------------------------------------------------------------------
    // Method Vec3d LP_moonposREC(double GPStime)
    //------------------------------------------------------------------------------
    /**
     * Return the position of the Moon in rectangular coordinates (ECI) from low-precision (LP) Lunar coordinates
     *
     * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
     *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X.
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     *
     * @return Position vector of the Moon in ECI Cartesian coordinates [m]
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::LP_moonposREC(double GPStime)
                            {
                            Vec3d moopos = Vec3d::Zero();
                            Mat3x3d Tx_eps = Mat3x3d::Zero();
                            Vec3d Vec_rL = Vec3d::Zero();
                            
                            double eps, JD, T;
                            double L_0_, L_0, l_, l, lp_, lp, D_, D, F_, F, dL, S, h, N, intpart;
                            double L_, L, B, R, cosB;
    
                            // Obliquity of J2000 ecliptic
                            eps = astro::EPS_J2000*DEG2RAD;
                            // Modified Julian date
                            JD = GPS2JD(GPStime);
                            // Julian centuries since J2000
                            T = (JD - timescales::JD_J2000)/36525.0;

                            // Mean elements of lunar orbit
                            L_0_ = 0.606433 + 1336.851344*T;
                            L_0 = modf(L_0_, &intpart); // Mean longitude [rev] w.r.t. J2000 equinox
                            
                            l_ = 0.374897 + 1325.552410*T;
                            l = PI2*modf(l_, &intpart); // Moon's mean anomaly [rad]
                            
                            lp_ = 0.993133 +   99.997361*T;
                            lp  = PI2*modf(lp_, &intpart); // Sun's mean anomaly [rad]
                            
                            D_ = 0.827361 + 1236.853086*T;
                            D   = PI2*modf(D_, &intpart); // Diff. long. Moon-Sun [rad]
                            
                            F_ = 0.259086 + 1342.227825*T;
                            F   = PI2*modf(F_, &intpart); // Argument of latitude

                            // Ecliptic longitude
                            dL = +22640*sin(l) - 4586*sin(l-2*D) + 2370*sin(2*D) +  769*sin(2*l) 
                                 -668*sin(lp) - 412*sin(2*F) - 212*sin(2*l-2*D) - 206*sin(l+lp-2*D)
                                 +192*sin(l+2*D) - 165*sin(lp-2*D) - 125*sin(D) - 110*sin(l+lp)
                                 +148*sin(l-lp) - 55*sin(2*F-2*D);
                                 
                            L_ = L_0 + dL/1296.0e3;

                            L = PI2*modf(L_, &intpart); // [rad] w.r.t. equinox of J2000
                            L = L + 1.3972*DEG2RAD*T; // w.r.t. equinox of Mjd_TT

                            // Ecliptic latitude
                            S  = F + (dL + 412*sin(2*F) + 541*sin(lp))/RAD2ARCS;
                            
                            h  = F - 2*D;
                            N  = -526*sin(h) + 44*sin(l + h) - 31*sin(-l + h) - 23*sin(lp + h) 
                                 +11*sin(-lp + h) - 25*sin(-2*l + F) + 21*sin(-l + F);
                        
                            B = ( 18520.0*sin(S) + N )/RAD2ARCS; // [rad]
                        
                            cosB = cos(B);

                            // Distance [m]
                            R = 385000e3 - 20905e3*cos(l) - 3699e3*cos(2*D - l) - 2956e3*cos(2*D)
                                -570e3*cos(2*l) + 246e3*cos(2*l - 2*D) - 205e3*cos(lp - 2*D) 
                                -171e3*cos(l + 2*D) - 152e3*cos(l + lp - 2*D);   

                            // Equatorial coordinates
                            Tx_eps = Rot_x(-eps);
                            Vec_rL << R*cos(L)*cosB, R*sin(L)*cosB, R*sin(B);
                        
                            moopos = Tx_eps*Vec_rL;
                            
                            return(moopos);
                            };   
    //------------------------------------------------------------------------------
    // Method Vec3d sunposRAD(double GPStime)
    //------------------------------------------------------------------------------
    /**
     * Return right ascension, declination and range of the sun (ECI)
     *
     * @note This function is based on the NASA SPICE library's function spkpos_c
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     *
     * @return Right ascension, declination [rad] and range [m] in ECI
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::sunposRAD(double GPStime)
                            {
                            Vec3d sunpos_rad = Vec3d::Zero();
                            Vec3d sunpos = Vec3d::Zero();
                            double sunpos_arr[3];
                            
                            double ra, dec, range;
                            
                            sunpos = sunposREC(GPStime);
                            for( int i = 0; i < 3 ; i++ ) sunpos_arr[i] = sunpos(i);
                            
                            // Convert in right ascension, declination and range, 
                            recrad_c(sunpos_arr, &range, &ra, &dec);
                            
                            sunpos_rad(0) = ra;
                            sunpos_rad(1) = dec;
                            sunpos_rad(2) = range;
                            
                            return(sunpos_rad);
                            };                               
    //------------------------------------------------------------------------------
    // Method double xi_angle(double GPStime, Vec3d& SC_pos)
    //------------------------------------------------------------------------------
    /**
     * Return the angle (xi) between the spacecraft position vector(ECI) and the Sun position vector (ECI)
     *
     * @note This function is based on the NASA SPICE library's function spkpos_c
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     * @param SC_pos      Spacecraft position vector
     *
     * @return xi angle [rad]
     */
    //------------------------------------------------------------------------------ 
    double SOLSYS::xi_angle(double GPStime,
                            const Vec3d& SC_pos)
                            {
                            Vec3d sunpos = Vec3d::Zero();
                            
                            sunpos = sunposREC(GPStime);
                            
                            double xi_angle = acos( sunpos.dot(SC_pos)/(sunpos.norm()*SC_pos.norm()) );
                            
                            return(xi_angle);
                            };
    //------------------------------------------------------------------------------
    // Method void eclipse(double GPStime, const Vec3d& SC_pos, bool& umbra, bool& penumbra)
    //------------------------------------------------------------------------------
    /**
     * Return the eclipse condition of the spacecraft at a certain time
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     * @param SC_pos      Spacecraft position vector
     *
     * @return umbra and penumbra conditions (true or false)
     */
    //------------------------------------------------------------------------------ 
    void SOLSYS::eclipse(double GPStime,
                         const Vec3d& SC_pos,
                         bool& umbra,
                         bool& penumbra)
                        {
                        umbra = false;
                        penumbra = false;
                        
                        Vec3d s_u = sunposREC(GPStime);
                        s_u = s_u.normalized();
                        double r_dot_s = SC_pos.dot(s_u);
                        Vec3d r_s = r_dot_s*s_u;
                        Vec3d delta = SC_pos - r_s;
                        double r_s_norm = r_s.norm();
                        double delta_norm = delta.norm();
                        
                        double D_Earth = 2.0*astro::R_EARTH;
                        double D_Sun = 2.0*astro::R_SUN;
                        
                        double Xi_p = D_Earth*astro::AU/(D_Sun + D_Earth);
                        double alpha_p = asin( D_Earth/(2.0*Xi_p) );
                        // Penumbra terminator
                        double p_t = (Xi_p + r_s_norm)*tan(alpha_p);
                        
                        double Xi_u = D_Earth*astro::AU/(D_Sun - D_Earth);
                        double alpha_u = asin( D_Earth/(2.0*Xi_u) );
                        // Umbra terminator
                        double u_t = (Xi_u + r_s_norm)*tan(alpha_u);
                        
                        if(r_dot_s < 0.0) // For r_dot_s > 0.0 the spacecraft is always in the Sun
                          {
                          if( (delta_norm > u_t) && (delta_norm <= p_t) )
                              {
                              umbra = false;
                              penumbra = true;    
                              }
                          if( delta_norm <= u_t)
                              {
                              umbra = true;
                              penumbra = true;    
                              }
                          }
                        };
    //------------------------------------------------------------------------------
    // Method Vec3d sundir(double GPStime, Vec3d& SC_pos)
    //------------------------------------------------------------------------------
    /**
     * Return the unit vector (ECI) of sun direction from the spacecraft's center of mass
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     * @param SC_pos      Spacecraft position vector
     *
     * @return Unit vector of spacecraft-to-sun direction in ECI frame
     */
    //------------------------------------------------------------------------------ 
    Vec3d SOLSYS::sundir_u(double GPStime,
                           const Vec3d& SC_pos)
                          {
                          Vec3d sunpos = Vec3d::Zero();
                          Vec3d sunsc = Vec3d::Zero();
                          Vec3d sunsc_u = Vec3d::Zero();
                          
                          sunpos = sunposREC(GPStime);
                        
                          sunsc = (sunpos - SC_pos);
                          sunsc_u = sunsc.normalized();
                        
                          return(sunsc_u);
                          };
    //------------------------------------------------------------------------------
    // Method double sun_angle(double GPStime, Vec3d& SC_pos, Vec3d& v3D_u)
    //------------------------------------------------------------------------------
    /**
     * Brief Return the dot product angle between the sun-spacecraft direction unit vector and a generic unit vector
     *
     * @param GPStime     GPS epoch (seconds) of the spacecraft position vector
     * @param SC_pos      Spacecraft position vector
     * @param v3D_u       Input unit vector
     *
     * @return Cosine of the angle between the sun-spacecraft direction unit vector and the input unit vector
     */
    //------------------------------------------------------------------------------ 
    double SOLSYS::sun_angle(double GPStime, Vec3d& SC_pos, Vec3d& v3D_u)
                        {
                        Vec3d sunsc_u = Vec3d::Zero();                        
                        
                        sunsc_u = sundir_u(GPStime, SC_pos);
                        
                        double cos_alpha = sunsc_u.dot(v3D_u);
                        
                        return(cos_alpha);
                        };

}; // End of namespace solarsystem



