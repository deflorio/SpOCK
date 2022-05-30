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

#include <Atmosphere.h>
#include <HarrisPriester.h>
#include <Constants.h>
#include <IO_utils.h>

using namespace std;
using namespace math;
using namespace constants;

extern "C"
        {
        // JB2008 Fortran
        extern void jb2008_(double*, double*, double*, double*, double*, double*,double*, double*, double*, double*, double*, int*, double*, double*);
        extern void solfsmy_(double*, double*, double*, double*, double*, double*, double*, double*, double*, char*);
        extern void dtcval_(double*, int*, char*, char*);
        // NRLMSISE-00 Fortran
        extern void gtd7d_(int*, float*, float*, float*, float*, float*, float*, float*, float*, int*, float*, float*);
        }

namespace atmosphere
    {
    double ATMO::f10 = 0.0;
    double ATMO::s10 = 0.0;
    double ATMO::y10 = 0.0;
    double ATMO::f10b = 0.0;
    double ATMO::s10b = 0.0;
    double ATMO::y10b = 0.0;
    double ATMO::m10 = 0.0;
    double ATMO::m10b = 0.0;
    int ATMO::dstdtc = 0.0;
    double ATMO::init_doy = 0.0;
    bool ATMO::idx_locked = false;
    int ATMO::swind = 0;
    char ATMO::c_dtc_file[200] = " ";
    char ATMO::c_dst_file[200] = " ";
    char ATMO::c_solfsmy_file[200] = " ";
    //SpW_row ATMO::SpaceWeather_idx_row = SpW_row::Zero();
    //------------------------------------------------------------------------------
    // ATMO implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    ATMO::~ATMO() {};
    
    //------------------------------------------------------------------------------
    // Method void AtmosphericDensity(double time, const Ref<const VectorXd>& orbstate)
    //------------------------------------------------------------------------------
    /**
     * Compute atmospheric density with a chosen model
     *
     * @see The JB2008 model's atmospheric density is computed by means of fortran routines.
     *      JB2008 model: http://sol.spacenvironment.net/~JB2008/
     * @see The implementation of the Harris-Priester model is based on
     *      Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
     *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X.
     *
     * @param time      GPS epoch (seconds) of the input state
     * @param SC_pos    Spacecraft state vector
     *
     */
    //------------------------------------------------------------------------------
    void ATMO::AtmosphericDensity(double time,
                                  const Ref<const VectorXd>& orbstate)
                                    {
                                    Vec3d posECI, posECEF, SCpos_RAD, SCpos_lonlath, sunpos;
                                    static double sat[3];
                                    posECI = orbstate.segment(0,3);
                                    posECEF = orbstate.segment(6,3);
                                    
                                    SCpos_RAD = ECI2RAD(posECI);
                                    SCpos_lonlath = ECEF2lonlath(posECEF);
                                        
                                    if( modelname.compare("JB2008") == 0 )
                                        {
                                        Vector6d UTCdate = Vector6d::Zero();
                                        double doy, amjd, T1950, D1950;
                                        static double sun[2];
                                        //static double sun[2], sat[3], temperatures[2];
                                        sunpos = Solar.sunposRAD(time);
                                    
                                        sun[0] = sunpos[0];
                                        sun[1] = sunpos[1];
                                        
                                        sat[0] = SCpos_RAD(0); // [rad]
                                        sat[1] = SCpos_lonlath(1); // [rad]
                                        sat[2] = SCpos_lonlath(2)/1e3; // [km]
                                        
                                        //static double amjd, f10, s10, y10, f10b, s10b, y10b, m10, m10b, dstdtc;
                                        UTCdate = GPS2UTCdate(time, "ISOD");
                                        doy = UTCdate(1);
                                        
                                        amjd = GPS2MJD(time);
                                        D1950 = amjd - 33281.0;
                                        // Use 1 day lag for f10 and s10 for jb2008
                                        T1950 = D1950;
                                        
                                        // Read again indices in case a new day has started
                                        if(doy > init_doy)
                                          {
                                          solfsmy_(&T1950, &f10, &f10b, &s10, &s10b, &m10, &m10b, &y10, &y10b, c_solfsmy_file);
                                          dtcval_(&T1950, &dstdtc, c_dtc_file, c_dst_file);
                                          
                                          //cout << doy << endl;
                                          //cout << fixed << f10 << "," << f10b << "," << s10 << "," << s10b << "," << m10 << "," << m10b << "," << y10 << "," << y10b << endl;
                                          //cout << "dstdtc: " << dstdtc << endl;
                                          //cout << fixed << "Sunpos C++ = " << sunpos(0) << "," << sunpos(1) << "," << sunpos(2) << endl;
                                          //cout << fixed << "Sunpos Fortran = " << sun[0] << "," << sun[1] << endl;
                                          init_doy = doy; // Update init_doy value
                                          }
                                        
                                        jb2008_(&amjd, sun, sat, &f10, &f10b, &s10, &s10b, &m10, &m10b, &y10, &y10b, &dstdtc, jb2008_temperatures, &rho_atm);
                                        
                                        //cout << setprecision(15) << "rho_atm = " << rho_atm << "\n" << endl;
                                        //cout << "rho_atm = " << rho_atm << "\n" << endl;
                                        
                                        //cout << fixed << "Altitude: " << sat[2] << endl;
                                        //cout << fixed << "Temp: " << jb2008_temperatures[0] << "," << jb2008_temperatures[1] << endl;
                                        }
                                        
                                    if( modelname.compare("NRLMSISE-00") == 0 )
                                        {
                                        Vec3d LST = Vec3d::Zero();
                                        static float ap[7];
                                        static float alt, glat, glong, stl, f107a, f107, UTsecs;
                                        int ap_now, doy, massnum;
                                        
                                        glong = mathconst::RAD2DEG*SCpos_lonlath(0); // [deg]
                                        glat = mathconst::RAD2DEG*SCpos_lonlath(1); // [deg]
                                        alt = SCpos_lonlath(2)/1E3; // [km]
                                        
                                        LST = GPS2LST(time, SCpos_lonlath(0));
                                        stl = LST(0);
                                    
                                        if(!idx_locked)
                                            {
                                            while(  ( time - SpaceWeather_idx(swind,0) ) > 1.0*timescales::JULIAN_DAY )
                                                {
                                                swind++;
                                                if( swind > (SpaceWeather_idx.rows() + 1) )
                                                    {
                                                    cerr << "Current epoch is more than 24 hours larger than the last available in space weather indices file" << endl;
                                                    //return(0);
                                                    return;
                                                    }
                                                }
                                                
                                            //SpaceWeather_idx_row = SpaceWeather_idx.row(swind);
                                            idx_locked = true;
                                            }
                                        else if( ( time - SpaceWeather_idx(swind,0) ) > 1.0*timescales::JULIAN_DAY )
                                                {
                                                swind++;
                                                //SpaceWeather_idx_row = SpaceWeather_idx.row(swind);    
                                                }
                                        else if( swind > (SpaceWeather_idx.rows() + 1) )
                                                {
                                                cerr << "Current epoch is more than 24 hours larger than the last available in space weather indices file" << endl;
                                                //return(0);
                                                return;
                                                }
                                                
                                        doy = (int)SpaceWeather_idx(swind,1);
                                        UTsecs = time - SpaceWeather_idx(swind,0); // Seconds of day
                                        f107a = SpaceWeather_idx(swind,27);
                                        f107 = SpaceWeather_idx(swind - 1, 29);
                                           
                                        Ap_4days.segment(0,3) = SpaceWeather_idx.block<1,3>(swind - 3,18);
                                        Ap_4days.segment(3,8) = SpaceWeather_idx.block<1,8>(swind - 2,13);
                                        Ap_4days.segment(11,8) = SpaceWeather_idx.block<1,8>(swind - 1,13);
                                        Ap_4days.segment(19,8) = SpaceWeather_idx.block<1,8>(swind,13);
                                        //cout << swind << endl;
                                        ap_now = 19 + floor(stl/3);
                                        ap[0] = SpaceWeather_idx(swind,21);                     // DAILY AP
                                        ap[1] = Ap_4days(ap_now);                               // 3 HR AP INDEX FOR CURRENT TIME
                                        ap[2] = Ap_4days(ap_now - 1);                           // 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
                                        ap[3] = Ap_4days(ap_now - 2);                           // 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
                                        ap[4] = Ap_4days(ap_now - 3);                           // 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
                                        ap[5] = ( Ap_4days.segment(ap_now - 11,8).sum() )/8.0;  // AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR TO CURRENT TIME       
                                        ap[6] = ( Ap_4days.segment(ap_now - 19,8).sum() )/8.0;  // AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR TO CURRENT TIME
                                        
                                        massnum = 48;
                                        //int UTsecs_int = (int)UTsecs;
                                        
                                        //cout << doy << "  " << UTsecs << "  " << alt << "  " << glat << "  " << glong << "  " << stl << endl;
                                        //cout << f107a << "  " << f107 << endl;
                                        //cout << ap[0] << "  " << ap[1] << "  " << ap[2] << "  " << ap[3] << "  " << ap[4] << "  " << ap[5] << "  " << ap[6] << "  " << endl;
                                        
                                        //gtd7d_(&doy, &UTsecs, &alt, &glat, &glong, &stl, &f107a, &f107, ap, &massnum, nrlmsise00_rho_vec, jb2008_temperatures);
                                        //cout << doy << endl;
                                        gtd7d_(&doy, &UTsecs, &alt, &glat, &glong, &stl, &f107a, &f107, ap, &massnum, nrlmsise00_rho_vec, nrlmsise00_temperatures);
                                        
                                        rho_atm = nrlmsise00_rho_vec[5]*1E3; // [kg/m3]
                                        //cout << setprecision(20) << "nrlmsise00_rho_vec[5] = " << nrlmsise00_rho_vec[5] << endl;
                                        //cout << setprecision(20) << "rho_atm = " << rho_atm << endl;
                                        //cout << nrlmsise00_rho_vec[0] << "  " << nrlmsise00_rho_vec[1] << "  " << nrlmsise00_rho_vec[2] << "  " << nrlmsise00_rho_vec[3] << "  " << nrlmsise00_rho_vec[4] << "  " << nrlmsise00_rho_vec[5] << "  " << nrlmsise00_rho_vec[6] << "  " << nrlmsise00_rho_vec[7] << "  " << nrlmsise00_rho_vec[8] << "  " << endl;
                                        //cout << nrlmsise00_rho_vec[0] << "  " << nrlmsise00_rho_vec[1] << "  " << nrlmsise00_rho_vec[2] << "  " << nrlmsise00_rho_vec[3] << "  " << nrlmsise00_rho_vec[4] << "  " << nrlmsise00_rho_vec[5] << "  " << nrlmsise00_rho_vec[6] << "  " << nrlmsise00_rho_vec[7] << "  " << endl;
                                        //cout << jb2008_temperatures[0] << "  " << jb2008_temperatures[1] << endl;
                                        }
                                        
                                    if( modelname.compare("Harris-Priester") == 0 )
                                        {
                                        double HP_alt_MIN, HP_alt_MAX, RA_lag, HP_prm;
                                        const int N_Coeff = 50;
                                        
                                        VectorNd<N_Coeff> Altitudes, rho_MIN, rho_MAX;
                                            
                                        /////// Get harcoded (HarrisPriester.h) Harris-Priester model ///////    
                                        get_HP(HP_alt_MIN, HP_alt_MAX, RA_lag, HP_prm, Altitudes, rho_MIN, rho_MAX);
                                        
                                        /////// Compute atmospheric density ///////
                                        
                                        // Spacecraft altitude
                                        double alt;
                                        // Sun declination, right asc.
                                        double Sun_dec, Sun_RA, c_dec;
                                        // Harris-Priester modification
                                        double c_psi2;
                                        // Altitude and density parameters
                                        double alt_MIN, alt_MAX, d_MIN, d_MAX;
                                        // Sun position
                                        Vec3d  sunpos;
                                        // Unit vector bulge_u towards the apex of the diurnal bulge ECI
                                        Vec3d bulge_u;
                                  
                                        // Spacecraft altitude
                                        alt = SCpos_lonlath(2)/1E3; // [km]
                                  
                                        if( alt <= HP_alt_MIN || alt >= HP_alt_MAX ) 
                                            {
                                            cerr << "Spacecraft altitude is outside Harris-Priester model altitude limits" << endl;
                                            return;
                                            }
                                  
                                        // Sun right ascension and declination
                                        sunpos = Solar.sunposRAD(time);
                                        Sun_RA  = sunpos(0);
                                        Sun_dec = sunpos(1);
                                  
                                        // Unit vector bulge_u towards the apex of the diurnal bulge ECI
                                        c_dec = cos(Sun_dec);
                                        bulge_u(0) = c_dec*cos(Sun_RA + RA_lag);
                                        bulge_u(1) = c_dec*sin(Sun_RA + RA_lag);
                                        bulge_u(2) = sin(Sun_dec);
                                  
                                        // Cosine of half angle between satellite position vector and apex of diurnal bulge                                  
                                        c_psi2 = 0.5 + 0.5*posECI.dot(bulge_u)/posECI.norm();
                                  
                                        // Height index search and exponential density interpolation
                                        int ih = 0; // Section index reset
                                        for (int i=0; i < N_Coeff - 1; i++) // Loop over N_Coeff altitude regimes
                                            {
                                            if( alt >= Altitudes(i) && alt < Altitudes(i+1) ) 
                                                {
                                                ih = i; // ih identifies altitude section
                                                break;
                                                }
                                            }
                                  
                                        alt_MIN = ( Altitudes(ih) - Altitudes(ih+1) )/log( rho_MIN(ih+1)/rho_MIN(ih) );
                                        alt_MAX = ( Altitudes(ih) - Altitudes(ih+1) )/log( rho_MAX(ih+1)/rho_MAX(ih) );
                                  
                                        d_MIN = rho_MIN(ih)*exp( (Altitudes(ih) - alt)/alt_MIN );
                                        d_MAX = rho_MAX(ih)*exp( (Altitudes(ih) - alt)/alt_MAX );
                                  
                                        // Density computation
                                        rho_atm = ( d_MIN + (d_MAX - d_MIN)*pow(c_psi2,HP_prm) )*1.0e-12; // [kg/m3]   
                                        }
                                    
                                    //cout << scientific <<setprecision(10) << rho_atm << endl;
                                    //cout << nrlmsise00_temperatures[0] << "  " << nrlmsise00_temperatures[1] << endl;
                                    
                                    //return(rho_atm);    
                                    };
    //------------------------------------------------------------------------------
    // Method void SetSurfaceParameters(double A, double C_D, Vec3d& n, const string mat)
    //------------------------------------------------------------------------------
    /**
     * Set up surface parameters and orientation required to compute the solar radiation pressure
     *
     * @param A       Surface area [m^2]
     * @param C_D     Solar radiation pressure coefficient
     * @param n       3-dimensional unit vector normal to the surface (ECI)
     * @param mat     Material composing the surface (@see file Constants.h for otpions)
     *
     */
    //------------------------------------------------------------------------------
    void ATMO::SetSurfaceParameters(SC::Face face_prms,
                                    Vec4d& q_attitude)
                                    {
                                    Area = face_prms.Area;
                                    Vec3d n_body;
                                    n_body = face_prms.n;
                                    
                                    Vec4d q_attitude_inv;
                                    q_attitude_inv = q_inv(q_attitude); // Attitude state inverted quaternion
                                    n = TransbyQ(n_body, q_attitude_inv); // Surface normal in ECI frame
                                    
                                    string mat = face_prms.Material;
                                    
                                    Vec3d aer_coefficients;
                                    aer_coefficients = materials::SC_aerodynamics; // See Constants.h
                                    
                                    rho1 = aer_coefficients(0); // Tangential momentum exchange coefficien
                                    rho2 = aer_coefficients(1); // Normal momentum exchange coefficient
                                    rho3 = aer_coefficients(2); // Mlecular speed ratio
                                    };
    //------------------------------------------------------------------------------
    // Method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate)
    //------------------------------------------------------------------------------
    /**
     * Implementation of the abstract method of class SPACEENV
     * Compute the force applied on a surface by the solar radiation pressure using a multiple-surfaces model
     *
     * @param time        GPS epoch (seconds) of the input state
     * @param orbstate    Spacecraft state vector (ECI)
     *
     * @return 3-dimensional solar radiation pressure force on the surface
     * considered
     *
     * @note In this function coefficient CF (@see method SetSurfaceParameters)
     * represents a scale factor to be estimated and should have a nominal value of 1
     */
    //------------------------------------------------------------------------------
    Vec3d ATMO::field_vec(double time,
                          const Ref<const VectorXd>& orbstate)
                            {
                            Vec3d F_ATM = Vec3d::Zero();
                            Vec3d v = Vec3d::Zero();
                            Vec3d u_rel = Vec3d::Zero();
                            Vec3d posECI = Vec3d::Zero();
                            Vec3d velECI = Vec3d::Zero();
                            double v_norm;
                            bool force_on = false;
                            
                            posECI = orbstate.segment(0,3);
                            velECI = orbstate.segment(3,3);
                            
                            v = velECI.normalized();// v is the spacecraft velocity unit vector
                            v_norm = velECI.norm();
                            // Compute atmospheric density
                            //rho_atm = AtmosphericDensity(time, orbstate);
                            AtmosphericDensity(time, orbstate);
                            
                            if( Drag_Model.compare("Panels") == 0 )
                                {
                                u_rel = -v;
                                
                                if( u_rel.dot(n) <= 0.0 ) force_on = true; // The surface is exposed to the flow
                                //cout << "u_rel.dot(n) = " << u_rel.dot(n) << " force_on = " << force_on << endl;
                                
                                //if(force_on) cout << "Area*( v.dot(n) ) = " << Area*( v.dot(n) ) << endl;
                                if(force_on) F_ATM = -CF*rho_atm*v_norm*v_norm*Area*( v.dot(n) )*( rho1*v + ( rho2*rho3 + (2.0 - rho1 - rho2)*(v.dot(n)) )*n );
                                //cout.precision(20);
                                //if(force_on) cout << "F_ATM = " << F_ATM << endl;
                                //cout << "F_ATM = " << F_ATM << endl;
                                //if( fabs(F_ATM(0)) > 1e-7 || fabs(F_ATM(1)) > 1e-7 || fabs(F_ATM(2)) > 1e-7) cout << "F_ATM = " << F_ATM << endl;
                                //cout << "rho1 = " << rho1 << " rho2 = " << rho2 << " rho3 = " << rho3 << endl;
                                }
                            else if( Drag_Model.compare("RefArea") == 0 )
                                {
                                F_ATM = -(1.0/2.0)*CF*rho_atm*v_norm*v_norm*Area*v;
                                }
                            
                            return(F_ATM);
                            };
    //------------------------------------------------------------------------------
    // Method getmodel_coeff()
    //------------------------------------------------------------------------------
    /**
     * Get surface optical coefficients
     *
     * @Read space weather indices file and put the content in matrix SpaceWeather_idx
     */
    //------------------------------------------------------------------------------
    void ATMO::getmodel_coeff()
                            {
                            if( modelname.compare("JB2008") == 0 )
                                {    
                                Vector6d UTCdate = Vector6d::Zero();
                                init_doy = 0.0;
                                
                                string dtc_file = modelfilepath + "/atmosphere/JB2008/DTCFILE.TXT";
                                string dst_file = modelfilepath + "/atmosphere/JB2008/DSTFILE.TXT";
                                string solfsmy_file = modelfilepath + "/atmosphere/JB2008/SOLFSMY.TXT";
                                
                                //char c_dtc_file[dtc_file.size() + 1];
                                strcpy(c_dtc_file, dtc_file.c_str());
                                
                                //char c_dst_file[dst_file.size() + 1];
                                strcpy(c_dst_file, dst_file.c_str());
                                
                                //char c_solfsmy_file[solfsmy_file.size() + 1];
                                strcpy(c_solfsmy_file, solfsmy_file.c_str());
                                
                                UTCdate = GPS2UTCdate(init_epoch, "ISOD");
                                init_doy = UTCdate(1);
                                
                                double amjd, T1950, D1950;
                                        
                                amjd = GPS2MJD(init_epoch);
                                D1950 = amjd - 33281.0;
                                // Use 1 day lag for f10 and s10 for jb2008
                                T1950 = D1950;
                                
                                solfsmy_(&T1950, &f10, &f10b, &s10, &s10b, &m10, &m10b, &y10, &y10b, c_solfsmy_file);
                                dtcval_(&T1950, &dstdtc, c_dtc_file, c_dst_file);
                                
                                //cout << init_doy << endl;
                                //cout << fixed << f10 << "," << f10b << "," << s10 << "," << s10b << "," << m10 << "," << m10b << "," << y10 << "," << y10b << endl;
                                //cout << "dstdtc: " << dstdtc << endl;
                                }  
                              
                            if( modelname.compare("NRLMSISE-00") == 0 )
                                {
                                string spaceweather_file = modelfilepath + "/spaceweather/CssiSpaceWeather_indices.txt";
                                SpaceWeather_idx = read_SpaceWeather(spaceweather_file.c_str(),init_epoch, simduration + 86400); // inittime variable of class PROP
                                }
                            
                            };
    

}; // End of namespace spaceenvironment



