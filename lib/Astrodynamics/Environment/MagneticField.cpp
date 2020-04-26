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
#include <fstream>
#include <string>
#include <math.h>

#include <MagneticField.h>
#include <Constants.h>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;

extern "C"
        {
        extern void igrf13syn_(int*, double*, int*, double*, double*, double*,double*, double*, double*, double*); // Fortran version of IGRF13 model
        
        extern void geomag_(int*, float*, char*); // Fortran version of WMM2020 model
        extern void geomg1_(float*, float*, float*, float*, float*, float*, float*, float*, float*); // Fortran version of WMM2020 model
        //INIT, MAXORD, PI, DTR, A, B, RE, A2, B2, C2, A4, B4, C4, OTIME, OALT, OLAT, C, CD, P, DP, SP, CP, FN, FM, PP, K;
        
        //struct init
        //        {
        //        int maxord;
        //        float pi;
        //        float dtr;
        //        float a;
        //        float b;
        //        float re;
        //        float a2;
        //        float b2;
        //        float c2;
        //        float a4;
        //        float b4;
        //        float c4;
        //        float otime;
        //        float oalt;
        //        float olat;
        //        float c;
        //        float cd;
        //        float p;
        //        float dp;
        //        float sp;
        //        float cp;
        //        float fn;
        //        float fm;
        //        float pp;
        //        float k;
        //            } init_;
        }

namespace magnetic
    {
    //------------------------------------------------------------------------------
    // MAGFIELD implementation
    //------------------------------------------------------------------------------
    Vec3d MAGFIELD::magfield = Vec3d::Zero();
    float MAGFIELD::init_year = 0.0;
    char MAGFIELD::c_wmm_file[200] = " ";
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    MAGFIELD::~MAGFIELD() {};
    //------------------------------------------------------------------------------
    // Method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate, const string aux)
    //------------------------------------------------------------------------------
    /**
     * Earth's magnetic field vector in ECEF or ECI frame
     *
     * @note This function is based on the free source C++ library GeographicLib
     * @see http://geographiclib.sourceforge.net
     *
     * @param time       GPS epoch (seconds) of the input state
     * @param orbstate   Spacecraft position vector (ECEF)
     *
     * @return Earth's magnetic field vector [nanotesla] in ECEF frame (default) or ECI if base
     * class' member 'refsys' is set to "ECI"
     */
    //------------------------------------------------------------------------------  
    Vec3d MAGFIELD::field_vec(double time,
                              const Ref<const VectorXd>& orbstate)
                            {
                            //Vec3d magfield = Vec3d::Zero();  
                            
                            Vector6d UTCdate = Vector6d::Zero();
                            Vec3d lonlath, SC_posECEF;
                            SC_posECEF = orbstate;
                            double year, lon, lat, h;
                            double Bx, By, Bz;//, H, F, D, I;
                              
                            // Find year
                            UTCdate = GPS2UTCdate(time, "ISOC");
                                    
                            year = UTCdate(0);
                            //double sec_J2000 = GPS2ET(time); // Seconds past 1 Jan 2000 11:58:55.816 UTC
                            //year = 2000.0 + sec_J2000/tyear_c();
                            
                            lonlath = ECEF2lonlath(SC_posECEF);
                            lon = lonlath(0);  lat = lonlath(1);  h = lonlath(2);
                            
                            if( modelname.compare("Dipole") == 0 )
                                {
                                double a, g10, g11, h11, r, r2, r5;
                                r = SC_posECEF.norm();
                                r2 = r*r;
                                r5 = r*r*r*r*r;
                                Vec3d m, B;
                                a = 6371.2E3; // [km]
                                g10 = -29496.57; // [nTesla]
                                g11 = -1586.42; // [nTesla]
                                h11 = 4944.26; // [nTesla]
                                
                                m << g11, h11, g10;
                                m = a*a*a*m;
                                
                                B = ( 3.0*( m.dot(SC_posECEF) )*SC_posECEF - r2*m )/r5;
                                
                                magfield = B;
                                }
                            
                            //if( modelname.compare("EMM") == 0 )
                            //    {
                            //    VectorNd<16> magfield_all;
                            //    magfield_all = EMM_GeoMagneticElements(time, lonlath);
                            //    
                            //    // Bx the northerly component, By the easterly component of the magnetic field and Bz the vertical (down) component of the magnetic field [nanotesla].
                            //    Bx = magfield_all(0);
                            //    By = magfield_all(1);
                            //    Bz = magfield_all(2);
                            //    
                            //    Vec3d magfield_SEZ(-Bx, By, -Bz);
                            //    
                            //    magfield = SEZ2ECEF(magfield_SEZ, lon, lat);
                            //    }
                            
                            //size_t found = modelname.find("IGRF");
                            //if (found!=std::string::npos)
                            //    {
                            //    string modelfile = "data/magneticfield/IGRF/" + modelname + ".COF";
                            //    
                            //    double igrf_magfield[14];
                            //    
                            //    Vector6d UTCdate = Vector6d::Zero();
                            //    double doy, longitude, latitude, alt, decyear;
                            //        
                            //    UTCdate = GPS2UTCdate(time, "ISOD");
                            //        
                            //    doy = UTCdate(1);
                            //        
                            //    decyear = year + doy/timescales::JULIAN_YEAR_DAYS;
                            //    
                            //    longitude = RAD2DEG*lon;   // Longitude [deg]
                            //    latitude = RAD2DEG*lat;   // Latitude [deg]
                            //    alt = h/1e3;
                            //    
                            //    //cout << longitude << "   " << latitude << endl;
                            //
                            //    igrf(modelfile.c_str(), decyear, longitude, latitude, alt, igrf_magfield);
                            //        
                            //    // Bx the northerly component, By the easterly component of the magnetic field and Bz the vertical (down) component of the magnetic field [nanotesla].
                            //    Bx = igrf_magfield[0];
                            //    By = igrf_magfield[1];
                            //    Bz = igrf_magfield[2];
                            //    
                            //    Vec3d magfield_SEZ(-Bx, By, -Bz);
                            //    
                            //    magfield = SEZ2ECEF(magfield_SEZ, lon, lat);
                            //    }
                            
                            size_t found = modelname.find("IGRF");
                            if(found!=std::string::npos)
                                {
                                //double igrf_magfield[4];
                                
                                Vector6d UTCdate = Vector6d::Zero();
                                double doy, elong, colat, alt, decyear, f;
                                int isv = 0;
                                int itype = 1;
                                    
                                UTCdate = GPS2UTCdate(time, "ISOD");
                                    
                                doy = UTCdate(1);
                                    
                                decyear = year + doy/timescales::JULIAN_YEAR_DAYS;
                                
                                elong = RAD2DEG*mod(lon,PI2);   // East longitude [deg] [0-360]
                                //latitude = RAD2DEG*lat;   // Latitude [deg]
                                colat = 90.0 - RAD2DEG*lat; // Colatitude [0-180]
                                alt = h/1E3;
                                
                                //cout << elong << "   " << latitude << endl;
                            
                                // Bx the northerly component, By the easterly component of the magnetic field and Bz the vertical (down) component of the magnetic field [nanotesla].
                                igrf13syn_(&isv,&decyear,&itype,&alt,&colat,&elong,&Bx,&By,&Bz,&f);
                                
                                Vec3d magfield_SEZ(-Bx, By, -Bz);
                                
                                magfield = SEZ2ECEF(magfield_SEZ, lon, lat);
                                }
                                
                            if( modelname.compare("WMM2020") == 0 )
                                {
                                Vector6d UTCdate = Vector6d::Zero();
                                float doy, longitude, latitude, alt, decyear, gv, Bx_f, By_f, Bz_f, year_f;
                                
                                year_f = float(year);
                                    
                                UTCdate = GPS2UTCdate(time, "ISOD");
                                    
                                doy = UTCdate(1);
                                    
                                decyear = float(year + doy/timescales::JULIAN_YEAR_DAYS);
                                
                                longitude = float(RAD2DEG*lon);   // Longitude [deg]
                                latitude = float(RAD2DEG*lat);   // Latitude [deg]
                                alt = float(h/1E3);
                                
                                if(year_f > init_year)
                                  {
                                  int maxdeg = 12;
                                  
                                  geomag_(&maxdeg, &year_f, c_wmm_file);
                                  
                                  init_year = year_f;
                                  }
                                
                                //cout << longitude << "  " << latitude << "  " << alt << endl;
                                // Bx the northerly component, By the easterly component of the magnetic field and Bz the vertical (down) component of the magnetic field [nanotesla].
                                //int maxdeg = 12;
                                //geomag_(&maxdeg, &year_f);
                                //cout << init_.maxord << endl;
                                geomg1_(&alt, &latitude, &longitude, &decyear, &Bx_f, &By_f, &Bz_f, &gv, &year_f);
                                
                                Bx = Bx_f;
                                By = By_f;
                                Bz = Bz_f;
                                
                                //cout << Bx << "  " << By << "  " << Bz << endl;
                                
                                Vec3d magfield_SEZ(-Bx, By, -Bz);
                                
                                magfield = SEZ2ECEF(magfield_SEZ, lon, lat);
                                }    
                           
                            if(refsys.compare("ECI") == 0)
                                {
                                Vec3d Bvec_ECI;
                                Bvec_ECI = v3D_transform(time, magfield, "ITRF93", "J2000");
                                magfield = Bvec_ECI;
                                }
                                
                            //cout << magfield[0] << "  " << magfield[1] << "  " << magfield[2] << endl;
                           
                            return(magfield);
                            };
                            
                            
    //------------------------------------------------------------------------------
    // Method getmodel_coeff()
    //------------------------------------------------------------------------------
    /**
     * Get magnetic field model coefficients
     *
     * @Read magnetiv field model coefficients file
     */
    //------------------------------------------------------------------------------
    void MAGFIELD::getmodel_coeff()
                                {
                                if( modelname.compare("WMM2020") == 0 )
                                    {
                                    Vector6d UTCdate = Vector6d::Zero();
                                    init_year = 0.0;
                                    int maxdeg = 12;
                                    
                                    string wmm_file = modelfilepath + "/magneticfield/WMM/" + modelname + ".COF";
                                    
                                    //char c_wmm_file[wmm_file.size() + 1];
                                    strcpy(c_wmm_file, wmm_file.c_str());
                                        
                                    UTCdate = GPS2UTCdate(init_epoch, "ISOC");
                                    
                                    init_year = float(UTCdate(0));
                                    
                                    if(init_year < 2020.0 || init_year > 2025.0)
                                        {
                                        cerr << "Epoch validity range for Earth's magnetic field model WMM2020 is years 2020 to 2025" << endl;
                                        exit(EXIT_FAILURE);
                                        }
                                    else
                                        {
                                        geomag_(&maxdeg, &init_year, c_wmm_file);
                                        } 
                                    }
                                };                        

}; // End of namespace spaceenvironment



