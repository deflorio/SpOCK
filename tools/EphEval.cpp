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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
//#include <windows.h>

#include <typeinfo>
#include <thread>
#include <string>

#include <Constants.h>
#include <VarTypes.h>
#include <Transformations.h>
#include <IO_utils.h>

// External libraries: Eigen
#include <Eigen/Core>

using namespace std;
using namespace constants;

using namespace boost;

Vector6d Compute_AN(VectorNd<2>& time_vec, VectorNd<2>& u_vec, VectorNd<2>& a_vec, VectorNd<2>& ex_vec, VectorNd<2>& ey_vec, VectorNd<2>& inc_vec, VectorNd<2>& Om_vec, bool& valid);

int main(int argc, char *argv[])
    {
    if(argc != 2)
      {
      // Tell the user how to run the program
      cerr << "Usage: " << argv[0] << " Path to file       OrbitEphemeris.csv" << endl;
      return 1;
      }
      
    ////////////////////////////// Output files //////////////////////////////
    string Orbel1_file_name, MeanOrbel_AN_file_name;
    
    Orbel1_file_name = "../output/OrbitEphemerides/Orbel.csv";
    ofstream Orbel1_file;
    Orbel1_file.open(Orbel1_file_name);
    
    MeanOrbel_AN_file_name = "../output/OrbitEphemerides/MeanOrbel_AN.csv";
    ofstream MeanOrbel_AN_file;
    MeanOrbel_AN_file.open(MeanOrbel_AN_file_name);
    
    ////////////////////////////// Variables //////////////////////////////
    string OrbEph1_file(argv[1]);
    
    MatrixXd loaded_OrbEph1;
    Vector6d ECIstate1, oscorbel, meanorbel1;
    int OrbEph1_rows;
    double timeGPS;
    VectorNi<2> IntplInt;
    
    bool header = false;
    bool valid = true;
    const int Eph_cols = 14;
    
    ////////////////////////////// Load orbit ephemeris 1 //////////////////////////////
    loaded_OrbEph1 = read_csvfile(OrbEph1_file.c_str(), Eph_cols, header);
    OrbEph1_rows = loaded_OrbEph1.rows();
    
    ////////////////////////////// Interpolation of ECI orbit state in ephemeris 1 on times of ephemeris 2 and computation of orbital elements //////////////////////////////
    //VectorNd<13> Orbel1 = VectorNd<13>::Zero();
    
    Vector6d MeanOrbel_AN = Vector6d::Zero();
    //int AN_idx = 0;
    
    double t_1 = 0.0;
    double a = 0.0;  // Previous step's osculating semi-major axis
    double ex = 0.0; // Previous step's osculating semi-major axis
    double ey = 0.0; // Previous step's osculating semi-major axis
    double inc = 0.0;    // Previous step's osculating semi-major axis
    double Om = 0.0; // Previous step's osculating semi-major axis
    double u = 0.0;  // Previous step's argument of latitude
    double a_1 = 0.0;    // Previous step's osculating semi-major axis
    double ex_1 = 0.0;   // Previous step's osculating semi-major axis
    double ey_1 = 0.0;   // Previous step's osculating semi-major axis
    double inc_1 = 0.0;  // Previous step's osculating semi-major axis
    double Om_1 = 0.0;   // Previous step's osculating semi-major axis
    double u_1 = 0.0;    // Previous step's argument of latitude
    
    VectorNd<2> time_vec = VectorNd<2>::Zero();
    VectorNd<2> u_vec = VectorNd<2>::Zero();
    VectorNd<2> a_vec = VectorNd<2>::Zero();
    VectorNd<2> ex_vec = VectorNd<2>::Zero();
    VectorNd<2> ey_vec = VectorNd<2>::Zero();
    VectorNd<2> inc_vec = VectorNd<2>::Zero();
    VectorNd<2> Om_vec = VectorNd<2>::Zero();
    
    
    
    for(int i = 0; i < OrbEph1_rows; i++)
        {
        timeGPS = loaded_OrbEph1(i,1);
        
        // Compute osculating and mean orbital elements of orbit 2
        ECIstate1 = loaded_OrbEph1.row(i).segment(2,6);
        
        // Compute osculating and mean orbital elements of orbit 1
        oscorbel = rv2oe(ECIstate1, valid);
        meanorbel1 = osc2mean(oscorbel, valid);
        
        //Orbel1(0) = timeGPS;
        //Orbel1.segment(1,6) = oscorbel;
        //Orbel1.segment(7,6) = meanorbel1;
        
        // Orbital elements ephemeris to csvfile
        Orbel1_file << fixed << timeGPS << "," << oscorbel(0) << "," << oscorbel(1) << "," << oscorbel(2) << "," << oscorbel(3) << "," << oscorbel(4) << "," << oscorbel(5) << "," << meanorbel1(0) << "," << meanorbel1(1) << "," << meanorbel1(2) << "," << meanorbel1(3) << "," << meanorbel1(4) << "," << meanorbel1(5) << endl;
        
        // Compute mean elements at ascending node (AN)
        a = meanorbel1(0);
        ex = meanorbel1(1);
        ey = meanorbel1(2);
        inc = meanorbel1(3);
        Om = meanorbel1(4);
        u = mod(meanorbel1(5),mathconst::PI2);
        
        time_vec << timeGPS, t_1;
        a_vec << a, a_1;
        ex_vec << ex, ex_1;
        ey_vec << ey, ey_1;
        inc_vec << inc, inc_1;
        Om_vec << Om, Om_1;
        u_vec << u, u_1;
        
        //////////// Detect ascending node ////////////
        if(u_1 > u)
            {
            MeanOrbel_AN = Compute_AN(time_vec, u_vec, a_vec, ex_vec, ey_vec, inc_vec, Om_vec,valid);
            
            // Mean orbital elements at ascending node to csv file
            MeanOrbel_AN_file << fixed << MeanOrbel_AN(0) << "," << MeanOrbel_AN(1) << "," << MeanOrbel_AN(2) << "," << MeanOrbel_AN(3) << "," << MeanOrbel_AN(4) << "," << MeanOrbel_AN(5) << "," << 0.0 << endl;
            }
        
        t_1 = timeGPS;
        
        a_1 = a;
        ex_1 = ex;
        ey_1 = ey;
        inc_1 = inc;
        Om_1 = Om;
        u_1 = u;
        }
        
    ////////////////////////////// Write csv files //////////////////////////////
    //for(int i = 0; i < OrbEph1_rows; i++)
    //    {
    //    // Orbital elements ephemeris
    //    Orbel1_file << fixed << Orbel1(i,0) << "," << Orbel1(i,1) << "," << Orbel1(i,2) << "," << Orbel1(i,3) << "," << Orbel1(i,4) << "," << Orbel1(i,5) << "," << Orbel1(i,6) << "," << Orbel1(i,7) << "," << Orbel1(i,8) << "," << Orbel1(i,9) << "," << Orbel1(i,10) << "," << Orbel1(i,11) << "," << Orbel1(i,12) << endl;
    //    
    //    // Mean orbital elements at ascending node
    //    MeanOrbel_AN_file << fixed << MeanOrbel_AN(i,0) << "," << MeanOrbel_AN(i,1) << "," << MeanOrbel_AN(i,2) << "," << MeanOrbel_AN(i,3) << "," << MeanOrbel_AN(i,4) << "," << MeanOrbel_AN(i,5) << endl;
    //    }
    
    // Close files        
    Orbel1_file.close();
    MeanOrbel_AN_file.close();
    
    return(0);
    } // End of main()
    
    //------------------------------------------------------------------------------
    // Vector6d Compute_AN(VectorNd<2>& time_vec, VectorNd<2>& u_vec, VectorNd<2>& a_vec, VectorNd<2>& ex_vec, VectorNd<2>& ey_vec, VectorNd<2>& inc_vec, VectorNd<2>& Om_vec, bool& valid))
    //------------------------------------------------------------------------------
    /**
     * Compute the real longitude of the ascending node by third order interpolation based on on-board estimated state
     *
     * @return 6D vector with t_AN, a_AN, ex_AN, ey_AN, i_AN, O_AN
     * 
     */
    //------------------------------------------------------------------------------
    Vector6d Compute_AN(VectorNd<2>& time_vec,
                        VectorNd<2>& u_vec,
                        VectorNd<2>& a_vec,
                        VectorNd<2>& ex_vec,
                        VectorNd<2>& ey_vec,
                        VectorNd<2>& inc_vec,
                        VectorNd<2>& Om_vec,
                        bool& valid)
                                {
                                Vector6d AN_vec = Vector6d::Zero();
                                
                                double time = 0.0;
                                double t_1 = 0.0;
                                //double t_2 = 0.0;
                                //double t_3 = 0.0;
                                
                                double u = 0.0;
                                double u_1 = 0.0;
                                double a = 0.0;
                                double a_1 = 0.0;
                                double ex = 0.0;
                                double ex_1 = 0.0;
                                double ey = 0.0;
                                double ey_1 = 0.0;
                                double inc = 0.0;
                                double inc_1 = 0.0;
                                double Om = 0.0;
                                double Om_1 = 0.0;
                                //double lon = 0.0;
                                //double lon_1 = 0.0;
                                //double lon_2 = 0.0;
                                //double lon_3 = 0.0;
                                
                                double  dt = 0.0;          // Time between current time and time of AN [s]
                                //double  l0 = 0.0;          // Arranged long. at curr. time
                                //double  l1 = 0.0;          // Arranged long. at curr. time-1
                                //double  l2 = 0.0;          // Arranged long. at curr. time-2
                                //double  l3 = 0.0;          // Arranged long. at curr. time-3
                                
                                double t_AN = 0.0;
                                double a_AN = 0.0;
                                double ex_AN = 0.0;
                                double ey_AN = 0.0;
                                double i_AN = 0.0;
                                double O_AN = 0.0;
                                //double LAN = 0.0;
                                
                                time = time_vec(0);
                                t_1 = time_vec(1);
                                //t_2 = time_vec(2);
                                //t_3 = time_vec(3);
                                
                                u = u_vec(0);
                                u_1 = u_vec(1);
                                
                                a = a_vec(0);
                                a_1 = a_vec(1);
                                
                                ex = ex_vec(0);
                                ex_1 = ex_vec(1);
                                
                                ey = ey_vec(0);
                                ey_1 = ey_vec(1);
                                
                                inc = inc_vec(0);
                                inc_1 = inc_vec(1);
                                
                                Om = Om_vec(0);
                                Om_1 = Om_vec(1);
                                //lon = lon_vec(0);
                                //lon_1 = lon_vec(1);
                                //lon_2 = lon_vec(2);
                                //lon_3 = lon_vec(3);
                            
                                ////////////////////////// Compute current LAN (LAN) //////////////////////////
                            
                                // The current state is not necessarily at the AN exactly (argument of
                                // latitude u equal 0). Therefore it is necessary to interpolate the
                                // longitude between the previous and the current observation.
                                // Compute the time between u = 0 and the current observation
                            
                                // Check first possible division by 0
                                double u_ = u + mathconst::PI2 - u_1;
                                if( u_ == 0.0 )
                                    {
                                    valid = false;
                                    return(AN_vec);
                                    }
                            
                                dt = (u/u_)*(time - t_1);
                            
                                // Compute time of ascending node (t_AN)
                                t_AN = time - dt;
                            
                                //// 3rd order polynomial interpolation of LAN, before interpolation the longitudes are arranged for
                                //// proper interpolation (management of the case in which there is the change of value 0 - mathconst::PI2)
                                //l3 = lon_3;
                                //l2 = ( lon_2 > l3 ? lon_2 - mathconst::PI2 : lon_2 );
                                //l1 = ( lon_1 > l2 ? lon_1 - mathconst::PI2 : lon_1 );
                                //l0 = ( lon > l1 ? lon - mathconst::PI2 : lon );
                                //
                                //// Check possible division by 0 in the computation of LAN
                                //if( t_3-time == 0.0 || t_3-t_1 == 0.0 || t_3-t_2 == 0.0 || 
                                //  t_2-time == 0.0 || t_2-t_1 == 0.0 || 
                                //  t_1-time == 0.0 )
                                //    {
                                //    valid = false;
                                //    return(AN_vec);
                                //    }
                                //// Third order interpolation of the LAN
                                //LAN = (t_AN-t_2)*(t_AN-t_1)*(t_AN-time  )/(t_3-t_2)/(t_3-t_1)/(t_3-time  )*l3  
                                //  + (t_AN-t_3)*(t_AN-t_1)*(t_AN-time  )/(t_2-t_3)/(t_2-t_1)/(t_2-time  )*l2  
                                //  + (t_AN-t_3)*(t_AN-t_2)*(t_AN-time  )/(t_1-t_3)/(t_1-t_2)/(t_1-time  )*l1  
                                //  + (t_AN-t_3)*(t_AN-t_2)*(t_AN-t_1)/(time  -t_3)/(time  -t_2)/(time  -t_1)*l0; 
                                //
                                //// Conversion to range [0,mathconst::PI2[ (LAN could be smaller than 0)
                                //LAN = mod(LAN,mathconst::PI2);
                                
                                // Semi-major axis at ascending node
                                a_AN = a_1 + ( (a - a_1)/(time - t_1) )*(t_AN - t_1);
                                
                                // Eccentricity vector x-component
                                ex_AN = ex_1 + ( (ex - ex_1)/(time - t_1) )*(t_AN - t_1);
                                
                                // Eccentricity vector y-component
                                ey_AN = ey_1 + ( (ey - ey_1)/(time - t_1) )*(t_AN - t_1);
                                
                                // Inclination
                                i_AN = inc_1 + ( (inc - inc_1)/(time - t_1) )*(t_AN - t_1);
                                
                                // Right ascension of ascending node
                                if( (fabs(Om - Om_1) > mathconst::PI) && (Om < Om_1) ) Om = Om + mathconst::PI2;
                                else if( (fabs(Om - Om_1) > mathconst::PI) && (Om_1 < Om) ) Om_1 = Om_1 + mathconst::PI2;
                                    
                                O_AN = Om_1 + ( (Om - Om_1)/(time - t_1) )*(t_AN - t_1);
                                O_AN = mod(O_AN,mathconst::PI2);
                                
                                // Output
                                AN_vec << t_AN, a_AN, ex_AN, ey_AN, i_AN, O_AN;
                            
                                return(AN_vec);
                                };// End of Compte_AN
      
